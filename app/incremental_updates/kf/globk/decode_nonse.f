CTITLE DECODE_NONSEC

      subroutine decode_nonsec( buffer, indx)

      implicit none 

*     Routine to decode the non-secular entries read from the globk
*     apriori files.  The entries must be unique and any replicated
*     one are replaced by the latest ones read.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* PASSED VARIABLES
* indx    -- Position in buffer

      integer*4 indx
 
* buffer  -- buffer read from apriori file

      character*(*) buffer

* LOCAL VARIABLES
* is  -- Site number
* it  -- Non-secular term type (1-4 for offset, periodic, exponential and
*        logarithmic)
* i, j -- Loop counters
* jerr -- Error reading buffer
* date(5) -- Date read from line 
* jndx    -- Pointer in string

      integer*4 is, it, i,j, jerr, date(5), jndx, cndx

* vals -- Dummy value for read_line
* sectag -- Seconds tag for dates
* vals_read(6) -- Upto to 6 arguments read from line 

      real*8 vals, sectag, vals_read(6)

* done   -- Logical set true when finished reading line
* duplicate -- Logical set true if there are duplicate entries
*     in the nonsecular terms
* different -- Set true if new extended values differ from old 
*     values

      logical done, duplicate, different

* type -- Character string with type read from line
* cval -- Dummy character value
* fatal -- Fatal message if too many terms

      character*8 type, cval
      character*64 fatal

c      external globk_cmd_bd
    

****  OK, Start decoding.  See if we can find station name
      call get_cmd(buffer, gsite_names, gnum_sites, is,
     .                 indx )

      cndx = index(buffer,'!')
      if( cndx.eq.0 ) cndx = index(buffer,'#')
      if( cndx.ne.0 ) buffer(cndx:) = ' '

      if( is.gt.0 .and. is.ne.999999 ) then
*         OK, found site name.  Now start working through the line.
*         There can be multiple entries for each site on a line.
          done = .false.
          do while ( .not.done )
             call read_line(buffer, indx, 'CH',jerr,vals, type)
             if( type(1:1).eq.'!' .or. type(1:1).eq. '#' ) jerr = -1
             if( jerr.eq.0 ) then
*                OK, see if can match type
                 call casefold(type)
                 jndx = 1
                 call get_cmd(type, nonsec_types,num_nonsec_types,
     .                            it, jndx) 

                 if( it.gt.0  ) then 

*                    Found another non-secular entry.  Increment the 
*                    number and make not too many
                     num_nonsec = num_nonsec + 1
                     if( num_nonsec.gt. max_nonsec) then
                         write(fatal,120) max_nonsec
 120                     format('Too many nonsecular terms. ',
     .                           'Max allowed is ',i6)
                         call report_stat('FATAL','GLOBK',
     .                           'decode_nonsec',' ',fatal,0)
                     end if

*                    Save the site and type for this term
                     param_nonsec(1,num_nonsec) = is
                     param_nonsec(2,num_nonsec) = it

*                    Read the date from buffer.
                     call multiread(buffer,indx,'I4',jerr,date,cval,5)
                     sectag = 0.d0
                     call ymdhms_to_jd(date,sectag,
     .                      apr_val_nonsec(1,num_nonsec))
*                    Get the "parameter" either period or decay
*                    time.  For offset and rate value should be
*                    set to zero.
                     call read_line(buffer,indx,'R8', jerr, vals, cval)
                     apr_val_nonsec(2,num_nonsec) = vals
*                    Now get the  parameters for the term.
*                    The number of terms depends on the type.
                     if( it.eq. 1 .or. it.eq.2 ) then
*                        Offset and Rate change and periodic
*                        terms.   For offset and rate, the
*                        first argument is ignored, for 
*                        periodic it is period.
                         call multiread(buffer,indx,'R8', jerr,
     .                        vals_read,cval,6)
*                        Now convert from NEU to XYZ
                         call nonsec_convert('TOXYZ',2,vals_read,
     .                        apr_val_nonsec(3,num_nonsec),
     .                        apr_val_site(1,1,is))
                     else if( it.eq.3 .or. it.eq.4 ) then
*                        Exponential and logarithm. In this case
*                        we have just the three NE and U amplitudes
                         call multiread(buffer,indx,'R8', jerr,
     .                        vals_read,cval,3)
*                        Clear the other terms in values read
                         do i = 4,6
                            vals_read(i) = 0.d0
                         enddo
*                        Now convert from NEU to XYZ
                         call nonsec_convert('TOXYZ',1,vals_read,
     .                        apr_val_nonsec(3,num_nonsec),
     .                        apr_val_site(1,1,is))
                     endif

*                    Now check if the entry just added is unique
                     do i = 1, num_nonsec-1
                        duplicate = .true.
                        do j = 1, 2
                           if( param_nonsec(j,i).ne.
     .                         param_nonsec(j,num_nonsec) ) 
     .                                        duplicate = .false.
                           if( apr_val_nonsec(j,i).ne.
     .                         apr_val_nonsec(j,num_nonsec) )
     .                                        duplicate = .false.
                        end do
                        if( duplicate ) then
*                          Print a message if the new values do not
*                          match the old values
                           different = .false.
                           do j = 3,8
* MOD TAH 210429: Replaced .ne. with test of difference magntitude to
*                 avoid reporting small differences.
C                             if( apr_val_nonsec(j,i).ne.
C    .                            apr_val_nonsec(j,num_nonsec) )
C    .                            different = .true.
                              if( abs(apr_val_nonsec(j,num_nonsec)-
     .                            apr_val_nonsec(j,i)).gt.1e-8 )
     .                            different = .true.
                           enddo
                           
                           if( different )
     .                     write(*,180) i,gsite_names(is),type,
     .                           num_nonsec 
 180                       format('**NOTE** Non-secular term ',i5,
     .                            ' at site ',a8,' Type ',a8,
     .                            ' being replaced by # ',i5)
                           do j = 3,8
                              apr_val_nonsec(j,i) = 
     .                            apr_val_nonsec(j,num_nonsec)
                           end do
                           num_nonsec = num_nonsec - 1
                           exit
                        end if
                     end do

                     indx = indx + jndx  ! Save where we are in the line
                else
C                     call report_error('DECODE_EXTENDED TYPE',it,
C    .                   'decod', buffer,0,'DECODE_Nosec')  ! Reports ends of lines
                     jerr = -1
                     done = .true.
                end if
             else
*               Nothing is left on line so we are done
                done = .true.
             end if
          end do
      end if

****  Thats all
      return
      end

CTITLE NONSEC_CONVERT

      subroutine nonsec_convert(direct, num, vals_in, vals_out, xyz) 

      implicit none 

*     Routine to convert the non-secular position components from
*     NEU to XYZ or visa versa.  

* PASSED VARIABLES
* num  -- Lead dimenstion of the values (ie. if num=2, then the vals
*      are paired (e.g. cos and sin terms in NE and U; if 1 then a
*      column vector of NEU)

      integer*4 num

* direct -- String containing direction of conversion (TOXYZ or
*     TONEU)

      character*(*) direct

* vals_in(num,3) -- Input values to be converted
* vals_out(num,3) -- Output values to be returned
* xyz(3)  -- Cartesian coordinates of site

      real*8 vals_in(num,3), vals_out(num,3), xyz(3)

* LOCAL VARIABLES
* i,j  -- Loop variable

      integer*4 i,j

* loc_in(3), loc_out(3) -- In and out values of terms
* loc_coord(3) -- Geodetic latitude, longtiude and height
* rot_mat(3,3) -- Rotation matrix from XYZ to NEU

      real*8 loc_in(3), loc_out(3), loc_coord(3), rot_mat(3,3)


****  OK, see which direction we are going
      if( direct(1:5).eq.'TOXYZ' ) then

         do i = 1, num
            do j = 1, 3
               loc_in(j) = vals_in(i,j)
            end do

*           Convert the vector over to the other frame 
            call rotate_geod(loc_in, loc_out, 'NEU', 'XYZ', xyz, 
     .                       loc_coord, rot_mat)

*           Save the result
            do j = 1,3
               vals_out(i,j) = loc_out(j)
            end do
         end do  
      else if( direct(1:5).eq.'TONEU' ) then

         do i = 1, num
            do j = 1, 3
               loc_in(j) = vals_in(i,j)
            end do

*           Convert the vector over to the other frame 
            call rotate_geod(loc_in, loc_out, 'XYZ', 'NEU', xyz, 
     .                       loc_coord, rot_mat)
*           Save the result
            do j = 1,3
               vals_out(i,j) = loc_out(j)
            end do
         end do  
      else 
*        Error in call
         call report_stat('FATAL','GLOBK','nonsec_convert',
     .                    direct, 'Bad direction argument',0)
      end if

****  Thats all
      return
      end

CTITLE GET_NONLOG

      subroutine get_nonlog(site,val_log)

      implicit none 

*     Rutine to scan the non-secular terms and return the 3 values for
*     siye number site

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* PASSED VARIABLES
      integer*4 site   ! site number
      real*8 val_log(3)   ! NEU log apriori values used.

* LOCAL VARIABLES
      integer*4 i  ! loop counter

*     Scan over the non-secular terms to see if we have a log
*     value for this site
      do i = 1,3
          val_log(i) = 0.d0
      end do


      do i = 1, num_nonsec
         if( param_nonsec(1,i).eq.site .and. 
     .       param_nonsec(2,i).eq. 4        ) then   ! type 4 is log
             call nonsec_convert('TONEU',1,apr_val_nonsec(3,i),
     .       	  val_log, apr_val_site(1,1,site))
         end if
      end do

      return
      end
