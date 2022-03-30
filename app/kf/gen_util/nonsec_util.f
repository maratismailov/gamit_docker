CTITLE DECODE_NONSUT

      subroutine decode_nonsut( buffer, indx, gsite_names,gnum_sites,
     .      ref_xyz )

      implicit none 

*     Routine to decode the non-secular entries read from the globk
*     apriori files.  The entries must be unique and any replicated
*     one are replaced by the latest ones read.

*     All sites names must already be read before calling.

      include '../gen_util/nonsec_util.h'

* PASSED VARIABLES
* indx    -- Position in buffer

      integer*4 indx
      integer*4 gnum_sites   ! Number of sites in name list

      real*8 ref_xyz(3,2,gnum_sites) ! Apriori site coordinates

      character*(*) buffer 
      character*(*) gsite_names(*)   ! Names of sites 

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
*     in the nonsutular terms
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
                 call get_cmd(type, nonsut_types,num_nonsut_types,
     .                            it, jndx) 

                 if( it.gt.0  ) then 

*                    Found another non-secular entry.  Increment the 
*                    number and make not too many
                     num_nonsut = num_nonsut + 1
                     if( num_nonsut.gt. max_nonsut) then
                         write(fatal,120) max_nonsut
 120                     format('Too many nonsutular terms. ',
     .                           'Max allowed is ',i6)
                         call report_stat('FATAL','GLOBK',
     .                           'decode_nonsut',' ',fatal,0)
                     end if

*                    Save the site and type for this term
                     param_nonsut(1,num_nonsut) = is
                     param_nonsut(2,num_nonsut) = it

*                    Read the date from buffer.
                     call multiread(buffer,indx,'I4',jerr,date,cval,5)
                     sectag = 0.d0
                     call ymdhms_to_jd(date,sectag,
     .                      apr_val_nonsut(1,num_nonsut))
*                    Get the "parameter" either period or decay
*                    time.  For offset and rate value should be
*                    set to zero.
                     call read_line(buffer,indx,'R8', jerr, vals, cval)
                     apr_val_nonsut(2,num_nonsut) = vals
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
                         call nonsut_convert('TOXYZ',2,vals_read,
     .                        apr_val_nonsut(3,num_nonsut),
     .                        ref_xyz(1,1,is))
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
                         call nonsut_convert('TOXYZ',1,vals_read,
     .                        apr_val_nonsut(3,num_nonsut),
     .                        ref_xyz(1,1,is))
                     endif

*                    Now check if the entry just added is unique
                     do i = 1, num_nonsut-1
                        duplicate = .true.
                        do j = 1, 2
                           if( param_nonsut(j,i).ne.
     .                         param_nonsut(j,num_nonsut) ) 
     .                                        duplicate = .false.
                           if( apr_val_nonsut(j,i).ne.
     .                         apr_val_nonsut(j,num_nonsut) )
     .                                        duplicate = .false.
                        end do
                        if( duplicate ) then
*                          Print a message if the new values do not
*                          match the old values
                           different = .false.
                           do j = 3,8
                              if( apr_val_nonsut(j,i).ne.
     .                            apr_val_nonsut(j,num_nonsut) )
     .                            different = .true.
                           enddo
                           
                           if( different )
     .                     write(*,180) i,gsite_names(is),type
 180                       format('**NOTE** Non-secular term ',i5,
     .                            ' at site ',a8,' Type ',a8,
     .                            ' being replaced')
                           do j = 3,8
                              apr_val_nonsut(j,i) = 
     .                            apr_val_nonsut(j,num_nonsut)
                           end do
                           num_nonsut = num_nonsut - 1
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

CTITLE NONSUT_CONVERT

      subroutine nonsut_convert(direct, num, vals_in, vals_out, xyz) 

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
         call report_stat('FATAL','GLOBK','nonsut_convert',
     .                    direct, 'Bad direction argument',0)
      end if

****  Thats all
      return
      end

CTITLE EVAL_NONSUT

      subroutine eval_nonsut(is, epoch, num_nonsut, param_nonsut,
     .                       apr_val_nonsut, dsol, tran_est )

*     Routine to evaluate the value of the non-secular contributions
*     to the apriori station position
* NOD TAH 031024: Passed tran_est into the routine.  If Bit 16 is
*     set, then log coefficients have been estimated anf therefore
*     should not be added to the station coordinates.  

      include '../includes/const_param.h'

* PASSED VARIABLES
* is  -- Site number
* num_nonsut -- Number of non-secular terms
* param_nonsut(2,num_nonsut) -- Parameters for non-secular terms
* tran_est  -- If Bit 16 is set then log terms are estimated and
*     therefore should not be added to coordinates.  Used in site_OC
*     so that we do not doubly correct the coordinates for the log term

      integer*4 is, num_nonsut, param_nonsut(2,num_nonsut), 
     .          tran_est

* epoch    -- Time at which postion is required
* apr_val_nonsut(8,num_nonsut)
* dsol(3)  -- The change in XYZ at this time

      real*8 epoch, apr_val_nonsut(8,num_nonsut), dsol(3)

* LOCAL VARIABLES
* i,j -- Loop variables

      integer*4 i, j

* dt  -- Change in time in years between solution epoch and offset
*        rate term
* arg -- Argument for the periodic, exponential or logarithm terms

      real*8 dt, arg 

* fatal -- String with any fatal error message written into it

      character*64 fatal

* kbit -- Check if bit is set
      logical kbit

****  Loop over the non-secular terms seeing which ones apply to this
*     site.
      do j = 1,3
         dsol(j) = 0.d0
      end do

      do i = 1, num_nonsut
         if( is.eq.param_nonsut(1,i) ) then
*           OK, this term applies.  Based on type see what we should
*           do.
            if ( param_nonsut(2,i).eq.1 ) then
*              offset and rate change.  See if we are past the time
*              of the change
               if( epoch.ge.apr_val_nonsut(1,i) ) then
*                  We are past the change point.  Compute the change
*                  in position
                   dt = (epoch-apr_val_nonsut(1,i))/365.25d0
                   do j = 1,3
                      dsol(j) = dsol(j) + apr_val_nonsut(2*j+1,i) +
     .                                    apr_val_nonsut(2*j+2,i)*dt
                   end do
               end if
            else if( param_nonsut(2,i).eq.2 ) then
*              Periodic term.  Get the argument and evaluate cose and sine
               arg = 2*pi*(epoch-apr_val_nonsut(1,i))/
     .                           apr_val_nonsut(2,i)
               do j = 1,3
                  dsol(j) = dsol(j) + cos(arg)*apr_val_nonsut(2*j+1,i)+
     .                                sin(arg)*apr_val_nonsut(2*j+2,i)
               end do
            else if( param_nonsut(2,i).eq.3 ) then
*              Exponential decay term
               if( epoch.ge. apr_val_nonsut(1,i)) then
                  arg = -(epoch-apr_val_nonsut(1,i))/apr_val_nonsut(2,i)
                  do j = 1,3
                     dsol(j) = dsol(j) + 
     .                        (1-exp(arg))*apr_val_nonsut(j+2,i)
                  end do
               end if
            else if( param_nonsut(2,i).eq.4 ) then
*              Logarithmic  term
               if( epoch.ge.apr_val_nonsut(1,i) .and.
     .             .not.kbit(tran_est,16) ) then
* MOD TAH 030124: Changed argument to match enfit definition
* MOD TAH 030610: Changed definition again to match tsview and enfit
* MOD TAH 031024: Added check on tran_est to see if log terms have 
*                 been estimated.
                  arg = 1+(epoch-apr_val_nonsut(1,i))/
     .                     apr_val_nonsut(2,i)
                  do j = 1,3
                     dsol(j) = dsol(j) + log(arg)*apr_val_nonsut(j+2,i)
                  end do
               end if
            else 
*              Strange term that should not exit.
               write(fatal,220) is, i, param_nonsut(2,i)
 220           format('Invalid non-secular term for site # ',i4,
     .                ' term ',i5,' Saved type ',i5)
               call report_stat('FATAL','GLOBK','eval_nonsut',' ',
     .                           fatal,0)
            end if
         end if
      end do

***** Thats all
      return
      end


CTTILE BLOCKDATA NONSUT_BD

      block data nonsut_db

      implicit none

*     Block data for extended apr types


      include '../gen_util/nonsec_util.h'

      data nonsut_types / 'OFFSET  ',
     .                    'PERIODIC',
     .                    'EXP     ',
     .                    'LOG     '  /

      end


