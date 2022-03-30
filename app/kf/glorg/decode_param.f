CTITLE DECODE_PARAM
 
      subroutine decode_param( line, indx, param_num,  num_pn, ierr, 
     .                         ref, apr, options)

      implicit none 
 
*     This routine will return the parameter number for the force
*     and equate commands.  The return is either by reading the
*     numerical value or by decoding the station name and parameter
*     type.  (Only for station position and velcoity at the moment).
*     If there is an error decoding the string it is returned in ierr.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
      
* PASSED VARIABLES
 
*   indx    - Current position in string (updated)
*   ierr    - Error status on decode (0 if OK, read error
*           - otherwise)
*   ref     - Reference site number.  Set to 0 on first call

* MOD TAH 150130: Changed param_num into an array so that mutiple sites
*   can be returned.  Limit max_eq.  Added num_pn for number of values
*   in array
*   param_num(max_eq)   - Parameter numbers decoded (0 if error)
*   num_pn  - Number of values in the param_num array
 
      integer*4 indx, ierr, ref, options
      integer*4 param_num(max_eq), num_pn

*   apr     - apriori value.  Is set with first call using ref = 0

      real*8 apr
 
*   line    - Line read from input
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   trimlen - Length of string
*   jndx    - Local positioning in string
*   jerr    - Local IOSTAT error on reading value
*   iel, jel, kel   - General element numbers from get_cmd
 
      integer*4 trimlen, jndx, jerr, i, iel, jel, kel, lel

      integer*4 iels(max_eq), ! List of station number returned when
                              ! mutiple sites
     .          len_next      ! Length of next word


*   next_word   - Next word in string.
*   types(15)   - Types of parameters for station names
*   misc(8)     - Miscellaneous parameters for translation and scale
 
      character*8 next_word, types(15), cdums, cdumg, 
     .            misc(8)

*   dapr  - Difference in apriori values
      real*8 dapr

*   iout  - Output unit number
      integer*4 iout

*   kbit  - Check bit status
      logical kbit

      common / progcon / iout
 
      data types / 'XPOS    ', 'YPOS    ', 'ZPOS    ',
     .            'XDOT    ', 'YDOT    ', 'ZDOT    ',
     .            'NPOS    ', 'EPOS    ', 'UPOS    ',
     .            'NDOT    ', 'EDOT    ', 'UDOT    ',
     .            'NLOG    ', 'ELOG    ', 'ULOG    '  /

      data misc  / 'VXTRAN' , 'VYTRAN', 'VZTRAN', 'VSCALE',
     .             'RXTRAN' , 'RYTRAN', 'RZTRAN', 'RSCALE' /

 
****  Get the next word from the string
      num_pn = 0
      ierr = 0
      call GetWord( line, next_word, indx )
      call casefold( next_word )
* MOD TAH 150130: Remove any wild cards
      if( next_word(trimlen(next_word):trimlen(next_word)).eq.'@' ) 
     .     next_word(trimlen(next_word):trimlen(next_word)) = ' '
      if( next_word(trimlen(next_word):trimlen(next_word)).eq.'*' ) 
     .     next_word(trimlen(next_word):trimlen(next_word)) = ' '
 

*     See if comment command used
      if( next_word(1:1).eq.'!' .or. next_word(1:1).eq.'#' ) then
          next_word = ' '
      end if

      if( trimlen(next_word).gt.0 ) then
*         See if we can find station name
          jndx = 1
          call get_cmd(next_word, gsite_names, gnum_sites, 
     .                     iel, jndx)

* MOD TAH 150130: If the error is duplicate names; then get the list
*         that matches
          if( iel.gt.0 ) then
* MOD TAH 152023: Make sure site is used.
              if( kbit(guse_site,iel) ) then
                 num_pn = 1
              else
                 write(*,90) next_word
 90              format('WARNING: Site ',a,
     .               ' used in equate but not in solution')
                 iel = 0
*                Get next word so no error on ndot/edot/.. etc
                 call GetWord( line, next_word, indx)
                 return
             endif
          endif
*
          if( iel.lt.0 ) then   ! Ambiguious command 
             len_next = trimlen(next_word)
             do i = 1, gnum_sites
                if( next_word(1:len_next).eq.
     .              gsite_names(i)(1:len_next)  .and.
     .              kbit(guse_site,i) ) then
                    num_pn = num_pn+1
                    iels(num_pn) = i
                endif
             enddo
             if( num_pn.gt.0 ) iel = iels(1)
          endif

*         See if we match a satellite name
          jndx = 1
          call get_cmd(next_word, gsvs_names, gnum_svs, 
     .                     kel, jndx)

*         See if matches misc name
          jndx = 1
          call get_cmd(next_word, misc,8, lel, jndx)
          
*                             ! Site name found, Get next argument
          if( iel.gt.0 ) then
              call GetWord( line, next_word, indx)
              call casefold( next_word )
              jndx = 1
              call get_cmd(next_word, types, 15, jel, jndx )
*                                 ! OK, compute parameter number
              if( Jel.gt.0 .and. Jel.le.12 ) then
*                                                 ! Set back to XYZ
                  if( jel.gt.6 ) jel = jel - 6
 
*                 Now see if position or rate and set type.
                  if( jel.gt.3 ) then
                      jel = jel - 3
                      kel = 2
                  else
                      kel = 1
                  end if
 
*                 Now get the parameter number
                  if( num_pn.eq.1 ) then
                      param_num(1) = parn_site(jel, kel, iel )
                  else
                      do i = 1, num_pn
                          param_num(i) = parn_site(jel, kel, iels(i))
                      enddo
                  endif

* MOD TAH 980929: Either set or check the apriori value.  Check
*                 if reference site number is set.
                  if( ref.eq. 0 ) then
                      do i = 1, num_pn
                         if( kbit(guse_site,iels(i)) ) then
                            apr = apr_val_site(jel, kel, iels(i))
                            ref = iels(i)
                            exit
                         endif
                      enddo
                  endif
*                 Check for consistency
                  do i = 1,num_pn 
                     dapr = apr - apr_val_site(jel, kel, iels(i))

                     if( abs(dapr).gt.1.d-4 .and. 
     .                   kbit(guse_site,iels(i)).and.
     .                   kbit(guse_site,ref)      ) then

                        call eqfixa(apr_val_site,site_epoch, 
     .                      parn_site, parm_change, 
     .                      gsite_names, jel, kel, ref,iels(i),
     .                      iout, options, gepoch_out)
                     end if
                  enddo 
              else if ( Jel.gt.0 .and. Jel.le.15 ) then
* MOD TAH 030610: Log estimates:
                  Jel = Jel - 12
                  do i = 1, num_pn
                     param_num(i) = parn_log(jel,iels(i))
                  enddo 
                  if( ref.eq.0 ) then
                      do i = 1, num_pn
                          if( kbit(guse_site,iels(i)) ) then                        
                             ref = iels(i)
                             apr = apr_val_log(jel,iels(i))
                             exit
                          endif
                      end do
                  endif 
*                 See if there is a difference
                  do i = 1,num_pn
                     dapr = apr - apr_val_log(jel,iels(i))
                     if( abs(dapr).gt.1.d-4 .and.
     .                   kbit(guse_site,iels(i)).and.
     .                   kbit(guse_site,ref)      ) then
                         call eqfixlog(jel, ref, iels(i), iout,
     .                                 options )
                     endif
                  enddo 
*                                 ! Unknown or ambiguous type
              else
                  write(*,100) next_word
100               format('** Error decoding parameter token ',a8,
     .                   '. Most likely non-existant station')
* MOD TAH 980415: Get next word on assumption that it is operation
*                 (i.e., npos, epos,....)
                  call GetWord( line, next_word, indx)

                  param_num(1) = 0
              end if
          else if ( kel.gt.0 ) then

*             Not a station name, but stallite name matches
              jndx = 1
              call GetWord( line, next_word, indx ) 
              call casefold( next_word )
              call svel_to_code(cdums, cdumg, next_word, jel, 'ATOC')
              if( jel.gt.0 ) then
                  param_num(1) = parn_svs(jel, kel )
                  num_pn = 1
              else
                  write(*,120) gsvs_names(kel), next_word
 120              format('**WARNING** No match found for ',a,
     .                    ' element type ',a)
                  param_num(1) = 0
              end if
              
          else if ( lel.gt.0 ) then
*             Match on miscellaneous parameters
              if( lel.ge. 1 .and. lel.le.3 ) 
     .                    param_num(1) = parn_tran(lel,1)
              if( lel.ge. 5 .and. lel.le.7 ) 
     .                    param_num(1) = parn_tran(lel-4,2)
              if( lel.eq. 4 ) param_num(1) = parn_scale(1)
              if( lel.eq. 8 ) param_num(1) = parn_scale(2)
              num_pn = 1
          else

*             Could not find station name or Satellite. Try a direct read.
* MOD TAH 980415: Make sure it is a number
              call check_num( next_word, jerr )
              if( jerr.eq.0 ) then
                 read(next_word,*, iostat=jerr ) param_num(1)
                 num_pn = 1
              end if
              if( jerr.ne.0 ) then
                  write(*,150) next_word
150               format('** Error decoding parameter token ',a8,
     .                   '. Either non-existant station or parameter')
                  param_num(1) = 0
                  num_pn = 0
              else
 
*                 See if value too large.
                  if( param_num(1).gt.num_glb_parn ) then
                      write(*,200) param_num(1), line(1:trimlen(line))
  200                 format(' ** ERROR: Parameter number ',i5,
     .                    ' from string:',/,a,/,
     .                    ' Too Large.  Entry being ignored')
                      param_num(1) = 0
                  end if
              end if
          end if
 
      else
*         Nothing seems to be in string.  Set param_num to zero, and
*         error to -1
          param_num(1) = 0
          ierr = -1
      end if
 
****  Thats all.
      return
      end
 
