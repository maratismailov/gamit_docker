CTITLE GET_APR_POSITIONS
 
      subroutine get_apr_positions
 
      implicit none

*     Routine to open the file containing the list of apriori site
*     and source positions and their velocities and to save these
*     values in the globk_common.
* MOD TAH 991110: Updated to read non-secular terms
* MOD TAH 141006: Updated so that EXTENDED PERIODIC terms will be ignored
*     Implemented by adding -PER to file name (case sensitive)
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

*   indx        - Indx for readline
*   is          - site or source number
*   fn          - File number
 
      integer*4 indx, is, i, fn, ierr, jerr, trimlen

*   kbit        - Check to see is bit is set

      logical kbit
      logical no_per  ! No periodic terms to be retained
      integer*4 ind_no_per   ! Index in string of end or no-per terms
 
*   cval        - Dummy characters for readline
 
      character*8 cval
 
*   buffer      - Line read from file
 
      character*120 buffer, cpbuff
 
*   pos_epoch   - Position epoch when velocities given
*   vals(9)     - Values read from the file
*   source_vals(4)  - Source position and rate converted to mas
 
      real*8 pos_epoch, vals(9), source_vals(4)

***** Clear the array that says we have updated apriori 
      do i = 1, gnum_sites
         call sbit(gapr_updated,i,0)
      end do
 
****  Loop over the apriori files that have been specified:
      do fn = 1, num_apr_files

* MOD TAH 141006: See if -PER added to apriori file name
          ind_no_per = index(glb_apr_file(fn),'-PER')
          if( ind_no_per.gt.0 ) then
               no_per = .true.
               ind_no_per = ind_no_per-1
          else
               ind_no_per = trimlen(glb_apr_file(fn))
               no_per = .false.
          end if
 
*****     Open the apriori file
          write(*,150)fn, glb_apr_file(fn)(1:ind_no_per)
 150      format('Reading Apriori File ',i3,1x,a)
 
          open(101, file=glb_apr_file(fn)(1:ind_no_per), 
     .              iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open',glb_apr_file(fn),0,
     .                  'GET_APR_POSITIONS')

*         Read through file finding site and source names
 
          do while ( ierr.eq.0 )
 
              read(101,'(a)', iostat=ierr) buffer

* MOD TAH 080724: See if reference frame name is passed
              if( (buffer(1:16).eq.'+REFERENCE_FRAME' .or.
     .             buffer(1:16).eq.'+REFERENCE FRAME')  .and.
     .             ierr.eq.0                            ) then
                  cpbuff = buffer(17:)
                  call trimlead(cpbuff)
                  reference_frame = cpbuff
              endif

              call casefold(buffer)

              if( ierr.ne.-1) then
                  call report_error('IOSTAT',ierr,'read',glb_apr_file,
     .                              0,'GET_APR_POSITIONS')
              end if
 
*             Decode if buffer is not a comment, and no error
              if( ierr.eq.0 .and. buffer(1:1).eq.' ' ) then
 
*                 See if site name
                  indx = 1
                  call get_cmd(buffer, gsite_names, gnum_sites, is,
     .                             indx )
 
*                                                     ! site name found
                  if( is.gt.0 .and. is.ne.999999 ) then
                      call multiread(buffer,indx,'R8',jerr, vals, cval,
     .                                7)

* MOD TAH 980417: Changed to jerr on multiread so that file will continue
*                 to be read.
                      if( jerr.eq.0 ) then
 
*                        Now assign values
                         call assign_ema_val8(vals, 6, apr_val_site, 
     .                                        is, 1)
 
*                        Convert epoch to JD
                         call decyrs_to_JD( vals(7), pos_epoch)
                         site_epoch(is) = pos_epoch
   
*                        Set bit to say we have updated site
                         call sbit(gapr_updated,is,1)
                      else 
                         write(*,200) gsite_names(is)
 200                     format('** New coordinates of ',a8,
     .                      ' not used due to currupt apr_file line')
                      end if
 
*                                                  ! site name found
                  end if
 
*****             See if source
                  indx = 1
                  call get_cmd(buffer, gsource_names, gnum_sources,
     .                             is, indx )
 
*                                                 ! source name found
                  if( is.gt.0 .and. is.ne.999999 ) then
                      call multiread(buffer,indx,'R8',jerr, vals, cval,
     .                                9)
*                                         ! RA, Dec and Ra dec dot
*                                     ! RA and Dec in h/d min sec
*                                     ! and epoch (decimal years)
*                                                       ! RA
                      call con_mas(vals(1), source_vals(1))
*                                                    ! convert to mas
                      source_vals(1) = source_vals(1)*15
*                                                    ! Dec
                      call con_mas(vals(4), source_vals(2))
*                                             ! RA rate mas/year
                      source_vals(3) = vals(7)
*                                             ! Dec rate mas/year
                      source_vals(4) = vals(8)
 
                      call assign_ema_val8(source_vals, 4, 
     .                                  apr_val_source, is, 1)
 
*                     Convert epoch to JD
                      call decyrs_to_JD( vals(9), pos_epoch)
                      source_epoch(is) = pos_epoch
 
                  end if

* MOD TAH 991110: See if extended non-secular term
                  indx = 1
                  call GetWord(buffer, cval, indx)
                  call casefold(cval) 
                  if( cval(1:8).eq.'EXTENDED' ) then
*                     OK, we have extended non-secular model entry
C                      print *,'Extended ',buffer(1:80)
                      if( no_per ) then   ! Check to see if PERIODIC 
                              ! before decoding
                          if ( index(buffer,'PERIODIC').eq.0 ) then
                               call decode_nonsec(buffer, indx)
                          endif
                      else   ! Regular decoding of EXTENDED entry
                          call decode_nonsec(buffer, indx)
                      endif
                  end if

*                             ! No file reading error
              end if
*                         ! Looping until EOF found
          end do
 
          close(101)

      end do 

***** Now scan the list to see if we have sites that have not
*     been updated.  If we do try to use pre-earthquake or rename
*     positions
      do i = 1, gnum_sites
         if( .not.kbit(gapr_updated,i) ) then
             call try_upd_apr(i)
         end if
      end do

***** Thats all
      return
      end
 
