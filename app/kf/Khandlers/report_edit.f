CTITLE REPORT_EDIT
 
      subroutine report_edit( iout )
 
      implicit none

 
*     Routine to report the various edting conditions which were
*     used during this run.
*                                 12:08 PM  WED., 11  MAR., 1987
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_names.h'
 
*   down_date(5,2)  - Times for down weighed periods
*   i,j,k   - Loop counters
*   iout    - Ouput LU number
*   lnp     - Position for writing in current line
*   num_out - Number of entries on a line.  (used to put carriage
*           - return at the ends of lines)
*   trimlen - HP function to get string length
 
      integer*4 down_date(5,2), i,j,k, iout, lnp, num_out, trimlen
 
*   elev_out    - Output minimum elevation angle (deg) NOT IMPLEMENTED
 
*     real*4 elev_out
 
*   down_epoch  - Julian date for down times
*   sectag      - seconds tag on Julian date to date conversion
 
      real*8 down_epoch, sectag
 
*   header_out  - Indicates that the title line for a section
*           - hase been written
*   kbit    - Function to check if bit turned on
 
      logical header_out, kbit
 
*   nl      - Causes a new line to started
 
      character*4 nl
 
*   buffer  - Buffer used to output some lines which need to be
*           - sequentially built up.
 
      character*80 buffer
 
      data nl / '(1x)' /
 
 
***** Start with the FRNGE quality, and maskes
 
      write(iout,100)
  100 format(/' EDIT CONDITIONS',/,
     .        ' ---------------' )
 
      write(iout,120) data_type
  120 format(' Data type        ',o6)

      write(iout,140) min_frnge_quality, data_mask, ion_mask
  140 format(' FRNGE quality    ', a,t25,
     .       ' Data Mask  ',o12,t50,
     .       ' Ion  Mask  ',o12)
 
      write(iout,160) wvr_code_limit, av_obs_dur
  160 format(' WVR code limit   ',i3,t25,
     .       ' AV_OBS_DUR (sec) ',i4)
 
      write(iout,180) n_sigma_limit, tol_close
  180 format(' Sigma limit      ',f5.1,t25,
     .       ' Closure tolerance',f5.1 )
 
      write(iout,200) correlator_noise
 200  format(' Correlator noise ',f5.1,' (ps), ',f5.1,' (fs/s)')
 
*     Now do the down weighted times
      header_out = .false.
 
      do i = 1, down_num
          if( .not. header_out ) then
              write(iout,300)
  300         format(' EDITED TIME INTERVALS',/
     .               ' BASELINE',t25,'START',t50,'STOP')
              header_out = .true.
          end if
 
*                         ! Get site names (if site number <0 then use
          do j = 1,2
*                         ! ALL)
              if( down_sites(j,i).gt.0 .and.
     .            down_sites(j,i).ne.999999 ) then
                  write(buffer((2+(j-1)*9):),310)
     .                site_names(down_sites(j,i))
 310              format(a)
              else
                  write(buffer((2+(j-1)*9):),310) 'ALL'
              end if
          end do
 
*         Now add start and stop times
          do j = 1,2
              down_epoch = start_epoch + down_times(j,i)
C             call epoc_8(down_date(2,j),down_date(3,j),down_date(1,j),
C    .                    down_date(4,j),down_date(5,j), down_epoch )
              call jd_to_ymdhms( down_epoch, down_date(1,j), sectag )
          end do
 
          write(buffer(25:),340) ((down_date(k,j),k=1,5),j=1,2)
  340     format(2(i4.4,'/',i2.2,'/',i2.2,1x,i2.2,':',i2.2,:,' -- '))
 
*         Now write the buffer
          write(iout,'(a)') buffer
      end do
 
***** Now do the sources which have been not used
 
*                     ! Clear the buffer, so that code below wont
      buffer = ' '
*                     ! output anything if no source deletes
      header_out = .false.
      num_out = 0
      do i = 1, num_sources
 
*                                          ! This source not used, add to
          if( kbit(down_source,i) ) then
*                                          ! line
 
              if( .not.header_out ) then
                  write(iout,400)
  400             format(' SOURCES NOT USED IN THE SOLUTION')
                  header_out = .true.
              end if
 
              num_out = num_out + 1
              if( mod(num_out,9).eq.0 ) then
                  write(iout,440) buffer(1:trimlen(buffer))
                  num_out = 1
                  buffer = ' '
              end if
 
              lnp = (num_out-1)*9 + 2
              write(buffer(lnp:),420) source_names(i)
  420         format(a8)
          end if
      end do
 
      if( trimlen(buffer).gt.0 ) then
          write(iout,440) buffer(1:trimlen(buffer))
  440     format(a)
      end if
 
      write(iout,nl)
 
      end
 
