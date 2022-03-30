      program test_RTr

      implicit none

*     Routine to read rinex file records and return the data
*     in the TrackRT data structure.  Normally this is done
*     with trackRT_saveobs routine in the realtime comm module.

      include 'trackRT.h'      ! Real-time data structures 
      include 'trackRTComm.h'  ! Common block
      
      integer*4 i,ep

      num_site = 3
      call getarg(1,rx_file(1))
      call getarg(2,rx_file(2))
      call getarg(3,rx_file(3))
      site_names(1) = 'SIT1'
      site_names(2) = 'SIT2'
      site_names(3) = 'SIT3'


      num_rtobs = 1
      rx_not_open = .true.
      ep = 0 
      do while( num_rtobs.gt.0 )
         call rx_to_obsInt
         print *,'num_rtobs ',ep, num_rtobs
         ep = ep + 1
         do i = 1, num_rtobs
            write(*,120) ep, i, RT_StatID(i),RT_satNum(i), 
     .                   RT_MJD_obs(i), RT_L1(i), RT_L2(i), 
     .                   RT_P1(i), RT_P2(i)
 120        format(i4,i3,1x,a4,1x,I3,1x,F15.5,1x,4F15.3)
         end do
      end do

      end

