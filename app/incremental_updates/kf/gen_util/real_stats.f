CTITLE REAL_STATS
      subroutine real_stats(times,res,rerr,num_ts, sig_scale,
     .                      taufin)

      implicit none 

*     Routine to compute realistic sigmas from time series 
*     residuals
* MOD TAH 081231: Added feature to use negative error to denote
*     deleted data.
* MOD MAF 20211109: Increased max_av to 300 from 150

* PARAMETERS 
      integer*4 maxdt     ! Maximum number of values in one interval
                          ! (Related to max_ts in include)
      parameter ( maxdt = 3000 )   ! Time less than total length
      integer*4 ntau      ! Number of trial exponential decay times
      parameter ( ntau = 15 )
      integer*4 max_av     ! Maximum number of averaging times
      parameter ( max_av = 300 )  ! Corresponds to 2000 days/7 days 


* PASSED VARIABLES
      integer*4 num_ts    ! Number of values in time series
      real*8 times(num_ts),  ! JDs of times series values
     .       res(num_ts),    ! residuals (m)
     .       rerr(num_ts)    ! Sigmas (m).  If negative data
                             ! not be used.

      real*8 sig_scale       ! Scale factor to multiply slope
                             ! sigma by
      real*8 taufin          ! Time constant 

* LOCAL VARIABLES
      real*8 stt, ent  ! Start and stop times
      real*8 minav, maxav  ! Min and Max averaging times (days)
      integer*4 numav      ! Number averaging intervals to use

      integer*4 i,j,it        ! Loop counters
      integer*4 num        ! Number if mean values in averaging interval
      integer*4 numdt      ! Number of residuals in current time range
      real*8 t, dt         ! Current time and Current averaging interval

      real*8 resdt(maxdt), errdt(maxdt) ! Residuals and sigmas in current
                           ! interval
      real*8 summ, sumw    ! Sums for mean and weight in current interval
      real*8 wmean, varm   ! Mean and variance in current interval
      real*8 sumavs        ! Sum for chi**2 of weighted means

      real*8 chi2(max_av), tims(max_av)  ! Chi**2 and averaging times for
                           ! each averaging interval
      real*8 ef(max_av)    ! Exponential function for fitting 
      real*8 taus(ntau)    ! Trial exponential decay times
      real*8 alpha(ntau)   ! Scaling coefficients to match variances
      real*8 alsum         ! Summation variable for getting alpha
      real*8 rmsum         ! RMS fit to chi**2 behavior 
      real*8 rmmin         ! Min RMS fir to chi**2

      data taus / 1.d0, 2.d0, 4.d0, 8.d0, 16.d0, 32.d0,
     .           64.d0, 128.d0, 256.d0, 512.d0, 1024.d0,
     .           2048.d0, 4096.d0, 8192.d0, 16384.d0 / 
      

****  Get the time range (assumes time sorted data)
      stt = times(1)
      ent = times(num_ts)
      minav = 7.d0
      maxav = (ent-stt)/10.d0
      if( maxav.lt.minav*3 ) then
          minav = int(maxav/4)+1
      end if
      numav = int(maxav/minav)
* MOD TAH 190618: Check numav 
      if( numav.gt.max_av ) then
         write(*,110) numav, max_av
 110     format('NUMAV value ',i4,' too large, setting to max ',I4)
         numav = maxav
      endif

*     Now loop over the different averaging times, get the nrms of the 
*     residuals
      do i = 1, numav
          dt = i*minav
* MOD TAH 190618: Check to make bounds are not exceeded
          if( dt.gt.maxdt )then
              write(*,120) i, dt, maxdt
 120          format('REAL_STATS Averging duration exceeded. ',
     .               'Count ',i4,' DT ',I5,' maxdt ',i5)
              exit
          endif
          num = 0
          sumavs = 0.d0

c         do t = stt,ent,dt
          do it = 0, nint((ent-stt)/dt)
              t = stt + it*dt
*             get the data in this interval
              call getdata(t,t+dt,times,res,rerr,num_ts, 
     .                     resdt, errdt, numdt)
              if( numdt.gt. maxdt ) then
                  write(*,140) numdt, maxdt
 140              format('Number returned residuals ',i6,
     .                   ' > maxdt ',i6)  
                  stop 'REAL_STATS Bounds exceeded' 
              endif
              if( numdt.gt. 0 ) then
*                 Compute mean in this section
                  summ = 0.d0
                  sumw = 0.d0
                  do j = 1, numdt
                     summ = summ + resdt(j)/errdt(j)**2
                     sumw = sumw + 1.d0/errdt(j)**2
                  end do
                  wmean = summ/sumw
                  varm  = 1.d0/sumw
*                 Sum the mean and weight into stats at this averaging
                  num = num + 1
                  sumavs = sumavs +  wmean**2/varm
              end if
          end do
*         Save the chi**2 at this interval
          chi2(i) = sumavs/num
          tims(i) = dt
      end do
c      write(*,220) numav,(i,tims(i),chi2(i),i=1,numav)
c 220  format('NumAV ',i4,3(1x,i4,2e14.4))

*     OK, we now variation of chi**2 with length of averaging interval
*     Now fit to an exponential model

      rmmin = 1.d30
      do i = 1, ntau
*        Generate the error function at this tau
         do j = 1, numav
            ef(j) = 1.d0 - exp(-tims(j)/taus(i))
         end do
*        Find coefficient needed for scaling
         alsum = 0.d0
         do j = 1, numav
            alsum = alsum + chi2(j)/ef(j)
         enddo
         alpha(i) = alsum/numav
*        See how well it fits
         rmsum = 0.d0
         do j = 1, numav
            rmsum = rmsum + (chi2(j)-ef(j)*alpha(i))**2
         end do
*        If this is minimum so far save value
         if( rmsum.lt.rmmin ) then
             sig_scale = alpha(i)
             taufin = taus(i)
             rmmin = rmsum
         end if
      end do

*     Return the square root of the value
      sig_scale = sqrt(sig_scale)

****  Thats all
      return
      end

CTITLE GETDATA

      subroutine getdata(ts,te,times,res,rerr,num_ts, 
     .                     resdt, errdt, numdt)

      implicit none 

*     Routine to return the data in the interval t to t+dt
* MOD TAH 081231: Added feature to use negative error to denote
*     deleted data.

          
* PASSED VARIABLES
      integer*4 num_ts    ! Number of values in time series
      real*8 ts, te       ! Start and stop times of interval
      real*8 times(num_ts),  ! JDs of times series values
     .       res(num_ts),    ! residuals (m)
     .       rerr(num_ts)     ! Sigmas (m)
      real*8 resdt(*), errdt(*) ! Residuals and sigmas in current
      integer*4 numdt      ! Number of residuals in current time range

* LOCAL VARIABLES
      integer*4 i

*     Loop over the data and add the values to the array
      numdt = 0
      do i = 1,num_ts
         if( times(i).ge.ts .and. times(i).lt.te .and.
     .       rerr(i).gt.0  ) then
            numdt = numdt + 1
            resdt(numdt) = res(i)
            errdt(numdt) = rerr(i)
         end if
      end do

****  Thats all
      return
      end

                          
