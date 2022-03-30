***** poly01 <weigh=yes/no>  <order=0/1>   <fake=yes/no>
*****
*****

****** comparison with svdlin. estimates are ok.
***** Uncertainties are not. Murray estimates uncertainties as follows:
***** The uncertainties are obtained from the squareroot of (multiplying the covariance-relevant element) with the chi**2)
***** dof a : N_obs - 1
***** dof b : N_obs - 2


* calculates (weighted)mean and fits  a 1st order polynomial(y=a+bx)
** input x(km) y(mm) s(mm)


* n_max : maximum data points
      integer*4 n_max
      parameter (n_max = 9999)

* x data, y data, s sigmas, covariance-weights
      real*8 x(n_max),y(n_max),s(n_max),w(n_max)

* loops
      integer*4 i
* number of observations
      integer*4 n_obs
*     account for weigths
      character*10 weigh
* intermediate parameters  sumyw,sumw,sums,sumxw,sumxyw,sumxxw,det
      real*8 sumyw,sumw,sums,sumxw,sumxyw,sumxxw,det
* arbitrary x-values for nicer plots
      real*8 fake_x
* (weighted) mean
      real*8 wmean,unc_mean
* a  + b X     constant :a     slope : b
      real*8 constant,slope,unc_constant,unc_slope
* order  0: 0th order polynomial
* order  1: 1th order polynomial
      character*1 order
*     chisquare
      real*8 chisqr
* normalized rms, weighted rms
      real*8 nrms, wrms
* write fake start and stop x
      character*10 fake
* number of command line arguments
      integer*4 n_arg, iargc
***** poly01 <weigh=yes/no>  <order=0/1>   <fake=yes/no>
*****
*****

      n_arg = iargc()
      weigh = "no"
      order = "0"
      fake  = "no"
      if ( n_arg .ge. 1 )  call  rcpar(1,weigh)
      if ( n_arg .ge. 2 )  call  rcpar(2,order)
      if ( n_arg .ge. 3 )  call  rcpar(3,fake)
c      write(20,*)  n_arg,weigh,order,fake

      do i=1,n_max
      w(i) = 1.d0
      s(i) = 1.d0
      if(weigh(1:3).eq."yes") then
        read(*,*, end = 900 ) x(i),y(i),s(i)
        if(s(i).ne.0.d0) w(i)=1.d0/(s(i)*s(i))
      else
        read(*,*, end = 900 ) x(i),y(i)
      endif
      enddo
 900  n_obs = i -1

c      write(*,'(i5,3f10.3)') (i, x(i),y(i),w(i),i=1,n_obs)

****** WMean **********
*      SUM y*w
*      SUM   w
      sumyw = 0.d0
      sumw  = 0.d0
      sums  = 0.d0
      do i=1,n_obs
        sumyw  = sumyw + y(i) * w(i)
        sumw   = sumw  +        w(i)
        sums   = sums  + s(i)*s(i)
      enddo
      wmean = sumyw / sumw
c where did the original theory for unc mean come from?
c current unc_mean from King/McClusky 950904.
c      unc_mean     = sqrt ( sums )
       unc_mean     = sqrt ( 1/sumw )

      if ( order .eq. "1" ) then
******* y = A + B*x
      sumxw  = 0.d0
      sumxyw = 0.d0
      sumxxw = 0.d0
      do i=1,n_obs
        sumxxw  = sumxxw + x(i)*x(i)*w(i)
        sumxyw  = sumxyw + x(i)*y(i)*w(i)
        sumxw   = sumxw  + x(i)     *w(i)
      enddo
      det = sumw* sumxxw - sumxw * sumxw
      constant = ( sumxxw * sumyw - sumxw * sumxyw ) / det
      slope    = (-sumxw  * sumyw + sumw  * sumxyw ) / det
      unc_constant = sqrt ( sumxxw / det )
      unc_slope    = sqrt ( sumw / det )
      endif

      chisqr  = 0.0

* order 0 :
      if ( order .eq. "0" ) then
        do i=1,n_obs
          chisqr = chisqr + (( y(i) - wmean ) / s(i))**2
        enddo

        if ( n_obs .gt. 1 ) then
          nrms = sqrt ( chisqr / (n_obs - 1. ) )
          wrms = sqrt ( ( n_obs /  (n_obs - 1. ) ) *  chisqr / sumw )
        else
          nrms = 999.99
          wrms = 999.99
        endif

      endif

      if ( order .eq. "1" ) then
* order 1 :
      do i=1,n_obs
        chisqr = chisqr + (( y(i) - constant - slope * x(i)) / s(i))**2
      enddo
c      write (*,*) "n_obs ", n_obs
      if ( n_obs .gt. 2 ) then
        nrms = sqrt ( chisqr / (n_obs - 2. ) )
        wrms = sqrt ( ( n_obs /  (n_obs - 2. ) ) *  chisqr / sumw )
      else
        nrms = 999.99
        wrms = 999.99
      endif
      endif



      if (weigh(1:3).eq."yes") then
        write(*,600) wmean,unc_mean,
     >       constant,unc_constant,slope,unc_slope, n_obs, nrms,wrms
 600    format("  wMEAN:",2x,2f10.3,
     >       " constant:",2x,2f10.3,
     >       "    slope:",2x,2f20.10,
     >       "      obs:",2x,i5,
     >       "     nrms:",2x,f10.3,
     >       "     wrms:",2x,f10.3)
      else
        write(*,601) wmean,unc_mean,
     >       constant,unc_constant,slope,unc_slope, n_obs, nrms,wrms
 601    format("  wMEAN:",2x,2f10.3,
     >       " constant:",2x,2f10.3,
     >       "    slope:",2x,2f20.10,
     >       "      obs:",2x,i5,
     >       "     nrms:",2x,f10.3,
     >       "     wrms:",2x,f10.3)
      endif

      if ( order .eq. "1" ) then
***   this is for repeatability plots
        if ( x(1) .le.  x(2) ) then
          fake_x = x(1) / 2.
        else
          fake_x = x(1) * 2.
        endif
        if(fake(1:3) .eq."yes") write(*,610)fake_x,constant+slope*fake_x
*     write increasing order
        do i=1,n_obs
          write(*,610) x(i),constant+slope*x(i)
        enddo
        if ( x(n_obs) .gt.  x(n_obs-1) ) then
          fake_x = x(1) * 2.
        else
          fake_x = x(1) / 2.
        endif
        if(fake(1:3) .eq."yes") write(*,610)fake_x,constant+slope*fake_x

      else  if ( order .eq. "0" ) then
***   this is for repeatability plots
        if(fake(1:3) .eq."yes")  then
          fake_x = 00000.d0
          write(*,610)fake_x,wmean
*     write increasing order
          do i=1,n_obs
            write(*,610) x(i),wmean
          enddo
          fake_x = 100000.d0
          write(*,610)fake_x,wmean
        endif
      endif
 610  format(f20.4,5x,f20.4)


      stop
      end
