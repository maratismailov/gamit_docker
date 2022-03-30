***** poly01 <weigh=yes/no>  <order=0/1> 
*****
*****

* calculates (weighted)mean, scatter and nrms 
** input y(km) s(mm)


* n_max : maximum data points 
      integer*4 n_max
      parameter (n_max = 9999)

*  y data, s sigmas, covariance-weights 
      real*8 y(n_max),s(n_max),w(n_max)
 
* scatterfrom the mean
      real*8 scatter(n_max)
* loops 
      integer*4 i
* number of observations
      integer*4 n_obs 
*     account for weigths 
      character*10 weigh,fake
* intermediate parameters  sumyw,sumw,sums
      real*8 sumyw,sumw,sums

* (weighted) mean
      real*8 wmean,unc_mean

*     chisquare
      real*8 chisqr
* normalized rms , weighted rms 
      real*8 nrms , wrms

* number of command line arguments 
      integer*4 n_arg,iargc,order
***** poly01 <weigh=yes/no>  <order=0/1> 
*****
*****

      n_arg = iargc()
      weigh = "no"
      order = 0
      fake  = "no"
      if ( n_arg .ge. 1 )  call  rcpar(1,weigh)

c      write(20,*)  n_arg,weigh

      do i=1,n_max
      w(i) = 1.d0
      if(weigh(1:3).eq."yes") then 
        read(*,*, end = 900 ) y(i),s(i)
        if(s(i).ne.0.d0) w(i)=1.d0/(s(i)*s(i))
      else 
        s(i) = 1.
        read(*,*, end = 900 ) y(i)
      endif
      enddo
 900  n_obs = i -1



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
      unc_mean     = sqrt ( sums )

      chisqr = 0.
      do i=1,n_obs
        scatter(i) =  y(i) - wmean
        chisqr = chisqr +  scatter(i)*scatter(i) / s(i)/ s(i) 
      enddo

      if ( n_obs .gt. 1 ) then
        nrms = sqrt ( chisqr / (n_obs - 1 ) )
        wrms = sqrt ( ( n_obs /  (n_obs - 1 ) ) *  chisqr / sumw )
      else 
          nrms = 999.99
          wrms = 999.99
      endif 


      if (weigh(1:3).eq."yes") then
        write(*,600) wmean,unc_mean,n_obs, nrms,wrms 
 600    format("  wMEAN:",2x,2f20.3,
     >       "      obs:",2x,i5,
     >       "     nrms:",2x,f10.3,
     >       "     wrms:",2x,f10.3)
      else
        write(*,601) wmean,unc_mean, n_obs, nrms,wrms
 601    format("  wMEAN:",2x,2f20.3,
     >       "      obs:",2x,i5,
     >       "     nrms:",2x,f10.3,
     >       "     wrms:",2x,f10.3)
      endif
      do i=1,n_obs
        write(*,610) i, scatter(i) 
 610    format(1x,i5,5x, f20.3)
      enddo
      
      
      stop
      end
