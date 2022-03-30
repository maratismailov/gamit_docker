      Subroutine xyzpole( icall, vec, sig, correl ) 

c     Convert Euler vectors from xyz components to position and rate 
c     and vice versa, including covariances

c     R. King 28 Jan 97                   

c        icall = positve  :  Input is rotation vector (wx, wy, wz), output lat, long, and rate
c              = negative :  Input is lat, long, and rate, output is wx, wy, wz
c              = +/- 1       Convert Euler pole only
c              = +/- 2       Do covariances also
c
c        vec(3)  :    Input/ouput values (wx,wy,wz or lat,lon,rate)
c                
c        sig(3)  :    Input/output sigmas (wx,wy,wz or lat,lon,rate)
c
c        cov(3,3):    Input/output correlation matrix    

c      Units of vector are rad/My;  of location, degrees
                  
      implicit none
  
      integer icall,ipivot(3),i,j
                   
      real*8 vec(3),sig(3),correl(3,3),wx,wy,wz,lat,lon,rate
     .     , dx(3,3),dxt(3,3),cov(3,3),cov_out(3,3),dum3(3)
     .     , work33(3,3),pi
               
c     debug   
c      real*8 dx1(3,3),check(3,3)
                
      pi = 4.d0*atan(1.d0)
         
c    convert rotation vector to pole position and rate
c    (need for covariance in either direction)

c      if( icall.gt. 0 .or. icall.eq.-2 ) then
      if( icall.gt. 0 ) then
         wx = vec(1)
         wy = vec(2)
         wz = vec(3)
         lon = atan2(wy,wx)
         lat = atan2(wz,sqrt(wx*wx+wy*wy))
         rate = sqrt(wx*wx+wy*wy+wz*wz) 
      endif
           
c    convert pole position and rate to rotation vector

      if( icall.lt.0 ) then
        lat = vec(1)*pi/180.d0
        lon = vec(2)*pi/180.d0   
        sig(1) = sig(1)*pi/180.d0
        sig(2) = sig(2)*pi/180.d0
        rate = vec(3) 
        vec(1) = rate*cos(lat)*cos(lon)
        vec(2) = rate*cos(lat)*sin(lon)
        vec(3) = rate*sin(lat)
      endif
           
c     convert the uncertainties if requested

      if( iabs(icall).ge.2 ) then

c       convert the sigmas and correlations into covariances
        do i=1,3
          do j=1,3
            cov(i,j) = correl(i,j)*sig(i)*sig(j)
          enddo
        enddo

c       compute the Jacobian [ d(w) / d(lon,lat,rat) ]
c       dw/d(lat)
        dx(1,1) = -rate*sin(lat)*cos(lon)
        dx(2,1) = -rate*sin(lat)*sin(lon)
        dx(3,1) = rate*cos(lat)
c       dw/d(lon)
        dx(1,2) = -rate*cos(lat)*sin(lon)
        dx(2,2) =  rate*cos(lat)*cos(lon)
        dx(3,2) = 0.d0
c       dw/d(rate)
        dx(1,3) = cos(lat)*cos(lon)
        dx(2,3) = cos(lat)*sin(lon)
        dx(3,3) = sin(lat)
            
c       if vector to lat/lon/rate, need to invert
        if( icall.gt.0 ) then     
c          do i=1,3
c            do j=1,3
c              dx1(i,j) = dx(i,j)
c            enddo
c          enddo
          call gauss_elim( dx, dum3, ipivot, 3, 0, 3 )    
c             print *,'dx inv',dx
c             call matmpy(dx,dx1,check,3,3,3)
c             print *,'check ',check
        endif

c       pre- and post-multiply the covariance by the Jacobian (dx=J))
c       C(out) = J * C(in) * J(transpose)
                    
c          print *,'cov ',cov
        call matmpy(dx,cov,work33,3,3,3)
c          print *,'work33 ',work33     
        call transp(dx,dxt,3,3)
        call matmpy(work33,dxt,cov_out,3,3,3)
c          print *,'cov_out ',cov_out

c       compute the sigmas and correlations from the covariances
        do i=1,3
          sig(i) = sqrt(cov_out(i,i)) 
        enddo
        do i=1,3
          do j=1,3
            correl(i,j) = cov_out(i,j)/sig(i)/sig(j)
          enddo
        enddo
       
      endif


c       if the output is pole position, convert lat long to degrees
               
        if( icall.gt.0 ) then
          vec(1) = lat * 180.d0/pi
          vec(2) = lon * 180.d0/pi
          vec(3) = rate
          if( icall.eq.2 ) then
            sig(1) = sig(1) * 180.d0/pi
            sig(2) = sig(2) * 180.d0/pi
           endif 
        endif

      return
      end

