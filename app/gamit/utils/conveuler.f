c     Convert Euler vectors
                   
      implicit none

      integer i,j              

      real*8 in(3),in_stats(6),wx,wy,wz,lat,lon,rate,pi
      real*8 dum3(3),dum33(3,3),mat(3,3),sig(3),sigma(3)
      real*8 correl(3,3),corr(3,3),vcv(3,3),cov(3,3),vec(3)
      real*8 r2d,d2r,a,b,azim,rate_sig,c2g,s2g

      character*1 in_type,units,units_stats,type_stats,frame_stats
      character*1 stats,answer
                
      pi = 4.*atan(1.d0)
      r2d = 180.d0/pi
      d2r = pi/180.d0
      
c     initialize correlation and vcv matricies
      do i=1,3
        do j=1,3
	   vcv(i,j) = 0.d0
	   cov(i,j) = 0.d0
	   if ( i .eq. j ) then
             correl(i,j) = 1.d0
	   else
	     correl(i,j) = 0.d0
	   endif
        enddo
      enddo
      
c     By default don't calculate statistics.
      stats = 'n'
      
c     Get pole location information from the user
      write(*,'(a,a)') 
     . 'Enter lat long (deg) w (rad or deg /My) + r or d (rad or deg)'
     . ,' or wx wy wz  + r or d (rad or deg)'
      read(*,*) in,units
   
c     Which pole type do we have?  assume that either lat or long > 1
      in_type = 'w'
      if( abs(in(1)).gt.1. .or. abs(in(2)).gt.1. ) in_type = 's'
                          
      if( in_type.eq.'w' ) then
         if( units.eq.'d' ) then
             do i=1,3
               in(i) = in(i)*d2r
             enddo  
         endif
         wx = in(1)
         wy = in(2)
         wz = in(3)
         call xyzpole( 1,in,dum3,dum33 )
         lat = in(1)
         lon = in(2)
         rate = in(3)

      else 
        if( units.eq.'r' ) then
             do i=1,2
               in(i) = in(i)*r2d
             enddo  
        endif        
        lat = in(1)
        lon = in(2) 
        if ( units.eq.'d' ) in(3) = in(3)*d2r   
        rate = in(3)
        call xyzpole( -1,in,dum3,dum33 )
        wx = in(1)
        wy = in(2)
        wz = in(3)
      endif
      
c     Find out if the user wants to enter pole statistics information
      write(*,'(a,a)')
     .     'Do you want to enter pole vcv/correl matrix (m) or '
     .     ,'error ellipse (e) (m/e/n)? [no statistics (n)]?'
      read(*,*) answer
      if (answer .eq. 'm' )then
        stats = 'y'
        write(*,'(a,a,a/,a,a)') 
     .  'Enter vcv/correl vector (upper triangle) + r or d (rad or deg)'
     .  ,' + v or p (vcv or correl) + g or c (geodetic or cartesian)'
     .  ,' diag values of correl matrix = sigma'
     .  ,'Eg: 3344e-10 2978e-10 1781e-10 2889e-10 1795e-10 1140e-10'
     .  ,' r v c'
        read(*,*) in_stats,units_stats,type_stats,frame_stats
	
        mat(1,1) = in_stats(1)
        mat(2,1) = in_stats(2)
        mat(3,1) = in_stats(3)
        mat(1,2) = in_stats(2)
        mat(2,2) = in_stats(4)
        mat(3,2) = in_stats(5)
        mat(1,3) = in_stats(3)
        mat(2,3) = in_stats(5)
        mat(3,3) = in_stats(6)
      
c     Compute the sigma and correlations from the covariance                                   
        if ( type_stats .eq. "v" ) then
          do i=1,3
	     sigma(i) = sqrt(mat(i,i)) 
	     if ( units_stats.eq.'d' ) sigma(i) = sigma(i)*d2r  
          enddo
          do i=1,3
            do j=1,3
	      vcv(i,j) = mat(i,j)
	      if ( units_stats.eq.'d' ) vcv(i,j) = vcv(i,j)*d2r**2  
              correl(i,j) = vcv(i,j)/sigma(i)/sigma(j)
            enddo
          enddo
        endif
      
c     Compute the sigma and covariances from the correlations       
        if ( type_stats .eq. "p" ) then
          do i=1,3
            sigma(i) = mat(i,i)
	    if ( units_stats.eq.'d' ) sigma(i) = sigma(i)*d2r   
          enddo
          do i=1,3
            do j=1,3
	      if ( i . ne. j ) then
                correl(i,j) = mat(i,j)
	        vcv(i,j) = correl(i,j)*sigma(i)*sigma(j)
	      else
	        correl(i,j) = 1.0d0
	        vcv(i,j) = mat(i,j)**2
		if ( units_stats.eq.'d' ) vcv(i,j) = vcv(i,j)*d2r**2  
	      endif
            enddo
          enddo
        endif

c     Statistics were given in error ellipse format.	
      else if ( answer .eq. "e" ) then
        stats = 'y'
	frame_stats = 'g'
	write(*,'(a,a)') 
     .  'Enter error ellipse vector a, b, azim, rotation rate sigma '
     .  ,'+ r or d (rad or deg) for azim and rotation rate sigma'
        read(*,*) a,b,azim,rate_sig,units_stats
        if ( units_stats.eq.'d' ) then
	  a = a*d2r
	  b = b*d2r
	  azim = azim*d2r
	  rate_sig = rate_sig*d2r
	endif	
        c2g = dcos(2.d0*azim)
        s2g = dsin(2.d0*azim)
        vcv(1,1) = (a**2+b**2+(a**2-b**2)*c2g)/2.d0
        vcv(2,2) = (a**2+b**2+(b**2-a**2)*c2g)/2.d0
        vcv(1,2) = (a**2-b**2)*s2g/2.d0
	vcv(2,1) = vcv(1,2)
	vcv(3,3) = rate_sig**2	
c       compute sigmas (radians)
        do i=1,3
           sigma(i) = sqrt(vcv(i,i))
        enddo
c       compute correlation matrix 
        do i=1,3
           do j=1,3
              correl(i,j) = vcv(i,j)/sigma(i)/sigma(j)
           enddo
        enddo
	
      endif
	
      if ( stats .eq. "y" .and. frame_stats .eq. "c" ) then
          vec(1) = wx
	  vec(2) = wy
          vec(3) = wz
	  sig(1) = sigma(1)
	  sig(2) = sigma(2)
	  sig(3) = sigma(3)
	  do i=1,3
            do j=1,3
              corr(i,j) = correl(i,j)
	    enddo
          enddo
	
          call xyzpole( 2,vec,sig,corr )
	
	  write(6,100)wx*r2d,wy*r2d,wz*r2d,rate*r2d
100       format(/,'Cartesian Pole',/
     .    ,'    Wx       Wy       Wz     rate     (deg/My) ',/,4f9.5)

          write(6,200)(sigma(i)*r2d,i=1,3)
     .    ,((correl(i,j),j=1,3),i=1,3),((vcv(i,j)*r2d**2,j=1,3),i=1,3)
200       format('Cartesian Statistics',/
     .    , 'Sigma    Wx      Wy       Wz    (deg/My)',/,'     ',3f9.5,/
     .    , 'Correlation matrix ',/,3(3f12.5/)
     .    , 'Covariance matrix',/,3(3f15.9/))

          write(6,250) wx*r2d,sigma(1)*r2d,wy*r2d,sigma(2)*r2d,wz*r2d,
     .                sigma(3)*r2d,correl(1,2),correl(1,3),correl(2,3)
250       format('Globk PLATE Cartesian format (XYZ)',/
     .    ,'  Wx (deg/my)   +-    Wy (deg/my)   +-    Wz (deg/My)'
     .    ,'   +-      RhoXY  RhoXZ  RhoYZ   ',/,6f10.6,3f10.6,' XYZ',/)	
  
          write(6,300) vec(1), vec(2), vec(3)*r2d
300       format('Geodetic Pole',/
     .    ,'  Lat       Long      rate(deg/My)',/,3f10.5 )

c       convert lat and long sigma's to radians
          do i=1,2
	    sig(i) = sig(i)*d2r
	  enddo
c       convert the sigmas and correlations into covariances (radians^2)
          do i=1,3
            do j=1,3
              cov(i,j) = corr(i,j)*sig(i)*sig(j)
            enddo
          enddo
     
          write(6,400)(sig(i)*r2d,i=1,3)
     .    ,((corr(i,j),j=1,3),i=1,3),((cov(i,j)*r2d**2,j=1,3),i=1,3)
400       format('Geodetic Statistics',/
     .    , 'Sigma   lat      long     rate(deg/My)',/,'     ',3f9.5,/
     .    , 'Correlation matrix ',/,3(3f12.5/)
     .    , 'Covariance matrix',/,3(3f15.9/))

          write(6,450)vec(1),sig(1)*r2d,vec(2),sig(2)*r2d,vec(3)*r2d,
     .                sig(3)*r2d,corr(1,2),corr(1,3),corr(2,3) 
450       format('Globk PLATE Geodetic format (LLM)',/
     .    ,'  Lat. (deg)    +-    Long (deg)    +-       Mag (deg/My)'
     .    ,'   +-   RhoLtLg RhoLtMg RhoLgMg',/,6f10.4,3f8.4,' LLM',/)	
     
c     Compute error ellipse components
           a = dsqrt(0.5d0 * (cov(2,2)+cov(1,1)
     .      + dsqrt((cov(2,2)-cov(1,1))**2 + 4*cov(1,2)**2)))*r2d 
           b = dsqrt(0.5d0*(cov(2,2)+cov(1,1)
     .      - dsqrt((cov(2,2)-cov(1,1))**2 + 4*cov(1,2)**2)))*r2d
           azim = 90 + (atan(2*cov(1,2)/(cov(1,1)-cov(2,2)))/2.d0)*r2d
  
      else if ( stats .eq. "y" .and. frame_stats .eq. "g" ) then
      
          vec(1) = lat
	  vec(2) = lon
          vec(3) = rate
	  sig(1) = sigma(1)*r2d
	  sig(2) = sigma(2)*r2d
	  sig(3) = sigma(3)
	  do i=1,3
            do j=1,3
              corr(i,j) = correl(i,j)
	    enddo
          enddo
	
          call xyzpole( -2,vec,sig,corr )     
      
 	  write(6,100)wx*r2d,wy*r2d,wz*r2d,rate*r2d
	 
c       convert the sigmas and correlations into covariances
          do i=1,3
            do j=1,3
              cov(i,j) = corr(i,j)*sig(i)*sig(j)*r2d**2
            enddo
          enddo

          write(6,200)(sig(i)*r2d,i=1,3)
     .    ,((corr(i,j),j=1,3),i=1,3),((cov(i,j),j=1,3),i=1,3)
	  
          write(6,250) wx*r2d,sig(1)*r2d,wy*r2d,sig(2)*r2d,wz*r2d,
     .                sig(3)*r2d,corr(1,2),corr(1,3),corr(2,3)

          write(6,300) lat,lon,rate*r2d
     
          write(6,400)(sigma(i)*r2d,i=1,3)
     .    ,((correl(i,j),j=1,3),i=1,3),((vcv(i,j)*r2d**2,j=1,3),i=1,3)

          write(6,450)lat,sigma(1)*r2d,lon,sigma(2)*r2d,rate*r2d,
     .                sigma(3)*r2d,correl(1,2),correl(1,3),correl(2,3) 
     
c     Compute error ellipse components
          a = dsqrt(0.5d0 * (vcv(2,2)+vcv(1,1)
     .      + dsqrt((vcv(2,2)-vcv(1,1))**2 + 4*vcv(1,2)**2)))*r2d 
          b = dsqrt(0.5d0*(vcv(2,2)+vcv(1,1)
     .      - dsqrt((vcv(2,2)-vcv(1,1))**2 + 4*vcv(1,2)**2)))*r2d
          azim = 90 + (atan(2*vcv(1,2)/(vcv(1,1)-vcv(2,2)))/2.d0)*r2d
     
      endif
      
      write(*,'(a)') 'Pole Position Summary'
      write(*,'(a,a)') '     Wx      Wy      Wz     Rate(rad)   Lat'    
     .  ,'     Lon     Wx     Wy      Wz     rate(deg)'
        write(*,'(4f9.5,2f8.2,4f8.4)') wx,wy,wz,rate,lat,lon
     .   ,wx*r2d,wy*r2d,wz*r2d,rate*r2d 
     
c     Print error ellipse 
      if ( stats .eq. 'y' ) then
        write(6,500)a,b,azim
500     format(/,'Error ellipse',/
     .      ,'     a        b     azim    (deg)',/,3f9.4)
      endif	
               
      stop
      end
