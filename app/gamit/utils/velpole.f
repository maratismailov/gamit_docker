
      Subroutine velpole( pos, wvec, wsig, wcorrel, svec, ssig, scorrel)
                     
c     Compute a station's velocity and uncertainty from plate rotation
c     angular velocity vector and covariances

c     R. King 29 Jan 97                   

c     Input
c     
c       pos(3)       :  Cartesian geocentric position vector of station (meters)
c                
c       wvec(3)      :  Cartesian angular velocity vector for plate motion (rad/My)
c                
c       wsig(3)      :  Uncertainties of angular velocity components (rad/My)
c
c       wcorrel(3,3) :  Correlation matrix for angular velocity components    
c
c     Output
c
c       svec(3)      :  Cartesian geocentric velocity vector of station  (m/yr)   
c  
c       ssig(3)      :  Uncertainties of station velocities (m/yr)
c
c       scorrel(3,3) :  Correlation matrix for station velocities 
                  
      implicit none
  
      integer i,j
                   
      real*8 wvec(3),wsig(3),wcorrel(3,3)
     .     , pos(3),svec(3),ssig(3),scorrel(3,3)
     .     , jac(3,3),jact(3,3),cov(3,3),cov_out(3,3)
     .     , work33(3,3)
               
                  
c   compute the station velocity   r x w

      svec(1) = pos(3)*wvec(2) -  pos(2)*wvec(3)
      svec(2) = pos(1)*wvec(3) -  pos(3)*wvec(1)
      svec(3) = pos(2)*wvec(1) -  pos(1)*wvec(2)

c   compute the Jacobian [ d svec / d wvec ]  meters / (rad/My)

c     d v / d wx
      jac(1,1) =  0.d0
      jac(2,1) = -pos(3)
      jac(3,1) =  pos(2) 
c     d v / d wy
      jac(1,2) =  pos(3)
      jac(2,2) =  0.d0
      jac(3,2) = -pos(1) 
c     d v /d wz
      jac(1,3) = -pos(2)
      jac(2,3) =  pos(1)
      jac(3,3) =  0.d0            
c       add a small effect for the radial (diagonal) elements to keep 
c       the propagation to local "up" finite.   Use 1.d-6 of the off-  
c       -diagonal partials (6.d6), which makes the radial partial
c       6 m/ (rad/My) 
      jac(1,1) = 6.d0
      jac(2,2) = 6.d0
      jac(3,3) = 6.d0
   
    
c    convert the plate rotation sigmas and correlations into covariances 

        do i=1,3
          do j=1,3
            cov(i,j) = wcorrel(i,j)*wsig(i)*wsig(j)
          enddo
        enddo

c     pre- and post-multiply the covariance by the Jacobian 
c       C(out) = J * C(in) * J(transpose)
                    
c          print *,'cov ',cov
        call matmpy(jac,cov,work33,3,3,3)
c          print *,'work33 ',work33     
        call transp(jac,jact,3,3)
        call matmpy(work33,jact,cov_out,3,3,3)
c          print *,'cov_out ',cov_out

c     compute the sigmas and correlations from the covariances

        do i=1,3
          ssig(i) = sqrt(cov_out(i,i)) 
        enddo
        do i=1,3
          do j=1,3
            scorrel(i,j) = cov_out(i,j)/ssig(i)/ssig(j)
          enddo
        enddo

c      change units from m/Myr to m/yr

        do i=1,3
          svec(i) = svec(i)/1.d6
          ssig(i) = ssig(i)/1.d6
        enddo

      return
      end
     

      subroutine gdetic( fin,out,a,finv )
 
      implicit none
c
c     given geodetic coordinates in fin return geocentric
c     coordinates in out 


c     formulae are closed form:
c
c       tan(lat')=((1-f)**2 + 2*f*h/a) * tan(lat)
c
c      S.A. Gourevitch 6/25/81
c      Shortened for velpole by R. King 1/29/97

      real*8 fin(3),out(3)
     .     , zero,one,two,piotwo,onefsq,h,lat,hilat
     .     , fden,finv,fnum,a,f,fac,del

      data zero,one,two/0.d0,1.d0,2.d0/
c
c high latitude cutoff
      data hilat/5.d-4/
c
c pi/2
      piotwo = two*datan(one)
c
      lat=fin(1)
      h=fin(3)
      f= one/finv
      onefsq = (one-f)**2
c
c latitude : take into account lat=90
      del=dabs(lat)-piotwo
      if(dabs(del).gt.hilat) then
         out(1)=datan((onefsq+two*f*h/a) * dtan(lat))
      else
         out(1)=dsign(piotwo + del/(onefsq+two*f*h/a),lat)
      endif
c
c longitude
      out(2)=fin(2)
c
c radius
      fnum=dcos(lat)**2 + onefsq**2 * dsin(lat)**2
      fden=dcos(lat)**2 + onefsq    * dsin(lat)**2
      fac=dsqrt(fnum/fden)
      out(3)= a*fac + h*(one-f*f*dsin(two*lat)/two)

      return
      end


      Subroutine xyzneu ( lat, lon, xyz, xyzsig, xyzcor
     .                  , neu, neusig, neucor )
                     
c     Map coordinates or velocities and uncertainties from geocentric to local 

c     R. King 29 Jan 97                   

c     Input             
c
c       lat          :  latitude of the station (deg)
c
c       lon          :  longitude of the station (deg)
c     
c       xyz(3)       :  geocentric baseline or velocity (right-handed xyz)
c                
c       xyzsig(3)    :  geocentric baseline or velocity sigma
c                
c       xyzcor(3,3)  :  geocentric correlation matrix
c                
c      Output
c
c       neu(3)       :  local baseline or velocity (left-handed neu)
c                
c       neusig(3)    :  local basline or velocity sigma
c                
c       neucor(3,3)  :  local correlation matrix

c
      implicit none
  
      integer i,j
                   
      real*8 lat,lon,xyz(3),xyzsig(3),xyzcor(3,3)
     .     , neu(3),neusig(3),neucor(3,3)
     .     , sinlat,coslat,sinlon,coslon
     .     , jac(3,3),jact(3,3),cov(3,3),cov_out(3,3)
     .     , work33(3,3),pi,rad2deg
               
      
      pi = 4.d0*atan(1.d0)
      rad2deg = 180.d0/pi


c   compute the coordinates 

      sinlat = sin(lat/rad2deg)
      coslat = cos(lat/rad2deg)
      sinlon = sin(lon/rad2deg)
      coslon = cos(lon/rad2deg)

c     north       
      neu(1) = -sinlat*coslon*xyz(1) - sinlat*sinlon*xyz(2) 
     .         + coslat*xyz(3) 
c     east                          
      neu(2) = -sinlon*xyz(1) + coslon*xyz(2)
c     up
      neu(3) = coslat*coslon*xyz(1) + coslat*sinlon*xyz(2) 
     .         + sinlat*xyz(3)


c   compute the Jacobian [ dxyz/dneu ]
          
c     d(neu) / dx
      jac(1,1) = -sinlat*coslon
      jac(2,1) = -sinlon
      jac(3,1) =  coslat*coslon  
c     d(neu) / dy
      jac(1,2) = -sinlat*sinlon
      jac(2,2) =  coslon
      jac(3,2) =  coslat*sinlon
c     d(neu) / dz
      jac(1,3) =  coslat
      jac(2,3) =  0.
      jac(3,3) =  sinlat

    
c    convert the xyz sigmas and correlations into covariances 

        do i=1,3
          do j=1,3
            cov(i,j) = xyzcor(i,j)*xyzsig(i)*xyzsig(j)
          enddo
        enddo

c     pre- and post-multiply the covariance by the Jacobian 
c       C(out) = J * C(in) * J(transpose)
c          print *,'jac ',jac          
c          print *,'cov ',cov
        call matmpy(jac,cov,work33,3,3,3)
c          print *,'work33 ',work33     
        call transp(jac,jact,3,3)
        call matmpy(work33,jact,cov_out,3,3,3)
c          print *,'cov_out ',cov_out

c     compute the local sigmas and correlations from the covariances

        do i=1,3
          neusig(i) = sqrt(cov_out(i,i)) 
        enddo
        do i=1,3
          do j=1,3
            neucor(i,j) = cov_out(i,j)/neusig(i)/neusig(j)
          enddo
        enddo

      return
      end

          
                     
      Subroutine getmag( vec, vsig, vcor, mag, msig )
                     
c     Compute the magnitude of a vector and its uncertainty           

c     R. King 31 Jan 97                   

c     Input
c     
c       vec(3)       :  Cartesian vector                    
c                
c       vsig( 3)     :  Cartesian uncertainties
c                
c       vcor(3,3)     :  Correlation matrix 
c
c     Output
c
c       mag         :  Magnitude of the vector   
c  
c       msig        :  Uncertainty of the magnitude
                  
      implicit none
  
      integer i,j
                
      real*8 vec(3),vsig(3),vcor(3,3),mag,msig
     .     , jac(3),cov(3,3),var_out
     .     , work3(3)
                                
c   compute the magnitude

      mag = sqrt( vec(1)**2 + vec(2)**2  + vec(3)**2 )

c   compute the Jacobian [ d mag / d vec ]  

c     d mag / dx
      jac(1) = vec(1)/mag
      jac(2) = vec(2)/mag
      jac(3) = vec(3)/mag

c    convert the vector sigmas and correlations into covariances 

        do i=1,3
          do j=1,3
            cov(i,j) = vcor(i,j)*vsig(i)*vsig(j)
          enddo
        enddo

c     pre- and post-multiply the covariance by the Jacobian 
c       C(out) = J * C(in) * J(transpose)
                    
c          print *,'cov ',cov
        call matmpy(cov,jac,work3,3,3,1)
c          print *,'work3 ',work3     
        call matmpy(work3,jac,var_out,1,3,1)
c          print *,'var_out ',var_out

c     compute the sigma of the magnitude

        msig = sqrt(var_out)

      return
      end
     
      Subroutine magdir( ne, nesig, necor, r, rsig, az, azsig, razcor )
                     
c     Compute the magnitude of a vector and its uncertainty           

c     R. King 31 Jan 97                   

c     Input
c       
c       ne(2)       :  North, east horizontal vector                    
c                
c       nesig(2)    :  Carteisian uncertainties
c                
c       necor       :  Correlation between north and east 
c
c     Output
c
c       r           :  Magnitude of the vector   
c        
c       rsig        :  Uncertainty of the magnitude
c
c       az          :  Azimuth of vector (degrees)
c
c       azsig       :  Uncertainty of azimuth (degrees)   
c
c       razcor      :  Correlation between magnitude and azimuth
                  
      implicit none
  
      real*8 ne(2),nesig(2),necor,r,rsig,az,azsig,razcor
     .     , n,e,jac(2,2),jact(2,2),cov(2,2),cov_out(2,2)
     .     , work22(2,2),pi,rad2deg
                          

      pi = 4.d0*atan(1.d0)
      rad2deg = 180.d0/pi

       
c   humanize the notation

      n = ne(1)
      e = ne(2)
                       
c   compute the magnitude and direction

      r = sqrt( n**2 + e**2 )
      az  = atan2( e,n )

c   compute the Jacobian [ d r,az / d n,e ]  

c     use dir = acos( n/r)  and
c         dir = asin( e/r ) to avoid indeterminancy of atan function and its derivatives  
              
c     d r / d n
      jac(1,1) = n/r
c     d theta / d n
      jac(2,1) = -sin(az)/r
c     d r / d e
      jac(1,2) = e/r
c     d theta / d e
      jac(2,2) = cos(az)/r

c    convert the vector sigmas and correlations into covariances 

      cov(1,1) = nesig(1)**2
      cov(2,2) = nesig(2)**2
      cov(1,2) = necor*nesig(1)*nesig(2)
      cov(2,1) = cov(1,2)
            
c     pre- and post-multiply the covariance by the Jacobian 
c       C(out) = J * C(in) * J(transpose)
                       
c          print *,'nesig necor ',nesig,necor
c          print *,'cov ',cov
      call matmpy(jac,cov,work22,2,2,2)
c          print *,'work22 ',work22    
      call transp(jac,jact,2,2)     
c          print *,'jac jact ',jac,jact
      call matmpy(work22,jact,cov_out,2,2,2)
c          print *,'cov_out ',cov_out
                             

c     compute the sigmas and correlations of magnitude and direction
          
      rsig = sqrt(cov_out(1,1))
      azsig = sqrt(cov_out(2,2))
      razcor = cov_out(1,2)/rsig/azsig

      return
      end
     


