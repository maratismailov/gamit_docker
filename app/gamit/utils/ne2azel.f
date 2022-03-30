      Program ne2azel

c     Utility to convert a north-east velocity and uncertainty to az-el
c     R. King 5 Feb 97

      implicit none
                         
      integer i

      real*8 n,nsig,e,esig,necor,ne(2),nesig(2),r,rsig,az,azsig,razcor
     .     , pi,rad2deg
       
      pi = 4.d0*atan(1.d0)
      rad2deg = 180.d0/pi

      write(*,'(/,a,/)') 
     . 'Enter n-vel nv-sig  e-vel ev-sig correlation  (mm/yr)'
      read(*,*) n,nsig,e,esig,necor

      ne(1) = n
      ne(2) = e
      nesig(1) = nsig
      nesig(2) = esig
      
      call magdir( ne, nesig, necor, r, rsig, az, azsig, razcor )
          
      write(*,'(/,a,/,24x,9(f8.2))')
     .'Local velocities (mm/yr)  North     +/-     East    +/-     Rate
     .    +/-    Az       +/-    Correl.'
     .  ,(ne(i),nesig(i),i=1,2)
     . , r,rsig,az*rad2deg,azsig*rad2deg,razcor  

      stop
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
     

