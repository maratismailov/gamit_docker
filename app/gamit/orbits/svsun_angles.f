      Subroutine svsun_angles(xsv,xsun,betar,sunlon,svbcos)
  
      implicit none

c  Input 
      real*8 xsv(6),xsun(6)
c       x     E-centered inertial coordinates of the satellite
c       xsun  E-centered inertial coordinates of the sun

c  Output
      real*8 betar,sunlon,svbcos
c       betar    Sun angle wrt to the orbital plane (radians)
c       long     Sun angle with respect to the orbital node (radians)
c       svbcos   Cosine of the angle between the SV position and the Sun
                                                       
        
c  Internal                                    
      integer*4 i 
      real*8 vnu(3),zv(3),sunsv(3),sundist,svdist,svrate
     .     , xsvhat(3),sunhat(3),vsvchat(3),orbnorm(3),h(3)
     .     , R1(3,3),R3(3,3),ROT(3,3),sunhatorb(3)
     .     , ombig,finc,twopi

c  Function
      real*8 dot

      data zv/0.d0,0.d0,1.d0/, twopi/6.283185307179586d0/

c  Sun wrt SV 
      do i=1,3
         sunsv(i) = xsun(i) - xsv(i)
      enddo  
  
c  Unit vectors of Sun, SV position, and SV velocity
      sundist = dsqrt(sunsv(1)**2+sunsv(2)**2+sunsv(3)**2)   
      svdist =  dsqrt(xsv(1)**2+xsv(2)**2+xsv(3)**2)
      svrate =  dsqrt(xsv(4)**2+xsv(5)**2+xsv(6)**2)
      do i=1,3
        sunhat(i) = sunsv(i)/sundist
        xsvhat(i) = xsv(i)/svdist 
        vsvchat(i) = xsv(i+3)/svrate       
      enddo
                             
c Orbit normal (Kouba convention) and angular momentum for ascending node (keplr convention)
      call cross( vsvchat,xsvhat,orbnorm )  
      do i=1,3
        h(i)= -orbnorm(i)
      enddo
  
c Angle between the Sun and the orbital plane
      betar = dacos( dot(orbnorm,sunhat) )                 

c Cosine of angle between the SV and the Sun normal   
      svbcos = dot(xsvhat,sunhat) 
              
c Ascending node and inclination of the orbit in the E-equatorial frame (from routine keplr.f)
      call cross(zv,h,vnu)      
      ombig = datan2(vnu(2),vnu(1))
      finc = datan2(dsqrt(h(1)**2+h(2)**2+h(3)**2),h(3))
      if(finc.lt.0.d0) finc = finc+twopi

c Transform Sun coordinates into the orbital frame  Xorb = R1(finc))*R3(ombig) * X 
      R3(1,1) = dcos(ombig)
      R3(1,2) = dsin(ombig)
      R3(1,3) = 0.d0
      R3(2,1) = -dsin(ombig)
      R3(2,2) = dcos(ombig)
      R3(2,3) = 0.d0
      R3(3,1) = 0.d0 
      R3(3,2) = 1.d0
      R3(3,3) = 0.d0
      R1(1,1) = 1.d0
      R1(1,2) = 0.d0
      R1(1,3) = 0.d0
      R1(2,1) = 0.d0
      R1(2,2) = dcos(finc)
      R1(2,3) = dsin(finc)
      R1(3,1) = 0.d0
      R1(3,2) = -dsin(finc)
      R1(3,3) = dcos(finc)
      call matmpy(R3,R1,ROT,3,3,3)
      call matmpy(ROT,sunhat,sunhatorb,3,3,1) 
      sunlon = datan2(sunhatorb(2),sunhatorb(1))
 
      return
      end

