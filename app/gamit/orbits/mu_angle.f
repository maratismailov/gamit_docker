      Subroutine mu_angle( xsv,vsvc,xsun,mu,debug,iprndb )

c     rwk 170520: Compute the orbit angle mu using vector products
                                                      
c         xsv   radius vector of the satellite
c         vsvc  velocity vector of the satellite
c         xsun  Sun wrt Earth
c         mu    angle between orbit midnight and the satellite
c       For debug only:
c         debug   T/F
c         iprndb  PRN number for debug printout 
c         

c     The positive orbit normal is the radius vector crossed with the velocity vector
c        ON = R x V
c     Crossing the ON with the Sun vector gives a vector pointing 90-deg from midnight
c        A = ON x S
c     The angle between A and R is mu + 90 degrees.


      real*8 xsv(3),vsvc(3),xsun(3),mu
     .     , svmag,svvmag,sunmag,er(3),ev(3),en(3),es(3),eg(3),ea(3)
     .     , ermag,enmag,egmag,eamag,dot,sinalpha,cosalpha 
     .     , pi,dtr 
                                                   
      logical debug
      integer*4 iprndb 
              
      pi=4.d0*datan(1.d0) 
      dtr=pi/180.D0    
                 
cd      print *,'MU_ANGLE xsv xsun debug iprndb '
cd     .        ,xsv,xsun,debug,iprndb 
c  get the radius, velocity, and sun  unit vectors
      svmag = dsqrt(dot(xsv,xsv))
      svvmag = dsqrt(dot(vsvc,vsvc))
      sunmag = dsqrt(dot(xsun,xsun))
      do i=1,3
        er(i) = xsv(i)/svmag 
        ev(i) = vsvc(i)/svvmag
        es(i) = xsun(i)/sunmag 
      enddo  

c  get the orbit normal unit vector
      call cross(er,ev,en)   
      enmag = dsqrt(dot(en,en))  
      do i=1,3
        en(i) = en(i)/enmag
      enddo

c  get the vector 90 degrees from midnight in the orbital plane 
      call cross(en,es,eg)                                 
      egmag = dsqrt(dot(eg,eg))
      do i=1,3
        eg(i) = eg(i)/egmag
      enddo   

      if( debug.and.iprndb.ne.0 ) then 
        print *,'er ',er
        print *,'ev ',ev 
        print *,'en ',en
        print *,'eg ',eg
      endif
                            
c  get the angle ('alpha') between the 90-degree vector and the SV vector
c  ( ea is also an orbit-normal direction but will have  zero at midnight )
      call cross(eg,er,ea)  
      eamag = dsqrt(dot(ea,ea)) 
      ermag = dsqrt(dot(er,er))
      sinalpha = eamag/(egmag*ermag)
      cosalpha = dot(eg,er)/(egmag*ermag)
      mu = datan2(sinalpha,cosalpha)/dtr - 90.d0  

      if( debug.and.iprndb.ne.0 ) 
     .   print *,'sinalpha cosalpha mu ',sinalpha,cosalpha,mu
                                     
      return
      end
