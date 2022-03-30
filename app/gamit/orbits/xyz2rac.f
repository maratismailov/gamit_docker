      Subroutine xyz2rac( posvel,dxyz,drac )

c     convert delta x,y,z to radial, along-track, cross-track

c       posvel(6)    position,velocity of SV
c       dxyz(3)      cartesian difference
c       drac(3)      rad,along,cross difference

      implicit none
                       
      real*8 posvel(6),dxyz(3),drac(3)
      real*8 radlen,sprodr,vellen,c(3),sproda
      real*8 dot
c      real*8 c1,c2,c3,crosslen,sprodc


c     length of position vector
      radlen = dsqrt( dot(posvel,posvel) )   
                                 
c     projection of dxyz onto position vector 
      sprodr = dot(dxyz,posvel)

c     radial difference
      drac(1) = sprodr/radlen

c     length of velocity vector
      vellen = dsqrt( dot(posvel(4),posvel(4)) ) 

c     projection of dxyz onto velocity vector
      sproda = dot(dxyz,posvel(4))

c     along-track difference
      drac(2) = sproda/vellen

c     cross-product of position and velocity
c      c1 = posvel(2)*posvel(6) -
c     .     posvel(3)*posvel(5)
c      c2 = posvel(3)*posvel(4) -
c     .     posvel(1)*posvel(6)
c      c3 = posvel(1)*posvel(5) -
c     .     posvel(2)*posvel(4)
                               
c      print *,'posvel dxyz ',posvel,dxyz 
c      print *,'c1 c2 c3 ',c1,c2,c3
c     length of cross-track vector
c      crosslen=dsqrt(c1**2+c2**2+c3**2)

c     projection of dxyz onto cross-track vector
c      sprodc = dxyz(1)*c1 + dxyz(2)*c2 + dxyz(3)*c3
c
c cross track differences given by the above scalar product divided by the
c length of the the pol/vel cross product vector.
c
c      drac(3) = sprodc/crosslen
c
c      print *,'c1 c2 c3 sprodc drac3 ',c1,c2,c3,sprodc,drac(3)

      call cross_unit(posvel(1),posvel(4),c)
      drac(3) = dot(dxyz,c)
c      print *,'c drac3 ',c,drac(3)
           
      return
      end
