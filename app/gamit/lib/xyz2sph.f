      Subroutine xyz2sph( pos,latr,lonr,rad )

c     Subroutine to convert Cartesian coordinates to spherical
c     Output lat lon are in radians

      implicit none
    
      real*8 pos(3),latr,lonr,rad

      rad = dsqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
      lonr = datan2(pos(2),pos(1)) 
      latr = dasin(pos(3)/rad)  

      return
      end
