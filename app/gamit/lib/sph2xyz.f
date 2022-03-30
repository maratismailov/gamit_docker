      Subroutine sph2xyz( latr,lonr,rad,pos )

c     Subroutine to convert spherical coordinates in radians to Cartesian

      implicit none
    
      real*8 pos(3),latr,lonr,rad
                          
      pos(1) = rad*dcos(latr)*dcos(lonr)
      pos(2) = rad*dcos(latr)*dsin(lonr)
      pos(3) = rad*dsin(latr)
 
      return
      end
