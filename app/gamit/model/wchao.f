      real*8 function wchao(elev)
c
c     JLD 870215 Wet Chao mapping function
c
c     ELEV: Elevation angle, radians
c
      real*8 elev
      real*8 c, d
c
      data c / 0.35D-03 /
      data d / 0.17D-01 /
c
      wchao = 1.0D+00 / ( sin(elev) + c / ( tan(elev) + d ) )
c
      end
