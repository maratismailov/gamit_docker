      real*8 function dchao(elev)
c
c     JLD 870215 Dry Chao mapping function
c
c     ELEV: Elevation angle, radians
c
      real*8 elev
      real*8 a, b
c
      data a / 0.143D-02 /
      data b / 0.455D-01 /
c
      dchao = 1.0D+00 / ( sin(elev) + a / ( tan(elev) + b ) )
c
      end
