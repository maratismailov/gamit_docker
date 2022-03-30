      subroutine geotab(frame,mode,radius,finv,tx,ty,tz)

      implicit real*8(a-h,o-z)

      character*6 frame
      integer mode

c     using WGS84 as the default frame
c      if (frame.eq.'WGS84 '.or.frame.eq.'wgs84 ') then
         radius = 6378137.0d0
         finv = 298.257222101d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0
c      endif

      if (frame.eq.'WGS72 '.or.frame.eq.'wgs72 ') then
         radius = 6378135.0d0
         finv = 298.26d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 4.5d0
      endif
      
      if (frame.eq.'NAD27 '.or.frame.eq.'nad27 ') then
c        these are for NAD27 (Clarke 1866 Ellipsoid)
         radius = 6378206.4d0
         finv = 294.9786982d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0
c      tx = -12.01d0
c      ty = 162.97d0
c      tz = 189.74d0
      endif

      if (frame.eq.'NAD83 '.or.frame.eq.'nad83 ') then
         radius = 6378137.0d0
         finv = 298.257222101d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0
      endif

      if (frame.eq.'AGD84 '.or.frame.eq.'agd84 ') then
         radius = 6378160.0d0
         finv = 298.2500d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0
      endif

      if (frame.eq.'KRASO '.or.frame.eq.'kraso ') then
         radius = 6378245.0d0
         finv = 298.3d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0
      endif
c
c     used by IGN (Clarke 1880)
      if (frame.eq.'CLK80 '.or.frame.eq.'clk80 ') then
         radius = 6378249.2d0
         finv = 293.4660208d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0
      endif
c
c     mak entered INT24 - International 1924 ellipsoid
      if (frame.eq.'INT24'.or.frame.eq.'int24') then
         radius = 6378388.0d0
         finv = 297.0d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0
      endif

      return
      end

