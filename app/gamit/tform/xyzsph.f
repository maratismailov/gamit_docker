      subroutine xyzsph(x,alat,along,rad)
c
c  converts cartesian x,y,z to spherical coordinates
c
      implicit real*8(a-h,o-z)
      dimension x(3)
c
      twopi= 8.0*atan(1.0)
c
      rad= sqrt(x(1)**2+x(2)**2+x(3)**2)
      along= atan2(x(2),x(1))
c      if( along.lt.0.d0 ) along=along+twopi
      alat= asin(x(3)/rad)
      along= along*360.d0/twopi
      alat=  alat*360.d0/twopi
c
      return
      end
