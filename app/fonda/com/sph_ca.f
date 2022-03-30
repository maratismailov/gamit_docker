      subroutine sph_ca(alat,alon,xin,yin,zin,
     .  xout,yout,zout,job)
c	
c     transform a vector from spherical coordinate to 
c     Cartesian coordinate or vice versa.
c
c     xin,yin,zin   : components of input vector
c     xout,yout,zout: components of output vector
c     job  :  1 = spherical to Cartesian
c             2 = Cartesian to spherical
c     alat,alon     : unit radian
c
cmk   This is a topocentric <-> Cartesian transformation
cmk   Topocentric is in the form ENU, not NEU!
cmk   Cartesian is in the form XYZ

      implicit real*8 (a-h,o-z)
      integer job

      s1 = dsin(alat)
      c1 = dcos(alat)
      s2 = dsin(alon)
      c2 = dcos(alon)

      if (job.eq.1) then
      xout = -xin*s2+(zin*c1-yin*s1)*c2
      yout =  xin*c2+(zin*c1-yin*s1)*s2
      zout =  yin*c1+zin*s1
      endif

      if (job.eq.2) then
      xout = -xin*s2+yin*c2
      yout = -(xin*c2+yin*s2)*s1+zin*c1
      zout = (xin*c2+yin*s2)*c1+zin*s1
      endif

      return
      end
