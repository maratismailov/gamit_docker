      real*8 function multip (el,r0,hi,wavel)

c     calculate path difference due to multipath

c     input
c        el       elevation angle in radians
c        r0       modeled range from antenna to sat in meters
c        hi       antenna height above ground in meters
c        wavel    wavelength in meters
c     output
c        multip   observed phase minus predicted multipath

c     notation
c        i is for the small triangle including the Instrument
c        r is for Received
c        0 is for the direct path
c        s is for the large triange including the Satelxite

      real*8 hi,xi,ri,x0,r0,hs,xs,rs,el,wavel,pi,phis,phi0
      real*8 phii
      complex yi,y0,yr,zi

c     function to return fractional part
      real*8 frac

c     function to return argument of complex number
      real*8 phase

      zi = csqrt((-1,0))
      pi = 4.0d0 * datan(1.0d0)

      if (hi .gt. 0.0d0) then

c        vertical distance of sat above reflector
         hs = hi + r0 * dsin(el)

c        horizontal distances
         x0 = r0 * dcos(el)
         xs = x0 * hs/(hi+hs)
         xi = x0 * hi/hs

c        hypoteneuses
         ri = dsqrt(hi*hi + xi*xi)
         rs = dsqrt(hs*hs + xs*xs)

c        waves
c        direct path
         phi0 = frac(r0/wavel)
         y0 = 1.0 * cexp(zi*cmplx(phi0))

c        sat to reflector
         phis = frac(rs/wavel)

c        reflector to antenna (include pi phase shift)
c        reflection coefficient is 0.2
         phii = frac(ri/wavel)
         yi = 0.2 * cexp(zi*cmplx(phis + pi + phii))

c        sum them
         yr = y0 + yi

c        find the phase shift due to multipath
         multip = phase(yr) - phi0

      else
         multip = 0.0d0
      endif


      return
      end

      real*8 function phase(z)
c     return argument of a complex number
      complex z
      real*8 realp,znorm

      phase = dacos(realp(z)/znorm(z))

      if (aimag(z) .lt. 0.) phase = -phase

      return

      end

      real*8 function znorm(z)
      complex z
      real*8 realp
c     return magnitude of complex number

      znorm = dsqrt(aimag(z)*aimag(z) + realp(z)*realp(z))

      return

      end
      
      real*8 function realp(z)
c     return real part of z
      complex z
      complex conjg

      realp = dble(z + conjg(z))/2.

      return

      end





