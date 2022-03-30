      subroutine azimuc(isit1,isit2,dt,coef,mode)
c
c     get coefficients for azimuth abservation
c     all final partial derivatives are related to global 
c     Cartesian coordinates and velocities
c
c     unit:
c        dt  : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer isit1,isit2,mode,i,j
      dimension coef(12),cjaco(9)
      dimension avec(3),cvec(3)
c
c     use geocentric coordinate, then transform to geodetic
      xd = x(isit2)-x(isit1)
      yd = y(isit2)-y(isit1)
      zd = z(isit2)-z(isit1)

c     put deflection correction
      sla = slat(isit1)+defn(isit1)
      slo = slon(isit1)+defe(isit1)/dcos(slat(isit1))

      s1 = dsin(slo)
      c1 = dcos(slo)
      s2 = dsin(sla)
      c2 = dcos(sla)
      dx = -xd*s1    + yd*c1
      dy = -xd*c1*s2 - yd*s1*s2 + zd*c2

c     horizontal baseline length squared
      bl2 = dx**2+dy**2

c     coefficients for local planar coordinate
      avec(1) = -dy/bl2
      avec(2) =  dx/bl2
      avec(3) =  0.0d0

c     transform from local to global Cartesian coordinate
      if (pmode.eq.2) then
         do i=1,3
            cvec(i) = avec(i)
         enddo
      else
         call getjac(slo,sla,srad(isit1),cjaco,5)
         call axb(1,3,3,avec,cjaco,cvec,1,0)
      endif

c     partials for first station
      do 15 j = 1,3
         coef(j)   = cvec(j)
         coef(j+3) = cvec(j)*dt
 15   continue

c     partials for second station
      do 17 i = 1,6
         coef(i+6) = -coef(i)
 17   continue

      if (pmode.eq.2) then
c        local coordinate case
         call getjac(x(isit1),y(isit1),z(isit1),cjaco,4)
         a1 = cvec(1)*cjaco(2)+cvec(2)*cjaco(5)+cvec(3)*cjaco(8)
         a2 = cvec(1)*cjaco(1)+cvec(2)*cjaco(4)+cvec(3)*cjaco(7)
         a3 = cvec(1)*cjaco(3)+cvec(2)*cjaco(6)+cvec(3)*cjaco(9)
         coef(2) = a1/rtod/3.6d3
         coef(1) = a2/rtod/3.6d3
         coef(3) = a3
         call getjac(x(isit2),y(isit2),z(isit2),cjaco,4)
         a1 = cvec(1)*cjaco(2)+cvec(2)*cjaco(5)+cvec(3)*cjaco(8)
         a2 = cvec(1)*cjaco(1)+cvec(2)*cjaco(4)+cvec(3)*cjaco(7)
         a3 = cvec(1)*cjaco(3)+cvec(2)*cjaco(6)+cvec(3)*cjaco(9)
         coef(8) = -a1/rtod/3.6d3
         coef(7) = -a2/rtod/3.6d3
         coef(9) = -a3
         do i = 1,3
            coef(3+i) = avec(i)*dt*1.0d-6
            coef(9+i) = -coef(3+i)
         enddo
      endif

c     print some stuff:
c     do i=1,3
c        print *,'AZIMUC: dt,coef ', dt,coef(i),coef(i+3)
c     enddo

      return
      end
c
