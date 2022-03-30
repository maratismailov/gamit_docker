      subroutine baslc(isit1,isit2,bl,dt,coef,mode)
c
c     get coefficients for baseline length abservation
c     all final partial derivatives are related to global
c     Cartesian coordinates and velocities
c
c     unit:
c        dt  : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer isit1,isit2,mode,i
      dimension coef(12),cjaco(12),tmp(6)
      fac = 1.0d0
c
c     coefficients (xyz coordinate)	
      tv = vx(isit2)-vx(isit1)
      coef(1) = -(x(isit2)-x(isit1)+tv*dt)/bl
      tv = vy(isit2)-vy(isit1)
      coef(2) = -(y(isit2)-y(isit1)+tv*dt)/bl
      tv = vz(isit2)-vz(isit1)
      coef(3) = -(z(isit2)-z(isit1)+tv*dt)/bl
      coef(4) = dt*coef(1)
      coef(5) = dt*coef(2)
      coef(6) = dt*coef(3)
c      print*,'bl,omc,wght ',bl,omc,wght
      do 25 i = 1,6
      coef(i+6) = -coef(i)
 25   continue
      if (pmode.le.2) goto 100

      tmp(1) = coef(1)
      tmp(2) = coef(2)
      tmp(3) = coef(3)
      call getjac(x(isit1),y(isit1),z(isit1),cjaco,4)
      a1 = tmp(1)*cjaco(2)+tmp(2)*cjaco(5)+tmp(3)*cjaco(8)
      a2 = tmp(1)*cjaco(1)+tmp(2)*cjaco(4)+tmp(3)*cjaco(7)
      a3 = tmp(1)*cjaco(3)+tmp(2)*cjaco(6)+tmp(3)*cjaco(9)
      coef(2) = a1/rtod/3.6d3
      coef(1) = a2/rtod/3.6d3
      coef(3) = a3
      call getjac(x(isit2),y(isit2),z(isit2),cjaco,4)
      a1 = tmp(1)*cjaco(2)+tmp(2)*cjaco(5)+tmp(3)*cjaco(8)
      a2 = tmp(1)*cjaco(1)+tmp(2)*cjaco(4)+tmp(3)*cjaco(7)
      a3 = tmp(1)*cjaco(3)+tmp(2)*cjaco(6)+tmp(3)*cjaco(9)
      coef(8) = -a1/rtod/3.6d3
      coef(7) = -a2/rtod/3.6d3
      coef(9) = -a3
      tmp(1) = coef(4)
      tmp(2) = coef(5)
      tmp(3) = coef(6)
      call getjac(x(isit1),y(isit1),z(isit1),cjaco,5)
      a1 = tmp(1)*cjaco(2)+tmp(2)*cjaco(5)+tmp(3)*cjaco(8)
      a2 = tmp(1)*cjaco(1)+tmp(2)*cjaco(4)+tmp(3)*cjaco(7)
      a3 = tmp(1)*cjaco(3)+tmp(2)*cjaco(6)+tmp(3)*cjaco(9)
      coef(5) = a1*fac
      coef(4) = a2*fac
      coef(6) = a3*fac
      call getjac(x(isit2),y(isit2),z(isit2),cjaco,5)
      a1 = tmp(1)*cjaco(2)+tmp(2)*cjaco(5)+tmp(3)*cjaco(8)
      a2 = tmp(1)*cjaco(1)+tmp(2)*cjaco(4)+tmp(3)*cjaco(7)
      a3 = tmp(1)*cjaco(3)+tmp(2)*cjaco(6)+tmp(3)*cjaco(9)
      coef(10) = -a1*fac
      coef(11) = -a2*fac
      coef(12) = -a3*fac
c
 100  continue
      return
      end
c
c   --------------------------------------------------------
      subroutine baslcv(isit1,isit2,bl,dt,coef,mode)
c
c     get coefficients for baseline length rate abservation
c     all final partial derivatives are related to global
c     Cartesian coordinates and velocities
c
c     unit:
c        dt  : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer isit1,isit2,mode,i
      dimension coef(12),cjaco(12),tmp(6)
      fac = 1.0d0
c
c     coefficients (xyz coordinate)	
      coef(1) = 0.0d0
      coef(2) = 0.0d0
      coef(3) = 0.0d0
      tv = vx(isit2)-vx(isit1)
      coef(4) = -(x(isit2)-x(isit1)+2.0d0*tv*dt)/bl
      tv = vy(isit2)-vy(isit1)
      coef(5) = -(y(isit2)-y(isit1)+2.0d0*tv*dt)/bl
      tv = vz(isit2)-vz(isit1)
      coef(6) = -(z(isit2)-z(isit1)+2.0d0*tv*dt)/bl
c    
      do 25 i = 1,6
      coef(i+6) = -coef(i)
 25   continue
      if (pmode.le.2) goto 100

      tmp(1) = coef(1)
      tmp(2) = coef(2)
      tmp(3) = coef(3)
      call getjac(x(isit1),y(isit1),z(isit1),cjaco,4)
      a1 = tmp(1)*cjaco(2)+tmp(2)*cjaco(5)+tmp(3)*cjaco(8)
      a2 = tmp(1)*cjaco(1)+tmp(2)*cjaco(4)+tmp(3)*cjaco(7)
      a3 = tmp(1)*cjaco(3)+tmp(2)*cjaco(6)+tmp(3)*cjaco(9)
      coef(2) = a1/rtod/3.6d3
      coef(1) = a2/rtod/3.6d3
      coef(3) = a3
      call getjac(x(isit2),y(isit2),z(isit2),cjaco,4)
      a1 = tmp(1)*cjaco(2)+tmp(2)*cjaco(5)+tmp(3)*cjaco(8)
      a2 = tmp(1)*cjaco(1)+tmp(2)*cjaco(4)+tmp(3)*cjaco(7)
      a3 = tmp(1)*cjaco(3)+tmp(2)*cjaco(6)+tmp(3)*cjaco(9)
      coef(8) = -a1/rtod/3.6d3
      coef(7) = -a2/rtod/3.6d3
      coef(9) = -a3
      tmp(1) = coef(4)
      tmp(2) = coef(5)
      tmp(3) = coef(6)
      call getjac(x(isit1),y(isit1),z(isit1),cjaco,5)
      a1 = tmp(1)*cjaco(2)+tmp(2)*cjaco(5)+tmp(3)*cjaco(8)
      a2 = tmp(1)*cjaco(1)+tmp(2)*cjaco(4)+tmp(3)*cjaco(7)
      a3 = tmp(1)*cjaco(3)+tmp(2)*cjaco(6)+tmp(3)*cjaco(9)
      coef(5) = a1*fac
      coef(4) = a2*fac
      coef(6) = a3*fac
      call getjac(x(isit2),y(isit2),z(isit2),cjaco,5)
      a1 = tmp(1)*cjaco(2)+tmp(2)*cjaco(5)+tmp(3)*cjaco(8)
      a2 = tmp(1)*cjaco(1)+tmp(2)*cjaco(4)+tmp(3)*cjaco(7)
      a3 = tmp(1)*cjaco(3)+tmp(2)*cjaco(6)+tmp(3)*cjaco(9)
      coef(10) = -a1*fac
      coef(11) = -a2*fac
      coef(12) = -a3*fac
c
 100  continue
      return
      end
c
