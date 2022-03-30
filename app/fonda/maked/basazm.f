      real*8 function basazm(isit1,isit2,dt,idim)
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

c     calculate azimuth from sit1 to sit2. unit=radian
c     azimuth is insensitive to the height change
  
      integer isit1,isit2,idim,iek

      factr = 1.0d-3
      dtf = dt*factr
c
c     geodetic coordinate
c     
      vxdt = (vx(isit2)-vx(isit1))*dtf
      vydt = (vy(isit2)-vy(isit1))*dtf
      vzdt = (vz(isit2)-vz(isit1))*dtf
      xd = x(isit2)-x(isit1)+vxdt
      yd = y(isit2)-y(isit1)+vydt
      zd = z(isit2)-z(isit1)+vzdt
      s1 = dsin(slo(isit1))
      c1 = dcos(slo(isit1))
      s2 = dsin(sla(isit1))
      c2 = dcos(sla(isit1))
      dx = -xd*s1+yd*c1
      dy = -xd*c1*s2-yd*s1*s2+zd*c2
  
c     episodic coordinate jump correction
      if (cmode.eq.2) then
         call chk_quake_list(isit1,dt,iek,fac)
         if (iek.gt.0) then
           dx = dx-quake_ce(iek)*fac
           dy = dy-quake_cn(iek)*fac
         endif
         call chk_quake_list(isit2,dt,iek,fac)
         if (iek.gt.0) then
           dx = dx+quake_ce(iek)*fac
           dy = dy+quake_cn(iek)*fac
         endif
      endif
      
      basazm = datan2(dx, dy)
      if (basazm .lt. 0.0d0) basazm = basazm+pi*2.0d0

      return
      end
c--------------------------------------------------------------
      real*8 function baslen(isit1,isit2,dt,idim,mode)
c
c     idim : 2-D or 3-D
c     mode = 1: geodetic coordinate
c     mode = 2: Cartesian coordinate
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
      logical jump
      integer isit1,isit2,idim,iek,mode

      factr = 1.0d-3
      dtf = dt*factr
      jump = .false.
c
c     geodetic coordinate
      if (mode.eq.1) then
         vxdt = (vx(isit2)-vx(isit1))*dtf
         vydt = (vy(isit2)-vy(isit1))*dtf
         vzdt = (vz(isit2)-vz(isit1))*dtf
         xd = x(isit2)-x(isit1)+vxdt
         yd = y(isit2)-y(isit1)+vydt
         zd = z(isit2)-z(isit1)+vzdt
         s1 = dsin(slo(isit1))
         c1 = dcos(slo(isit1))
         s2 = dsin(sla(isit1))
         c2 = dcos(sla(isit1))
         dx = -xd*s1+yd*c1
         dy = -xd*c1*s2-yd*s1*s2+zd*c2
         dz = 0.0d0
         if (idim.eq.3) dz = xd*c1*c2+yd*s1*c2+zd*s2
      endif
c
c     geocentric coordinate
      if (mode.eq.2) then
         x1 = x(isit1)+vx(isit1)*dtf
         y1 = y(isit1)+vy(isit1)*dtf
         x2 = x(isit2)+vx(isit2)*dtf
         y2 = y(isit2)+vy(isit2)*dtf
         dx = x2-x1
         dy = y2-y1
         dz = 0.0d0
         if (idim.eq.3) then
            z1 = z(isit1)+vz(isit1)*dtf
            z2 = z(isit2)+vz(isit2)*dtf
            dz = z2-z1
         endif
      endif
  
c     episodic coordinate jump correction
      if (cmode.eq.2) then
         s1 = dsin(slo(isit1))
         c1 = dcos(slo(isit1))
         s2 = dsin(sla(isit1))
         c2 = dcos(sla(isit1))
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         call chk_quake_list(isit1,dt,iek,fac)
         if (iek.gt.0) then
           xd = xd-quake_ce(iek)*fac
           yd = yd-quake_cn(iek)*fac
           zd = zd-quake_cu(iek)*fac
           jump = .true.
         endif
         call chk_quake_list(isit2,dt,iek,fac)
         if (iek.gt.0) then
           xd = xd+quake_ce(iek)*fac
           yd = yd+quake_cn(iek)*fac
           zd = zd+quake_cu(iek)*fac
           jump = .true.
         endif
         if (jump) then
            if (mode.eq.1) then
               dx = dx+xd
               dy = dy+yd
               dz = dz+zd
            else
               dx = dx-xd*s1-yd*c1*s2+zd*c1*c2
               dy = dy+xd*c1-yd*s1*s2+zd*c2*s1
               if (idim.eq.3) dz = dz+yd*c2+zd*s2
            endif
         endif
      endif
c
      baslen = dsqrt(dx**2+dy**2+dz**2)
      
      return
      end
c     

