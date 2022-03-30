      real*8 function basazm(isit1,isit2,dt,mode)
c
c     mode = 1:  azimuth observable in local astronomical frame
c     mode = 2:  azimuth observable in geodetic frame
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer isit1,isit2,mode
      real*8 lat1,lon1,lat2,lon2

c     calculate azimuth from sit1 to sit2. unit = radian
c     

      factr = 1.0d0
      dtf = dt*factr

c
c     first, get vector components in Cartesian coordinate
c     vxdt = (vx(isit2)-vx(isit1))*dtf
c     vydt = (vy(isit2)-vy(isit1))*dtf
c     vzdt = (vz(isit2)-vz(isit1))*dtf
c     xd = x(isit2)-x(isit1)+vxdt
c     yd = y(isit2)-y(isit1)+vydt
c     zd = z(isit2)-z(isit1)+vzdt
      
c     following does the same thing, but need coords at measurement epoch below
      x1 = x(isit1)+vx(isit1)*dtf
      y1 = y(isit1)+vy(isit1)*dtf
      z1 = z(isit1)+vz(isit1)*dtf
      x2 = x(isit2)+vx(isit2)*dtf
      y2 = y(isit2)+vy(isit2)*dtf
      z2 = z(isit2)+vz(isit2)*dtf
      xd = x2-x1
      yd = y2-y1
      zd = z2-z1

c     print*,'x,y,z : ',x1,y1,z1,x2,y2,z2


c
c     second, rotate to local coordinate
c     consider deflection correction
cmk   this seems to be computed at the point prior to dt being applied

      call geoxyz(radius,finv,tx,ty,tz,lat1,lon1,dh1,x1,y1,z1,2,he1)

      sla = lat1+defn(isit1)
      slo = lon1+defe(isit1)/dcos(lat1)
      s1 = dsin(slo)
      c1 = dcos(slo)
      s2 = dsin(sla)
      c2 = dcos(sla)
      dx = -xd*s1+yd*c1
      dy = -xd*c1*s2-yd*s1*s2+zd*c2
      dl1 = dsqrt(dx*dx+dy*dy)

      if (mode.eq.2) then
         dtf = dtf/radius
         c1 = dcos(lat1)
         vedt = (ve(isit2)/dcos(lat2)-ve(isit1)/c1)*dtf
         vndt1 = vn(isit1)*dtf
         vndt2 = vn(isit2)*dtf
         del = lon2-lon1+vedt
         dx = dsin(del)
         c1 = dcos(del)
         s1 = dsin(slat(isit1)+vndt1)
         e2 = 2.0d0/finv-1.0d0/finv**2
         c3 = (1.0d0-e2)*dtan(lat2)+
     .      e2*s1/dcos(lat2)
         c2 = dtan(lat2)+vndt2
         c3 = c3*dsqrt((1.d0-e2*dsin(lat2)**2)/(1.d0-e2*s1**2))
         dy = dcos(lat1+vndt1)*c2-s1*c1

         basazm = datan2(dx, dy)
         
         ceta = defn(isit1)*dsin(basazm)-defe(isit1)*dcos(basazm)
         dazm = -ceta*dtor*(srad(isit2)-srad(isit1))/dl1
         basazm = basazm+dazm

         if (basazm.lt.0.d0) basazm = basazm+pi*2.0d0

         return
      endif

      basazm = datan2(dx, dy)
      if (basazm.lt.0.d0) basazm = basazm+pi*2.0d0

c     print*,' basazm: ',basazm

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
      include 'solvem.fti'
      integer isit1,isit2,idim,mode
      real*8 lat1,lon1

      factr = 1.0d0
      dtf = dt*factr
c
c     geodetic coordinate
      if (mode.eq.1) then
c        first, get vector components in Cartesian coordinate
c        vxdt = (vx(isit2)-vx(isit1))*dtf
c        vydt = (vy(isit2)-vy(isit1))*dtf
c        vzdt = (vz(isit2)-vz(isit1))*dtf
c        xd = x(isit2)-x(isit1)+vxdt
c        yd = y(isit2)-y(isit1)+vydt
c        zd = z(isit2)-z(isit1)+vzdt

c     following does the same thing, but need coords at measurement epoch below
         x1 = x(isit1)+vx(isit1)*dtf
         y1 = y(isit1)+vy(isit1)*dtf
         z1 = z(isit1)+vz(isit1)*dtf
         x2 = x(isit2)+vx(isit2)*dtf
         y2 = y(isit2)+vy(isit2)*dtf
         z2 = z(isit2)+vz(isit2)*dtf
         xd = x2-x1
         yd = y2-y1
         zd = z2-z1

         call geoxyz(radius,finv,tx,ty,tz,lat1,lon1,dh1,x1,y1,z1,2,he1)

c        second, rotate to local coordinate
         sla = lat1        
         slo = lon1 
c        put deflection correction
         if (idim.eq.2) then
         sla = lat1+defn(isit1)
         slo = lon1+defe(isit1)/dcos(lat1)
         endif
         s1 = dsin(slo)
         c1 = dcos(slo)
         s2 = dsin(sla)
         c2 = dcos(sla)
         dx = -xd*s1+yd*c1
         dy = -xd*c1*s2-yd*s1*s2+zd*c2
         dz = 0.0d0
         if (idim.eq.3) dz = xd*c1*c2+yd*s1*c2+zd*s2
      endif

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
c
c     This is mark-to-mark chord distance 
      baslen = dsqrt(dx**2+dy**2+dz**2)
      

      return
      end
c
