      subroutine vmodel(ierr)
c
c     derive velocity field by specified velocity model
c     Here we adopt geodetic reference frame.
c
c     velocity mode : (vmode)
c         1. directly read geocentric velocity values
c         2. directly read geodetic velocity values
c         3. using model parameters (trainsition, rotation, gradient, ...)
c         4. using dislocation model
c         5. ......
c
c     mopt:
c         1 = transition
c         2 = rotation
c         3 = gradient
c         4 - 10 = ???
c
c     unit:
c         x, y, z : meter
c        vx,vy,vz : mm/year
c        ve,vn,vu : mm/year
c        wrot     : degree/year
c        dvx,dvy,dvz : mm/year/km
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
      
      integer ierr,i0,ig,i,j,isit,i2,i1
      dimension vd(4)

      ierr = 0
      pi2 = 2.0d0*pi
c
c     vmode = 3
      if (vmode.eq.3) then
      open (17,file=modfil,status='old',err=1000)
c     get group number
      read (17,*) ngp
      i0 = 0
      do 10 ig = 1,ngp
c     get site number of every group
      read (17,*) ngsit
      print*, ' group =',ig,'  site # =',ngsit
c     get model structure option
      read (17,*) (mopt(i),i=1,10)
      write (6,*) (mopt(i),i=1,10)
c     get site ID
      read (17,*) (idgs(i),i=1,ngsit)
      print*, ' sit ID:', (idgs(i),i=1,ngsit)
      if (mopt(1).gt.0) read (17,*) vxt, vyt
      if (mopt(2).gt.0) then
         read (17,*) xrc, yrc, omega
c         print*, ' rotation parameter:',xrc, yrc, omega
         xrc = xrc*dtor
         yrc = yrc*dtor
         omega = omega*dtor
      endif
      if (mopt(3).gt.0) then 
         read (17,*) xori,yori,vxori,vyori,vazi,(vd(j),j=1,4)
         write(6,'(9f6.2)') xori,yori,vxori,vyori,vazi,(vd(j),j=1,4)
         xori = xori*dtor
         yori = yori*dtor
         vazi = vazi*dtor
c        rotate to local x-y coordinate
         s1 = dsin (vazi)
         c1 = dcos (vazi)
         s2 = s1*s1
         c2 = c1*c1
         sc = s1*c1
         dvxx = vd(1)*s2-vd(2)*sc-vd(3)*sc+vd(4)*c2
         dvxy = vd(1)*sc+vd(2)*s2-vd(3)*c2-vd(4)*sc
         dvyx = vd(1)*sc-vd(2)*c2+vd(3)*s2-vd(4)*sc
         dvyy = vd(1)*c2+vd(2)*sc+vd(3)*sc+vd(4)*s2
         f = 1.0d0/finv
         e2 = 2.0*f-f*f
         rn = radius/(1.0d0-e2*dsin(yori)**2)
c         print*,'vx/dx,vx/dy,vy/dx,vy/dy: ',dvxx,dvxy,dvyx,dvyy
         wa = 0.5d0*(dvxy-dvyx)
         e12 = 0.5d0*(dvxy+dvyx)
         gam1 = dvxx-dvyy
         gam2 = dvxy+dvyx
         cita = 0.5d0*datan(-gam2/gam1)*rtod
c         print*,'e12,omega,cita: ',e12,wa,cita
      endif
c     combine all velocity factors
      do 20 isit = 1,ngsit
      i2 = idgs(isit)
      ve(i2) = 0.0d0
      vn(i2) = 0.0d0
      vu(i2) = 0.0d0
      if (mopt(1).gt.0) ve(i2) = ve(i2)+vxt
      if (mopt(1).gt.0) vn(i2) = vn(i2)+vyt
      if (mopt(2).gt.0) then
         dlo = slo(i2)-xrc
         dla = sla(i2)-yrc
         dx = radius*dcos(yrc)*dlo
         dy = radius*dla
         ve(i2) = ve(i2)-dy*omega
         vn(i2) = vn(i2)+dx*omega
      endif
      if (mopt(3).gt.0) then
         clo = slo(i2)-xori
         cla = sla(i2)-yori
         s1 = dsin(vazi)
         c1 = dcos(vazi)
         f1 = dcos(yori)
         dlo = -cla*s1*c1/f1+clo*c1*c1
         dla = cla*s1*s1-clo*s1*c1*f1
         dx = rn*dcos(yori)*dlo
         dy = rn*dla
         ve(i2) = ve(i2)+vxori+(dvxx*dx+dvxy*dy)*1.0d-3
         vn(i2) = vn(i2)+vyori+(dvyx*dx+dvyy*dy)*1.0d-3
      endif
      a1 = sla(i2)
      a2 = slo(i2)
      call sph_ca(a1,a2,ve(i2),vn(i2),vu(i2),v1,v2,v3,1)
      vx(i2) = v1
      vy(i2) = v2
      vz(i2) = v3
 20   continue
      i0 = i0+i1
 10   continue
      close (17)
      goto 500
      endif
c
c     vmode = 4
c
 500  continue
 1000 continue
      return
      end
