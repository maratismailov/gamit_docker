c--------------------------------------------------------------
      subroutine calc_omc_block(it,isit1,isit2,dt,ib,ibd,obs,cal,
     .                          omc,oc,ill)
c
c     data type: (it)
c         1. astrometric azimuth
c         2. horizontal angle
c         3. horizontal direction
c         4. baseline length
c         5. zenith height
c         6. leveling
c
c        11. astrometric azimuth rate
c        12. horizontal angle rate
c        13. horizontal direction rate
c        14. baseline length rate
c        15. zenith height rate
c        16. leveling rate
c
c        21. 3-D geocentric coordinate
c        22. 3-D geocentric velocity
c        23. 3-D geodetic coordinate
c        24. 3-D geodetic velocity
c        25. 3-D geocentric baseline vector
c        26. 3-D geocentric baseline rate vector
c        27. 3-D spherical baseline vector
c        28. 3-D spherical baseline rate vector
c
c        31. 3-D spherical coordinate
c        32. 3-D spherical frame velocity
c        33. 3-D Cartesian coordinate
c        34. 3-D Cartesian frame velocity
c
c     variables:
c        obs:  observation data
c        cal:  calculated velues
c        omc:  postfit
c        oc :  prefit
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer it,isit1,isit2,ib,ibd
      integer i1,i2,id,isit3,i,is1,is2,iterm
      dimension tempc(18),omc(3),oc(3),obs(3),cal(3)
      dimension coef(30),atmp(3),btmp(3)
      logical ill
c
c     from m to mm
      fac = 1.0d3
      small = 1.0d-15
      ill = .false.
c
c     3-D geocentric coordinate
      if (it.eq.21) then
c        i2 = ibd-ibd/3*3
c        if (i2.eq.0) i2 = 3
       do 10 i2= 1,3
         i1 = (isit1-1)*6+i2
         id = map(i1)
         adj = 0.0d0
         if (id.gt.0) adj = bnorm(id)
         if (i2.eq.1) cal(1) = x(isit1)+adj
         if (i2.eq.2) cal(2) = y(isit1)+adj
         if (i2.eq.3) cal(3) = z(isit1)+adj
         omc(i2) = (obs(i2)-cal(i2))*fac
         oc(i2) = omc(i2)+adj*fac
 10    continue 
      endif
c
c     3-D geocentric velocity
      if (it.eq.22) then
         is1 = jtoi(isit1)
         if (is1.gt.0) cal(1) = slnvx(is1)*fac
         if (is1.gt.0) cal(2) = slnvy(is1)*fac
         if (is1.gt.0) cal(3) = slnvz(is1)*fac
         if (is1.le.0) cal(1) = vx(isit1)*fac
         if (is1.le.0) cal(2) = vy(isit1)*fac
         if (is1.le.0) cal(3) = vz(isit1)*fac
         oc(1) = obs(1)-vx(isit1)*fac
         oc(2) = obs(2)-vy(isit1)*fac
         oc(3) = obs(3)-vz(isit1)*fac
         omc(1) = obs(1)-cal(1)
         omc(2) = obs(2)-cal(2)
         omc(3) = obs(3)-cal(3)
      endif
c
c     3-D geodetic velocity
      if (it.eq.24) then
         is1 = jtoi(isit1)
         if (is1.gt.0) cal(1) = slnve(is1)*fac
         if (is1.gt.0) cal(2) = slnvn(is1)*fac
         if (is1.gt.0) cal(3) = slnvu(is1)*fac
         if (is1.le.0) cal(1) = ve(isit1)*fac
         if (is1.le.0) cal(2) = vn(isit1)*fac
         if (is1.le.0) cal(3) = vu(isit1)*fac
         oc(1) = obs(1)-ve(isit1)*fac
         oc(2) = obs(2)-vn(isit1)*fac
         oc(3) = obs(3)-vu(isit1)*fac
         omc(1) = obs(1)-cal(1)
         omc(2) = obs(2)-cal(2)
         omc(3) = obs(3)-cal(3)
      endif
c
c     3-D Cartesian baseline vector
      if (it.eq.25) then
         do 50 i = 1,3
            call zero1d(1,18,coef)
            coef(6+i) = 1.0d0
            coef(i) = -1.0d0
            coef(9+i) = dt
            coef(3+i) = -dt
            iterm = 12
            call calc_res(2,isit1,isit2,isit3,iterm,dt,coef,cal(i))
            if (i.eq.1) tmp = x(isit2)-x(isit1)+
     .         (vx(isit2)-vx(isit1))*dt
            if (i.eq.2) tmp = y(isit2)-y(isit1)+
     .         (vy(isit2)-vy(isit1))*dt
            if (i.eq.3) tmp = z(isit2)-z(isit1)+
     .         (vz(isit2)-vz(isit1))*dt
            omc(i) = (obs(i)-tmp-cal(i))*fac
            oc(i)  = (obs(i)-tmp)*fac
            cal(i) = obs(i)-omc(i)/fac
 50      continue
      endif
c
c     3-D spherical baseline vector
      if (it.eq.27) then
         call getadj(2,isit1,isit2,isit3,tempc)
         x1 = x(isit1)
         y1 = y(isit1)
         z1 = z(isit1)
         call sphxyz(a2,a1,a3,x1,y1,z1,2)
         call getjac(a1,a2,a3,coef,5)
         omc(1) = (x(isit2)-x(isit1)+(vx(isit2)-vx(isit1))*dt)*fac
         omc(2) = (y(isit2)-y(isit1)+(vy(isit2)-vy(isit1))*dt)*fac
         omc(3) = (z(isit2)-z(isit1)+(vz(isit2)-vz(isit1))*dt)*fac
         call axb(3,3,1,coef,omc,oc,1,0)
         do i = 1,3
            oc(i) = obs(i)*fac-oc(i)
         enddo
         atmp(1) = tempc(7)-tempc(1)+(tempc(10)-tempc(4))*dt
         atmp(2) = tempc(8)-tempc(2)+(tempc(11)-tempc(5))*dt
         atmp(3) = tempc(9)-tempc(3)+(tempc(12)-tempc(6))*dt
         if (iomode(3).gt.0.and.iq_sit.gt.0) then
            call get_qk_adj(isit1,i1,dt,tempc)
            if (i1.gt.0) then
               do i = 1,3
                  atmp(i) = atmp(i)-tempc(i)
               enddo
            endif
            call get_qk_adj(isit2,i1,dt,tempc)
            if (i1.gt.0) then
               do i = 1,3
                  atmp(i) = atmp(i)+tempc(i)
               enddo
            endif
         endif
         call axb(3,3,1,coef,atmp,btmp,1,0)
         do i = 1,3
            omc(i) = oc(i)-btmp(i)*fac
            cal(i) = obs(i)-omc(i)/fac
         enddo
      endif
c
c     3-D geodetic baseline rate velocity
      if (it.eq.28) then
         is2 = jtoi(isit2)
         if (is2.gt.0) cal(1) = slnvn(is2)*fac
         if (is2.gt.0) cal(2) = slnve(is2)*fac
         if (is2.gt.0) cal(3) = slnvu(is2)*fac
         if (is2.le.0) cal(1) = vn(isit2)*fac
         if (is2.le.0) cal(2) = ve(isit2)*fac
         if (is2.le.0) cal(3) = vu(isit2)*fac
         is1 = jtoi(isit1)
         if (is1.gt.0) cal(1) = cal(1)-slnvn(is1)*fac
         if (is1.gt.0) cal(2) = cal(2)-slnve(is1)*fac
         if (is1.gt.0) cal(3) = cal(3)-slnvu(is1)*fac
         if (is1.le.0) cal(1) = cal(1)-vn(isit1)*fac
         if (is1.le.0) cal(2) = cal(2)-ve(isit1)*fac
         if (is1.le.0) cal(3) = cal(3)-vu(isit1)*fac
         oc(1) = (obs(1)-vn(isit2)+vn(isit1))*fac
         oc(2) = (obs(2)-ve(isit2)+ve(isit1))*fac
         oc(3) = (obs(3)-vu(isit2)+vu(isit1))*fac
         obs(1) = obs(1)*fac
         obs(2) = obs(2)*fac
         obs(3) = obs(3)*fac
         omc(1) = obs(1)-cal(1)
         omc(2) = obs(2)-cal(2)
         omc(3) = obs(3)-cal(3)
      endif
c
 100  continue
      return
      end
c
