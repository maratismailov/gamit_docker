c

      subroutine calc_omc_full(it,isit,dt,obs,cal,omc,oc,ill)
c
c     calculate residuals with full covariance data
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
c        27. 3-D geodetic baseline vector
c        28. 3-D geodetic baseline rate vector
c
c        31. 3-D spherical coordinate
c        32. 3-D spherical frame velocity
c        33. 3-D Cartesian coordinate
c        34. 3-D Cartesian frame velocity
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : m/year
c        err      : mm
c        time     : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer i1,id,is1,isit,it,iterm,isit2,isit3,ib
      real*8  obs,omc,oc,cal,cjaco(9),coef(18)
      dimension  obs(3),omc(3),oc(3),cal(3),tempd(6)
      logical ill
      integer ifrm
      common /iframe/ifrm
c
c     from m to mm
      fac = 1.0d3
      small = 1.0d-15
      ill = .false.
      tfac = dt
c
c     3-D spherical coordinate
      if (it.eq.31) then
         is1 = jtoi(isit)
         x1 = x(isit)+vx(isit)*tfac
         y1 = y(isit)+vy(isit)*tfac
         z1 = z(isit)+vz(isit)*tfac
         call sphxyz(a2,a1,a3,x1,y1,z1,2)
         call getjac(a1,a2,a3,cjaco,1)
         iterm = 6
         ib = nlive-iaux-jaux+(ifrm-1)*6
         if (ifrm.gt.0) then
            do i1 = 1,6
               tempd(i1) = bnorm(ib+i1)
            enddo
         endif
         do i1 = 1,3
            coef(i1) = cjaco(i1+3)
            coef(i1+3) = coef(i1)*dt
         enddo
         call calc_res(1,isit,isit2,isit3,iterm,dt,coef,adj)
         obs(1) = obs(1)*rtod
         cal(1) = (a2+adj)
         if (ifrm.gt.0) then
            do i1 = 1,3
               cal(1) = cal(1)+tempd(i1)*coef(i1)
            enddo
            cal(1) = cal(1)
     .         +tempd(4)*(coef(2)*z1-coef(3)*y1)
     .         +tempd(5)*(coef(3)*x1-coef(1)*z1)
     .         +tempd(6)*(coef(1)*y1-coef(2)*x1)
         endif
         cal(1) = cal(1)*rtod
         if (obs(1).lt.0.0d0.and.cal(1).gt.0.0d0) then
            a2 = pi*0.5d0-a2
            cal(1) = 90.0d0-cal(1)
         endif
         omc(1) = (obs(1)-cal(1))*3.6d3
         oc(1) = (obs(1)-a2*rtod)*3.6d3
         iterm = 6
         do i1 = 1,3
            coef(i1) = cjaco(i1)
            coef(i1+3) = coef(i1)*dt
         enddo
         call calc_res(1,isit,isit2,isit3,iterm,dt,coef,adj)
         obs(2) = obs(2)*rtod
         cal(2) = (a1+adj)
         if (ifrm.gt.0) then
            do i1 = 1,3
               cal(2) = cal(2)+tempd(i1)*coef(i1)
            enddo
            cal(2) = cal(2)
     .         +tempd(4)*(coef(2)*z1-coef(3)*y1)
     .         +tempd(5)*(coef(3)*x1-coef(1)*z1)
     .         +tempd(6)*(coef(1)*y1-coef(2)*x1)
         endif
         cal(2) = cal(2)*rtod
         if (obs(2).lt.0.0d0.and.cal(2).gt.0.0d0) then
            a1 = a1-pi*2.0d0
            cal(2) = cal(2)-360.0d0
         endif
         omc(2) = (obs(2)-cal(2))*3.6d3
         oc(2) = (obs(2)-a1*rtod)*3.6d3
         iterm = 6
         do i1 = 1,3
            coef(i1) = cjaco(i1+6)
            coef(i1+3) = coef(i1)*dt
         enddo
         call calc_res(1,isit,isit2,isit3,iterm,dt,coef,adj)
         cal(3) = (a3+adj)
         if (ifrm.gt.0) then
            do i1 = 1,3
               cal(3) = cal(3)+tempd(i1)*coef(i1)
            enddo
c            cal(3) = cal(3)
c     .         +tempd(4)*(coef(2)*z1-coef(3)*y1)
c     .         +tempd(5)*(coef(3)*x1-coef(1)*z1)
c     .         +tempd(6)*(coef(1)*y1-coef(2)*x1)
         endif
         omc(3) = (obs(3)-cal(3))*fac
         oc(3) = (obs(3)-a3)*fac
      endif
c
c     3-D spherical velocity
c
c     3-D geocentric coordinate
      if (it.eq.33) then
         is1 = jtoi(isit)
         i1 = (isit-1)*6+1
         id = map(i1)
         adj = 0.0d0
         if (id.gt.0) adj = bnorm(id)
         if (is1.gt.0) then
            adj = adj+slnvx(is1)*tfac
         else
            adj = adj+vx(isit)*tfac
         endif
         cal(1) = x(isit)+adj
         omc(1) = (obs(1)-cal(1))*fac
         oc(1) = (obs(1)-(x(isit)+vx(isit)*tfac))*fac
         i1 = i1+1
         id = map(i1)
         adj = 0.0d0
         if (id.gt.0) adj = bnorm(id)
         if (is1.gt.0) then
            adj = adj+slnvy(is1)*tfac
         else
            adj = adj+vy(isit)*tfac
         endif
         cal(2) = y(isit)+adj
         omc(2) = (obs(2)-cal(2))*fac
         oc(2) = (obs(2)-(y(isit)+vy(isit)*tfac))*fac
         i1 = i1+1
         id = map(i1)
         adj = 0.0d0
         if (id.gt.0) adj = bnorm(id)
         if (is1.gt.0) then
            adj = adj+slnvz(is1)*tfac
         else
            adj = adj+vz(isit)*tfac
         endif
         cal(3) = z(isit)+adj
         omc(3) = (obs(3)-cal(3))*fac
         oc(3) = (obs(3)-(z(isit)+vz(isit)*tfac))*fac
      endif
c
c     3-D geocentric velocity
      if (it.eq.34) then
         is1 = jtoi(isit)
         if (is1.gt.0) then
            cal(1) = slnvx(is1)
            cal(2) = slnvy(is1)
            cal(3) = slnvz(is1)
         else
            cal(1) = vx(isit)
            cal(2) = vy(isit)
            cal(3) = vz(isit)
         endif
c        get prefit o-c
         oc(1) = obs(1)-vx(isit)
         oc(2) = obs(2)-vy(isit)
         oc(3) = obs(3)-vz(isit)
c        get postfit o-c
         omc(1) = obs(1)-cal(1)
         omc(2) = obs(2)-cal(2)
         omc(3) = obs(3)-cal(3)
         do i1 = 1,3
            obs(i1) = obs(i1)*fac
            cal(i1) = cal(i1)*fac
            oc(i1) = oc(i1)*fac
            omc(i1) = omc(i1)*fac
         enddo
      endif
c
      return
      end
