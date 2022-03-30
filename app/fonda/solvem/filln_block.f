c
      subroutine filln_block(it,isit1,isit2,l_obs,dt,ib,ibd,
     .           obs,er,rho,omc,ill)
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
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer it,isit1,isit2,ib,ibd,i,i1,i2
      integer indx,l_obs,ier,iterm,j,j1,j3,j6
      real*8  er,rho,obs,omc,covo,cjaco
      dimension  er(3),rho(3),obs(3),omc(3),covo(9)
      dimension indx(24),coef(24),cjaco(9),temp(9),cef1(24)
      logical ill
c
      ill = .false.
      small = 1.0d-15
      wght = 1.0d0
      if (dabs(er(1)).gt.small) wght = 1.0d0/er(1)/er(1)
c
c     determin the index of parameters (later)
      i1 = (isit1-1)*6
      i2 = (isit2-1)*6
      do 10 i = 1,6
         indx(i) = i1+i
         indx(i+6) = i2+i
 10   continue
c
c     3-D geocentric Cartesian coordinate
      if (it.eq.21) then
c        from mm/year to m/year
         fac = 1.0d-3
         omc(1) = obs(1)-x(isit1)-vx(isit1)*dt
         omc(2) = obs(2)-y(isit1)-vy(isit1)*dt
         omc(3) = obs(3)-z(isit1)-vz(isit1)*dt
         covo(1) = er(1)**2
         covo(3) = er(2)**2
         covo(6) = er(3)**2
         covo(2) = rho(1)*er(1)*er(2)
         covo(4) = rho(2)*er(1)*er(3)
         covo(5) = rho(3)*er(2)*er(3)
         call cholsk(covo,cjaco,1,3,ier)
         cjaco(1) = covo(1)*dt
         cjaco(2) = covo(2)*dt
         cjaco(3) = covo(4)*dt
         cjaco(4) = covo(2)*dt
         cjaco(5) = covo(3)*dt
         cjaco(6) = covo(5)*dt
         cjaco(7) = covo(4)*dt
         cjaco(8) = covo(5)*dt
         cjaco(9) = covo(6)*dt
         call latwa(1,3,omc,covo,chi2,temp,1,2)
c
c        right hand terms
         call axb(3,3,1,covo,omc,coef,2,0)
         do i = 1,3
            i1 = indx(i)
            bnorm(i1) = bnorm(i1)+coef(i)
            bnorm(i1+3) = bnorm(i1+3)+coef(i)*dt
         enddo
         i1 = indx(1)
         call submtx(i1,i1,3,3,anorm,covo,1,3)
         call submtx(i1+3,i1,3,3,anorm,cjaco,2,3)
         do i = 1,6
            covo(i) = covo(i)*dt*dt
         enddo
         call submtx(i1+3,i1+3,3,3,anorm,covo,1,3)
c        unit: mm
         omc(1) = omc(1)*1.0d3
         omc(2) = omc(2)*1.0d3
         omc(3) = omc(3)*1.0d3
         goto 100
      endif
c
c     3-D geocentric velocity
      if (it.eq.22) then
c        from mm/year to m/year
         fac = 1.0d-3
         omc(1) = obs(1)-vx(isit1)/fac
         omc(2) = obs(2)-vy(isit1)/fac
         omc(3) = obs(3)-vz(isit1)/fac
         covo(1) = er(1)**2
         covo(3) = er(2)**2
         covo(6) = er(3)**2
         covo(2) = rho(1)*er(1)*er(2)
         covo(4) = rho(2)*er(1)*er(3)
         covo(5) = rho(3)*er(2)*er(3)
         call cholsk(covo,cjaco,1,3,ier)
         call latwa(1,3,omc,covo,chi2,temp,1,2)
         call axb(3,3,1,covo,omc,coef,2,0)
         do i = 1,3
            i1 = indx(i+3)
            bnorm(i1) = bnorm(i1)+coef(i)*fac
         enddo
         i1 = indx(4)
         call submtx(i1,i1,3,3,anorm,covo,1,3)
         goto 100
      endif
c
c     3-D geodetic velocity
      if (it.eq.24) then
c        from mm/year to m/year 
         fac = 1.0d-3
         omc(1) = obs(1)*fac-ve(isit1)
         omc(2) = obs(2)*fac-vn(isit1)
         omc(3) = obs(3)*fac-vu(isit1)
         covo(1) = er(1)**2
         covo(3) = er(2)**2
         covo(6) = er(3)**2
         covo(2) = rho(1)*er(1)*er(2)
         covo(4) = rho(2)*er(1)*er(3)
         covo(5) = rho(3)*er(2)*er(3)
         call cholsk(covo,coef,1,3,ier)
         call latwa(1,3,omc,covo,chi2,temp,1,2)
         a1 = slon(isit1)
         a2 = slat(isit1)
         a3 = radius
         call getjac(a1,a2,a3,coef,6)
c
c        right hand terms
         call axb(3,3,3,coef,covo,cjaco,3,0)
         call axb(3,3,1,cjaco,omc,temp,1,0)
         do i = 1,3
            i1 = indx(i+3)
            bnorm(i1) = bnorm(i1)+temp(i)
         enddo
         call atwa(3,3,coef,covo,cjaco,temp,1)
         i1 = indx(4)
         call submtx(i1,i1,3,3,anorm,cjaco,1,3)
         omc(1) = omc(1)*1.0d3
         omc(2) = omc(2)*1.0d3
         omc(3) = omc(3)*1.0d3
         goto 100
      endif
c
c     3-D Cartesian baseline vector
      if (it.eq.25) then
         fac = 1.0d0*dt
c        get o-c
         omc(1) = obs(1)-(x(isit2)-x(isit1)+
     .            (vx(isit2)-vx(isit1))*fac)
         omc(2) = obs(2)-(y(isit2)-y(isit1)+
     .            (vy(isit2)-vy(isit1))*fac)
         omc(3) = obs(3)-(z(isit2)-z(isit1)+
     .            (vz(isit2)-vz(isit1))*fac)
c        print* ,' in filln_block'
c        print* ,'obs = ',obs
c        print*, sname(isit1),'--',sname(isit2)
c        print* ,'x1',isit1,' ',x(isit1)
c        print* ,'x2',isit2,' ',x(isit2)
c        print* ,'y1',isit1,' ',y(isit1)
c        print* ,'y2',isit2,' ',y(isit2)
c        print* ,'z1',isit1,' ',z(isit1)
c        print* ,'z2',isit2,' ',z(isit2)
c        print* ,'vx1',isit1,' ',vx(isit1)
c        print* ,'vx2',isit2,' ',vx(isit2)
c        print* ,'vy1',isit1,' ',vy(isit1)
c        print* ,'vy2',isit2,' ',vy(isit2)
c        print* ,'vz1',isit1,' ',vz(isit1)
c        print* ,'vz2',isit2,' ',vz(isit2)
c        print* , 'fac', fac
c        print* ,'omc ', omc
         covo(1) = er(1)**2
         covo(3) = er(2)**2
         covo(6) = er(3)**2
         covo(2) = rho(1)*er(1)*er(2)
         covo(4) = rho(2)*er(1)*er(3)
         covo(5) = rho(3)*er(2)*er(3)
         call cholsk(covo,cjaco,1,3,ier)
         call latwa(1,3,omc,covo,chi2,temp,1,2)
         cjaco(1) = -covo(1)
         cjaco(2) = -covo(2)
         cjaco(3) = -covo(4)
         cjaco(4) = -covo(2)
         cjaco(5) = -covo(3)
         cjaco(6) = -covo(5)
         cjaco(7) = -covo(4)
         cjaco(8) = -covo(5)
         cjaco(9) = -covo(6)
         iterm = 12
         do i = 1,3
            cef1(i) = -1.0
            cef1(i+3) = -dt
            cef1(i+6) = 1.0
            cef1(i+9) = dt
         enddo
         if (iomode(3).gt.0.and.iq_sit.gt.0)
     .      call add_quake_c(isit1,isit2,iterm,dt,cef1,indx)
c        right hand terms
         call axb(3,3,1,cjaco,omc,coef,1,0)
         do i = 1,3
            i1 = indx(i)
            i2 = indx(i+6)
            bnorm(i1) = bnorm(i1)+coef(i)
            bnorm(i1+3) = bnorm(i1+3)+coef(i)*dt
            bnorm(i2) = bnorm(i2)-coef(i)
            bnorm(i2+3) = bnorm(i2+3)-coef(i)*dt
            if (iterm.gt.12) then
               i1 = indx(i+12)
               bnorm(i1) = bnorm(i1)-coef(i)*cef1(i+12)
            endif
            if (iterm.gt.15) then
               i1 = indx(i+15)
               bnorm(i1) = bnorm(i1)-coef(i)*cef1(i+15)
            endif
         enddo
         i1 = indx(1)
         i2 = indx(7)
         call submtx(i1,i1,3,3,anorm,covo,1,3)
         call submtx(i2,i2,3,3,anorm,covo,1,3)
         call submtx(i2,i1,3,3,anorm,cjaco,2,3)
         if (iterm.gt.12) then
            do i = 1,3
               do j = 1,3
                  j1 = (i-1)*3+j
                  cef1(j1) = -cjaco(j1)*cef1(12+i)*cef1(12+j)
                  if (j.le.i) temp(i*(i-1)/2+j) = cef1(j1)
               enddo
            enddo
            j3 = indx(13)
            call submtx(j3,j3,3,3,anorm,temp,1,3)
         endif
         if (iterm.gt.15) then
            j6 = indx(16)
            call submtx(j6,j6,3,3,anorm,temp,1,3)
            do i = 1,9
               cef1(i) = -cef1(i)
            enddo
            if (j6.gt.j3)
     .      call submtx(j6,j3,3,3,anorm,cef1,2,3)
            if (j3.gt.j6)
     .      call submtx(j3,j6,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.12) then
            do i = 1,3
               do j = 1,3
                  j1 = (i-1)*3+j
                  cef1(j1) = cjaco(j1)*cef1(12+j)
               enddo
            enddo
            call submtx(j3,i1,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.15) then
            call submtx(j6,i2,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.12) then
            do i = 1,9
               cef1(i) = -cef1(i)
            enddo
            call submtx(j3,i2,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.15) then
            call submtx(j6,i1,3,3,anorm,cef1,2,3)
         endif
         do i = 1,9
            cjaco(i) = cjaco(i)*dt
         enddo
         call submtx(i2+3,i1,3,3,anorm,cjaco,2,3)
         call submtx(i2,i1+3,3,3,anorm,cjaco,2,3)
         do i = 1,9
            cjaco(i) = -cjaco(i)
         enddo
         call submtx(i1+3,i1,3,3,anorm,cjaco,2,3)
         call submtx(i2+3,i2,3,3,anorm,cjaco,2,3)
         if (iterm.gt.12) then
            do i = 1,9
               cef1(i) = -cef1(i)*dt
            enddo
            call submtx(j3,i1+3,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.15) then
            call submtx(j6,i2+3,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.12) then
            do i = 1,9
               cef1(i) = -cef1(i)
            enddo
            call submtx(j3,i2+3,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.15) then
            call submtx(j6,i1+3,3,3,anorm,cef1,2,3)
         endif
         do i = 1,6
            covo(i) = covo(i)*dt*dt
         enddo
         call submtx(i1+3,i1+3,3,3,anorm,covo,1,3)
         call submtx(i2+3,i2+3,3,3,anorm,covo,1,3)
         do i = 1,9
            cjaco(i) = -cjaco(i)*dt
         enddo
         call submtx(i2+3,i1+3,3,3,anorm,cjaco,2,3)
         goto 100
      endif        
c
c     3-D spherical baseline vector
      if (it.eq.27) then
         fac = 1.0d0*dt
         x1 = x(isit1)
         y1 = y(isit1)
         z1 = z(isit1)
         call sphxyz(a2,a1,a3,x1,y1,z1,2)
         call getjac(a1,a2,a3,cjaco,5)
         omc(1) = x(isit2)-x(isit1)+(vx(isit2)-vx(isit1))*fac
         omc(2) = y(isit2)-y(isit1)+(vy(isit2)-vy(isit1))*fac
         omc(3) = z(isit2)-z(isit1)+(vz(isit2)-vz(isit1))*fac
         call axb(3,3,1,cjaco,omc,covo,1,0)
         omc(1) = obs(1)-covo(1)
         omc(2) = obs(2)-covo(2)
         omc(3) = obs(3)-covo(3)
         covo(1) = er(1)**2
         covo(3) = er(2)**2
         covo(6) = er(3)**2
         covo(2) = rho(1)*er(1)*er(2)
         covo(4) = rho(2)*er(1)*er(3)
         covo(5) = rho(3)*er(2)*er(3)
         call cholsk(covo,coef,1,3,ier)
         call latwa(1,3,omc,covo,chi2,temp,1,2)
         iterm = 12
         do i = 1,3
            cef1(i) = -1.0
            cef1(i+3) = -dt
            cef1(i+6) = 1.0
            cef1(i+9) = dt
         enddo
         if (iomode(3).gt.0.and.iq_sit.gt.0)
     .      call add_quake_c(isit1,isit2,iterm,dt,cef1,indx)
         call getjac(a1,a2,a3,coef,6)
         call axb(3,3,3,coef,covo,cjaco,3,0)
         call axb(3,3,1,cjaco,omc,temp,1,0)
         do i = 1,3
            i1 = indx(i)
            i2 = indx(i+6)
            bnorm(i1) = bnorm(i1)-temp(i)
            bnorm(i1+3) = bnorm(i1+3)-temp(i)*dt
            bnorm(i2) = bnorm(i2)+temp(i)
            bnorm(i2+3) = bnorm(i2+3)+temp(i)*dt
            if (iterm.gt.12) then
               i1 = indx(i+12)
               bnorm(i1) = bnorm(i1)+temp(i)*cef1(i+12)
            endif
            if (iterm.gt.15) then
               i1 = indx(i+15)
               bnorm(i1) = bnorm(i1)+temp(i)*cef1(i+15)
            endif
         enddo
         call atwa(3,3,coef,covo,cjaco,temp,1)
         do i = 1,6
            covo(i) = cjaco(i)
         enddo
         cjaco(1) = covo(1)
         cjaco(2) = covo(2)
         cjaco(3) = covo(4)
         cjaco(4) = covo(2)
         cjaco(5) = covo(3)
         cjaco(6) = covo(5)
         cjaco(7) = covo(4)
         cjaco(8) = covo(5)
         cjaco(9) = covo(6)
         i1 = indx(1)
         i2 = indx(7)
         call submtx(i1,i1,3,3,anorm,covo,1,3)
         call submtx(i2,i2,3,3,anorm,covo,1,3)
         if (iterm.gt.12) then
            do i = 1,3
               do j = 1,3
                  j1 = (i-1)*3+j
                  cef1(j1) = cjaco(j1)*cef1(12+i)*cef1(12+j)
               enddo
            enddo
            j = indx(13)
            call submtx(j,j,3,3,anorm,cef1,1,3)
         endif
         if (iterm.gt.15) then
            j1 = indx(16)
            call submtx(j1,j1,3,3,anorm,cef1,1,3)
            do i = 1,9
               cef1(i) = -cef1(i)
            enddo
            if (j1.gt.j)
     .      call submtx(j1,j,3,3,anorm,cef1,2,3)
            if (j.gt.j1)
     .      call submtx(j,j1,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.12) then
            do i = 1,3
               do j = 1,3
                  j1 = (i-1)*3+j
                  cef1(j1) = -cjaco(j1)*cef1(12+j)
               enddo
            enddo
            call submtx(j,i1,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.15) then
            call submtx(j1,i2,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.12) then
            do i = 1,9
               cef1(i) = -cef1(i)
            enddo
            call submtx(j,i2,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.15) then
            call submtx(j1,i1,3,3,anorm,cef1,2,3)
         endif
         do i = 1,9
            temp(i) = -cjaco(i)
         enddo
         call submtx(i2,i1,3,3,anorm,temp,2,3)
         do i = 1,9
            temp(i) = temp(i)*dt
         enddo
         call submtx(i2+3,i1,3,3,anorm,temp,2,3)
         call submtx(i2,i1+3,3,3,anorm,temp,2,3)
         do i = 1,9
            cjaco(i) = -temp(i)
         enddo
         call submtx(i1+3,i1,3,3,anorm,cjaco,2,3)
         call submtx(i2+3,i2,3,3,anorm,cjaco,2,3)
         if (iterm.gt.12) then
            do i = 1,9
               cef1(i) = -cef1(i)*dt
            enddo
            call submtx(j,i1+3,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.15) then
            call submtx(j1,i2+3,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.12) then
            do i = 1,9
               cef1(i) = -cef1(i)
            enddo
            call submtx(j,i2+3,3,3,anorm,cef1,2,3)
         endif
         if (iterm.gt.15) then
            call submtx(j1,i1+3,3,3,anorm,cef1,2,3)
         endif
         do i = 1,6
            covo(i) = covo(i)*dt*dt
         enddo
         call submtx(i1+3,i1+3,3,3,anorm,covo,1,3)
         call submtx(i2+3,i2+3,3,3,anorm,covo,1,3)
         do i = 1,9
            cjaco(i) = -cjaco(i)*dt
         enddo
         call submtx(i2+3,i1+3,3,3,anorm,cjaco,2,3)
         do i = 1,3
            omc(i) = omc(i)*1.0d3
         enddo
         goto 100
      endif        
c
c     3-D spherical baseline rate vector
      if (it.eq.28) then
         fac = 1.0d3
         x1 = x(isit1)
         y1 = y(isit1)
         z1 = z(isit1)
         a1 = slon(isit1)
         a2 = slat(isit1)
         a3 = radius
         call getjac(a1,a2,a3,coef,6)
         omc(2) = obs(2)-(ve(isit2)-ve(isit1))
         omc(1) = obs(1)-(vn(isit2)-vn(isit1))
         omc(3) = obs(3)-(vu(isit2)-vu(isit1))
         covo(1) = er(1)**2
         covo(3) = er(2)**2
         covo(6) = er(3)**2
         covo(2) = rho(1)*er(1)*er(2)
         covo(4) = rho(2)*er(1)*er(3)
         covo(5) = rho(3)*er(2)*er(3)
         call cholsk(covo,coef,1,3,ier)
         call latwa(1,3,omc,covo,chi2,temp,1,2)
         call getjac(a1,a2,a3,coef,6)
         do i = 1,3
            tt = coef(i*3-2)
            coef(i*3-2) = coef(i*3-1)
            coef(i*3-1) = tt
         enddo
         call axb(3,3,3,coef,covo,cjaco,3,0)
         call axb(3,3,1,cjaco,omc,temp,1,0)
         do i = 1,3
            i1 = indx(i)
            i2 = indx(i+6)
            bnorm(i1+3) = bnorm(i1+3)-temp(i)
            bnorm(i2+3) = bnorm(i2+3)+temp(i)
         enddo
         call atwa(3,3,coef,covo,cjaco,temp,1)
         do i = 1,6
            covo(i) = cjaco(i)
         enddo
         i1 = indx(1)
         i2 = indx(7)
         call submtx(i1+3,i1+3,3,3,anorm,covo,1,3)
         call submtx(i2+3,i2+3,3,3,anorm,covo,1,3)
         cjaco(1) = -covo(1)
         cjaco(2) = -covo(2)
         cjaco(3) = -covo(4)
         cjaco(4) = -covo(2)
         cjaco(5) = -covo(3)
         cjaco(6) = -covo(5)
         cjaco(7) = -covo(4)
         cjaco(8) = -covo(5)
         cjaco(9) = -covo(6)
         call submtx(i2+3,i1+3,3,3,anorm,cjaco,2,3)
         do i = 1,3
            omc(i) = omc(i)*1.0d3
         enddo
         goto 100
      endif
c
c     kick out outliers
 120  if (iomode(10).gt.0) then
         continue
      else
         
         print 123,ib,it,sname(isit1),isit1,sname(isit2),
     .   isit2,obs(1),omc(1)
 123  format (' FILLN: outlier at obs. # ',i4,' type ',i3,1x,
     . a8,'(',i4,') to ',a8,'(',i4,')',' obs = ',
     . f20.4,' omc = ',1pg12.4)
      endif
      ill = .true.
 122  format(2x,i5,4x,i3,2i5,f10.3,3x,f12.4)
c
 100  continue
      return
      end
c
