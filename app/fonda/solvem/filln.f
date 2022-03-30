c
      subroutine filln(it,isit1,isit2,dt,ib,ibd,
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
c        27. 3-D geodetic baseline vector
c        28. 3-D geodetic baseline rate vector
c
c        31. 3-D spherical coordinate
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer it,isit1,isit2,ib,ibd,iterm,i,i1,i2,istsav,indsav
      integer indx
      real*8  er,rho,obs,omc
      dimension indx(30),coef(30)
      dimension indsav(12),coesav(12)
      logical ill
      common/saved/azisav,coesav,obssav,esave,tsave,istsav,indsav
c
      ill = .false.
      small = 1.0d-15
      wght = 1.0d0
      if (dabs(er).gt.small) wght = 1.0d0/er/er
c
c     determin the index of parameters (later)
      iterm = 12+jaux
      i1 = (isit1-1)*6
      i2 = (isit2-1)*6
      do 10 i = 1,6
         indx(i) = i1+i
         indx(i+6) = i2+i
 10   continue
      if (pmode.eq.3) then
         do i = 1,jaux
            indx(12+i) = nsit*6+i
         enddo
         dphi1 = slat(isit1)-slat(2)
         dphi2 = slat(isit2)-slat(2)
         dlam1 = slon(isit1)-slon(2)
         dlam2 = slon(isit2)-slon(2)
         dtime = dt*365.2422d0*24.0d0*6.0d1
         dfac = dtime*rtod*3.6d3
      endif
c
c     independent observation: add an auxiliary parameter
c     differnce angle: no additional parameter, correlated obs.
      if (it.eq.2) then
         iterm = 13+jaux
         if (ibd.eq.1.or.isit1.ne.istsav.or.dabs(dt-tsave)
     .       .gt.1.0d-3.or.dabs(obs).lt.1.0d-5) iaux = iaux+1
         indx(iterm) = nsit*6+jaux+iq_sit*3+iaux
c        unit in second for the auxiliary parameter
         coef(iterm) = -dtor/3.6d3
      endif
c
c     astrometric azimuth
      if (it.eq.1) then
         ba = basazm(isit1,isit2,dt,1)
         omc = obs-ba
cmk      Subtract 360 degrees if omc is in LH quadrants
         if (omc.gt.pi) omc = omc-pi*2.0d0
         if (omc.lt.-pi) omc = omc+pi*2.0d0
         omcs = dabs(omc)*rtod*3.6d3
         if (omcs.gt.cria) goto 120
c         print*,'azimuth o-c:',isit1,isit2,obs,ba,omc
c        get coefficients
         call azimuc(isit1,isit2,dt,coef,1)
c        get coefficients for CD parameters
         if (pmode.eq.3)
     .   call addaux(iterm,coef,dphi1,dphi2,dlam1,dlam2,dfac,1)
c        earthquake correction coefficients
         if (iomode(3).gt.0.and.iq_sit.gt.0) 
     .      call add_quake_c(isit1,isit2,iterm,dt,coef,indx)
c        fill the normal matrix
         call norms(iterm,indx,coef,wght,omc,1)
c        cumulative chi squares
         chi2 = chi2+omc*wght*omc
c        unit: second
         omc = omc*rtod*3600.0d0
         goto 100
      endif
c
c     horizontal angle
c     oriented direction approach as suggested by Kurt

      if (it.eq.2) then

         if (ibd.eq.1.or.isit1.ne.istsav.or.dabs(dt-tsave)
     .       .gt.1.0d-3.or.dabs(obs).lt.1.0d-5) then
            azisav = basazm(isit1,isit2,dt,1)
            obssav = azisav-obs
            istsav = isit1
            tsave = dt
            esave = er
            call azimuc(isit1,isit2,dt,coef,1)

c           get coefficients for CD parameters
            if (pmode.eq.3)
     .      call addaux(iterm,coef,dphi1,dphi2,dlam1,dlam2,dfac,1)

            omc = 0.0
            call norms(iterm,indx,coef,wght,omc,1)
            goto 100
         endif

         ba = basazm(isit1,isit2,dt,1)
         temp = ba-obssav

         if (temp.lt.0.0d0) temp = pi*2.0d0+temp

         omc = obs - temp

cmk      Subtract 360 degrees if omc is in LH quadrants
         if (omc.gt.pi) omc = omc-2.0d0*pi
         if (omc.lt.-pi) omc = omc+20.d0*pi

         omcs = dabs(omc)*rtod*3.6d3

         if (omcs.gt.cria) goto 120

c        print*,'horiz. angle o-c:',isit1,isit2,obs,temp,omc
         call azimuc(isit1,isit2,dt,coef,1)

c        get coefficients for CD parameters
         if (pmode.eq.3)
     .     call addaux(iterm,coef,dphi1,dphi2,dlam1,dlam2,dfac,1)

c        earthquake correction coefficients
         if (iomode(3).gt.0.and.iq_sit.gt.0) 
     .      call add_quake_c(isit1,isit2,iterm,dt,coef,indx)

c        fill the normal matrix
         wght = 1.0d0/(er**2)
         call norms(iterm,indx,coef,wght,omc,1)
         chi2 = chi2+omc*wght*omc
c        unit: second
         omc = omc*rtod*3600.0d0
         goto 100
      endif
c
c     horizontal direction (using differenced data and diagonal weight)
      if (it.eq.3) then

         if (ibd.eq.1.or.isit1.ne.istsav.or.dabs(dt-tsave)
     .       .gt.1.0d-3.or.dabs(obs).lt.1.0d-5) then
            azisav = basazm(isit1,isit2,dt,1)
            obssav = obs
            istsav = isit1
            tsave = dt
            esave = er

            do i = 1,12
               indsav(i) = indx(i)
            enddo

            call azimuc(isit1,isit2,dt,coesav,1)
c           get coefficients for CD parameters

            if (pmode.eq.3)
     .      call addaux(iterm,coef,dphi1,dphi2,dlam1,dlam2,dfac,1)

            goto 100
         endif

         ba = basazm(isit1,isit2,dt,1)
         temp = ba-azisav
         if (temp.lt.0.0d0) temp = pi*2.0d0+temp
         omc = obs-obssav
c        if (omc.lt.0.0d0) omc = pi*2.0d0+omc
         omc = omc-temp

cmk      Subtract 360 degrees if omc is in LH quadrants 
         if (omc.gt.pi) omc = omc-pi*2.0d0
         if (omc.lt.-pi) omc = omc+pi*2.0d0
         omcs = dabs(omc)*rtod*3.6d3

c        print*,'horiz. direction o-c:',isit1,isit2,obs,temp,omc,
c    .   omcs

         if (omcs.gt.cria) goto 120
         call azimuc(isit1,isit2,dt,coef,1)

c        get coefficients for CD parameters
         if (pmode.eq.3)
     .     call addaux(iterm,coef,dphi1,dphi2,dlam1,dlam2,dfac,1)

c        fill the normal matrix
         do i = 1,6
            coef(i) = coef(i)-coesav(i)
            coef(i+12) = -coesav(i+6)
            indx(i+12) = indsav(i+6)
         enddo

         iterm = 18
         wght = 1.0d0/(er**2+esave**2)
         call norms(iterm,indx,coef,wght,omc,1)
         chi2 = chi2+omc*wght*omc
         omc = omc*rtod*3600.0d0
         goto 100
      endif
c
c     baseline length
      if (it.eq.4) then
         bl = baslen(isit1,isit2,dt,3,2)
         omc = obs-bl
         if (dabs(omc).gt.cril) goto 120

c        get coefficients (xyz coordinate)	
         call baslc(isit1,isit2,bl,dt,coef,1)

c        get coefficients for CD parameters
         if (pmode.eq.3)
     .   call addaux(iterm,coef,dphi1,dphi2,dlam1,dlam2,dfac,1)

c        earthquake correction coefficients
         if (iomode(3).gt.0.and.iq_sit.gt.0) 
     .      call add_quake_c(isit1,isit2,iterm,dt,coef,indx)
c
c        fill the normal matrix
         call norms(iterm,indx,coef,wght,omc,1)
         chi2 = chi2+omc*wght*omc
         omc = omc*1.0d3
         goto 100
      endif
c
c     baseline length rate
      if (it.eq.14) then
         bl = baslen(isit1,isit2,dt,3,2)
         tv = vx(isit2)-vx(isit1)
         xd = (x(isit2)-x(isit1)+tv*dt)*tv
         tv = vy(isit2)-vy(isit1)
         yd = (y(isit2)-y(isit1)+tv*dt)*tv
         tv = vz(isit2)-vz(isit1)
         zd = (z(isit2)-z(isit1)+tv*dt)*tv
         cal = (xd+yd+zd)/bl
c        unit: mm/yr
         omc = obs-cal*1.0d3
         omckm = omc*1.0d-3
         if (dabs(omc).gt.criv*1.0d3) goto 120
c        get coefficients (xyz coordinate)	
         call baslcv(isit1,isit2,bl,dt,coef,1)
c        get coefficients for CD parameters
         if (pmode.eq.3)
     .   call addaux(iterm,coef,dphi1,dphi2,dlam1,dlam2,dfac,1)
c
c        fill the normal matrix
         call norms(iterm,indx,coef,wght,omckm,1)
         chi2 = chi2+omckm*wght*omckm
         goto 100
      endif
 
c     kick out outliers
 120  if (iomode(10).gt.0) then
         continue
      else
         
         print 123,ib,it,sname(isit1),isit1,sname(isit2),
     .   isit2,obs,omc
 123  format (' FILLN: outlier at obs. # ',i4,' type ',i3,1x,
     . a8,'(',i4,') to ',a8,'(',i4,')',' obs = ',
     . f20.4,' omc = ',1pg12.4, 'rad or m(?)')
      endif

      ill = .true.

 122  format(2x,i5,4x,i3,2i5,f10.3,3x,f12.4)
c
 100  continue
      return
      end
c
