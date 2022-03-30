c--------------------------------------------------------------
      subroutine getomc(it,l_obs,isit1,isit2,ibd,obs,omc,time,ill)
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

      logical ill,idr
      integer it,isit1,isit2,ibd,istsav,indsav,i2sit,i2,l_obs
      dimension coesav(12),indsav(12),obs(l_obs),tmp(3),tmp1(3)
      common/saved/azisav,coesav,obssav,esave,tsave,istsav,indsav
      common/angle/i2sit,idr
c
      ill = .false.
      idr = .false.
      small = 1.0d-15
c     from m to mm
      facm = 1.0d3
c
c     astrometric azimuth
      if (it.eq.1) then
         ba = basazm(isit1,isit2,time,1)
         omc = obs(1)-ba

         if (omc.gt.pi) omc = omc-pi*2.0d0
         if (omc.lt.-pi) omc = omc+pi*2.0d0

         omc = omc*rtod*3600.0d0

         if (dabs(omc).gt.cria) goto 120

         goto 100
      endif
c
c     horizontal angle or direction
      if (it.eq.2.or.it.eq.3) then

c        if this is the first observation in the file, 
c        the first of a pair or,
c        0.365 days (= 8.74 hrs) have passed since the last obs
cmk      Note this time limit! 
         if ((ibd.eq.1).or.(isit1.ne.istsav).or.(dabs(time-tsave)
     .       .gt.1.0d-3)) then
            azisav = basazm(isit1,isit2,time,1)
            obssav = obs(1)
            istsav = isit1
            tsave = time
            i2sit = isit2
            goto 100
         endif

         ba = basazm(isit1,isit2,time,1)
         temp = ba-azisav

         if (temp.lt.0.0d0) temp = pi*2.0d0+temp

         omc = obs(1)-obssav-temp

         if (omc.gt.pi) omc = omc-pi*2.0d0
         if (omc.lt.-pi) omc = omc+pi*2.0d0

c        unit: second
         omc = omc*rtod*3600.0d0

c        check to see if the omc is greater than the 
c        user specified criteria
         if (dabs(omc).gt.cria) goto 120

c        this is effectively a new occupation
         idr = .true.
         goto 100
      endif
c
c     baseline length
      if (it.eq.4) then
         bl = baslen(isit1,isit2,time,3,2)
         omc = obs(1)-bl

         if (dabs(omc).gt.cril) goto 120

c        unit: mm
         omc = omc*facm
         goto 100
      endif
c
c     baseline length rate
      if (it.eq.14) then
         bl = baslen(isit1,isit2,0.0d0,3,2)
         xd = (x(isit2)-x(isit1))*(vx(isit2)-vx(isit1))
         yd = (y(isit2)-y(isit1))*(vy(isit2)-vy(isit1))
         zd = (z(isit2)-z(isit1))*(vz(isit2)-vz(isit1))
         cal = (xd+yd+zd)/bl
c        unit: mm/yr
         omc = obs(1)-cal*1.0d3
         if (dabs(omc).gt.criv*1.0d3) goto 120
         goto 100
      endif
c
c     3-D geocentric coordinate
      if (it.eq.21) then
c        determine x, y or z
         i2 = ibd-ibd/3*3
         if (i2.eq.0) i2 = 3
         if (i2.eq.1) omc = obs(1)-x(isit1)
         if (i2.eq.2) omc = obs(1)-y(isit1)
         if (i2.eq.3) omc = obs(1)-z(isit1)
         if (dabs(omc).gt.cric) goto 120
c        unit: mm
         omc = omc*facm
         goto 100
      endif
c
c     3-D geocentric velocity
      if (it.eq.22) then
c        determine vx or vy
         i2 = ibd-ibd/3*3
         if (i2.eq.0) i2 = 3
         if (i2.eq.1) omc = (obs(1)-vx(isit1)*1.0d3)
         if (i2.eq.2) omc = (obs(1)-vy(isit1)*1.0d3)
         if (i2.eq.3) omc = (obs(1)-vz(isit1)*1.0d3)
         if (dabs(omc).gt.criv*3.0d3) goto 120
         goto 100
      endif
c
c     3-D geodetic velocity
      if (it.eq.24) then
         omc = dabs(obs(1)-ve(isit1)*1.0d3)
         omc = omc+dabs(obs(2)-vn(isit1)*1.0d3)
         omc = omc+dabs(obs(3)-vu(isit1)*1.0d3)
         if (dabs(omc).gt.criv*3.0d3) goto 120
         goto 100
      endif
c
c     3-D Cartesian baseline vector
      if (it.eq.25) then
         omc = 0.0d0
         fac = 1.0d0*time
         omc = omc+obs(1)-x(isit2)+x(isit1)
     .         -(vx(isit2)-vx(isit1))*fac
         omc = omc+obs(2)-y(isit2)+y(isit1)
     .         -(vy(isit2)-vy(isit1))*fac
         omc = omc+obs(3)-z(isit2)+z(isit1)
     .         -(vz(isit2)-vz(isit1))*fac
         if (dabs(omc).gt.cric) goto 120
c        unit: mm
         omc = omc*facm
         goto 100
      endif
c
c     3-D geodetic baseline vector
      if (it.eq.27) then
         fac = 1.0d0*time
         omc = 0.0d0
         x1 = x(isit1)
         y1 = y(isit1)
         z1 = z(isit1)
         call sphxyz(a2,a1,a3,x1,y1,z1,2)
         call getjac(a1,a2,a3,coesav,5)
         tmp(1) = x(isit2)-x(isit1)+(vx(isit2)-vx(isit1))*fac
         tmp(2) = y(isit2)-y(isit1)+(vy(isit2)-vy(isit1))*fac
         tmp(3) = z(isit2)-z(isit1)+(vz(isit2)-vz(isit1))*fac
         call axb(3,3,1,coesav,tmp,tmp1,1,0)
         omc = omc+obs(1)-tmp1(1)
         omc = omc+obs(2)-tmp1(2)
         omc = omc+obs(3)-tmp1(3)
         omc = omc*facm
         goto 100
      endif
c
c     3-D geodetic baseline rate vector
      if (it.eq.28) then
         omc = dabs(obs(1)-(vn(isit2)-vn(isit1))*1.0d3)
         omc = omc+dabs(obs(2)-(ve(isit2)-ve(isit1))*1.0d3)
         omc = omc+dabs(obs(3)-(vu(isit2)-vu(isit1))*1.0d3)   
         if (dabs(omc).gt.criv*3.0d3) goto 120
         goto 100
      endif
c
c     kick out outliers
c120  print*,' outlier: #',ibd,' type=',it,
c    .   ' sit =',isit1,isit2,' o-c =',omc,obs
 120  if (iomode(10).gt.0) write(24,130) ibd,it,sname(isit1),
     .   sname(isit2),omc,obs
 130  format(2x,i5,i5,5x,2(a8,1x),f16.6,f16.6)
      ill = .true.
c
 100  continue
      return
      end
c
