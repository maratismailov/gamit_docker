      subroutine frmobs(nobs)
c
c     construct simulated observation based on itype.
c     for the time being, only consider full combination.
c
c     itype: (it)
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
c     unit:
c         x, y, z : meter
c        vx,vy,vz : mm/year
c        aerr     : mm      berr : mm/km
c        itime(*,1) : year
c        itime(*,2) : day of year
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      integer*4 julday
      integer nobs,i1,isit,iobs,ie,iyear,iday,month,igp,jobs 
      integer it,ig,isit0,iobs0
c     common/saved/azibas,isbas
c
c     experiment loop
      i1 = 0
      isit = 0
      iobs = 0
      time0 = rtime
      do 50 ie = 1,iexp
         iyear = itime(ie,1)
         iday = itime(ie,2)
         time1 = 1900.0d0+julday(month,iday,iyear,2)/365.2422d0
         igp = igroup(ie)
         jobs = 0
         it = itype(ie)
c     
c        site group loop
         do 10 ig = 1,igp
            isit0 = isit
            iobs0 = iobs
            isit = isit+iesit(i1+ig)
            ea = aerr(i1+ig)
            eb = berr(i1+ig)
c
c           get observation combination 
            call combin(it,isit,isit0,time1,time0,ea,eb,iobs)
            igobs(i1+ig) = iobs-iobs0
            jobs = jobs+igobs(i1+ig)
 10      continue
         i1 = i1+igp
         inobs(ie) = jobs
 50   continue
      nobs = iobs
c
      return
      end
c--------------------------------------------------------------
      subroutine combin(it,isit,isit0,time1,time0,ea,eb,iobs)
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      integer it,isit,isit0,iobs,jsit,isit1,isit2,is1,iobs0  
      integer i1,is2,i2

      common/sava/azisav
c
c     from mm to meter
      factr = 1.0d-3
c
c     time interval
      dt = time1-time0
c
c     site number in this group
      jsit = isit-isit0
c
 10   format (1x,'No observable in this group. site =',i2)
c
c     conventional survey
      if (it.le.6) then
         if (jsit.le.1) then
            write (6,10) jsit
            return
         endif
         isit1 = isit0+1
         isit2 = isit-1
         if (it.eq.2.or.it.eq.3) isit2 = isit1
         do 30 is1 = isit1,isit2
            iobs0 = iobs
            i1 = lblsit(is1)
            do 20 is2 = is1+1,isit
               iobs = iobs + 1
               i2 = lblsit(is2) 
               idobs(iobs,1) = i1
               idobs(iobs,2) = i2
c              astrometric azimuth (unit: degree)
               if (it.eq.1.or.it.eq.3) then
                  temp = basazm(i1,i2,dt,3)
                  data(iobs) = temp*rtod
                  erd(iobs) = ea 
               endif
c              horizontal angle
               if (it.eq.2) then
                  if (iobs-iobs0.eq.1) then
                  azisav = basazm(i1,i2,dt,3)
                  data(iobs) = 0.0d0
                  erd(iobs) = ea
                  else
                  temp = basazm(i1,i2,dt,3)-azisav
                  if (temp.lt.0.0d0) temp = pi*2.0d0+temp
                  data(iobs) = temp*rtod
                  erd(iobs) = ea
                  endif
               endif
c              baseline length (meter)
               if (it.eq.4) then
                  bsl = baslen(i1,i2,dt,3,2)
                  data(iobs) = bsl
                  erd(iobs) = ea+eb*bsl
               endif
 20         continue
 30      continue
      endif
c
c     3-D geocentric coordinate
      if (it.eq.21) then
         do 40 is1 = isit0+1,isit
            iobs = iobs + 1
            i1 = lblsit(is1)
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = x(i1)+vx(i1)*dt*factr
            erd(iobs) = ea
            iobs = iobs + 1
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = y(i1)+vy(i1)*dt*factr
            erd(iobs) = ea
            iobs = iobs + 1
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = z(i1)+vz(i1)*dt*factr
            erd(iobs) = ea
 40      continue
      endif
c
c     3-D geocentric velocity
      if (it.eq.22) then
         do 60 is1 = isit0+1,isit
            iobs = iobs + 1
            i1 = lblsit(is1)
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = vx(i1)
            erd(iobs) = ea
            iobs = iobs + 1
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = vy(i1)
            erd(iobs) = ea
            iobs = iobs + 1
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = vz(i1)
            erd(iobs) = ea
 60      continue
      endif
c
c     3-D geodetic velocity
      if (it.eq.24) then
         do 80 is1 = isit0+1,isit
            iobs = iobs + 1
            i1 = lblsit(is1)
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = ve(i1)
            erd(iobs) = ea
            iobs = iobs + 1
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = vn(i1)
            erd(iobs) = ea
            iobs = iobs + 1
            idobs(iobs,1) = i1
            idobs(iobs,2) = 0
            data(iobs) = vu(i1)
            erd(iobs) = ea
 80      continue
      endif
c
c     baseline vector in spherical frame
      if (it.eq.27) then
         if (jsit.le.1) then
            write (6,10) jsit
            return
         endif
         isit1 = isit0+1
         isit2 = isit-1
         do 130 is1 = isit1,isit2
            iobs0 = iobs
            i1 = lblsit(is1)
            do 140 is2 = is1+1,isit
               iobs = iobs + 1
               i2 = lblsit(is2) 
               idobs(iobs,1) = i1
               idobs(iobs,2) = i2
c              baseline vector
 140        continue
 130     continue
      endif
c
 100  continue
      return
      end
c     

