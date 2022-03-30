      subroutine poscor(it,ifil,time1,isit1,isit2,dt,robs)
c
c     site position correction due to earthquake
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
c        41. deflection
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer it,ifil,isit1,isit2,i,ncor,myear,mday,msit,julday
c
c     read number of corrections
      rewind (ifil)
      read (ifil,*) ncor
c
      do 20 i = 1,ncor
         read (ifil,*) myear,mday,msit,dx,dy
         if (msit.ne.isit1.and.msit.ne.isit2) goto 20
         time2 = 1900.0d0+julday(1,mday,myear,2)/365.2422d0
         if (time1.lt.time2) goto 20   
         if (msit.eq.isit1) then
            azm = basazm(isit1,isit2,dt,1)
            aobs = robs+(dx*dsin(azm)+dy*dcos(azm))*1.d-2
            if (iomode(10).gt.0) write(24,*)
     .        time1,isit1,isit2,robs,aobs
            robs = aobs
         endif
         if (msit.eq.isit2) then
            azm = basazm(isit2,isit1,dt,1)
            aobs = robs+(dx*dsin(azm)+dy*dcos(azm))*1.d-2
            if (iomode(10).gt.0) write(24,*)
     .        time1,isit1,isit2,robs,aobs
            robs = aobs
         endif
c 
 20   continue     
c
      return
      end
c
