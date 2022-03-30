      subroutine obsred(idatf,it,dt,robs,err,vrho,
     .            l_obs,isit1,isit2)
c
c
c     read one observation data from data file
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
c     unit:
c         x, y, z : m
c        vx,vy,vz : mm/year
c        err      : mm for it <= 24
c        err      :  m for it > 24 
c        time     : year
c        error    : mm for trilatulation and 3-D survey
c                   second for triangulation
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer idatf,it,isit1,isit2,ide,imi,nlen,nblen
      integer match_name,l_obs,lift_arg
      character*8 name1, name2
      character*128  buf128
c     real*8  err,robs,vrho
      dimension  err(3),vrho(3),robs(3)
c
c     3-D data format is different from conventional survey
      l_obs = 1
      if (it.le.3) then
         read (idatf,'(a)') buf128
         read (buf128,*) time1,ide,imi,sec,err(1)
         nlen = lift_arg(buf128,name1,6)
         nlen = lift_arg(buf128,name2,7)
         nlen = nblen(name1)
         isit1 = match_name(nsit,nlen,sname,name1)
         if (isit1 .le. 0) goto 100
         nlen = nblen(name2)
         isit2 = match_name(nsit,nlen,sname,name2)
         if (isit2 .le. 0) goto 100
         robs(1) = dble(ide)+dble(imi)/6.0d1+sec/3.6d3
         robs(1) = robs(1)*dtor
c        convert err(1) to radians
         err(1) = err(1)/3.6d3/rtod
      else
         if (it.gt.3.and.it.lt.21) then
            read (idatf,'(a)') buf128
            read (buf128,*) time1,robs(1),err(1)
            nlen = lift_arg(buf128,name1,4)
            nlen = lift_arg(buf128,name2,5)
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            if (isit1 .le. 0) goto 100
            nlen = nblen(name2)
            isit2 = match_name(nsit,nlen,sname,name2)
            if (isit2 .le. 0) goto 100
c           error units unchanged
         endif
         if (it.gt.20.and.it.le.24) then
            l_obs = 3
            read (idatf,'(a)') buf128
c           print*, 'buf128 ',buf128
            read (buf128,*)
     .         time1,robs(1),err(1),robs(2),err(2),
     .         robs(3),err(3),vrho(1),vrho(2),vrho(3)
            nlen = lift_arg(buf128,name1,11)
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            isit2 = 0
            if (isit1 .le. 0) goto 100
         endif
         if (it.ge.25.and.it.le.30) then
c
c           errors are in meters here
c
            l_obs = 3
            read (idatf,'(a)') buf128
            read (buf128,*)
     .         time1,robs(1),err(1),robs(2),err(2),
     .         robs(3),err(3),vrho(1),vrho(2),vrho(3)
cc            print *,
cc     .      time1,robs(1),err(1),robs(2),err(2),
cc     .      robs(3),err(3),vrho(1),vrho(2),vrho(3)
c
            nlen = lift_arg(buf128,name1,11)
            nlen = lift_arg(buf128,name2,12)
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            if (isit1 .le. 0) goto 100
            nlen = nblen(name2)
            isit2 = match_name(nsit,nlen,sname,name2)
            if (isit2 .le. 0) goto 100
c           empirical baseline length related factor
cmk         dt is undefined here on the first iteration!
            bl = baslen(isit1,isit2,dt,3,2)
            bfac = dsqrt(1.0d0+bl*bl*1.0d-10)
            err(1) = err(1)*bfac
            err(2) = err(2)*bfac
            err(3) = err(3)*bfac
         endif
      endif
      dt = time1-rtime
c
      if (it.ge.4.and.it.le.24) then
c        convert errors to metres
         err(1) = err(1)*1.0d-3
         err(2) = err(2)*1.0d-3
         err(3) = err(3)*1.0d-3
      endif
c
      goto 200

 100  print*,' mismatch site name at OBSRED: ',
     .   name1,' ',name2,nlen,isit1,isit2
 200  continue
      return
      end

