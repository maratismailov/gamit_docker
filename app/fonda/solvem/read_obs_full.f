      subroutine read_obs_full(idatf,it,dt,iobs,sit_lst)
c
c     read full observation data from data file
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
c
c        31. 3-D spherical coordinate
c        32. 3-D spherical frame velocity
c        33. 3-D Cartesian coordinate
c        34. 3-D Cartesian frame velocity
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : mm/year
c        erd      : mm
c        time     : year
c        error    : mm for trilatulation and 3-D survey
c                   second for triangulation
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer idatf,it,isit1,k,i1,i,nlen,nblen
      integer match_name,iobs,sit_lst(iobs),lift_arg
      real*8  obs
      dimension obs(maxprm)
      character*8 name1
      character*128 line
      integer ista, iend
      equivalence (raw_obs,obs)
c
c     3-D spherical coordinates
      if (it.eq.31) then
         i1 = 0
         do 30 i = 1,iobs
            read (idatf,'(a128)') line
            read (line,*) time1,stla,stlo,stra
            k = lift_arg(line,name1,5)
            if (k.le.0) goto 30
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            if (isit1 .le. 0) goto 100
            obs(i1+1) = stla
            obs(i1+2) = stlo
            obs(i1+3) = stra
            i1 = i1+3
            sit_lst(i) = isit1
            dt = time1-rtime
 30      continue
         do 40 i = 1,iobs*3
            ista = i*(i-1)/2
            iend = ista+i
            read (idatf,110) i1,(cova(k),k=ista+1,iend)
c           print*, i1,(cova(k),k=ista+1,iend)
 40      continue
      endif
c
c     3-D Cartesian coordinates
      if (it.eq.33) then
         i1 = 0
         do 50 i = 1,iobs
            read (idatf,'(a)') line
            read (line,*) time1,stx,sty,stz
            k = lift_arg(line,name1,5)
            if (k.le.0) goto 50
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            if (isit1 .le. 0) goto 100
            obs(i1+1) = stx
            obs(i1+2) = sty
            obs(i1+3) = stz
            i1 = i1+3
            sit_lst(i) = isit1
            dt = time1-rtime
 50      continue
         do 60 i = 1,iobs*3
            ista = i*(i-1)/2
            iend = ista+i
            read (idatf,110) i1,(cova(k),k=ista+1,iend)
 60      continue
      endif
c
c     3-D Cartesian velocity
      if (it.eq.34) then
         i1 = 0
         do 70 i = 1,iobs
            read (idatf,'(a)') line
            read (line,*) time1,stx,sty,stz
            k = lift_arg(line,name1,5)
            if (k.le.0) goto 70
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            if (isit1 .le. 0) goto 100
            obs(i1+1) = stx
            obs(i1+2) = sty
            obs(i1+3) = stz
            i1 = i1+3
            sit_lst(i) = isit1
            dt = time1-rtime
 70      continue
         do 80 i = 1,iobs*3
            ista = i*(i-1)/2
            iend = ista+i
c           read (idatf,130) i1,(cova(k),k=ista+1,iend)
            read (idatf,110) i1,(cova(k),k=ista+1,iend)
 80      continue
      endif
c
 110  format (1x,i3,'. ',50(5(1x,d23.16),:,/,6x))
 130  format(i5,'. ', 10(1x,d11.4),:/,200(7x,10(1x,d11.4),:/) )
      goto 200

 100  print*,' mismatch site name at READ_OBS_FULL: ',name1
 200  continue
      return
      end

