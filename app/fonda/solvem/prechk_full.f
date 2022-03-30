      subroutine prechk_full(idatf)
c
c     check for bogus observations and "ill" sites before forming
c     normal equations
c
c     basic logics:
c       1. coordinate correction needs 1 set of coordinates input or
c          2 baseline length observations or 2 azimuth observations
c          or 3 direction observations for one epoch.
c       2. velocity needs 1 set of velocity input or
c          2 baseline length observations or 2 azimuth observations
c          or 3 direction observations for two epochs.
c       3. deflection is excluded from the counting.
c       4. velocity needs 1 year interval resolution.
c
c     solution combination mode: (smode)
c       1. Gauss-Markov model
c       2. Gauss-Helmert model
c       3. model coordinate approach
c       4. sequential updating model
c       5. Kalman filtering model
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
c        vx,vy,vz : mm/year
c        erd      : mm
c        time     : year
c        error    : mm for trilatulation and 3-D survey
c                   second for triangulation
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer idatf,iobs,i,j,ifile,isit1,ie1,it
      integer illsit,len,nblen,n1,nlen,match_name
      integer iuse(maxsit),vuse(maxsit)
      character*8 name1
      character*64 subfil
      equivalence (iuse,map)
      equivalence (vuse,map(maxsit+1))
c
c     print head of event record file
      if (iomode(10).gt.0) then
         write(24,'(15x,a)') 'check outliers of the data'
         call pline(24,64,'-',1)
         write(24,'(6x,a1,2x,3(1x,a4),6x,a3,9x,a3)')
     .      '#','type','sit1','sit2','o-c','obs'
      endif
c
      jaux = 0
c
c     read number of sequential combining files
      read (idatf,*) ifile
      do 50 i = 1,ifile
c
c        get sub-file name
         read (idatf,'(a)',err=200,end=200) subfil
         len = nblen(subfil)
         open (31,file=subfil(1:len),status='old',err=200)
         read (31,*) n1
         do 20 j = 1,n1
            read (31,'(a8)') name1
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            if (isit1 .le. 0) goto 20
            if (iuse(isit1).lt.minic) iuse(isit1) = iuse(isit1)+1
            if (iuse(isit1).ge.minic) vuse(isit1) = vuse(isit1)+1
 20      continue
         read (31,*) ie1
         read (31,*) ie1,it,iobs
ccc         if (it.eq.31.or.it.eq.32) jaux = jaux+6
         close(31)
 50   continue
c
c     identify "ill" sites
      illsit = 0
      do 80 i = 1,nsit
ccc      if (iuse(i).lt.minic) illsit = illsit+1
         if (iuse(i).lt.minic) then 
             illsit = illsit+1
             print*, 'ill site :isit(',i,')=',sname(i)
         endif
 80   continue
      if (illsit.gt.0) then
         if (iomode(10).gt.0) 
     .      write (24,'(2x,i4,a)') illsit,
     .      ' ill sites have been identified!' 
         print*,illsit,' ill sites have been identified!'
         print*,' in prechk_full'
      endif

 200  continue

      return
      end
c
