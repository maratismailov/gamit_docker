      subroutine prechk(idatf,timemin,timemax)
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
      logical ill,idr
      integer idatf,iobs,ie,iobs0,iobs1,isit1,isit2,ib
      integer it,i2sit,ie1,illsit,l_obs
c     integer record(4000,4)
      real*8  record(4000,4)
      real*8  timemax,timemin
      real*8  err,vrho,robs,omc1
      dimension  err(3),vrho(3),robs(3)
      common/angle/i2sit,idr

      timemax = 0.
      timemin = 3000.
c
c     print head of event record file
      if (iomode(10).gt.0) then
         write(24,'(15x,a)') 'check outliers of the data'
         call pline(24,64,'-',1)
         write(24,'(6x,a1,2x,3(2x,a4,3x),7x,a3,14x,a3)')
     .      '#','type','sit1','sit2','o-c','obs'
      endif
c
c     temporary assuming all data are fine in list mode
      if (idatf.eq.14) goto 80
c
c     read number of experiment
      read (idatf,*) iexp
      iobs = 0
      jobs = 0
      do 50 ie = 1,iexp
c        read experiment index, obs. type and number
         read (idatf,*) ie1,it,iobs1
         iobs0 = iobs
         if (it.le.20) iobs = iobs+iobs1
         if (it.gt.20.and.it.le.40) iobs = iobs+iobs1*3
         isit2 = 0
c     read observation data     
      do 20 ib = 1,iobs1
cmk      dt is uninitialised at this point in the first cycle!
         call obsred(idatf,it,dt,robs,err,vrho,l_obs,isit1,isit2)
         time1 = dt+rtime
         timemax = max(time1,timemax)
         timemin = min(time1,timemin)

c        check to make sure the from and to sites are 
c        not the same
         if (isit1.eq.isit2) then
             print *, 'The from-site is the same as the to-site :',
     .          sname(isit1)
         endif
c
c        3-D data as model coordinate
         if (smode.eq.3.and.it.gt.20) goto 20
         if (it.gt.30) goto 20
c
c        check o-c, filter out all unusual values
         call getomc(it,l_obs,isit1,isit2,ib,robs,omc,dt,ill)
c        
c        compress direction and az observations
c        if observation is a direction or az and a new occupation
c         (idr=true) 
cmk      why do we do this???
         if (it.eq.2.or.it.eq.3.or.it.eq.12.or.it.eq.13) then
cmk         if (idr) isit1 = i2sit
         endif
c
c        record all effective observation.
c        1. set up time resolution threshold as 0.5 year (0.5 days?)
c        2. sorting record table by site index and time

         if (.not.ill) call recod(isit1,isit2,it,dt,omc,record)      
 
 20   continue     
c
 50   continue
c
c     identify "ill" sites
 80   call remedy(idatf,record,illsit)
      if (illsit.gt.0) then
         if (iomode(10).gt.0) 
     .      write (24,'(2x,i4,a)') illsit,
     .      ' ill sites have been identified!' 
         print*,illsit,' ill sites have been identified!'
         print*,' in prechk'
      endif

      return
      end
c
