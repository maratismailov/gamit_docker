      subroutine redinf(ifile,noisl)
c
c     read input file, arrange all informations.
c     basic logics:
c     1. due to the existance of velocity, all observations
c        are time-dependent
c     2. different type of observations could be mixed
c
c     coordinate system : (comode)
c         1. geocentric spherical coordinate (lat,lon,radius)
c         2. geodetic (ellipsoidal) coordinate (lat,lon,height)
c         3. geocentric Cartesian coordinate (x,y,z)
c         4. topocentric coordinate (x(east),y(north),z(height))
c         Usually, we name the geocentric coordinate as 'global' frame.
c         Meanwhile, the topocentric coordinate is named as 'local' frame.
c
c     coordinate mode : (cmode)
c         1. using network file
c         2. add earthquake jump
c         3. ......
c
c     velocity mode : (vmode)
c         1. directly read geocentric velocity values
c         2. directly read geodetic velocity values
c         3. using model parameters (trainsition, rotation, gradient, ...)
c         4. using dislocation model
c         5. ......
c
c     itype: (it)
c         1. astrometric azimuth
c         2. horizontal angle
c         3. horizontal direction
c         4. baseline length
c         5. zenith height
c         6. leveling
c
c        11. 3-D geocentric coordinate
c        12. 3-D geocentric velocity
c        13. 3-D geodetic coordinate
c        14. 3-D geodetic velocity
c
c     unit:
c         x, y, z : km
c        vx,vy,vz : mm/year
c        vxdx,vxdy,vydx,vydy : mm/km/year
c        aerr     : second    (for triangulation)
c        aerr     : mm        berr  : mm/km   (for trilatulation)
c        itime(*,1) : year
c        itime(*,2) : day of year
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
c
      character*80 buf80
      integer i1,i2,ie,is,i
      integer ifile,noisl,ierr,isit,it,igp,ig,j,is1
      integer lift_arg
c
c     read velocity and noise mode
      read (ifile,*) cmode,vmode,nmode
      print*,' cmode, vmode, nmode =',cmode, vmode, nmode
c    
c     earthquake induced coordinate jumps
      if (cmode.eq.2) then
         read (ifile,'(a)') buf80
         i = lift_arg(buf80,quakfil,1)
         call cmodel(ierr)
      endif
c
c     velocity model
c
      if (vmode.le.2) then
c      rewind (ifile)
      do 20 isit = 1,nnet
      if (vmode.eq.1) read (ifile,22) vx(isit),vy(isit),vz(isit)
      if (vmode.eq.2) then
         read (ifile,22) ve(isit),vn(isit),vu(isit)
         s1 = dsin(sla(isit))
         c1 = dcos(sla(isit))
         s2 = dsin(slo(isit))
         c2 = dcos(slo(isit))
         vx(isit) = -ve(isit)*s2-vn(isit)*s1*c2+vu(isit)*c1*c2
         vy(isit) = ve(isit)*c2-vn(isit)*s1*s2+vu(isit)*c1*s2
         vz(isit) = vn(isit)*c1+vu(isit)*s1
      endif
 20   continue
      endif
 22   format (64x,3f8.3)
c
c     we always use geodetic coordinate to construct velocity model
c    
      if (vmode.eq.3) then
         read (ifile,'(a)') buf80
         i = lift_arg(buf80,modfil,1)
      endif
      if (vmode.gt.2) call vmodel(ierr)
c
c     read number of experiment
      read (ifile,*) iexp
      i1 = 0
      i2 = 0
      do 50 it = 1,iexp
c
c        read time, observation type and group number
         read (ifile,*) itime(it,1),itime(it,2),itype(it),igp
         igroup(it) = igp
         do 30 ig = 1,igp
            i2 = i2+1
c
c           read site number and error model
            read (ifile,*) isit,aerr(i2),berr(i2)
            read (ifile,*) (lblsit(is),is = i1+1,i1+isit)
            iesit(i2) = isit
            i1 = i1+isit
 30      continue
 50   continue
c
c     get the effective site number
      ie = 0
      do 70 i = 1,i1
         is = lblsit(i)	
         if (ie.eq.0) go to 90
         do 80 j = 1,ie
            is1 = itoj(j)
            if (is.eq.is1) goto 70
 80      continue
 90      ie = ie+1
         itoj(ie) = is
         jtoi(is) = ie
 70   continue
      nsit = ie
c
c     read noise level
      read (ifile,*) noisl
c
c     read output mode
      read (ifile,*) outmod
c
      return
      end
