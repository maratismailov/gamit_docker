      subroutine frmnor(idatf,nobs)
c
c     read observation data and construct normal matrix
c     All partial derivatives are related to geocentric coordinates.
c     Transformation among different coordinates is neccesary.
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
      logical ill,new
      
      integer idatf,n1,i,igo,iobs,iv,iout
      integer ie,k,iobs0,iobs1,ib,it,ibd,isit1,isit2
      integer igp,nobs,ikt,ie1,l_obs
      real*8  err,vrho,robs,omc
      dimension  err(3),vrho(3),robs(3),omc(3)
c
c     initialization
      iaux = 0
      rewind (idatf)
      read (idatf,*) n1
      do i = 1,n1
         read (idatf,'(2x)')
      enddo
c
c     print head of event record file
      if (iomode(10).gt.0) then
         call pline(24,64,'-',2)
         if (iomode(3).gt.0.and.iq_optn.eq.2) then
         write(24,'(10x,a)') 'coseismic deformation correction'
         call pline(24,64,'.',1)
         write(24,'(2x,a4,2x,a3,2(1x,a4),4x,a8,12x,a8)')
     .      'year','day','sit1','sit2','original','modified'
         else
            write (24,'(10x,a)') 'list of unused observations'
            call pline(24,64,'.',1)
            write(24,'(4x,a3,4x,a4,2x,a4,2(5x,a4),2(5x,a6))')
     .         'No.','time','type','sit1','sit2','reason','misfit'
         endif
      endif
c
c     read number of experiments
      read (idatf,*) iexp
      igo = 0
      iobs = 0
      jobs = 0
      iv = 0
      ntype = 0
      iout = 0
      do 50 ie = 1,iexp
c
c        read experiment index, obs. type and number
c        ntype is counter of the number of data types
         read (idatf,*) ie1,it,iobs1
         if (ntype.eq.0) then
            ntype = 1
            listyp(1) = it
         endif
         new = .true.
         do 15 k = 1,ntype
            ikt = listyp(k)
            if (it.eq.ikt) new = .false.
 15      continue
         if (new) then
            ntype = ntype+1
            listyp(ntype) = it
         endif
         iobs0 = iobs
         if (it.le.20) iobs = iobs+iobs1
         if (it.gt.20.and.it.le.40) iobs = iobs+iobs1*3
c
c     read observation data
      do 20 ib = iobs0+1,iobs0+iobs1
         ibd = ib-iobs0
         call obsred(idatf,it,dt,robs,err,vrho,l_obs,isit1,isit2)
         time1 = dt+rtime
         obs = robs(1)
         er  = err(1)
         rho = vrho(1)
         if (it.lt.21) then
            if (jtoi(isit1).eq.0.or.jtoi(isit2).eq.0) then
               iout = iout+1
               if (iomode(10).gt.0) 
     .         write (24,60) iout,time1,it,
     .         sname(isit1),sname(isit2),'ill site'
               goto 20
            endif
         endif
c
c        position correction due to earthquake
         if (iomode(3).gt.0.and.it.lt.11.and.iq_optn.eq.2) 
     .      call poscor(it,15,time1,isit1,isit2,dt,robs(1))
c
c        get coefficients and fill normal matrix
         if (it.le.20) call filln(it,isit1,isit2,dt,ib,ibd,obs,er,
     .                      rho,omc(1),ill)
         if (it.gt.20) call filln_block(it,isit1,isit2,l_obs,dt,ib,ibd,
     .                      robs,err,vrho,omc,ill)
c
         if (ill) then
            do i = 1,l_obs
               iout = iout+1
               if (iomode(10).gt.0) 
     .         write (24,60) iout,time1,
     .         it,sname(isit1),sname(isit2),' big o-c',omc(i)
            enddo
         endif
c 
 20   continue     
c
      if (iomode(10).gt.0) write (24,'(1x,i6,a,i5)') 
     .   iout,' obs. have been removed from exp. ',ie
c
 10   continue
      igo = igo+igp
 50   continue
      nobs = iobs
      write(*,'(2a,f8.2)') ' cumulated prefit residual',
     .                       ' weighted squares = ',chi2
c     print*,' cumulated prefit residual weighted squares = ',chi2
 60   format (1x,i5,f10.3,i3,4x,2(1x,a8),2x,a8,2x,f10.3)
c
      return
      end
c
