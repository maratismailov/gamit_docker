      subroutine read_obs_block(idatf,it,dt,iobs0,iobs1)
c
c     read observation data and construct normal matrix
c     All partial derivatives are related to geocentric coordinates.
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
      logical ill
      integer idatf,i,iobs1,iv,iout
      integer ie,iobs0,ib,it,ibd,isit1,isit2
      integer l_obs
      real*8  err,vrho,robs,omc
      dimension  err(3),vrho(3),robs(3),omc(3)
c
      iv = 0
      iout = 0
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
         if (iomode(3).gt.0.and.it.lt.11.and.iq_optn.eq.2) 
     .      call poscor(it,15,time1,isit1,isit2,dt,robs(1))
c
c        get coefficients and fill normal matrix
         if (it.le.20) call filln(it,isit1,isit2,dt,ib,ibd,obs,er,
     .                      rho,omc(1),ill)
         if (it.gt.20) call filln_block(it,isit1,isit2,l_obs,dt,ib,ibd,
     .                      robs,err,vrho,omc,ill)
c
         if (it.eq.14.and..not.ill) then
           write(19,'(4(1x,f10.4))') slon(isit1)*rtod,slat(isit1)*rtod,
     .         slon(isit2)*rtod,slat(isit2)*rtod
         endif
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
 60   format (1x,i5,f10.3,i3,4x,2(1x,a8),2x,a8,2x,f10.3)
c
      return
      end
c
