      subroutine frmnor_seqe(idatf,nobs)
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
      logical ill
      integer idatf,n1,i,ifile,iobs,iout
      integer ie,iobs1,it,ixp,iobs0
      integer nobs,ie1,ie2,sit_lst(maxsit)
      integer len,nblen,ifrm
      character*64 subfil
      common /iframe/ifrm
c
c     initialization
      iaux = 0
      ifrm = 0
      rewind (idatf)
      read (idatf,*) n1
      do i = 1,n1
         read (idatf,'(2x)')
      enddo
c
c     print head of event record file
      if (iomode(10).gt.0.and.iq_optn.eq.2) then
         call pline(24,64,'-',2)
         if (iomode(3).gt.0) then
         write(24,'(10x,a)') 'coseismic deformation correction'
         call pline(24,64,'.',1)
         write(24,'(2x,a4,2x,a3,2(1x,a4),4x,a8,12x,a8)')
     .      'year','day','sit1','sit2','original','modified'
         else
            write (24,'(10x,a)') 'list of unused observations'
            call pline(24,64,'.',1)
            write(24,'(4x,a3,3x,a4,3x,a4,2(5x,a4),2(5x,a6))')
     .         'No.','time','type','sit1','sit2','reason','misfit'
         endif
      endif
c
c     read number of sequential combining files
      read (idatf,*) ifile
      do 50 ie = 1,ifile
         read (idatf,'(a)',err=200,end=200) subfil
         len = nblen(subfil)
         open (31,file=subfil(1:len),status='old',err=200)
         print*,'   Begin to process <<< ',subfil(1:len)
         read (31,*) n1

         do i = 1,n1
            read (31,'(2x)')
         enddo

         read (31,*) ie1

         do 20 ixp = 1,ie1

         iobs = 0
         read (31,*) ie2,it,iobs1
         iobs0 = iobs

         if (it.le.20) iobs = iobs+iobs1
         if (it.gt.20) iobs = iobs+iobs1*3
         if (jaux.gt.0.and.(it.eq.31.or.it.eq.32)) ifrm = ifrm+1
         if (it.le.30) then
            call read_obs_block(31,it,dt,iobs0,iobs1)
         else
            call read_obs_full(31,it,dt,iobs1,sit_lst)
         endif

         time1 = dt+rtime
         if (it.le.30) goto 20
         call filln_full(it,iobs1,dt,sit_lst,ill)

         if (ill) then
            continue
         endif

 20      continue     
         close (31)
c
      if (iomode(10).gt.0) write (24,'(1x,i6,a,i5)') 
     .   iout,' obs. have been removed from exp. ',ie
c
      write(*,'(2a,f8.2)') ' cumulated prefit residual weighted '
     .                    ,'squares = ',chi2
c     print*,' cumulated prefit residual weighted squares = ',chi2
 50   continue
      nobs = iobs
c
 200  continue
c
      return
      end
c
