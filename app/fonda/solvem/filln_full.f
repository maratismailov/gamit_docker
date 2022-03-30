c
      subroutine filln_full(it,iobs,dt,sit_lst,ill)
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
c        27. 3-D spherical baseline vector
c        28. 3-D spherical baseline rate vector
c
c        31. 3-D spherical coordinate
c        32. 3-D spherical frame velocity
c        33. 3-D Cartesian coordinate
c        34. 3-D Cartesian frame velocity
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer it,iobs,ib,ibd,i,i1,i2,isit,i3,isav,iele
      integer sit_lst(iobs),ier,indx_tmp(maxprm)
      real*8  obs,cjaco
      dimension  obs(maxprm), cjaco(9),fac(10)
      logical ill
      integer irow,icol,iek(10),imatch,j,k,k1,j1,ibase,k2
      integer ifrm,ipoi
      common /iframe/ifrm
      equivalence (raw_obs,obs)
c
      ill = .false.
      p2 = 2.0d0*pi
c     unit of velocity is m/year
      fact = 1.0d0*dt
c
c     initialize
      irow = iobs*3
      icol = nsit*6
      ibase = 0
      k2 = 6
      if (it.eq.31.and.ifrm.gt.0) k2 = 12 
      call zero1d(1,irow*k2,aprm)
c
c     3-D spherical coordinate
      if (it.eq.31) then
c
c     form the design martrix
      i1 = 0
      i3 = 0
      imatch = 0
      ipoi = nsit*6+iq_sit*3+(ifrm-1)*6
      do 10 ib = 1,iobs
         isit = sit_lst(ib)
         x1 = x(isit)+vx(isit)*fact
         y1 = y(isit)+vy(isit)*fact
         z1 = z(isit)+vz(isit)*fact
         call sphxyz(a2,a1,a3,x1,y1,z1,2)
         call getjac(a1,a2,a3,cjaco,1)
         ibd = (isit-1)*6
         if (iomode(3).gt.0.and.iq_sit.gt.0) then
            call chk_quake_list(isit,dt,iek,fac,imatch)
         endif
         isav = i3
         iele = 6+imatch*3
         if (ifrm.gt.0) iele = iele+6
c        observation order = n,e,u
         do i = 1,3
            indx_tmp(i1+i) = i3+iele
            if (i.eq.1) i2 = 3
            if (i.eq.2) i2 = 0
            if (i.eq.3) i2 = 6
            do j = 1,3
               indx_ele(i3+j) = ibd+j
               indx_ele(i3+j+3) = ibd+j+3
               aprm(i3+j) = cjaco(i2+j)
               aprm(i3+j+3) = cjaco(i2+j)*dt
               if (ifrm.gt.0) then
                  indx_ele(i3+j+6) = ipoi+j
                  indx_ele(i3+j+9) = ipoi+j+3
                  aprm(i3+j+6) = cjaco(i2+j)
               endif
            enddo
            if (ifrm.gt.0) then
               aprm(i3+10) = cjaco(i2+2)*z1-cjaco(i2+3)*y1
               aprm(i3+11) = cjaco(i2+3)*x1-cjaco(i2+1)*z1
               aprm(i3+12) = cjaco(i2+1)*y1-cjaco(i2+2)*x1
            endif
            i3 = i3+iele
         enddo
         if (imatch.gt.0) then
            do j = 1,imatch
               j1 = iek(j)
               quake_use(j1) = quake_use(j1)+1
               do i = 1,3
                  i2 = isav+(i-1)*iele+(j-1)*3+6
                  if (ifrm.gt.0) i2 = i2+6
                  if (i.eq.1) k1 = 3
                  if (i.eq.2) k1 = 0
                  if (i.eq.3) k1 = 6
                  do k = 1,3
                     indx_ele(i2+k) = icol+(j1-1)*3+k
                     aprm(i2+k) = cjaco(k1+k)*fac(j)
                  enddo
               enddo
               if (j1.gt.ibase) ibase = j1
            enddo
         endif
         if (obs(i1+2).lt.0.0d0) a1 = a1-p2
         obs(i1+1) = obs(i1+1)-a2
         obs(i1+2) = obs(i1+2)-a1
         obs(i1+3) = obs(i1+3)-a3
         i1 = i1+3
 10   continue
      endif
c
c     3-D Cartesian coordinate
      if (it.eq.33) then
c     form the transpose of design martrix
      i1 = 0
      i3 = 0
      imatch = 0
      do 30 ib = 1,iobs
         isit = sit_lst(ib)
         x1 = x(isit)+vx(isit)*fact
         y1 = y(isit)+vy(isit)*fact
         z1 = z(isit)+vz(isit)*fact
         ibd = (isit-1)*6
         if (iomode(3).gt.0.and.iq_sit.gt.0) then
            call chk_quake_list(isit,dt,iek,fac,imatch)
         endif
         isav = i3
         iele = 2+imatch
c        observation order = x,y,z
         do i = 1,3
            indx_tmp(i1+i) = i3+iele
            indx_ele(i3+1) = ibd+i
            indx_ele(i3+2) = ibd+i+3
            aprm(i3+1) = 1.0d0
            aprm(i3+2) = dt
            i3 = i3+iele
         enddo
         if (imatch.gt.0) then
            do j = 1,imatch
               j1 = iek(j)
               quake_use(j1) = quake_use(j1)+1
               do i = 1,3
                  indx_ele(isav+j+(i-1)*iele+2) = icol+(j1-1)*3+i
                  aprm(isav+j+(i-1)*iele+2) = fac(j)
               enddo
               if (j1.gt.ibase) ibase = j1
            enddo
         endif
c ------ velocity considered
         obs(i1+1) = obs(i1+1)-x1
         obs(i1+2) = obs(i1+2)-y1
         obs(i1+3) = obs(i1+3)-z1
         i1 = i1+3
 30   continue
      endif
c
c     3-D Cartesian velocities
      if (it.eq.34) then
      i1 = 0
      i3 = 0
      iele = 1
      do 50 ib = 1,iobs
         isit = sit_lst(ib)
         vx1 = vx(isit)
         vy1 = vy(isit)
         vz1 = vz(isit)
         ibd = (isit-1)*6+3
c        observation order = Vx,Vy,Vz
         do i = 1,3
            indx_tmp(i1+i) = i3+iele
            indx_ele(i3+1) = ibd+i
            aprm(i3+1) = 1.0d0
            i3 = i3+iele
         enddo
c ------ consider velocity later
         obs(i1+1) = obs(i1+1)-vx1
         obs(i1+2) = obs(i1+2)-vy1
         obs(i1+3) = obs(i1+3)-vz1
         i1 = i1+3
 50   continue
      endif
      ibase = ibase*3 
      if (ifrm.gt.0) ibase = iq_sit*3+ifrm*6 
      if (iaux.gt.0) ibase = iq_sit*3+jaux+iaux
c
c     Reorder the index array 
      icol = icol+ibase
      call reordr_indx(irow,icol,indx_tmp)
c
c     get the weighting matrix
      call nrmscl(cova,obs,scale,irow,1,0)
      call cholsk(cova,obs,1,irow,ier)
c     print*, 'ierr=',ier,'FILLN_FULL'
      call nrmscl(cova,obs,scale,irow,2,0)
c
      call latwa(1,irow,obs,cova,chi2,scale,1,2)
c
      call axb(irow,irow,1,cova,obs,scale,2,0)
      call laxb_indx(icol,irow,1,iele,aprm,scale,bnorm,
     .               indx_row,indx_ele,1)
c
c     form the normal matrix
      call latwa_indx(icol,irow,iele,aprm,cova,anorm,scale,
     .                indx_row,indx_ele,2)
c      do i = nsit*6+1,icol
c         j = i*(i-1)/2+nsit*6
c         print*,(anorm(j+i2),i2=1,i-nsit*6)
c      enddo
c
c     clean array aprm
      call zero1d(1,irow*k2,aprm)
      
      goto 100
c
c     kick out outliers (not finished)
 120  if (iomode(10).gt.0) then
         continue
      else
         continue
      endif
      ill = .true.
 122  format(2x,i5,4x,i3,2i5,f10.3,3x,f12.4)
c
 100  continue
      return
      end
c
