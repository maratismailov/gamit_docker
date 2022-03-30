      subroutine residl(ifile,idatf,nobs)
c
c     calculate residuals
c
c     modified by FERHAT Gilbert June 28, 1995
c     problem with horizontal angle equals to zero
c     version 1.00
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
c        err      : mm
c        time     : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      character*6 unit(3)
      character*25 desit(41)
      integer ifile,idatf,nobs,indsav,iazsav,i,igo,iobs,ie
      integer iobs0,iobs1,it,imsit,iefob,ib,ibd,ieff,ie1,k,isit1
      integer isit2,igp,istsav,iauxsav,l_obs,il
c     real*8  err,vrho,robs,omc,oc,cal
      dimension  err(3),vrho(3),robs(3),omc(3),oc(3),cal(3)
      dimension indsav(12),coesav(12),covo(6),temp(6)
      character*1 symb
      logical ill
      common/saved/azisav,coesav,obssav,esave,tsave,iazsav,indsav
      common/save2/istsav,iauxsav
c
      data (desit(i),i=1,6)/'astrometric azimuth      ',
     . 'horizontal angle         ','horizontal direction     ',
     . 'baseline length          ','zenith height            ',
     . 'leveling                 '/
      data (desit(i),i=11,16)/'astrometric azimuth rate ',
     . 'horizontal angle rate    ','horizontal direction rate',
     . 'baseline length rate     ','zenith height rate       ',
     . 'leveling rate            '/
      data (desit(i),i=21,24)/
     . '3-D geocentric coordinate','3-D geocentric velocity  ',
     . '3-D geodetic coordinate  ','3-D geodetic velocity    '/
      data (desit(i),i=25,28)/
     . '3-D geoc baseline vector ','3-D geoc bsln rate vector',
     . '3-D geod baseline vector ','3-D geod bsln rate vector'/
      data (desit(i),i=31,34)/
     . '3-D spherical coordinate ','3-D spherical velocity   ',
     . '3-D Cartesian coordinate ','3-D Cartesian velocity    '/
c
c     print title
      print*,' ------  Residual calculation begins ......'
      print*,'   Observation Type    # Obs.  Pre-rms  Post-rms'
      write (ifile,*) '*  ----------  Residual list  ----------'
c
      rewind (idatf)
      read (idatf,*) nsit
      do i = 1,nsit
         read (idatf,'(2x)') 
      enddo
    
      iauxsav = 0
c     form coefficients for every observation
      read (idatf,*) iexp
      igo = 0
      iobs = 0
      rms1 = 0.0d0
      rms2 = 0.0d0
      do 50 ie = 1,iexp
         if (it.ne.2) iazsav = 0
c
c        read experiment index, obs. type and number
         read (idatf,*) ie1,it,iobs1
         iobs0 = iobs
         if (it.le.20) k = iobs1
         if (it.gt.20.and.it.le.40) k = iobs1*3
         iobs = iobs+k
         write (ifile,'(''* exp. : '',i3,5x,''type:'',
     .          a25,5x,''No:'',i4)')
     .          ie,desit(it),k
         if (it.le.3) then
            unit(1) = '(deg.)'
            unit(2) = '(deg.)'
            unit(3) = '(sec.)'
         else
            if (it.ge.11.and.it.le.13) then
               unit(1) = '(sc/y)'
               unit(2) = '(sc/y)'
               unit(3) = '(sc/y)'
            else
            if (it.eq.22.or.it.eq.24.or.(it.ge.14.and.it.le.20)) then
               unit(1) = '(mm/y)'
               unit(2) = '(mm/y)'
               unit(3) = '(mm/y)'
            else
               unit(1) = '(  m )'
               unit(2) = '(  m )'
               unit(3) = '( mm )'
            endif
            if (it.eq.26.or.it.eq.28) then
               unit(1) = '(mm/y)'
               unit(2) = '(mm/y)'
               unit(3) = '(mm/y)'
            endif
            endif
         endif
         write (ifile,70) (unit(k),k=1,3), unit(3),unit(3)
         prerms = 0.0d0
         pstrms = 0.0d0
         wrms  = 0.0d0
         imsit = 0
         iefob = 0
      do 20 ib = iobs0+1,iobs0+iobs1
         ibd = ib-iobs0
         call obsred(idatf,it,dt,robs,err,vrho,l_obs,isit1,isit2)
         time1 = dt+rtime
         if (it.lt.21.and.(jtoi(isit1).eq.0.or.jtoi(isit2).eq.0))
     .      goto 20
c
c        3-D data as model coordinate
         if (smode.eq.3.and.it.gt.20) goto 20
         if (it.gt.30) goto 20
         if (it.gt.20.and.jtoi(isit1).eq.0) goto 20
c
         dt = time1-rtime
         iefob = iefob+l_obs

c
c        get calculated value
         if (it.gt.20.and.it.le.30) then
            call calc_omc_block(it,isit1,isit2,dt,ib,ibd,robs,
     .           cal,omc,oc,ill)
         else
            call calval(it,isit1,isit2,dt,ib,ibd,robs(1),
     .           cal(1),omc(1),oc(1),ill)
         endif
c        if (it.eq.3.and.dabs(oc).le.1.0d-15)
         if (it.eq.3.and.dabs(oc(1)).le.1.0d-15)
     .      imsit = imsit+1
c
         fac = 1.0d0
         if (it.le.3) fac = (3.6d3*rtod)**2
         if (it.ge.4.and.it.le.30) fac = 1.0d6
         if (.not.ill.and.l_obs.eq.3) then
            covo(1) = fac*err(1)**2
            covo(3) = fac*err(2)**2
            covo(6) = fac*err(3)**2
            covo(2) = fac*vrho(1)*err(1)*err(2)
            covo(4) = fac*vrho(2)*err(1)*err(3)
            covo(5) = fac*vrho(3)*err(2)*err(3)
c           get weighting matrix
            call nrmscl(covo,robs,scale,3,1,0)
            call cholsk(covo,temp,1,3,ie1)
            call nrmscl(covo,robs,scale,3,2,0)
            call atwa(1,3,omc,covo,wr,temp,1)
            wrms = wrms+wr
         endif
         if (.not.ill.and.l_obs.eq.1) then
            wr = omc(1)**2/fac/(err(1)**2)
            wrms = wrms+wr
         endif
c
c        rms
         fac = dsqrt(fac)
         do il = 1,l_obs
            if (ill) then
               symb = '*'
               iefob = iefob-1
            else 
               prerms = prerms+oc(il)**2
               pstrms = pstrms+omc(il)**2
               symb = ' '
            endif
            sensibilite = robs(il)*err(il)*fac
c
c           output residuals
c           if (it.eq.2) then   ! it=2=horizontal angle
c              write (ifile,85) symb,sname(isit1),sname(isit2),
c    .         time1,robs(il),cal(il),err(il)*fac,omc(il),oc(il)
c    .         sensibilite
c85            format (a1,1x,a8,2x,a8,1x,f10.3,2(1x,f15.5,1x),
c    .                 f12.5,2f15.5,f12.5)
c           else 
               write (ifile,80) symb,sname(isit1),sname(isit2),
     .         time1,robs(il),cal(il),err(il)*fac,omc(il),oc(il),
     .         (omc(il)/(err(il)*fac))
c           endif
         enddo
 20   continue     
         if (smode.eq.3.and.it.gt.20) goto 10
         if (it.gt.30) goto 10
         ib = iefob-imsit
         rms1 = rms1+prerms
         rms2 = rms2+pstrms
         ieff = ieff+ib
         if (ib.gt.0) then
            prerms = dsqrt(prerms/ib)
            pstrms = dsqrt(pstrms/ib)
            wrms = dsqrt(wrms/ib)
         else
            prerms = 0.0d0
            pstrms = 0.0d0
            wrms = 0.0d0
         endif
         
         print '(a25,i4,2f10.2)',desit(it),ib,prerms,pstrms
         write (ifile,*) '*.............. statistics ................'
         write (ifile,*) '*  observations:',ib
         write (ifile,*) '*  parameters:',' <need help from Tom>'
         write (ifile,*) '*  prefit rms :',unit(3),prerms
         write (ifile,*) '*  postfit rms:',unit(3),pstrms
         write (ifile,*) '*  weighted rms:   ',wrms
         write (ifile,*) '*..........................................'
         write (ifile,*) '*   '
 10   continue
c
 70   format ('*   Site1     Site2      Time   ',
     .        4x,' Obs.',a6,4x,'Adjusted ',a6,2x,'Sigma',
     .        a6,2x,'O-Calc',a6,4x,'Prefit ',a6,2x,'v/s')
 80   format (a1,1x,a8,2x,a8,1x,f10.4,2(1x,f15.5,1x),f12.5,2f15.5,f9.1)
      igo = igo+igp
 50   continue
      nobs = iobs
      if (ieff.le.0) goto 90
      rms1 = dsqrt(rms1/ieff)
      rms2 = dsqrt(rms2/ieff)
      write (ifile,*) '*   '
      write (ifile,*) '*---------------- statistics -----------------'
      write (ifile,*) '*  observations:',ieff
      write (ifile,*) '*  parameters:',nlive-ibnum-icnum-ncht 
      write (ifile,*) '*  prefit rms :','(unit??)',rms1
      write (ifile,*) '*  postfit rms:','(unit??)',rms2
      write (ifile,*) '*---------------------------------------------'
 90   continue
c
      return
      end
c--------------------------------------------------------------
      subroutine calval(it,isit1,isit2,dt,ib,ibd,obs,cal,omc,oc,ill)
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
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer it,isit1,isit2,ib,ibd,indsav,istsav,iazsav
      integer i,i1,i2,id,isit3,is1,iauxsav,iterm
      dimension tempc(18)
      dimension indsav(12),coesav(12),coef(30)
      logical ill
      common/saved/azisav,coesav,obssav,esave,tsave,iazsav,indsav
      common/save2/istsav,iauxsav
c
c     from m to mm
      fac = 1.0d3
      small = 1.0d-15
      ill = .false.
c
c     count the index of extra parameters
      if (it.eq.2) then
         if (ibd.eq.1.or.isit1.ne.istsav.or.dabs(dt-tsave)
     .       .gt.1.0d-3.or.dabs(obs).lt.1.0d-5) iauxsav = iauxsav+1
c        unit in second for the auxiliary parameter (angle)
         coef(13) = -dtor/3.6d3
      endif
c
c     astrometric azimuth
      if (it.eq.1) then
         iterm = 12
c        ba = basazm(isit1,isit2,0.0d0,1)
cmk      this needs to be computed at the observation epoch
cmk      not the reference epoch!
         ba = basazm(isit1,isit2,dt,1)
c        get coefficients
         call azimuc(isit1,isit2,dt,coef,1)
c        get residual adjustment
         call calc_res(2,isit1,isit2,isit3,iterm,dt,coef,adj)
         cal = ba+adj
c        conventional unit (obs,cal: degree, omc: second)
         obs = obs*rtod
         cal = cal*rtod
         if (cal.lt.0.0d0) cal = cal+360.0d0
c        
         omc = (obs-cal)*3600.0d0
c        oc values close to 360 degrees should be oc=oc-360
         oc = obs-temp*rtod
         if (oc.gt.180.0d0) oc = oc-360.0d0
         if (oc.le.-180.0d0) oc = oc+360.0d0
         oc = (obs-ba*rtod)*3600.0d0
         if (dabs(oc).gt.cria) ill = .true.
      endif
c
c     horizontal angle
      if (it.eq.2) then
         iterm = 13

         if (ibd.eq.1.or.isit1.ne.istsav.or.dabs(dt-tsave)
     .       .gt.1.0d-3.or.dabs(obs).lt.1.0d-5) then 

            azisav = basazm(isit1,isit2,dt,1)
            obssav = azisav-obs
            istsav = isit1
            tsave = dt

            call azimuc(isit1,isit2,dt,coef,1)
            call calc_res(2,isit1,isit2,iazsav,iterm,dt,coef,adj)

            cal =  dmod( (obs+adj)*rtod,360.0d0) 
            obs = obs*rtod

            if (cal.lt.0.0d0) cal = cal + 360.0d0 

            omc= obs-cal
c
            if (omc.ge.180.0d0)  omc = omc - 360.0d0
            if (omc.le.-180.0d0) omc = omc + 360.0d0                 
c           
            omc = omc * 3600.0d0
            oc = 0.0d0
            goto 100
         endif

         ba = basazm(isit1,isit2,dt,1)
         temp = ba-obssav

         if (temp.lt.0.0d0) temp = pi*2.0d0+temp

         call azimuc(isit1,isit2,dt,coef,1)
         call calc_res(2,isit1,isit2,iazsav,iterm,dt,coef,adj)

         cal = temp+adj
c        conventional unit (obs,cal: degree, omc: second)
         obs = obs*rtod
         cal = dmod(cal*rtod,360.0d0)
         if (cal.lt.0.0d0) cal = cal + 360.0d0 
         omc =  obs-cal
c
         if (omc.ge.180.0d0)  omc = omc - 360.0d0
         if (omc.le.-180.0d0) omc = omc + 360.0d0      
c
         omc = omc * 3600.0d0
         oc = dmod(obs-temp*rtod,360.0d0)
         if (oc.ge.180.0d0)  oc = oc - 360.0d0
         if (oc.le.-180.0d0) oc = oc + 360.0d0
         oc = oc * 3600.0d0
c        oc = dmod(oc,360.0d0) * 3600.0d0
         if (dabs(oc).gt.cria) ill = .true.
      endif
c
c     horizontal direction
      if (it.eq.3) then
         iterm = 18

         if (ibd.eq.1.or.isit1.ne.istsav.or.dabs(dt-tsave)
     .       .gt.1.0d-3.or.dabs(obs).lt.1.0d-5) then 

c           azisav = basazm(isit1,isit2,0.0d0,1)
cmk         this needs to be calculated at the observation epoch
cmk          not the reference epoch
            azisav = basazm(isit1,isit2,dt,1)
            obssav = obs
            istsav = isit1
            tsave = dt

            call azimuc(isit1,isit2,dt,coesav,1)

            iazsav = isit2
            obs = obs*rtod
            cal = obs
            omc = 0.0d0
            oc = 0.0d0
            goto 100
         endif

c        ba = basazm(isit1,isit2,0.0d0,1)
cmk         this needs to be calculated at the observation epoch
cmk          not the reference epoch
         ba = basazm(isit1,isit2,dt,1)
         temp = ba-azisav
         if (temp.lt.0.0d0) temp = pi*2.0d0+temp
         call azimuc(isit1,isit2,dt,tempc,1)

         do i = 1,6
            coef(i) = tempc(i)-coesav(i)
            coef(i+6) = tempc(i+6)
            coef(i+12) = -coesav(i+6)
         enddo

         call calc_res(3,isit1,isit2,iazsav,iterm,dt,coef,adj)

         cal = temp+adj
c        conventional unit (obs,cal: degree, omc: second)
         obs = (obs-obssav)*rtod

         if (obs.lt.0.0d0) obs = obs+360.0d0
         cal = cal*rtod
         if (cal.lt.0.0d0) cal = cal+360.0d0
         omc = (obs-cal)*3600.0d0
c        oc values close to 360 degrees should be oc=oc-360
         oc = obs-temp*rtod
         if (oc.gt.180.0d0) oc = oc-360.0d0
         if (oc.le.-180.0d0) oc = oc+360.0d0
         oc = oc*3600.0d0
         if (dabs(oc).gt.cria) ill = .true.

      endif
c
c     baseline length
      if (it.eq.4) then
         iterm = 12
         bl = baslen(isit1,isit2,dt,3,2)
         call baslc(isit1,isit2,bl,dt,coef,1)
         call calc_res(2,isit1,isit2,isit3,iterm,dt,coef,adj)
         cal = bl+adj
c        conventional unit (mm)
         omc = (obs-cal)*fac
         oc = (obs-bl)*fac
         if (dabs(oc).gt.cril*fac) ill = .true.
      endif
c
c     baseline length rate
      if (it.eq.14) then
         bl = baslen(isit1,isit2,dt,3,2)
         call baslcv(isit1,isit2,bl,dt,coef,1)
         call getadj(2,isit1,isit2,isit3,tempc)
         call axb(1,12,1,coef,tempc,adj,1,0)
         tv = vx(isit2)-vx(isit1)
         xd = (x(isit2)-x(isit1)+tv*dt)*tv
         tv = vy(isit2)-vy(isit1)
         yd = (y(isit2)-y(isit1)+tv*dt)*tv
         tv = vz(isit2)-vz(isit1)
         zd = (z(isit2)-z(isit1)+tv*dt)*tv
         cal = ((xd+yd+zd)/bl+adj)*fac
c        unit: mm/yr
         omc = obs-cal
         oc = omc+adj*fac
         if (dabs(oc).gt.criv*fac) ill = .true.
      endif
c
c     3-D geocentric coordinate
      if (it.eq.21) then
c        determine x, y or z
         print* , 'In subroutine RESIDL.F'
         i2 = ibd-ibd/3*3
         if (i2.eq.0) i2 = 3
         i1 = (isit1-1)*6+i2
         id = map(i1)
         adj = 0.0d0
         if (id.gt.0) adj = bnorm(id)
         if (i2.eq.1) cal = x(isit1)+adj
         if (i2.eq.2) cal = y(isit1)+adj
         if (i2.eq.3) cal = z(isit1)+adj
         omc = (obs-cal)*fac
         oc = omc+adj*fac
         if (dabs(oc).gt.cric) ill = .true.
      endif
c
c     3-D geocentric velocity
      if (it.eq.22) then
c        determine vx or vy
         i2 = ibd-ibd/3*3
         if (i2.eq.0) i2 = 3
         is1 = jtoi(isit1)
         if (i2.eq.1.and.is1.gt.0) cal = slnvx(is1)
         if (i2.eq.2.and.is1.gt.0) cal = slnvy(is1)
         if (i2.eq.3.and.is1.gt.0) cal = slnvz(is1)
         if (i2.eq.1.and.is1.le.0) cal = vx(isit1)
         if (i2.eq.2.and.is1.le.0) cal = vy(isit1)
         if (i2.eq.3.and.is1.le.0) cal = vz(isit1)
         if (i2.eq.1) oc = obs-vx(isit1)
         if (i2.eq.2) oc = obs-vy(isit1)
         if (i2.eq.3) oc = obs-vz(isit1)
         omc = obs-cal
         if (dabs(oc).gt.criv*fac) ill = .true.
      endif
c
c     3-D geodetic velocity
      if (it.eq.24) then
c        determine ve or vn
         i2 = ibd-ibd/3*3
         if (i2.eq.0) i2 = 3
         is1 = jtoi(isit1)
         if (i2.eq.1.and.is1.gt.0) cal = slnve(is1)
         if (i2.eq.2.and.is1.gt.0) cal = slnvn(is1)
         if (i2.eq.3.and.is1.gt.0) cal = slnvu(is1)
         if (i2.eq.1.and.is1.le.0) cal = ve(isit1)
         if (i2.eq.2.and.is1.le.0) cal = vn(isit1)
         if (i2.eq.3.and.is1.le.0) cal = vu(isit1)
         if (i2.eq.1) oc = obs-ve(isit1)
         if (i2.eq.2) oc = obs-vn(isit1)
         if (i2.eq.3) oc = obs-vu(isit1)
         omc = obs-cal
         if (dabs(oc).gt.criv*fac) ill = .true.
      endif
c
 100  continue
      return
      end
c
c--------------------------------------------------------------
      subroutine getadj(isit,isit1,isit2,isit3,tempc)
c
c     get adjustment in Cartesian coordinate
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer isit,isit1,isit2,isit3,i,is,is1,j,j1,j2,id
      dimension tempc(isit*6)
c
      do 20 i = 1,isit
         if (i.eq.1) is = isit1
         if (i.eq.2) is = isit2
         if (i.eq.3) is = isit3
         is1 = jtoi(is)

         do j = 1,3
            j1 = (is-1)*6+j
            j2 = (i-1)*6+j
            tempc(j2) = 0.0d0
            id = map(j1)
            if (id.gt.0) tempc(j2) = bnorm(id)
         enddo

         do j = 1,3
            j1 = (is1-1)*3+j
            j2 = (i-1)*6+3+j
            tempc(j2) = solutn(j1)
         enddo

 20   continue
c 
      return
      end
c
