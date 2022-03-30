      subroutine get_res_sgl(ifile,idatf,iobs1,it,prerms,pstrms,
     .                        unit,wrms,iefob,iobs,imsit)
c
c
c*aut modified by FERHAT Gilbert 
c*aut July 26, 1995
c
c*rol to have the sensibility
c*rol i.e the displacement due to a a small horizontal angle
c*rol here the small angle is 2 decimilligrad = 0.648 seconds
c 
c     get it < 30 residuals
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
      integer ifile,idatf,indsav,iazsav,iobs
      integer iobs0,iobs1,it,imsit,iefob,ib,ibd,ie1,k,isit1
      integer isit2,istsav,iauxsav,l_obs,il
c     real*8  err,vrho,robs,omc,oc,cal
      dimension  err(3),vrho(3),robs(3),omc(3),oc(3),cal(3)
      dimension indsav(12),coesav(12),covo(6),temp(6)
      character*1 symb
      logical ill
      common/saved/azisav,coesav,obssav,esave,tsave,iazsav,indsav
      common/save2/istsav,iauxsav
c
      if (it.ne.2) iazsav = 0
      iobs0 = iobs
      if (it.le.20) k = iobs1
      if (it.gt.20.and.it.le.40) k = iobs1*3
      iobs = iobs+k
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
c
c     iomode (17) is for file # 30 LAMBERT file
c     it = 2 = horizontal angle 
cc    print*, 'it= ',it
cc    print*, 'iomode(17)= ',iomode(17) 
      if ((it.eq.2.or.it.eq.3).and.iomode(17).eq.1) then
      write (ifile,75) (unit(k),k=1,3), unit(3),unit(3),'cm','(m)' 
     .   ,'(m)'
cc    print*, '2 it= ',it
      else
        write (ifile,70) (unit(k),k=1,3), unit(3),unit(3)
      endif
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
c
c           output residuals
            if ((it.eq.2.or.it.eq.3.or.it.eq.4)
     .         .and.iomode(17).eq.1) then   ! it=2=horizontal angle
cc             print*,' isit1= ',isit1
               call sensibility(isit1,isit2,sens,dd,dh)
cc             print*,' isit1= ',isit1
               write (ifile,85) symb,sname(isit1),sname(isit2),
     .         time1,robs(il),cal(il),err(il)*fac,omc(il),oc(il)
     .         ,sens,dd,dh
               write (*,85) symb,sname(isit1),sname(isit2),
     .         time1,robs(il),cal(il),err(il)*fac,omc(il),oc(il)
     .         ,sens,dd,dh
 85            format (a1,1x,a8,2x,a8,1x,f10.3,2(1x,f15.5,1x)
     .         ,f12.5,2f15.5,f8.1,1x,f10.3,f10.3)
            else 
               write (ifile,80) symb,sname(isit1),sname(isit2),
     .         time1,robs(il),cal(il),err(il)*fac,omc(il),oc(il)
            endif
cc            write (ifile,80) symb,sname(isit1),sname(isit2),
cc     .       time1,robs(il),cal(il),err(il)*fac,omc(il),oc(il)
         enddo
      if (iomode(17).eq.1) then
         close(30)
         close(31)
      endif
 20   continue     
c
 70   format ('*   Site1     Site2      Time   ',
     .        4x,' Obs.',a6,6x,'Calc.',a6,3x,'Sigma',
     .        a6,4x,'Adj.',a6,6x,'O-C ',a6) 
 75   format ('*   Site1     Site2      Time   ',
     .        4x,' Obs.',a6,6x,'Calc.',a6,3x,'Sigma',3x,
     .        a6,5x,'O-C',a6,4x,'Prefit',a6,2x,
     .        a2,'  Dist',a3,'    dh',a3)
 80   format (a1,1x,a8,2x,a8,1x,f10.3,2(1x,f15.5,1x),f12.5,2f15.5)
 
      return
      end
c
c
c *******************************************************************
c
      subroutine sensibility (isit1,isit2,sens,dd,dh)
c
c*aut FERHAT Gilbert
c*aut GRGS
c
c*ver July, 26 1995
c
c*par         input parameters 
c
c*par isit1 : site number 1
c     isit2 : site number 2
c
c*par        output parameters
c
c     sens : sensibilite
c     dd   : horizontal distance
c     dh   : denivele
c
c*************************************************************************
c   
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      real*8  x_lambert(maxsit)
      real*8  y_lambert(maxsit)
      real*8  h(maxsit)
      real*8  dd,dh,sens,an,obs3,sigma_angle
      real*8  xq(maxsit),yq(maxsit)
c     dd   :  horizontal distance
      integer izone(100),nb_site
      character*8  numero(maxsit)
      character*9  askipp
      character*3  code_alti
      character*22 skipp,line
      character*64 sub_file
c
c     character*8 sname2(maxsit)
      character*8 sitnam
      character*128 buffer
c
c    
      integer isite2b(maxsit)
c     integer idatf,iobs,ifile,illsit
      integer i,j,ie1,it
      integer nblen,n1,nlen,match_name
c     integer iuse(maxsit),vuse(maxsit)
      integer isit1,isit2,ios,icontrol1,icontrol2,i1,i2
      integer index,ijk,iilen,iii,nb_exp,iobs1,iobs2
      integer nlen1,im,jj,k,isite1,lift_arg
      character*8  name1,sitea,point_viseb(100)
c     character*64 subfil
ccc      equivalence (iuse,map)
ccc      equivalence (vuse,map(maxsit+1))
c
      i = 0
c  
c*************************************************************************
c      read the file # 30, Lambert file or cartographic coordinates xyh
c*************************************************************************
c
      if (iomode(17).eq.1)  then
         open (30,file=lambert,status='old',err=60)
         rewind(30)
         read(30,'(a)',iostat=ios,end=40,err=50) line
c        print*, line
         read(30,'(a)',iostat=ios,end=40,err=50) line
c        print*, line
 10      continue
         i = i + 1
         read(30,20,iostat=ios,end=40,err=50) askipp
     .    ,skipp,x_lambert(i),y_lambert(i),h(i),
     .    code_alti,izone(i)
         numero(i) = askipp(1:6)//askipp(8:9)
 20      format(a9,1x,a22,1x,3(f14.3),2x,a3,7x,i1)
c        write(*,25) numero(1)
c    .    ,skipp,x_lambert(1),y_lambert(1),h,
c    .    code_alti,izone(1)
c25      format(a8,1x,a22,1x,3(f14.3),2x,a3,7x,i1)
         goto 10 
      else
         print*, 'Suicide due to  LAMBERT file ...'
         stop
      endif
c
 40   continue
      icontrol1=0
      icontrol2=0
      nb_site=i
c     print*, 'nb_site= ',nb_site
c
      do i = 1,nb_site
          if (sname(isit1).eq.numero(i)) then
            x1 = x_lambert(i)
            y1 = y_lambert(i)
            h1 = h(i)
            i1 =i
            icontrol1=1
cc          print*, 'In subroutine sensiblity'
cc          print*, 'isit1= ', isit1
cc          print*, 'sname',sname(isit1),numero(i) 
           endif
c         print*, 'sname',sname(isit1),'  ',numero(i), 'i= ',i
      enddo      
c    
      if (icontrol1.eq.1) then
       do j = 1,nb_site
        if (sname(isit2).eq.numero(j)) then
            x2 = x_lambert(j)
            y2 = y_lambert(j)
            h2 = h(j)
            i2 = j
            dd= dsqrt( (x2-x1)**2 + (y2-y1)**2 )
            dh = h2 - h1
c           dd is in meters, sens in cm
c           sens for a angle variation of 2 dmgr = 0.648 seconds
c           this is for a first order network
            sens = 100.0d0 * dd * pi / (180.0d0*3600.0d0) * 0.648d0
            if (izone(i1).ne.izone(i2)) then
               print*, '2 sites in different Lambert zones'
               print*, sname(i1),sname(i2)
            endif
c           print*, sname(i1),sname(i2),' sens=',sens
         endif
       enddo
      endif
c
c**************************************************************
c     read sequential file list
c***************************************************************
c
      if (index.eq.1) then
          goto 46
      endif
      rewind(14)
c     read site number
      read (14,*) nsit
c     print*, 'Sens : nsit = ',nsit
c
c     read site name(id)
      do i = 1,nsit
         if (fcode.le.2) then
            read (14,'(a8)') sitnam
c          print*, 'Sens : sitnam', sitnam       
         endif
      enddo
c     read number of subfile to be read
      read (14,*) ijk
c
      read(14,'(a)',err=140,end=144) sub_file
 144  continue
      rewind(14)
c
      print *,' The subfile is = ',sub_file
      iilen = nblen(sub_file)
      print *,' The subfile is = ',sub_file(1:iilen)
c     close(31)
c      open(54,file='pro52_1.obs.12july',status='old',err=130)
      iii = 31
      open (iii,file=sub_file(1:iilen),status='old',err=130)
c
      read (iii,*) n1
      print*, 'n1 = ',n1
      do 21 j = 1,n1
         read (iii,'(a8)') name1
c        print*, 'name1= ',name1
ccc         nlen = nblen(name1)
ccc         isit1 = match_name(nsit,nlen,sname,name1)
ccc         if (isit1 .le. 0) goto 21
c?         if (iuse(isit1).lt.minic) iuse(isit1) = iuse(isit1)+1
c?         if (iuse(isit1).ge.minic) vuse(isit1) = vuse(isit1)+1
 21   continue
      read (iii,*) ie1
      print*, ie1
 22   continue
      read (iii,*,end=46) ie1,it,nb_exp
      print*, '***************************************************'
      print*, ie1
c     print*, ie1,' ',it,' ',nb_exp
ccc   if (it.eq.31.or.it.eq.32) jaux = jaux+6
c
c     it = 2 = horizontal angle
c     it = 3 = horizontal direction
c     for a relevement, nb_exp is at least 3
c
      if ((it.eq.2.or.it.eq.3).and.nb_exp.ge.3) then
          goto 1   
      else
          do i = 1 ,nb_exp
             read(iii,'(a)',end=46) buffer
          enddo
          goto 22
      endif
c
 1    continue
      do i = 1, nb_exp
cc       read(iii,120,end=200) an,iobs1,iobs2,obs3
         read(iii,'(a)',end=46) buffer
         print*,'in get_res_sgl.f'
c        print*, buffer
         read(buffer,*) an,iobs1,iobs2,obs3,
     &                  sigma_angle
         nlen1 = lift_arg(buffer,sitea,6)
         nlen1 = lift_arg(buffer,point_viseb(i),7)
         nlen1 = nblen(sitea)
         isite1 = match_name(nsit,nlen1,sname,sitea)
         if (isite1 .le. 0) goto 300                    
         nlen1 = nblen(point_viseb)
         isite2b(i)=match_name(nsit,nlen1,sname,point_viseb(i))
         if (isite2b(i) .le. 0) goto 300                    
c 120    format(2x,f9.4,4x,f9.3,5x,f5.1,9x,a8,2x,a8)
 120     format(2x,f9.4,3x,i3,3x,i2,3x,f11.8,3x,
     &          f6.4,2x,a8,2x,a8)
c        print*, 'sitea = ',sitea
c        print*, 'point_viseb i = ',point_viseb(i),' i= ',i
      enddo
c
      do i = 1,nb_site
         if (sname(isite1).eq.numero(i)) then
             xm = x_lambert(i)
             ym = y_lambert(i)
             hm = h(i)
             im =i
             icontrol1=1
c            print*, 'sname ',sname(isite1),numero(i), 'xyzm= ',
c    .               xm,' ',ym
          endif
c         print*, 'sname',sname(isite1),'  ',numero(i), 'i= ',i
      enddo      
c    
      if (icontrol1.eq.1) then
       do jj = 1 , nb_exp
        do j = 1,nb_site
         if (sname(isite2b(jj)).eq.numero(j)) then
            xq(jj) = x_lambert(j)
            yq(jj) = y_lambert(j)
c           print*, 'sname b ',sname(isite2b(jj)),numero(j),
c    .      'xyz q= ',xq(jj),' ',yq(jj),'jj= ',jj
cc          iq(jj) = j
c
c            if (izone(i1).ne.izone(i2)) then
c               print*, '2 sites in different Lambert zones'
c               print*, sname(i1),sname(i2)
c            endif
c           print*, sname(i1),sname(i2),' sens=',sens
          endif
         enddo
       enddo
      endif
c      iqq = 2 
c      if (iqq.eq.2) stop     
c
      do i = 1 , nb_exp
         do j = i+1 , nb_exp
            do k = i+2 , nb_exp
c               ii = i
c               ij = j
c               ik = k
c               if (ii.eq.ij) then
c                  goto 44
c               endif
c               if (ik.eq.ij) then
c                  goto 45
c               endif
c               if (i.eq.k) then
c                  goto 44
c               endif
               xa = xq(i)
               ya = yq(i)
               xb = xq(j)
               yb = yq(j)
               xc = xq(k)
               yc = yq(k)
               if (xa.eq.xb.and.ya.eq.yb) then
                  goto 44
               endif
               if (xa.eq.xc.and.ya.eq.yc) then
                  goto 45
               endif
               if (xb.eq.xc.and.yb.eq.yc) then
                  goto 45
               endif
ccc            print*, 'ijk=  ',i,' ',j,' ',k
ccc            print*, 'xm ym ',xm,' ',ym
ccc            print*, 'xa ya ',xa,' ',ya
ccc            print*, 'xb yb ',xb,' ',yb
ccc            print*, 'xc yc ',xc,' ',yc              
c              Gisements des lieux in degree
               gam = gisement(xm,ym,xa,ya)
               gbm = gisement(xm,ym,xb,yb)
               gcm = gisement(xm,ym,xc,yc)
               gab = gisement(xb,yb,xa,ya)
               gbc = gisement(xc,yc,xb,yb)
c
               gma = gisement(xa,ya,xm,ym)
               gmb = gisement(xb,yb,xm,ym)
               gmc = gisement(xc,yc,xm,ym)
ccc            print*, 'gma= ',gma
ccc            print*, 'gmb= ',gmb
ccc            print*, 'gmc= ',gmc
ccc            print*, ' '
c
c   points A, B and C in topographic sens
c
               bb = max(gma,gmb,gmc)
               cc = min(gma,gmb,gmc)
c
c   bb is is the biggest azimut
c
               if (bb.eq.gma) then
                   xc = xq(i)
                   yc = yq(i)
               endif
               if (bb.eq.gmb) then
                   xc = xq(j)
                   yc = yq(j)
               endif
               if (bb.eq.gmc) then
                   xc = xq(k)
                   yc = yq(k)
               endif
c
c    cc is the lowest azimut
c             
               if (cc.eq.gma) then
                   xa = xq(i)
                   ya = yq(i)
               endif
               if (cc.eq.gmb) then
                   xa = xq(j)
                   ya = yq(j)
               endif
               if (cc.eq.gmc) then
                   xa = xq(k)
                   ya = yq(k)
               endif
c
c search for the value, cc < middle < bb
c
               if (gma.gt.cc.and.gma.lt.bb) then
                   xb = xq(i)
                   yb = yq(i)
               endif
c
               if (gmb.gt.cc.and.gmb.lt.bb) then
                   xb = xq(j)
                   yb = yq(j)
               endif
c
               if (gmc.gt.cc.and.gmc.lt.bb) then
                   xb = xq(k)
                   yb = yq(k)
               endif
c
c recompute the azimut. Now the points A, B, C are sorted by azimut
c
ccc            print*, 'Now the points A B C are sorted by azimut'
ccc            print*, 'xm ym ',xm,' ',ym
ccc            print*, 'xa ya ',xa,' ',ya
ccc            print*, 'xb yb ',xb,' ',yb
ccc            print*, 'xc yc ',xc,' ',yc              
               dam = dsqrt ( (xa-xm)**2 + (ya-ym)**2 )
               dbm = dsqrt ( (xb-xm)**2 + (yb-ym)**2 )
               dcm = dsqrt ( (xc-xm)**2 + (yc-ym)**2 )
               dab = dsqrt ( (xb-xa)**2 + (yb-ya)**2 )
               dac = dsqrt ( (xc-xa)**2 + (yc-ya)**2 )
               dbc = dsqrt ( (xc-xb)**2 + (yc-yb)**2 )
ccc            print*,'dam = ',dam
ccc            print*,'dbm = ',dcm
ccc            print*,'dcm = ',dcm
ccc            print*,'dab = ',dab
ccc            print*,'dac = ',dac
ccc            print*,'dbc = ',dbc
c
c              s : sensibilite des lieux
c              deltaphi = 2 dmgr 
               deltaphi = 2.0d0
               sab = 0.157d0 * deltaphi * dam*dbm * 0.001d0 / dab
               sac = 0.157d0 * deltaphi * dam*dcm * 0.001d0 / dac
               sbc = 0.157d0 * deltaphi * dbm*dcm * 0.001d0 / dbc
ccc            print*,'sab = ',sab
ccc            print*,'sac = ',sac
ccc            print*,'sbc = ',sbc
c
c              Gisements des lieux in degree
               gam = gisement(xm,ym,xa,ya)
               gbm = gisement(xm,ym,xb,yb)
               gcm = gisement(xm,ym,xc,yc)
               gab = gisement(xb,yb,xa,ya)
               gbc = gisement(xc,yc,xb,yb)
c
ccc            print*, 'gam= ',gam
ccc            print*, 'gbm= ',gbm
ccc            print*, 'gcm= ',gcm
ccc            print*, 'gab= ',gab
ccc            print*, 'gbc= ',gbc
c
               garc_ab = gam + gbm - gab
               garc_bc = gbm + gcm - gbc
               garcab = dmod(garc_ab,360.0d0)
               garcbc = dmod(garc_bc,360.0d0)
ccc            print*, 'garcab= ',garcab
ccc            print*, 'garcbc= ',garcbc 
               if (garcab.lt.0.0d0) garcab = garcab + 360.0d0
               if (garcbc.lt.0.0d0) garcbc = garcbc + 360.0d0 
c             
ccc             print*, 'garcab= ',garcab
ccc             print*, 'garcbc= ',garcbc 
               a  = garcbc - garcab
c               a1 = 360.0d0 - a
c               a2 = 360.0d0 + a
                if (a.le.0.0d0)  a = a + 360.0d0
c               if (a1.le.180.0d0) then a = a2
c               
               if (a.le.180.0d0) then
ccc               print*, 'a= ',a,' OK below 180 degres'
               else
ccc              print*, ' pb a above 180 degres a= ', a
c                 stop
               endif
ccc            print*, 'a in degree',a
c              a in radian
               arad = a * dtor
c              a never equals zero because for a relevement we never 
c              chose points alined !!!!
c
               sm1 = 
     .  dsqrt ( sab**2 + sbc**2 + 2*sab*sbc*dcos(arad) ) / dsin(arad)
               sm2 = sbc / dsin(arad)
               sm3 = sab / dsin(arad)
c
ccc            print*, 'site = ',sitea,' sm1= ',
ccc  .          sm1,' ',sm2,'  ',sm3,' nb_exp= ',nb_exp
ccc            print*, 'smax= ',max( dabs(sm1),dabs(sm2),dabs(sm3) )
c
ccc            print*, '-----------------------------------------------'
 44            continue
 45            continue
          enddo
         enddo
      enddo
      goto 22
 46   continue
      index = 1
c
c     close(iii)
c      close(30)
      return
 50   print *, 'Error in reading LAMBERT file'
      stop
 60   print *, 'Error in opening LAMBERT file'
      stop
 130  continue
      print *, 'Stop in sens : Cannot read the sub filename'
      stop
 140  continue
      print *, 'Stop in sens'
      print *, 'Cannot read the name of the sequential'
      print *, 'specified in the file list'
      stop
 300  print *,' mismatch site name at SENS: ',
     .   sitea,' ',point_viseb(i),nlen,isite1,isite2b(i)
      stop
c
      end
c
c*****************************************************************
c
      real*8 function gisement(x1,y1,x2,y2)
c
c
c*par input parameters
c
c*par x1,y1 
c*par x2,y2
c
c*par output parameters
c 
c*par gisement : azimut in degree g21 : 2 --> 1
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      real*8 gi
c
ccc   print*, 'x y 12',x1,y1,x2,y2
      dx = x1 - x2
      dy = y1 - y2
ccc   print*, 'dx dy ',dx,dy
c
       if (dy.eq.0) then
         if (dx.gt.0.0d0) gi= pi/2.0d0
         if (dx.lt.0.0d0) gi=3.0d0*pi/2.0d0  
       endif
c
       if (dx.gt.0.0d0.and.dy.gt.0.0d0) then 
c          print*, '  dx>0         dy>0'
           gi = datan( dx /dy )
       endif
       if (dx.gt.0.0d0.and.dy.lt.0.0d0) then
c          print*, '  dx>0         dy < 0'
            gi = datan( dx /dy ) + pi 
       endif
       if (dx.lt.0.0d0.and.dy.lt.0.0d0) then
c          print*, '  dx < 0 dy < 0'
             gi = datan( dx /dy ) + pi
       endif
c
       if (dx.lt.0.0d0.and.dy.gt.0.0d0) then
c          print*, '  dx  <  0         dy  > 0'
             gi = datan( dx /dy ) + 2.0d0*pi
       endif
c
       gisement = gi  * 180.0d0/pi
ccc    print*, 'gi =', gi
ccc    print*, 'arctan2x y  =', datan2(dx,dy)*rtod 
ccc    print*, 'arctanx/y =', datan(dx/dy)*rtod 
       return

       end






