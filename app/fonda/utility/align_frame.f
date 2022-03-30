      program align_frame
c
c     align reference frame using 7-parameter transformation
c
c     options:
c        1 = 3d translation+rotation, scale
c        2 = update all 2d
c        3 = update 3d coordinate only
c        4 = update 2d coordinate only
c        5 = update 3d velocity only
c        6 = update 2d velocity only
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : mm/year
c        slat,slon: radian
c        srad     : m
c        ve,vn,vu : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)
c
      character*1 sym1,sym2,sym3,sym4
      character*8 sit_nam, tmpn
      character*12 fullnam
      character*64 infile,outfile,modfile
      integer i,id1,im1,id2,im2,is,l1,l2,l3,l4
      integer ios,lcmd,ioptn,iyr,imon,iday,ihr,imn,isec
      integer isit
      integer*4 iclarg
      character*128 line
      integer nblen
      dimension alat(100),alon(100),blat(100),blon(100)
c
c     check if the runstring is complete.
      l1 = iclarg(1,infile)
      l2 = iclarg(2,outfile)
      l3 = iclarg(3,modfile)
      if (l1.le.0.or.l2.le.0.or.l3.le.0) then
         print*,' Runstring: align_frame <infile> <outfile>',
     .      ' <modfile> <option>'
         print*,'  option: 1=7_para transf     2=3-D transl+rota'
         print*,'          3=2-D transl+rota   4=3-D translation'
         print*,'          5=2-D translation   6=2-D rotation'
         stop
      endif
      i = iclarg(4,tmpn)
      if (i.le.0) then
         ioptn = 3
      else
         read (tmpn,*) ioptn
         if (ioptn.le.0.or.ioptn.gt.6) ioptn = 3
      endif
c
c     open file
      open (11,file=infile,status='old',err=1000)
      open (12,file=outfile,status='unknown',err=1000)
      open (13,file=modfile,status='old',err=1000)
      print*,' Begin to calculate transformation parameter ...'
c
c     get core sites
      rewind (13)
      isit = 0
      dtor = datan(1.0d0)/45.0d0
      a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      b1 = 0.0d0
      b2 = 0.0d0
      do 20 i = 1,100
c        read site coordinte 
         ios = 0
         read (13,'(a)',iostat=ios,end=40,err=20) line
         lcmd = nblen(line)
c        skip comment lines
         if (line(1:1).ne.' '.or.lcmd.le.1) goto 20
         read (line(1:lcmd),24,iostat=ios,err=20) 
     .      sit_nam,sym1,id1,im1,sec1,
     .      sym2,id2,im2,sec2,arad,ve,vn,vu,core_rtime,
     .      sigx,sigy,sigz
 14         format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
     .              a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
         slat = dble(id1)+dble(im1)/60.0d0+sec1/3600.0d0
         slon = dble(id2)+dble(im2)/60.0d0+sec2/3600.0d0
c        westward longitude
         if (sym2.eq.'W') slon = 360.0d0-slon
c        southern latitutde
         if (sym1.eq.'S') slat = -slat
         alat(i) = slat*dtor
         alon(i) = slon*dtor
c        check if the name in the input file
         rewind (11)
         do 30 is = 1,1000
            ios = 0
c           read site coordinte 
            read (11,'(a)',iostat=ios,end=20,err=30) line
            lcmd = nblen(line)
c           skip comment lines
            if (line(1:1).ne.' '.or.lcmd.le.1) goto 30

            read (line(1:lcmd),14,iostat=ios,err=30) 
     .      tmpn,fullnam,sym3,l1,l2,s1,sym4,l3,l4,s2,brad,
     .      ve2,vn2,vu2,site_rtime,sx2,sy2,sz2
 24         format (1x,a8,14x,a1,2(i2,1x),f8.5,1x,
     .              a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)

            if (tmpn .eq. sit_nam) then
               isit = isit+1
               slat = dble(l1)+dble(l2)/60.0d0+s1/3600.0d0
               slon = dble(l3)+dble(l4)/60.0d0+s2/3600.0d0
c              westward longitude
               if (sym4.eq.'W') slon = 360.0d0-slon
c              southern latitutde
               if (sym3.eq.'S') slat = -slat
               blat(isit) = slat*dtor
               blon(isit) = slon*dtor
               alat(isit) = alat(i)
               alon(isit) = alon(i)
               a1 = a1+alat(isit)
               a2 = a2+alon(isit)
               b1 = b1+blat(isit)
               b2 = b2+blon(isit)
               a3 = a3+brad-arad
               goto 20
            endif
 30      continue
 20   continue
c
 40   if (isit.le.1) goto 2000
      c1 = dsqrt(1.0d0/dble(isit))
      x0 = (b2-a2)/dble(isit)
      y0 = (b1-a1)/dble(isit)
      x1 = a2/dble(isit) 
      y1 = a1/dble(isit) 
      c2 = 0.0d0
      do i = 1,isit
         tmp = (blon(i)-alon(i)-x0)*dcos(y1)
         c2 = c2+tmp**2
         c2 = c2+(blat(i)-alat(i)-y0)**2
      enddo
      if (c2.le.0.0d0) then 
         c2 = 0.0d0
      else
         c2 = dsqrt(1.0d0/c2)
      endif
      t1 = (b2-a2)*c1*dcos(y1)
      t2 = (b1-a1)*c1
      t3 = a3*c1
      t4 = 0.0d0
      do i = 1,isit
         t4 = t4-(blat(i)-alat(i)-y0)*(blon(i)-alon(i))
         t4 = t4+(blat(i)-alat(i))*(blon(i)-alon(i)-x0)
      enddo
      t4 = t4*dcos(y1)*c2
      print*,' t1,t2,t3,t4:',t1,t2,t3,t4
c
c     set up BIG DO loop to read until end of file...
      rewind (11)
      do 50 i = 1,1000
c        read site coordinte and velocity
         ios = 0
         read (11,'(a)',iostat=ios,end=100,err=50) line
         lcmd = nblen(line)
         if (lcmd.le.0) goto 50
c        copy comment lines
         if (line(1:1).ne.' ') then
            write (12,'(a)') line(1:lcmd)
            if (line(3:7).eq.'histo') then
               call getdat(iyr,imon,iday)
               call gettim(ihr,imn,isec,ios)
               write (12,'(a15,4(2x,d12.4),a8,i4,2(a1,i2),
     .            3x,i2,2(a1,i2))')
     .            '* aligned with:',t1,t2,t3,t4,'  time: ',
     .            iyr,'/',imon,'/',iday,ihr,':',imn,':',isec
            endif 
            goto 50
         endif
         read (line(1:lcmd),14,iostat=ios,err=50) 
     .      sit_nam,fullnam,sym1,id1,im1,sec1,
     .      sym2,id2,im2,sec2,arad,ve,vn,vu,
     .      site_rtime,sigx,sigy,sigz
         slat = dble(id1)+dble(im1)/60.0d0+sec1/3600.0d0
         slon = dble(id2)+dble(im2)/60.0d0+sec2/3600.0d0
c        westward longitude
         if (sym2.eq.'W') slon = 360.0d0-slon
c        southern latitutde
         if (sym1.eq.'S') slat = -slat
         clat = slat*dtor
         slon = slon*dtor
         slat = clat-(t2*c1+t4*c2*(slon-x1)*dcos(y1))
         slon = slon-(t1*c1-t4*c2*(clat-y1))/dcos(y1)
         arad = arad-t3*c1
c        restore original convention
         if (sym2.eq.'W') slon = datan(1.0d0)*8.0d0-slon
         if (sym1.eq.'S') slat = -slat
         call rtodms(slat,id1,im1,sec1,1)
         call rtodms(slon,id2,im2,sec2,1)
c
         write (12,14)
     .   sit_nam,fullnam,sym1,id1,im1,sec1,
     .   sym2,id2,im2,sec2,arad,ve,vn,vu,
     .   site_rtime,sigx,sigy,sigz

 50   continue  

      goto 100

 1000 print *,' Can not open the file :'
      stop 

 2000 print *,' Too few core sites.'
      stop 

 100  continue
      close (11)
      close (12)
      close (13)
c
      stop
      end
