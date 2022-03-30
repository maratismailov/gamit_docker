      program compare_coor
c
c     compare two coordinate sets
c
c     options:
c        1 = compare all 3d
c        2 = compare all 2d
c        3 = compare 3d coordinate only
c        4 = compare 2d coordinate only
c        5 = compare 3d velocity only
c        6 = compare 2d velocity only
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
      integer ios,lcmd,ioptn,len,iyr,imon,iday,ihr,imn,isec
      character*128 line
      integer nblen
      integer*4 iclarg
      logical same_lon,same_lat
c
c     check if the runstring is complete.
      len = iclarg(1,infile)
      l2 = iclarg(2,outfile)
      l1 = iclarg(3,modfile)
      if (l1.le.0.or.l2.le.0.or.len.le.0) then
         print*,' Runstring: compare_coor <infile> <outfile>',
     .      ' <modfile> <option>'
         print*,'  option: 1=3-D coor+velo   2=2-D coor+velo'
         print*,'          3=3-D coor        4=2-D coor'
         print*,'          5=3-D velo        6=2-D velo'
         stop
      endif
      i = iclarg(4,tmpn)
      if (i.le.0) then
         ioptn = 1
      else
         read (tmpn,*) ioptn
         if (ioptn.le.0.or.ioptn.gt.6) ioptn = 1
      endif
c
c     open file
      open (11,file=infile,status='old',err=1000)
      open (12,file=outfile,status='unknown',err=1000)
      open (13,file=modfile,status='old',err=1000)
      print*,' Begin to compare the coordinate files ...'
c
c     set up BIG DO loop to read until end of file...
      rewind (11)
      do 20 i = 1,1000
         rewind (13)
c        read site coordinte and velosity
         ios = 0
         read (11,'(a)',iostat=ios,end=100,err=20) line
         lcmd = nblen(line)
         if (lcmd.le.0) goto 20
c        copy comment lines
         if (line(1:1).ne.' ') then
            write (12,'(a)') line(1:lcmd)
            if (line(3:7).eq.'histo') then
               call getdat(iyr,imon,iday)
               call gettim(ihr,imn,isec,ios)
               write (12,'(a16,a,a5,a,a8,i4,2(a1,i2),3x,i2,2(a1,i2))')
     .            '* difference between: ',infile(1:len),' and ',
     .            modfile(1:l1),'  time: ',
     .            iyr,'/',imon,'/',iday,ihr,':',imn,':',isec
            endif
            goto 20
         endif
         read (line(1:lcmd),14,iostat=ios,err=20) 
     .      sit_nam,fullnam,sym1,id1,im1,sec1,
     .      sym2,id2,im2,sec2,arad,ve,vn,vu,
     .      site_rtime,sigx,sigy,sigz
 14         format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
     .              a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
         same_lon = .true.
         same_lat = .true.
c        check if the name in the modified network file
         do 30 is = 1,1000
            ios = 0
c           read site coordinte and velosity
            read (13,'(a)',iostat=ios,end=40,err=30) line
            lcmd = nblen(line)
c           skip comment lines
            if (line(1:1).ne.' '.or.lcmd.le.1) goto 30

            read (line(1:lcmd),24,iostat=ios,err=30) 
     .      tmpn,sym3,l1,l2,s1,sym4,l3,l4,s2,brad,ve2,vn2,vu2,
     .      site_rtim2,sigx2,sigy2,sigz2
 24         format (1x,a8,14x,a1,2(i2,1x),f8.5,1x,
     .              a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)

         if (tmpn .eq. sit_nam) then
            if (ioptn.gt.4) goto 25
c
c           using uncorrelated time
            delt = site_rtim2-site_rtime
            if (dabs(delt).lt.5.0d2) site_rtime = site_rtim2
c
c           different longitude definition
            if ((sym2.eq.'W'.and.sym4.eq.'E').or.
     .      (sym2.eq.'E'.and.sym4.eq.'W')) then
               id2 = id2-359+l3
               im2 = im2-59+l4
               sec2 = sec2-60.0d0+s2
            else 
               id2 = id2-l3
               im2 = im2-l4
               sec2 = sec2-s2
            endif
            if (im2.lt.0.and.sec2.gt.0.0d0) then
               im2 = im2+1
               sec2 = 60.0d0-sec2
               same_lon = .false.
            endif
            if (sec2.lt.0.0d0) then
               if (im2.le.0) then
                  same_lon = .false.
                  sec2 = -sec2
               else
                  im2 = im2-1
                  sec2 = 60.0d0+sec2
               endif
            endif
c
c           different laitude definition simply using modified one
            if (sym1.eq.sym3) then
               id1 = id1-l1
               im1 = im1-l2
               sec1 = sec1-s1
            else
               sym1 = sym3
               id1 = 89-id1-l1
               im1 = 60-im1-l2
               sec1 = 60.0d0-sec1-s1
            endif
            if (im1.lt.0.and.sec1.gt.0.0d0) then
               im1 = im1+1
               sec1 = 60.0d0-sec1
               same_lat = .false.
            endif
            if (sec1.lt.0.0d0) then
               if (im1.le.0) then
                  same_lat = .false.
                  sec1 = -sec1
               else
                  im1 = im1-1
                  sec1 = 60.0d0+sec1
               endif
            endif
            if (.not.same_lon) then
               if (sym2.eq.'W') then
                  sym2 = 'E'
               else
                  sym2 = 'W'
               endif
            endif
            if (.not.same_lat) then
               if (sym1.eq.'N') then
                  sym1 = 'S'
               else
                  sym1 = 'N'
               endif
            endif
            sigx = sigx2
            sigy = sigy2
            if (ioptn.eq.1.or.ioptn.eq.3) arad = arad-brad
            if (ioptn.eq.1.or.ioptn.eq.3) sigz = sigz2
            if (ioptn.eq.3.or.ioptn.eq.4) goto 35
      
 25         ve = ve-ve2
            vn = vn-vn2
            if (ioptn.eq.1.or.ioptn.eq.5) vu = vu-vu2
c
 35         write (12,14)
     .      sit_nam,fullnam,sym1,id1,im1,sec1,sym2,id2,im2,sec2,
     .      arad,ve,vn,vu,site_rtime,sigx,sigy,sigz
            goto 20
         endif

 30      continue
 40      write (12,14)
     .      sit_nam,fullnam,sym1,id1,im1,sec1,sym2,id2,im2,sec2,
     .      arad,ve,vn,vu,site_rtime,sigx,sigy,sigz

 20   continue  

      goto 100

 1000 print *,' Can not open the file :'
      stop 'in NET_UPDATE. '

 100  continue
      close (11)
      close (12)
      close (13)
c
      stop
      end
