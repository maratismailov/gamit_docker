      program calc_bsl_azm
c
c     calculate baseline length and azimuth between two sites
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
      character*1 sym3,sym4
      character*8 tmpn,sitnm1,sitnm2
      character*12 fullnam
      character*64 infile,outfile,netfile
      integer i,l1,l2,l3,l4,ioptn,j
      integer ios,lcmd,ifil1,ifil2,i1
      integer iclarg,mchkey,lift_arg
      character*128 line
      integer nblen
      logical same1,same2,d_t
      data pi/3.14159265d0/
c
c     check if the runstring is complete.
      l1 = iclarg(1,netfile)
      l2 = iclarg(2,infile)
      l3 = iclarg(3,outfile)
      if (l1.le.0.or.l2.le.0.or.l3.le.0) then
         print*,' calc_bsl_azm: calculate baseling length and'
         print*,'               azimuth between two sites'
         print*,' Runstring: calc_bsl_azm <netfile> <infile>',
     .      '<outfile> <rtime> <option>'
         print*,' option: 1=length+azimuth 2=dx, dy, dz '
         print*,' option: 3=de, dn, du     4=(undecided)'
         stop
      endif
c
c     set up reference time
      i = iclarg(4,fullnam) 
      if (i.le.0) then
         rtime = 0.0d0
      else
         read (fullnam,*) rtime
         d_t = .true.
         if (rtime.le.1.0d3.or.rtime.gt.3.0d3) d_t = .false.
      endif
c
c     option
      i = iclarg(5,fullnam) 
      if (i.le.0) then
         ioptn = 1
      else
         read (fullnam,*) ioptn
         if (ioptn.lt.1.or.ioptn.gt.3) ioptn = 1
      endif
c
c     open file
      open (10,file=netfile,status='old',err=1000)
      if (l2.le.2) then
         ifil1 = 5
      else
         ifil1 = 11
         open (ifil1,file=infile,status='old',err=1000)
      endif
      if (l3.le.2) then
         ifil2 = 6
      else
         ifil2 = 12
         open (ifil2,file=outfile,status='unknown',err=1000)
      endif
c
      rtod = 180.0/pi
      do 20 i = 1,1000
c        read site names
         ios = 0
         if (ifil1.eq.5) print*,' input names of sit1 and sit2:'
         read (ifil1,'(a)',iostat=ios,end=100,err=20) line
         lcmd = nblen(line)
         if (lcmd.le.0) goto 20
c        skip comment lines or exit
         if (ifil1.ne.5.and.line(1:1).ne.' ') then
            i1 = mchkey(line,'exit',lcmd,4)
            if (i1.gt.0) goto 100
            i1 = mchkey(line,'quit',lcmd,4)
            if (i1.gt.0) goto 100
            write (ifil2,'(a)') line(1:lcmd)
            goto 20
         endif
c        read two site coordinates
         j = lift_arg(line,sitnm1,1)
         j = lift_arg(line,sitnm2,2)
         if (ifil1.eq.5) then
            if (sitnm1.eq.'quit'.or.sitnm1.eq.'exit') goto 100
         endif
c
         same1 = .false.
         same2 = .false.
         rewind (10)
         do 40 j = 1,1000
            read (10,'(a)',iostat=ios,end=20,err=40) line
            lcmd = nblen(line)
            if (lcmd.le.0) goto 40
            if (line(1:1).ne.' ') goto 40
c           read site names fromnetwork file
            read (line(1:lcmd),14,iostat=ios,err=40)
     .      tmpn,sym3,l1,l2,s1,sym4,l3,l4,s2,brad,ve,vn,vu,
     .      site_rtime,sigx,sigy,sigz
 14         format (1x,a8,14x,a1,2(i2,1x),f8.5,1x,
     .              a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
            if (tmpn.eq.sitnm1) then
               same1 = .true.
               slat1 = dble(l1)+dble(l2)/60.0+s1/3.6d3
               if (sym3.eq.'S') slat1 = -slat1
               slon1 = dble(l3)+dble(l4)/60.0+s2/3.6d3
               if (sym4.eq.'W') slon1 = 360.0d0-slon1
               srad1 = brad
               dt1 = rtime-site_rtime
               ve1 = ve
               vn1 = vn
               vu1 = vu
            endif
            if (tmpn.eq.sitnm2) then
               same2 = .true.
               slat2 = dble(l1)+dble(l2)/60.0+s1/3.6d3
               if (sym3.eq.'S') slat2 = -slat2
               slon2 = dble(l3)+dble(l4)/60.0+s2/3.6d3
               if (sym4.eq.'W') slon2 = 360.0d0-slon2
               srad2 = brad
               dt2 = rtime-site_rtime
               ve2 = ve
               vn2 = vn
               vu2 = vu
            endif
            if (same1.and.same2) goto 50
 40      continue
         if (.not.same1.or..not.same2) goto 20
 50      slat1 = slat1/rtod
         slon1 = slon1/rtod
         slat2 = slat2/rtod
         slon2 = slon2/rtod
         call geotab('WGS84',1,radius,finv,tx,ty,tz)
cc         if (d_t) then
cc            sft_e = ve2*dt2-ve1*dt1
cc            sft_n = vn2*dt2-vn1*dt1
cc            sft_u = vu2*dt2-vu1*dt1
cc            call sphxyz(slat1,slon1,radius,sft_e,sft_n,sft_u,1)
cc         endif
         if (ioptn.eq.1) then
         call azmlen(slat1,slon1,srad1,slat2,slon2,srad2,azm,alen,1)
         azm = azm*rtod
         write (ifil2,16) sitnm1,sitnm2,azm,alen
         endif
         if (ioptn.ge.2) then
         call geoxyz(radius,finv,tx,ty,tz,slat1,slon1,srad1,
     .            x1,y1,z1,1,ht)
         call geoxyz(radius,finv,tx,ty,tz,slat2,slon2,srad2,
     .            x2,y2,z2,1,ht)
         xd = x2-x1
         yd = y2-y1
         zd = z2-z1
         if (ioptn.eq.3) then
         s1 = dsin(slon1)
         c1 = dcos(slon1)
         s2 = dsin(slat1)
         c2 = dcos(slat1)
         dx = -xd*s1+yd*c1
         dy = -xd*c1*s2-yd*s1*s2+zd*c2
         dz = xd*c1*c2+yd*s1*c2+zd*s2
         endif
         if (ioptn.eq.2) write (ifil2,17) sitnm1,sitnm2,xd,yd,zd
         if (ioptn.eq.3) write (ifil2,18) sitnm1,sitnm2,dx,dy,dz
         endif
 20   continue
 16   format (2x,'from ',a8,' to ',a8,':  azm =',f7.2,' len =',f12.3)
 17   format (2x,'from ',a8,' to ',a8,':  dx =',f11.2,' dy =',f11.2
     .        ' dz =',f11.2)
 18   format (2x,'from ',a8,' to ',a8,':  de =',f11.2,' dn =',f11.2
     .        ' du =',f11.2)
 100  continue
 1000 continue
      stop
      end
 
c-----------------------------------------------------------------------
      subroutine azmlen(sla1,slo1,sra1,sla2,slo2,sra2,azm,alen,mode)
c
c     mode = 1:  azimuth observable in local astronomical frame
c     mode = 2:  azimuth observable in geodetic frame
c
      implicit real*8(a-h,o-z)
c
      integer mode
      real*8 pi
      data pi/3.14159265d0/

c     calculate azimuth from sit1 to sit2. unit = radian
      call geotab('WGS84',1,radius,finv,tx,ty,tz)
c
c     first, get vector components in Cartesian coordinate
      call geoxyz(radius,finv,tx,ty,tz,sla1,slo1,sra1,
     .            x1,y1,z1,1,ht)
      call geoxyz(radius,finv,tx,ty,tz,sla2,slo2,sra2,
     .            x2,y2,z2,1,ht)
      xd = x2-x1
      yd = y2-y1
      zd = z2-z1
c
      s1 = dsin(slo1)
      c1 = dcos(slo1)
      s2 = dsin(sla1)
      c2 = dcos(sla1)
      dx = -xd*s1+yd*c1
      dy = -xd*c1*s2-yd*s1*s2+zd*c2
c
      if (mode.eq.2) then
         c1 = dcos(sla1)
         del = slo2-slo1
         dx = dsin(del)
         c1 = dcos(del)
         s1 = dsin(sla1)
         e2 = 2.0d0/finv-1.0d0/finv**2
         c3 = (1.0d0-e2)*dtan(sla2)+
     .      e2*s1/dcos(sla2)
         c2 = dtan(sla2)
         c3 = c3*dsqrt((1.d0-e2*dsin(sla2)**2)/(1.d0-e2*s1**2))
         dy = dcos(sla1)*c2-s1*c1

         azm = datan2(dx, dy)
         if (azm.lt.0.d0) azm = azm+pi*2.0d0
         alen = dsqrt(dx**2+dy**2+dz**2)
         return
      endif

      azm = datan2(dx, dy)
      dx = x2-x1
      dy = y2-y1
      dz = z2-z1
      if (azm.lt.0.d0) azm = azm+pi*2.0d0
      alen = dsqrt(dx**2+dy**2+dz**2)

      return
      end
