      subroutine getpri(iprif,sname2)
c
c     read priori coordinate and velocity values
c
c     coordinate system : (comode)
c         1. geocentric spherical coordinate (lat,lon,radius)
c         2. geodetic (ellipsoidal) coordinate (lat,lon,height)
c         3. geocentric Cartesian coordinate (x,y,z)
c         4. topocentric coordinate (x(east),y(north),z(height))
c         Usually, we name the geocentric coordinate as 'global' frame.
c         Meanwhile, the topocentric coordinate is named as 'local' frame.
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : m/year
c        slat,slon: radian
c        srad     : m
c        ve,vn,vu : m/year
c        time     : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      character*1 sym1,sym2
      character*8 sname2(maxsit)
      character*8 tmpn,tmpn2
      integer i,iprif,is,l1,l2,l3,l4,i2
      integer ios
      character*128 line
      logical found
      integer nfound
      integer nblen

      nfound = 0
      found = .true.
c
c     set up BIG DO loop to read until end of file...
      do 20 i = 1,nsit
      rewind (iprif)
      if (i.gt.1.and..not.found) 
     .   print*,' Missing site: ', sname(i-1)
 
c     have to read all the way to the end of the file.
c     there may be more stations in the apriori file 
c     than we know or care about. KLF
      do 30 is = 1,1000
      ios = 0

c     read site coordinate and velocity
      read (iprif,'(a)',iostat=ios,end=20,err=30) line

c     skip comment line
      i2 = nblen(line)
      if (line(1:1).ne.' '.or.i2.le.1) goto 30
       

c     Read coordinates depending on specified system (comode)
      if (comode.eq.1 .or. comode .eq. 2) then
         if (fcode.eq.1) then
            read (line(1:i2),14,iostat=ios,err=30) 
     .      tmpn,fullnm(i),sym1,l1,l2,s1,sym2,l3,l4,s2,arad,
     .      vxt,vyt,vzt,site_rtime,sigx,sigy,sigz
c14         format (1x,a8,14x,a1,2(i2,1x),f8.5,1x,
c    .              a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
 14         format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
     .          a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
         else if (fcode.eq.2) then 
            read (line(1:i2),24,iostat=ios,err=30) 
     .      tmpn2,fullnm(i),sym1,l1,l2,s1,sym2,l3,l4,s2,arad,
     .      vxt,vyt,vzt,site_rtime,sigx,sigy,sigz
c24         format (1x,a8,14x,a1,2(i2,1x),f8.5,1x,
c    .              a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
 24         format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
     .          a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
c           print*, '   '
c           print*, 'IN GETPRI'
c           print*, 'fcode,comode =',fcode,comode,
c    .      'char_8,geodetic'
c           write(*,24) tmpn2,sym1,l1,l2,s1,sym2,l3,l4,s2,arad,
c    .      vxt,vyt,vzt,site_rtime,sigx,sigy,sigz
         else if (fcode.eq.3) then 
            read (line(1:i2),34,iostat=ios,err=30) 
     .      tmpn,tmpn2,sym1,l1,l2,s1,sym2,l3,l4,s2,arad,vxt,vyt,vzt,
     .      site_rtime,sigx,sigy,sigz
c34         format (1x,a8,1x,a8,5x,a1,2(i2,1x),f8.5,1x,
c    .              a1,i3,1x,i2,1x,f8.5,f13.4,3f9.5,f9.3,3f8.4)
 34         format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
     .          a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
            fullnm(i) = tmpn2
         else
            print *,'GETPRI: unknown fcode = ',fcode
            stop 'in GETPRI -1'
         endif
      else
         read (line(1:i2),10,iostat=ios,err=30) 
     .      tmpn,xt,yt,zt,vxt,vyt,vzt,site_rtime,sigx,sigy,sigz
 10      format (1x,a8,5x,2f14.8,f13.4,3f9.5,f9.3,3f8.4)
      endif


c     Convert velocity Units to m/yr if needed
c      (specified in the driver file - 'Velocity Units:')
      if (iomode(18).ne.1) then
          vxt=vxt/1.0d3
          vyt=vyt/1.0d3
          vzt=vzt/1.0d3
      endif

      found = .false.
      if (fcode.eq.1 .and. tmpn  .eq. sname2(i)) found = .true.
      if (fcode.eq.2 .and. tmpn2 .eq. sname(i))  found = .true.
      if (fcode.eq.3 .and. 
     .    tmpn.eq.sname(i) .and. 
     .    tmpn2(1:4).eq.sname2(i)(1:4)) found = .true.
c     print*, 'found',found
c
      if (found) then
         nfound = nfound + 1

c        Check and see if site_rtime is reasonable
         if (site_rtime.lt.1.0d3.or.site_rtime.gt.3.0d3)
     .      site_rtime = rtime
         fact = (rtime-site_rtime)

      if (comode.eq.1) then
         slat(i) = dble(l1)+dble(l2)/60.0d0+s1/3600.0d0
         slon(i) = dble(l3)+dble(l4)/60.0d0+s2/3600.0d0
c        westward longitude
         if (sym2.eq.'W') slon(i) = 360.0d0-slon(i)
c        southern latitutde
         if (sym1.eq.'S') slat(i) = -slat(i)
         alat = slat(i)*dtor
         alon = slon(i)*dtor
         arad = arad
c        transform to Cartesian coordinate
         call sphxyz(alat,alon,arad,x1,y1,z1,1)
         call sph_ca(alat,alon,vxt,vyt,vzt,v1,v2,v3,1)

c        Bring coords up to epoch rtime
         x1 = x1+v1*fact
         y1 = y1+v2*fact
         z1 = z1+v3*fact
         x(i) = x1
         y(i) = y1
         z(i) = z1
         vx(i) = v1
         vy(i) = v2
         vz(i) = v3
c        transform to geodetic coordinate
         call geoxyz
     .   (radius,finv,tx,ty,tz,alat,alon,dh,x1,y1,z1,2,he)
         slat(i) = alat
         slon(i) = alon
         srad(i) = dh
         call sph_ca(alat,alon,v1,v2,v3,vxt,vyt,vzt,2)
         ve(i) = vxt
         vn(i) = vyt
         vu(i) = vzt

      else if (comode.eq.2) then
         slat(i) = dble(l1)+dble(l2)/60.0d0+s1/3600.0d0
         slon(i) = dble(l3)+dble(l4)/60.0d0+s2/3600.0d0
c        print*, 'slati,i ', slat(i),i
c        print*, 'sloni,i ', slon(i),i 
c        southern latitutde
         if (sym1.eq.'S') slat(i) = -slat(i)
c        westward longitude
         if (sym2.eq.'W') slon(i) = 360.0d0-slon(i)
         alat = slat(i)*dtor
         alon = slon(i)*dtor
         dh = arad
c        print*, 'dh = arad', arad

c        transform to Cartesian coordinate
         call geoxyz(radius,finv,tx,ty,tz,alat,alon,dh,x1,y1,z1,1,he)
c        print*, 'radius,finv,tx,ty,tz,alat,alon,dh,x1,y1,z1,1,he'
c        print*, radius,finv,tx,ty,tz,alat,alon,dh,x1,y1,z1,1,he

         ve(i) = vxt
         vn(i) = vyt
         vu(i) = vzt

c        Convert spherical to Cartesian
         call sph_ca(alat,alon,vxt,vyt,vzt,v1,v2,v3,1)

c        Bring coords up to epoch rtime
         x1 = x1+v1*fact
         y1 = y1+v2*fact
         z1 = z1+v3*fact
         x(i) = x1
         y(i) = y1
         z(i) = z1
         vx(i) = v1
         vy(i) = v2
         vz(i) = v3
c        transform back to geodetic coordinate
         call geoxyz(radius,finv,tx,ty,tz,alat,alon,dh,x1,y1,z1,2,he)
         slat(i) = alat
         slon(i) = alon
         srad(i) = dh

      else if (comode.eq.3 .or. comode .eq. 4) then
c        Bring coords up to epoch rtime
         xt = xt+vxt*fact
         yt = yt+vyt*fact
         zt = zt+vzt*fact
         x(i) = xt
         y(i) = yt
         z(i) = zt
         vx(i) = vxt
         vy(i) = vyt
         vz(i) = vzt
         dh = 0.0d0
         call geoxyz(radius,finv,tx,ty,tz,alat,alon,dh,xt,yt,zt,2,he)
         slat(i) = alat
         slon(i) = alon
         srad(i) = dh
         call sph_ca(alat,alon,vxt,vyt,vzt,v1,v2,v3,2)
         ve(i) = v1
         vn(i) = v2
         vu(i) = v3
      else
         print *,'GETPRI: unknown comode = ',comode
      endif
      if (found) goto 20
      endif

 30   continue

      if (ios .gt. 0) then
         print *,'GETPRI: file read error.  Format? Fcode = ',fcode
         print *,line
         call ferror (ios,6)
      endif

 20   continue  

      if (nfound .ne. nsit) then
         print *,'GETPRI: missing a station. nfound, nsit = ',
     .            nfound,nsit
         stop 'in GETPRI 1'
      endif
c     print*, 'x= ',x
c     print*, 'y= ',y
c     print*, 'z= ',z

c
      return
      end
