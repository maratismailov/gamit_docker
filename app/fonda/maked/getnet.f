      subroutine getnet(ifile)
c
c     get network site coordinates
c
c     coordinate system : (comode)
c         1. geocentric spherical coordinate (lat,lon,radius)
c         2. geodetic (ellipsoidal) coordinate (lat,lon,height)
c         3. geocentric Cartesian coordinate (x,y,z)
c         4. topocentric coordinate (x(east),y(north),z(height))
c         Usually, we name the geocentric coordinate as 'global' frame.
c         Meanwhile, the topocentric coordinate is named as 'local' frame.
c
c     network file style : (net_style)
c         1. GAMIT l-file 
c         2. GLOBK 
c         3. FONDA 
c         4. BLUE BOOK
c         5. USGS statab.lis
c         6  UCSD
c         7. IPGP 
c         8. (waiting list)
c
c     velocity mode : (vmode)
c         1. directly read geocentric velocity values
c         2. directly read geodetic velocity values
c         3. group movement (geodetic)
c         4. rotation (geodetic)
c         5. dislocation
c         6. ......
c
c     unit:
c         x, y, z : meter
c        vx,vy,vz : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
      character*1 ns,ew
      integer i,ifile,isit,ios,i2,l1,l2,l3,l4,nblen
      character*128 line

      print*,' net_style =',net_style
c
      rewind (ifile)
c     
c     get site name,coordinate and velocity
      isit = 0
      do 10 i = 1,1000
      ios = 0
c     read site coordinte and velosity
      read (ifile,'(a)',iostat=ios,end=50,err=10) line
c     skip comment line
      i2 = nblen(line)
      if (net_style.eq.3.and.(line(1:1).ne.' '.or.i2.le.1)) goto 10
      if (net_style.ne.3.and.i.eq.1) goto 10
      isit = isit+1
      jtoi(isit) = 0

      if (comode.eq.1) then
      if (net_style.eq.3) then
         read (line(1:i2),15) sname(isit),ns,l1,l2,s1,ew,l3,l4,s2,r1,
     .      vxt,vyt,vzt,site_rtime,sigx,sigy,sigz
      else
         read (line(1:i2),14) sname(isit)(1:4),ns,l1,l2,s1,ew,l3,l4,
     .      s2,r1
         sname(isit)(5:8) = '____'
      endif
      sla(isit) = dble(l1)+dble(l2)/60.0d0+s1/3600.0d0
c     southern latitude
      if (ns.eq.'S'.or.ns.eq.'s') sla(isit) = -sla(isit)
      slo(isit) = dble(l3)+dble(l4)/60.0d0+s2/3600.0d0
c     westward longitude
      if (ew.eq.'W'.or.ew.eq.'w') slo(isit) = 360.0d0-slo(isit)
      alat = sla(isit)*dtor
      alon = slo(isit)*dtor
      arad = r1
c     transform to Cartesian coordinate
      call sphxyz(alat,alon,arad,x1,y1,z1,1)
      x(isit) = x1
      y(isit) = y1
      z(isit) = z1
c     transform to gepodetic coordinate
      call geoxyz(radius,finv,tx,ty,tz,alat,alon,dh,x1,y1,z1,2,he)
      sla(isit) = alat
      slo(isit) = alon
      srad(isit) = dh
      endif
      if (comode.eq.2) then
      if (net_style.eq.3) then
         read (line(1:i2),15) sname(isit),ns,l1,l2,s1,ew,l3,l4,s2,r1,
     .      vxt,vyt,vzt,site_rtime,sigx,sigy,sigz
      else
         read (line(1:i2),14) sname(isit)(1:4),ns,l1,l2,s1,ew,l3,l4,
     .      s2,r1
         sname(isit)(5:8) = '____'
      endif
c     write (6,15) sname(isit),ns,l1,l2,s1,ew,l3,l4,s2,r1
      sla(isit) = dble(l1)+dble(l2)/60.0d0+s1/3600.0d0
c     southern latitude
      if (ns.eq.'S'.or.ns.eq.'s') sla(isit) = -sla(isit)
      slo(isit) = dble(l3)+dble(l4)/60.0d0+s2/3600.0d0
c     westward longitude
      if (ew.eq.'W'.or.ew.eq.'w') slo(isit) = 360.0d0-slo(isit)
      alat = sla(isit)*dtor
      alon = slo(isit)*dtor
      dh = r1
c     transform to Cartesian coordinate
      call geoxyz(radius,finv,tx,ty,tz,alat,alon,dh,x1,y1,z1,1,he)
      sla(isit) = alat
      slo(isit) = alon
      srad(isit) = dh
      x(isit) = x1
      y(isit) = y1
      z(isit) = z1
      endif
      if (comode.gt.2) then 
      read (line(1:i2),12) sname(isit),x1,y1,z1
      call geoxyz(radius,finv,tx,ty,tz,alat,alon,hght,x1,y1,z1,2,he)
      x(isit) = x1
      y(isit) = y1
      z(isit) = z1
      sla(isit) = alat
      slo(isit) = alon
      srad(isit) = he
      endif
 10   continue
 12   format (a8,12x,3(1x,f15.4))
 14   format (a4,13x,a1,2(i2,1x),f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)
 15   format (1x,a8,14x,a1,2(i2,1x),f8.5,1x,
     .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
c
 50   continue
c
      return
      end
