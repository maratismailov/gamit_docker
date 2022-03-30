      subroutine read_free_net(ifil1,ifil2,frame)
c
c     read apriori file with free format and create FONDA format site-table file
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
c        vx,vy,vz : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*8 name,stnm1
      character*6 frame
      character*128 line
      character*1 sym1,sym2
      integer id1,id2,im1,im2,ifil1,ifil2,ifil3,i,i1,k,len
      integer nblen
      dimension ax1(10)
c     
      write (ifil3,'(a)') ' Network distribution'
c
c     default format: free format (geocentric coordinate)
c      if (netfmt(1:1).eq.'*')
c     . netfmt(1:40) = '(a1,f7.3,1x,f7.3,6f8.1,f7.3,3f8.1,2x,a8)'
c
c     write FONDA style network file header
      call net_head(ifil2)
      call geotab(frame,1,semi,finv,tx,ty,tz)
c     screen display to avoid sleep
      print*,'   Processing  <<<  ',netfil
      stnm1 = 'glbk_apr'
      nsit = 0
      sx = 0.0
      sy = 0.0
      sz = 0.0
      do 20 i = 1,10000
c        shift format line by line
         read (ifil1,'(a)',end=50) line
c        skip comment line
         if (line(1:1).ne.' ') goto 20
         i1 = nblen(line)
         if (i1.lt.8) goto 20
c        free format reading does not allow name with number
         read(line(1:10),'(1x,a8)') name
         read (line(10:i1),*,err=20) 
     .        (ax1(k),k=1,7)
c        currently, only accept geocentric format
         bx = ax1(1)
         by = ax1(2)
         bz = ax1(3)
         vxt = ax1(4)
         vyt = ax1(5)
         vzt = ax1(6)
         time1 = ax1(7)
         tmp = bx**2+by**2+bz**2
         if (tmp.lt.1.0d13) goto 20
c        check site name list
c         if (ilst.gt.0) then
c            len = nblen(name)
c            k1 = match_name(ilst,len,st_list,name)
c         endif
         nsit = nsit+1
         len = nblen(name)
         if (len.lt.8) then
            do k = len+1,8
               name(k:k) = '_'
            enddo
         endif
c        from geocentric to geodetic frame
         call geoxyz(semi,finv,tx,ty,tz,alat,alon,ht,
     .      bx,by,bz,2,heit)
         call rtodms(alat,id1,im1,sela,1)
         sym1 = 'N'
         if (alat.lt.0.0d0) then
            sym1 = 'S'
            id1 = iabs(id1)
            im1 = iabs(im1)
            sela = dabs(sela)
         endif
         call rtodms(alon,id2,im2,selo,1)
c        use western longitude (follow the tradition)
         sym2 = 'W'
         if (alon.gt.0.0d0) then
            id2 = 359-id2
            im2 = 59-im2
            selo = 60.0d0-selo
            if (selo.ge.60.0d0) then
               selo = selo-60.0d0
               im2 = im2+1
            endif
            if (im2.ge.60) then
               im2 = im2-60
               id2 = id2+1
            endif
         endif
         call sph_ca(alat,alon,vxt,vyt,vzt,v1,v2,v3,2)
         write (ifil2,30) name,stnm1,
     .      sym1,id1,im1,sela,sym2,id2,im2,selo,ht,v1,v2,v3,
     .      time1,sx,sy,sz
c         if (iomode(5).gt.0.and.k1.gt.0) 
c     .      write (imfil,40) alon,alat,ve1,vn1,st1
 
 20   continue
c
 30   format (1x,a8,2x,a8,4x,a1,2(i2,1x),f8.5,1x,
c    .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
     .   a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
 
 50   continue
      return
      end

