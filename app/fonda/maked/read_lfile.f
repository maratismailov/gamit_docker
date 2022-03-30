      subroutine read_lfile(ifil1,ifil2,ifil3,nnet,netfmt,rtime)
c
c     read GAMIT l-file and create FONDA format site-table file
c
c     coordinate system : (comode)
c         1. geocentric spherical coordinate (lat,lon,radius)
c         2. geodetic (ellipsoidal) coordinate (lat,lon,height)
c         3. geocentric Cartesian coordinate (x,y,z)
c         4. topocentric coordinate (x(east),y(north),z(height))
c         Usually, we name the geocentric coordinate as 'global' frame.
c         Meanwhile, the topocentric coordinate is named as 'local' frame.
c
c     ifil1 :  input net coordinate file
c     ifil2 :  output net coordinate file
c     ifil3 :  output net coordinate mapping file
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)

      character*8 name
      character*80 netfmt
      character*12 name12,desc
      character*1 id1,id2
      logical     default
      integer ifil1,ifil2,ifil3,nnet,isit,i,ld,lm,md,mm,j
c     
      write (ifil3,'(a)') ' Network distribution'
      ve = 0.0d0
      vn = 0.0d0
      vu = 0.0d0
      ht = 0.0d0
      sigx = 0.0d0
      sigy = 0.0d0
      sigz = 0.0d0
      isit = 0
c
c     default format
      default = .false.
      if (netfmt(1:1).eq.'*') default = .true.
      if (default) netfmt(1:58) =
     . '(a4,1x,a12,a1,2(i2,1x),f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'
c     get reference epoch
      read (ifil1,'(6x,f8.3)') rnet_time
c
c     write FONDA style network file header
      call net_head(ifil2)
      do 20 i = 1,10000
c        shift format line by line
         if (default) then
            read (ifil1,fmt=netfmt,err=20,end=50) 
     .      name(1:4),name12,id1,ld,lm,sla,id2,md,mm,slo,ht
         else
            read (ifil1,fmt=netfmt,err=20,end=50) 
     .      name,name12,id1,ld,lm,sla,id2,md,mm,slo,ht
         endif
         s1 = dabs(dble(ld))+dble(lm)/60.0d0+sla/3600.0d0
         s2 = dabs(dble(md))+dble(mm)/60.0d0+slo/3600.0d0
         if (id1.eq.'S') s1 = -s1
         if (id2.eq.'W') s2 = -s2
c     fill 8 character site name with '_GPS'
         if (default) name(5:8) = '_GPS'
         j = remedy_space(name12,desc,12,'_',2)
c     count number of sites
         isit = isit+1
         write (ifil2,30) name,desc,
     .         id1,ld,lm,sla,id2,md,mm,slo,ht,ve,vn,vu,
     .         rnet_time,sigx,sigy,sigz
         write (ifil3,40) s2,s1,ve,vn,name
 20   continue
c
c     update nnet
      nnet = isit
c
 30   format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
c    .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
     .   a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
c
 50   continue
      return
      end
