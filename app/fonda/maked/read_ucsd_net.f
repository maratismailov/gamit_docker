      subroutine read_ucsd_net(ifil1,ifil2,ifil3,nnet,netfmt,rtime)
c
c     read UCSD site-table file and create FONDA format site-table file
c     ( SV6 format).
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
c         x, y, z : km
c        vx,vy,vz : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)

      character*8 name,desc
      character*80 netfmt
      character*1 id1,id2
      integer*4 ierr,ifil1,ifil2,ifil3,nnet,i,ld,lm,md,mm
c      
      write (ifil3,'(a)') ' Network distribution'
      ve = 0.0d0
      vn = 0.0d0
      vu = 0.0d0
      sigx = 0.0d0
      sigy = 0.0d0
      sigz = 0.0d0
c
c     default format
      if (netfmt(1:1).eq.'*')
     . netfmt(1:46) = '(6x,a8,3x,a1,2i3,f9.5,i4,i3,f9.5,3x,a1,f8.2)'
c
c     write FONDA style network file header
      call net_head(ifil2)
      do 20 i = 1,nnet
c        shift format line by line
         read (ifil1,fmt=netfmt,iostat=ierr, end=50)
     .      name,id1,ld,lm,sla,md,mm,slo,id2,ht
         print*,name,id1,ld,lm,sla,md,mm,slo,id2,ht
         s1 = dabs(dble(ld))+dble(lm)/60.0d0+sla/3600.0d0
c        westward longitude
         s2 = dabs(dble(md))+dble(mm)/60.0d0+slo/3600.0d0
         id1 = 'N'
         if (ld.lt.0) then
            id1 = 'S'
            ld = -ld
            s1 = -s1
         endif
         id2 = 'E'
         if (md.lt.0) then
            id2 = 'W'
            md = -md
            s2 = -s2
         endif
         write (desc,'(i3)') i
         desc(4:8) = '_ucsd'
         write (ifil2,30) name,desc,
     .         id1,ld,lm,sla,id2,md,mm,slo,ht,ve,vn,vu,
     .         rtime,sigx,sigy,sigz
         write (ifil3,40) s2,s1,ve,vn,name
         if(i.lt.nnet) read (ifil1,'(2x)')
 20   continue
c
c10   format (6x,a8,3x,a1,2i3,f9.5,i4,i3,f9.5,3x,a1,f8.2)
 30   format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
c    .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
     .   a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
c
 50   continue
      return
      end
