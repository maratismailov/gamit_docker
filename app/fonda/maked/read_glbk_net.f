      subroutine read_glbk_net(ifil1,ifil2,ifil3,nnet,netfmt,rtime)
c
c     read GLOBK output file and create FONDA format site-table file
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
c         x, y, z : m
c        vx,vy,vz : mm/year
c        time     : year
c
      implicit real*8(a-h,o-z)

      character*8 name,desc
      character*80 netfmt
      character*1 id1,id2
      integer ifil1,ifil2,ifil3,nnet,i,k,ld,lm,md,mm
      dimension ax1(10)
c     
      write (ifil3,'(a)') ' Network distribution'
      ve = 0.0d0
      vn = 0.0d0
      vu = 0.0d0
      ht = 0.0d0
      sigx = 0.0d0
      sigy = 0.0d0
      sigz = 0.0d0
c
c     default format
      if (netfmt(1:1).eq.'*')
     . netfmt(1:40) = '(a1,f7.3,1x,f7.3,6f8.1,f7.3,3f8.1,2x,a8)'
c
      print*,' Creating coordinate file ...'
c     write FONDA style network file header
      call net_head(ifil2)
      do 20 i = 1,10000
c        shift format line by line
         read (ifil1,fmt=netfmt,err=20,end=50) 
     .        id1,s2,s1,(ax1(k),k=1,10),name
c        skip comment line
         if (id1.eq.'*') goto 20
         id1 = 'N'
         if (s1.lt.0.0d0) then
            id1 = 'S'
            s1 = -s1
         endif
         id2 = 'E'
         if (s2.lt.0.0d0) then
            id2 = 'W'
            s2 = -s2
         endif
c
         ld = int(s1)
         lm = int((s1-ld)*60.0d0)
         sla = (s1-ld-dble(lm)/60.0d0)*3.6d3
         md = int(s2)
         mm = int((s2-md)*60.0d0)
         slo = (s2-md-dble(mm)/60.0d0)*3.6d3
         write (desc,'(i3)') i
         desc(4:8) = '_glbk'
         write (ifil2,30) name,desc,
     .         id1,ld,lm,sla,id2,md,mm,slo,ht,ve,vn,vu,
     .         rtime,sigx,sigy,sigz
         write (ifil3,40) s2,s1,ve,vn,name
 20   continue
c
c10   format (a1,f7.3,1x,f7.3,6f7.1,f8.3,3f7.1,1x,a8)
c10   format (a1,f7.3,1x,f7.3,6f8.1,f7.3,3f8.1,2x,a8)
c10   format (a1,f8.3,1x,f8.3,1x,6f7.1,f7.3,2x,3f7.1,1x,a8)
 30   format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
c    .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
     .   a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
c
 50   continue
      return
      end
