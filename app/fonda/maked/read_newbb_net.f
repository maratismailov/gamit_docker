      subroutine read_newbb_net(ifil1,ifil2,ifil3,nnet,netfmt,
     .                          iesit,sname,rtime)
c
c     read BLUE BOOK file and create FONDA format site-table file
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

      character*16 name,desc
      character*8  sname(200)
      character*80 netfmt
      character*120 line
      character*12  name2
      character*1 id1,id2
      integer   i,j,type,lad,lam,las,lod,lom,los,iht
      integer   ii,remedy_space,iesit(200),unique_nam1
      integer   ifil1,ifil2,ifil3,nnet,k
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
     . netfmt(1:39) = '(7x,i2,1x,i4,a16,13x,2(i3,i2,i7,a1),i6)'
   
c
c     write FONDA style network file header
      call net_head(ifil2)
      nnet = 0
      do 20 i = 1,10000
c        shift format line by line
         read (ifil1,'(a120)', end=50, err=20)
     .     line 
c        pick up the line with site coordinate
         read (line,'(7x,i2)',err=20) type
         if (type.ne.80) goto 20
         read (line,fmt=netfmt, err=20)
     .      type,ii,name,lad,lam,las,id1,lod,lom,los,id2,iht
c        modify site name from BLUE BOOK format to FONDA format
         j = remedy_space(name,desc,16,'_',2)
         nnet = nnet+1
         if (j.ge.8) then
            sname(nnet) = desc(1:8)
         else
            sname(nnet)(1:j) = desc(1:j)
            do k = j+1,8
               sname(nnet)(k:k) = '_'
            enddo
         endif
c        make sure the _NGS is added correctly
         if (sname(nnet)(5:8).eq.'____') then
            sname(nnet)(5:8) = '_NGS'
         endif
         name2 = desc(1:12)
         iesit(nnet) = ii
         
c        restore real number second and height
         sla = dble(las)*1.0d-5
         slo = dble(los)*1.0d-5
         ht = dble(iht)*1.0d-2

         s1 = dabs(dble(lad))+dble(lam)/60.0d0+sla/3600.0d0
         s2 = dabs(dble(lod))+dble(lom)/60.0d0+slo/3600.0d0
         if (id1.eq.'S') s1 = -s1
         if (id2.eq.'W') s2 = -s2
         write (ifil2,30) sname(nnet),name2,
     .         id1,lad,lam,sla,id2,lod,lom,slo,ht,ve,vn,vu,
     .         rtime,sigx,sigy,sigz
         write (ifil3,40) s2,s1,ve,vn,sname(nnet)
 20   continue
c
c10   format (6x,a8,3x,a1,2i3,f9.5,i4,i3,f9.5,3x,a1,f8.2)
 30   format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
c    .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
     .   a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
c
 50   continue
c
c     make all site names unique
      i = unique_nam1(nnet,iesit,sname)
      if (i.gt.0) then
         print*,'         ============'
         print*,i,' site names have been changed in your data file.'
         print*,' But the net_file keeps them untouched.'
         print*,' Please change these site names in net_file manually.'
         print*,'         ============'
      endif
c
      return
      end
