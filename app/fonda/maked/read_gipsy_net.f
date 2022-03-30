      subroutine read_gipsy_net(ifil1,ifil2,ifil3,nnet,netfmt,
     .                          iesit,sname,rtime)
c
c     read GIPSY statistics.gc  file and create FONDA format site-table file
c     ( SV6 format).
c     modify to fit new GIPSY statistics.gc or .gd file format
c     meanwhile still accept the old GIPSY format
c                            Danan Dong 06/16/95
c
c     coordinate system : (comode)
c         1. geocentric spherical coordinate (lat,lon,radius)  ***
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

*      character*16 name,desc
      character*8  sname(200)
      character*4 aname,bname
      character*80 netfmt
      character*120 line
      character*1 ns,ew 
      character*4 tail1,tail2,cord
      character*3 month,adum
*      integer type, ii,remedy_space, 
      integer   i,lad,lam,lod,lom
      integer   iesit(200),unique_nam1
      integer*4 julday
      integer imonth, iyear, iday, ifil1,ifil2,ifil3,nnet,idum1,idum2
      real*8 lamin, lomin, ht, las, los, flat, flon
      real*8 rtime
      logical ierr
c      
      write (ifil3,'(a)') '* Network distribution'
      ve = 0.0d0
      vn = 0.0d0
      vu = 0.0d0
      sigx = 0.0d0
      sigy = 0.0d0
      sigz = 0.0d0
c
c     default format
      if (netfmt(1:1).eq.'*')
     . netfmt(1:36) = '(1x,a4,8x,i2,a3,i2,2x,f15.9,6x,f9.4)'
   
c
c     write FONDA style network file header
      call net_head(ifil2)
      nnet = 0
      aname = '    '
      bname = '    '
c defining tails
      tail2 = '____'
      tail1 = '_GPS'
      do 20 i = 1,30000
         ierr = .false.
c        Reads 3 lines of a *.gc file to get sitename lat, lon, epoch, sig1 
c                            -bawden 3-94
c        new GIPSY statistics file has 6 lines for each station
         read (ifil1,'(a120)', end=50, err=20) line
         read (line,'(1x,a4,2x,a4)',err=20) bname,cord
         if (cord.eq.'    ') goto 20
         if (aname.eq.'    ') aname = bname
c        new site, output old site information
         if (aname.ne.bname) then
c           modify site name from GIPSY format to FONDA format
            nnet = nnet+1
            sname(nnet) = aname // tail1

            write (ifil2,30) sname(nnet),aname,tail2, 
     .         ns,lad,lam,las,ew,lod,lom,los,ht,ve,vn,vu,
     .         rtime,sigx,sigy,sigz
            write (ifil3,40) flon,flat,ve,vn,sname(nnet)
            aname = bname
         endif
      
         if (cord(1:3).eq.'RAD') then
            read (line,'(22x,f15.4,6x,f9.4)',err=24) dum1,dum2
            goto 22
 24         ierr = .true.
 22         if (cord.ne.'RADV') then
               ht = dum1
               if (ierr) then
                  sigz = 99.999
               else
                  sigz = dum2
               endif
            else
               if (ierr) then
                  vu = 0.0d0
                  goto 20
               endif
c              from mm/yr to m/yr
               vu = dum1*0.001d0
c              remove unrealistic estimates
               if (dum2.ge.2.0d2.or.dabs(vu).gt.1.0d0) vu = 0.0d0
            endif
         else 
            read (line,fmt=netfmt,err=26)
     .         bname,idum1,adum,idum2,dum1,dum2
            goto 28
 26         ierr = .true.
 28         if (cord.eq.'LAT ') then
               iyear = idum1
               month = adum 
               iday = idum2
               flat = dum1
               if (ierr) then
                  sigx = 99.999
               else
                  sigx = dum2
               endif
c modify north-south and east-west for lat & lon  --Bawden 3-94
               if (flat.lt.0) then 
                  ns='S'
               else
                  ns='N'
               end if
c setting neg lat and lon to positive numbers 
               temp = dabs(flat)
               lad = int(temp)
               lamin = temp-dble(lad)
c modify degrees.fraction to deg min sec.fraction 
               call deg2sec(lamin,lam,las)
c modify yymondd to epoch(year.partofyear)
               if (month.eq.'JAN') then
                  imonth=1
               else if (month.eq.'FEB') then
                  imonth=2
               else if (month.eq.'MAR') then
                  imonth=3
               else if (month.eq.'APR') then
                  imonth=4
               else if (month.eq.'MAY') then
                  imonth=5
               else if (month.eq.'JUN') then
                  imonth=6
               else if (month.eq.'JUL') then
                  imonth=7
               else if (month.eq.'AUG') then
                  imonth=8
               else if (month.eq.'SEP') then
                  imonth=9
               else if (month.eq.'OCT') then
                  imonth=10
               else if (month.eq.'NOV') then
                  imonth=11
               else
                  imonth=12
               end if
c  calling function julday  ~/fonda/com/. to do the conversion
               rtime = 1900.0d0+julday(imonth,iday,iyear,1)/365.25d0
            else if (cord.eq.'LON ') then
               flon = dum1
               if (ierr) then
                  sigy = 99.999
               else
                  sigy = dum2
               endif
               if (flon.lt.0) then
                  ew='W'
               else
                  ew='E'
               end if
               temp = dabs(flon)
               lod = int(temp)
               lomin = temp-dble(lod)
               call deg2sec(lomin,lom,los)
c           from mm/yr to m/yr
            else if (cord.eq.'LATV') then
               if (ierr) then
                  vn = 0.0d0
                  goto 20
               endif
               vn = dum1*0.001d0
c              remove unrealistic estimates
               if (dum2.ge.2.0d2.or.dabs(vn).gt.1.0d0) vn = 0.0d0
            else if (cord.eq.'LONV') then
               if (ierr) then
                  ve = 0.0d0
                  goto 20
               endif
               ve = dum1*0.001d0
c              remove unrealistic estimates
               if (dum2.ge.2.0d2.or.dabs(ve).gt.1.0d0) ve = 0.0d0
            endif
         endif
c
 20   continue
c
c10   format (6x,a8,3x,a1,2i3,f9.5,i4,i3,f9.5,3x,a1,f8.2)
 30   format (1x,a8,1x,a4,a4,5x,a1,2(i2,1x),f8.5,1x,
c    .        a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
     .   a1,i3,1x,i2,1x,f8.5,f13.4,2f8.2,f7.3,f9.3,3f8.4)
 40   format (1x,2f11.5,2f9.3,2x,a8)
c
 50   continue
c
c     last site
      if (cord.ne.'    ') then
c        modify site name from GIPSY format to FONDA format
         nnet = nnet+1
         sname(nnet) = bname // tail1
         write (ifil2,30) sname(nnet),bname,tail2, 
     .      ns,lad,lam,las,ew,lod,lom,los,ht,ve,vn,vu,
     .      rtime,sigx,sigy,sigz
         write (ifil3,40) flon,flat,ve,vn,sname(nnet)
      endif
c
c     make all site names unique
      i = unique_nam1(nnet,iesit,sname)
      if (i.gt.0) then
         print*,'         ============'
         print*,i,' site names have been changed in your data file.'
         print*,' But the net_file keeps untouched.'
         print*,' Please change these site names in net_file manually.'
         print*,'         ============'
      endif
c
      return
      end

      subroutine deg2sec(lmin,imins,secs)

c this subroutine will convert degrees to minutes and seconds

      double precision lmin, modsec, secs, min
      integer imins
                  
      min=lmin*60.0d0
      imins = int(min)
                              
      modsec = dmod(min,1.0d0)
                                    
      secs= modsec*60.0d0
      
      return
      end

