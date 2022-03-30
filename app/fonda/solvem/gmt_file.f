      subroutine gmt_file
c
c**************************************************************************
c
c*aut Gilbert FERHAT
c*aut GRGS 
c
c*ver July, 3 1995 v0_00
c
c*par          input parameters
c
c*par        : input data file (#13)
c
c*par          output parameters
c
c*par igmt   : GMT file (#29)
c
c*rol read an a priori coordinates and observation data files in order to
c*rol create a GMT instructions file
c
c**************************************************************************
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
      character*12 nom_du_site(900)
      character*1  sym1,sym2
      character*8  tempn2(900)
c
      integer igmt,kk,ij,ii,ios,l1,l2,l3,l4,nblen
      integer mn,mlength,nlength
      character*128 line
      dimension xlat_deg(900), xlong_deg(900)
c
      igmt = 29
      kk = 0
      ij = 0
      ii = 0
      xlong_max = -360.0d0 
      xlong_min =  360.0d0 
      xlat_max = -90.0d0
      xlat_min =  90.0d0 
c
c
c*******************************************************************
c     
c    read input data file
c
c*******************************************************************
      rewind (13)
c     print*, 'in gmt_file.f'
c
 10   continue
c     read site coordinate only
      read (13,'(a)',iostat=ios,end=20,err=30) line
c     print*, line
      goto 15
 30   goto 10 
 15   continue
c     skip comment line
      if (line(1:1).ne.' ') goto 10
c
c        comode = 2 = geodetic
      if (comode .eq. 2) then
c        fcode = 2 = char_8
         if (fcode.eq.2) then 
            ii = ii + 1
            read (line,24,iostat=ios)
     .      tempn2(ii),nom_du_site(ii),sym1,l1,l2,s1,
     .      sym2,l3,l4,s2,
     .      arad,vxt,vyt,vzt,
     .      site_rtime,sigx,sigy,sigz
 24         format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
     .              a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
c           write(*,24) 
c    .      tmpn2,nom_du_site(ii),sym1,l1,l2,s1,sym2,l3,l4,s2,
c    .      arad,vxt,vyt,vzt,
c    .      site_rtime,sigx,sigy,sigz
c
            xlong_deg(ii) = dble(l3)+dble(l4)/60.0d0+s2/3600.0d0
            xlat_deg(ii)  = dble(l1)+dble(l2)/60.0d0+s1/3600.0d0      
c           westward longitude
            if (sym2.eq.'W') xlong_deg(ii) = 360.0d0-xlong_deg(ii)
c           southern latitude
            if (sym1.eq.'S') xlat_deg(ii) = - xlat_deg(ii)
c
            xlong_max = max(xlong_deg(ii),xlong_max)
            xlat_max  = max(xlat_deg(ii) ,xlat_max)
            xlong_min = min(xlong_deg(ii),xlong_min)
            xlat_min  = min(xlat_deg(ii) ,xlat_min)
         endif
      endif
      goto 10
 20   continue
      close(13)
c
c**************************************************************
c     get the name of the sketch file name # 28
c         and a priori coordinate file name #13
c**************************************************************
c
      if (iomode(15).gt.0) then
         open(28,file=gmt_sketch_file,status='old',err=1000)
         mlength = nblen(gmt_sketch_file)
c        print*,' sketch file',gmt_sketch_file(1:mlength)
         close(28)
       else
         print*,'You have not specified the GMT sketch file'
         stop
      endif
c
      if (iomode(1).gt.0) then
         open(13,file=priori_coord,status='old',err=1000)
         nlength = nblen(priori_coord)
c        print*,' a priori file',priori_coord(1:nlength)
         close(13)
       else
         print*,'You have not specified the a priori value file'
         stop
      endif
c
c**************************************************************
c    create the GMT script
c**************************************************************
c
      write(igmt,100) '#!/bin/csh'
      write(igmt,100) '# GMT commands to draw a triangulation network'
      write(igmt,100) ' '
      write(igmt,100) '#MEASURE_UNIT = INCH '
      write(igmt,100) ' '
      write(igmt,100) 'psbasemap -JX6/8 -R\\'
c     write(*,*) ' *************************'
c     write(*,*) xlong_min-0.5,'/',xlong_max+0.5,'/'
c    &              ,xlat_min-0.5,'/',xlat_max+0.5
c     write(*,*) ' *************************'
      write(igmt,40) xlong_min-0.5,'/',xlong_max+0.5,'/'
     &              ,xlat_min-0.5,'/',xlat_max+0.5,' ',char(92)
 40   format(4(f5.1,a1))
c40   format(2(f5.1,a1),f4.1,a1,f4.1,a1,a1)
      write(igmt,100) ' -P -Bf0.1a1/WeSn -K -Y1.5 >! $1.ps'
      write(igmt,100) ' '
      write(igmt,100) '#scale bar '
      write(igmt,100) 'mapproject <<! -R -F -I -H0 -JX >!\\',
     & ' tmp.scale'
      write(igmt,100) ' 55.e4  5.0e4'
      write(igmt,100) ' 55.e4  7.0e4'
      write(igmt,100) ' 55.e4  5.0e4'
      write(igmt,100) ' 65.e4  5.0e4'
      write(igmt,100) ' 65.e4  7.0e4'
      write(igmt,100) ' 65.e4  5.0e4'
      write(igmt,100) '!'
      write(igmt,100) ' '
c
      write(igmt,100) ' '
      write(igmt,100) 'pscoast -JX -R -P -W1 -B -O -K -N1ta >> $1.ps '
      write(igmt,100) ' '
      write(igmt,100) '# add stations locations '
      write(igmt,100) 
     &' awk ''{print $2,$3}'' ',priori_coord(1:nlength),' | psxy -R -JX
     &-Ss0.1 -G0 -O -K >> $1.ps'
      write(igmt,100) ' '
      write(igmt,100) '# draw directions'
      write(igmt,100) 'psxy ',gmt_sketch_file(1:mlength),' ',char(92)
      write(igmt,100) ' -R -JX -H0 -M -W1/0/0/0ta200_20_1700_20:20
     &-P -K -O -: >> $1.ps'
      write(igmt,100) ' '
c
      write(igmt,100) '# add station numbers'
      write(igmt,100) '#station    Long_deg Lat_deg'
      write(igmt,100) ' pstext <<! -R -JX -H0  -K -O -Y0.18 >> $1.ps'
c
      do mn = 1, ii
      write(igmt,100) xlong_deg(mn),xlat_deg(mn),
     .     ' 7 0 29 10 ',tempn2(mn)
      enddo 
      write(igmt,100) '!' 
c
      write(igmt,100) ' '
      write(igmt,100) '# add station names'
      write(igmt,100) '#station    Long_deg Lat_deg'
      write(igmt,100) ' pstext <<! -R -JX -H0  -K -O -Y0.08
     . >> $1.ps '
      do ij = 1, ii
          write(igmt,100) xlong_deg(ij),xlat_deg(ij),' 7 0 29 1 ',
     &    nom_du_site(ij)
      enddo
      write(igmt,100) '!' 
      write(igmt,100) ' '
      write(igmt,100) ' '
      write(igmt,100) '# justify center of text on mid point of scale'
      write(igmt,100) 'mapproject <<! -R -F -I -H0 -JX  >! tmp.labels'
      write(igmt,100) ' 60.e4   5.4e4 12 0 1 2 100 km'
      write(igmt,100) '!'
      write(igmt,100) ' '
      write(igmt,100) '# Note use of -A to supress great circle'
      write(igmt,100)  'psxy tmp.scale -JX -R -H0 -K -O -L -W4 -A 
     . -Y-0.1 >> $1.ps'
      write(igmt,100) ' '
      write(igmt,100) 'pstext tmp.labels -H0 -R -JX -O   -K >> $1.ps'
      write(igmt,100) ' ' 
      write(igmt,100) '# must be last item'
      write(igmt,100) 'psxy <<! -H0 -R -JX -O   -X+1.1 -U" " >> $1.ps'
      write(igmt,100) ' 0. 0. '
      write(igmt,100) '!'
      write(igmt,100) ' '
      write(igmt,100) 'gs $1.ps'
c
 100  format (a)
      open(13,file=priori_coord,status='old',err=1000)
      rewind(13)
      return
c     troubled stop
 1000 print*,' Suicide due to gmt file'
      stop
      end

