      program displace
C
      implicit none
C
C     PROGRAM seismic displacement predictor        1/26/94
C	Tracadas
C       original source code is in 
C       /data9/brad1/projects/gps/disloc/displace.f
C
C     Modified by B. Souter                         2/1/94
C
C     Program to calculate coseismic displacements
C     using the elastic half space models of Okada
C
C MOD TAH 190104: Enabled the tensile fault component (was coded
C     but value was set to zero).
C
C***********************************************************************
C
C234567890         21234567890         41234567890         61234567890
C        11234567890         31234567890         51234567890         712
C
C     MAIN 

      include 'param.displace.h'

      real*8 dlat   ! Latitude range based on grid cells
      
      integer m, n
      
C     read the input file and find out whether to generate a grid of 
C     points, or whether specific stations are input.
C MOD TAH 190104: Return u3 as well
      call readinput(flong, flat, fstrike, depth, dip, length, width, 
     +     u1, u2, u3, scale, makegrid, numsta, numflts,
     +     xmin, xinc, xnum, ymin, yinc,ynum,slong, slat, 
     +     sname, origlong, origlat,infile,outfile)
C     set up an xy grid or transform lat and long into xy coordinates
      if (makegrid.eq.0) then
         call creategrid(xmin, xinc, xnum, ymin, yinc, ynum, x, y)
C        MOD TAH Estimate the latitude range and compute parallels
C        and approximte min/max longitude and latitude.
         dlat = (ynum*yinc/6371)*57.2958  ! degrees (only approximate needed)
         parallel1 = origlat-dlat/2 + 0.135*dlat
         parallel2 = origlat+dlat/2 - 0.135*dlat 
         minlon = origlong + xmin/6371*57.2958/cos(origlat/57.2958)
         maxlon = origlong + (xinc*xnum)/6371*57.2958/
     .                            cos(origlat/57.2958)
         minlat = origlat + ymin/6371*57.2958
         maxlat = origlat + (yinc*ynum)/6371*57.2958
      elseif (makegrid.eq.1) then
c     set up translation parameters for lambert conformable 
c     conical projection
         call calctrans(origlong, origlat, maxlon, minlon, maxlat, 
     +        minlat, parallel1, parallel2, slong(1), slat(1),
     +        numsta)
c     do the translation to xy coordinates in km
         call lamberttrans(origlong, origlat, maxlon, minlon,
     +        maxlat, minlat, parallel1, parallel2,
     +        numsta, 0, slong, slat, x, y)
      else 
         write(*,*) 'makegrid must equal 0 or 1'
         stop 'makegrid must equal 0 or 1'
      endif
      write(*,120) makegrid, minlon, maxlon, minlat, maxlat,  
     .             parallel1, parallel2
 120  format('DISPLACE G ',i2,' Long Rng ',2F10.3,' Lat Rng ',2F10.3,
     .       ' Parallels ',2F10.3,' deg')
      
C     Initialize displacement arrays.
      call initdispl(displn, disple, displu, numsta)
         
C     Find displacement for each fault
      do 180 m = 1, numflts

C     find xy coordinates of fault origin using lambert
C     conformable conical projection
         call lamberttrans(origlong, origlat, maxlon, minlon,
     +        maxlat, minlat, parallel1, parallel2,
     +        one, zero, flong(m), flat(m), fx, fy)

C     translate xy coordinates to fault origin
         call translatexy(fx,fy,x,y,tx,ty,numsta)

C     rotate xy coordinates 
         call rotatexy(fstrike(m),tx,ty,rx,ry,numsta)
 

C     find displacement at all rotated xy points
* MOD TAH 190104: Pass u3(m) as well.
         call getdispl(rx, ry, depth(m), dip(m), length(m), width(m), 
     +        u1(m), u2(m), u3(m), deltan, deltae, deltau, numsta)

C     remove the rotation of the displacements
         call unrotate(fstrike(m), deltan, deltae, numsta)
         
C     Add in this fault's displacements at each station
         call adddispl(deltan, deltae, deltau, 
     +        displn,disple, displu, numsta)

 180  continue
      
C     create an output file for use in gmt plotting routines
C MOD TAH 190104: Pass u3 as well.
      call writeoutput(numflts, flong, flat, fstrike, depth, dip,
     +     length, width, u1, u2, u3, scale, makegrid, numsta,  
     +     x, y, slong, slat, sname, displn, disple, displu,
     +     origlong, origlat,outfile)
      
      call writegmt( flong, flat, fstrike, depth, dip, length, 
     +     width, u1, u2, u3, scale, makegrid, numflts, numsta, x, y, 
     +     slong, slat, sname, displn, disple, displu,infile)
      
      stop
      end
      
************************************************************************

      subroutine readinput(flong, flat, fstrike, depth, dip, length, 
     +       width, u1, u2, u3, scale, makegrid, numsta, numflts,
     +       xmin, xinc, xnum, ymin, yinc, ynum, slong, slat, 
     +       sname, origlong, origlat,infile,outfile)

      implicit none

* MOD TAH Read U3 as well.  Line can be
* depth(m), dip(m), length(m), width(m), u1(m), u2(m), <u3(m)> <scale(m)>
* wherr u3 and scale string are optional.

      include 'param.displace.h'
      character*256 line
      character*8 word

      integer  m, n, ierr, jerr, indx
      integer*4 trimlen   ! Length of string

      logical done

      real*8 geod(3)   ! Geodetic coordinates
      real*8 dglat, dglong  ! Change from edge to center of fault plane

      character*20 u3str, scalestr  ! Strings into which u3 and scale may
                                    ! read.

      call get_runstring(infile, outfile)

C     Open up the input file.  Kill program if there is an error
      open(100, file=infile, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open', infile,
     +                  1,'readinput')

* MOD TAH 031026: Get the Globk hypocenter information:
*     Latutude, longitude and depth (note radius comes before depth in
*     globk eq entry.  Generates fatal if not all present.
* MOD TAH 170415: Changed input line to allow displace to compute values
*     and to added the code and date to eq_def lines
      read(100,'(a)',iostat=ierr) line
      call report_error('IOSTAT',ierr,'read',infile,1,
     .     'First line of input file')
      indx = 1
      call GetWord(line, word, indx )
      call casefold(word)
      if( word(1:4).eq.'EQ_D' ) then
*         New eq_def line passed, so decode
          call GetWord(line,eq_code,indx)
          call multiread(line,indx,'I4',jerr, eq_date, word, 5)
          call report_error('IOSTAT',jerr,'decod',line,0,'Date')
          glat = -99
      else    ! Read old style globk line      
          read(line,*, iostat=ierr) glat, glong, gdepth
          call report_error('IOSTAT',ierr,'read',infile,1,
     .     'GLOBK coordinate line -- check input')
          eq_code = '<CODE>'
          eq_date = 0
      end if
      
C     Read how many faults.
      read(100,'(A)', iostat=ierr) line
      read(line,*, iostat=ierr) numflts

      if ((numflts.gt.MAXFLT).or.(numflts.le.0)) then
         write(*,*) 'Displace: number of faults, ',numflts,
     +        ' exceeds maximum allowed, ',MAXFLT
         stop
      endif

      do 150 m = 1,numflts
         
C     Get fault characteristics
         read(100,'(A)', iostat=ierr) line
         read(line,*, iostat=ierr) flong(m), flat(m), fstrike(m)
         
         read(100,'(A)', iostat=ierr) line
* MOD TAH 190104: Test for option u3str and scalestr
         u3str = ''
         scalestr = ''
         read(line,*, iostat=ierr) depth(m), dip(m), length(m), 
     +        width(m), u1(m), u2(m), u3str, scalestr


* MOD TAH 190104: Additional logic in decoding line added
         if (ierr.eq.-1) then  ! not all arguments passed; 
                               ! see what we have
            if( u3str(1:1).eq.' ' ) then  ! Nothing after u2,
                                          ! Set u3=0 and default label
               scale(m) = 'mm'
               u3(m) = 0 
            else       ! Try to read u3 string as number
               call check_num(u3str,jerr) 
               if( jerr.eq.0 ) then
                  read(u3str,*) u3(m)
                  if( scalestr(1:1).eq.' ' ) then
                     scale(m) = 'mm'
                  else
                     scale(m) = scalestr
                  endif
               else    ! Assume scale is in u3str
                  scale(m) = u3str
               endif
            endif 
            ierr = 0
         else   ! All arguments passed.
            read(u3str,*,iostat=jerr) u3(m)
            scale(m) = scalestr
         endif
         
 150  continue

* MOD TAH 710415: Finish up globk position calcations 
      if( glat.eq.-99 ) then
*         Compute lat/long and depth of center of first fault
*         plane
          m = 1
 
          gdepth = depth(m) - width(m)*sin(dip(m)/rad_to_deg)/2
          dglat  = (length(m)*cos(fstrike(m)/rad_to_deg) +
     .              width(m)* cos(dip(m)/rad_to_deg)*
     .                        sin(fstrike(m)/rad_to_deg))/
     .              (2*earth_rad/1000)   ! km to angle rads

          dglong  = (length(m)*sin(fstrike(m)/rad_to_deg) -
     .              width(m)* cos(dip(m)/rad_to_deg)*
     .                        cos(fstrike(m)/rad_to_deg))/
     .              (2*cos(flat(m)/rad_to_deg)*earth_rad/1000) ! km to angle
          glat  = flat(m) + dglat*rad_to_deg
          glong = flong(m) + dglong*rad_to_deg 
 
      endif
 

*     Convert geodetic coordinates to XYZ
      geod(1) = pi/2 - glat*pi/180
      geod(2) = glong*pi/180
      geod(3) = -gdepth*1000.d0
      call GEOD_to_XYZ( geod, eq_pos )
      
      
C     Read makegrid, see whether grid is to be generated
      read(100,'(A)', iostat=ierr) line
      read(line,*, iostat=ierr) makegrid
      
C     If makegrid is 0 read in grid generating parameters grid
C     If makegrid is 1 read in numsta and station coordinates and names
      
      if (makegrid.eq.0) then 
         read(100,'(A)', iostat=ierr) line
         read(line,*, iostat=ierr) origlong, origlat
         read(100,'(A)', iostat=ierr) line
         read(line,*, iostat=ierr) xmin, xinc, xnum, ymin, yinc, ynum
         numsta = xnum * ynum
         if ((numsta.gt.MAXSTA).OR.(NUMSTA.LE.0)) then
            write(*,*) 'Displace: number of stations, ',numsta,
     +           ' exceeds maximum allowed, ',MAXSTA
            stop
         endif
      else 
c     set up origin of grid by calculating center later
         read(line,*, iostat=ierr) makegrid

         done = .false.
         numsta = 0
         do while ( .not. done )
            read(100, '(A)', iostat=ierr) line
            if( ierr.eq.0 ) then
                numsta = numsta + 1
                if( numsta.gt.MAXSTA ) then
                    write(*,*) 'Displace: number of stations, ',
     .                          numsta, ' exceeds maximum allowed, ',
     +                          MAXSTA
                    stop 'DISPLACE: Too many stations'
                end if
                n = numsta
                read(line,*, iostat=ierr) slong(n), slat(n), sname(n)
                if( ierr.ne.0 ) then
                   if( trimlen(line).gt.0 )    ! Blank line at end ignore
     .             write(*,110) ierr, line(1:80)
 110               format('IOSTAT error ',i4,' decoding ',a)
                   numsta = numsta - 1
                end if
            else
                done = .true.
            end if
         end do
      end if
      
C     Check for any read errors in the above statements.
      if( ierr.ne.-1 ) 
     .call report_error('IOSTAT',ierr,'read', infile,
     +     1,'readinput')
      write(*,*) 'Number of stations read ',numsta
      
      close(100)
      
      call checkinput(flong, flat, fstrike, depth, dip, length, 
     +     width, u1, u2, u3, scale, makegrid, numsta, numflts,
     +     xmin, xinc, xnum, ymin, yinc, ynum, slong, slat, 
     +     sname,infile)
      
      return
      end

************************************************************************

      subroutine creategrid(xmin, xinc, xnum, ymin, yinc, ynum, x, y)

      implicit none
       
      include 'param.displace.h'
      integer
     +     indx, n, m
      
      indx = 0
      do 200 n=1,xnum
         do 250 m=1,ynum
            indx = indx + 1
            x(indx) = xmin + ((n-1) * xinc)
            y(indx) = ymin + ((m-1) * yinc)
 250     continue
 200  continue
      
      return
      end
      
************************************************************************
      
      subroutine getdispl(rx, ry, de, di, len, wid, us, ud, ut,
     +     deltan, deltae, deltau, numsta)
      
C     compute displacement for each site using Okada method
C MOD TAH 190104: Added passing of Ut (tensile); u1 and u2 are us and ud ger

      implicit none
       
      include 'param.displace.h'
      integer
     +     n
      
      real
     +     de, di, len, wid, us, ud, ut

      double precision
     +     tempde, tempdi, templen, tempwid, tempus, tempud, temput,
     +     temprx, tempry, tempdele, tempdeln, tempdelu

      tempde = de
      tempdi = di
      templen = len
      tempwid = wid
      tempus = us
      tempud = ud
      temput = ut

      do 500 n = 1, numsta
         temprx = rx(n)
         tempry = ry(n)

* MOD TAH 190104: Added temput to calling arguments (tensile component)
         call okada(temprx,tempry,tempde,tempdi,templen,tempwid,
     +        tempus,tempud,temput, tempdele, tempdeln, tempdelu)
         deltae(n) = tempdele
         deltan(n) = tempdeln
         deltau(n) = tempdelu
 500  continue

      return
      end
      
************************************************************************
      
      subroutine writegmt( flong, flat, fstrike, depth, dip, length, 
     +     width, u1, u2, u3, scale, makegrid, numflts, numsta, x, y, 
     +     slong, slat, sname, displn, disple, displu,infile)

* MID TAH 190104: Incldued u3 (maybe zero).      

      implicit none

      include 'param.displace.h'
      
      integer numnodes
      parameter (numnodes = 4)
      
      integer ierr, n, m
      
C     Write the file with the fault nodes to fault.out 
      call writefault(numflts, flong, flat, fstrike, depth, dip,
     +     length, width)
C     Write the file with the fault traces 
      call writetrace(flong, flat, fstrike, dip, depth, length,
     +     width, numflts)
      
C     Write the file with the vectors 
      call writevectors(displn, disple, slong, slat, sname, 
     +     numsta, makegrid,infile)

C     Write the file with the heights 
      call writeheights(displu, slong, slat, sname, 
     +     numsta, makegrid,infile)
      
      return
      end
************************************************************************

      subroutine writeoutput(numflts, flong, flat, fstrike, depth, dip,
     +     length, width, u1, u2, u3, scale, makegrid, numsta, x, y, 
     +     slong, slat, sname, displn, disple, displu,
     +     origlong, origlat,outfile)

* MOD TAH 190104: Added U3 to output

      implicit none
       
      include 'param.displace.h'
c      character*256 infile, outfile
      
      integer
     +     ierr, n, m

      real*8 edist    ! Distance of point from earthqauke
      real*8 hamp, vamp  ! Horizontal and vertical amplitudes
                      ! and distance (m) from hypocenter
      real*8 eq_sh, eq_sv  ! Sigmas for distance dependent terms in globk
                      ! earthquake model
      character*18 eq_date_string  ! String with date or with labels,

****  Clear solution variables
      grad = 0.d0
      nume = 0
      neq = 0
      bh = 0
      bv = 0
      
C     Open the output file.
      open(300, file=outfile, iostat=ierr)
      call report_error('IOSTAT',ierr,'open', outfile,
     .     1,'writeoutput')
      
C     Write the fault characteristics to the file.
      write(300,440,iostat=ierr) 
 440  format('*     Long       Lat      Stike     Depth       Dip',
     .       '    Length     Width     U1 (SS)     U2 (DS)    U3 ',
     .       '(Tens) Unit')
      do 460 m = 1, numflts
         write(300,450,iostat=ierr) flong(m), flat(m), fstrike(m), 
     +        depth(m), dip(m), length(m), 
     +        width(m), u1(m), u2(m), u3(m), scale(m)
460   continue
      write(300,'(a)',iostat=ierr) '*' 
450   format('* ',7(f9.3,1x),3(F11.2,1x),1x, a8)

      if (makegrid.eq.0) then
        call transcoord(origlong, origlat, slong, slat, x, y, numsta)
        write(300, 405, iostat=ierr)
 405    format('*       x         y    longitude  latitude     east',
     .         '     north        up Grd      AbsNE     AbsdU',
     .         '     Dist',/, 
     .         '*     (km)      (km)     (deg)      (deg)      (mm)'
     .         '      (mm)      (mm) Grd       (mm)      (mm)',
     .         '      (km)')
        do 420 n = 1, numsta
           call geqmod(slong(n), slat(n), 
     .                  disple(n),displn(n),displu(n), 
     .                  hamp, vamp, edist)
           write(300, 400, iostat=ierr) x(n), y(n), slong(n), slat(n), 
     +          disple(n), displn(n), displu(n), hamp, vamp, 
     .          edist/1000.d0
400       format(1x,2f10.2,2f10.4,3f10.2, ' Grd ',2F10.2,1x,f10.2)
          call report_error('IOSTAT',ierr,'write', outfile,
     +                      1,'writeoutput')

420     continue

      else 
          write(300, 415, iostat=ierr)
          do 430 n = 1, numsta
            call geqmod(slong(n), slat(n), 
     .                  disple(n),displn(n),displu(n), 
     .                  hamp, vamp, edist)

            write(300, 410, iostat=ierr) slong(n), slat(n), sname(n), 
     +          disple(n), displn(n), displu(n), hamp, vamp, 
     .          edist/1000.d0
415         format('*  Lonitude   Latitude Site          dEast',
     .             '     dNorth    dUp   Mod     AbsdNE',
     .             '      AbsdU     Dist',/,
     .             '*  (deg)        (deg)                (mm)',
     .             '       (mm)     (mm)           (mm)       (mm)',
     .             '      (km)')

410         format(2f11.5,1x,a8,1x,3f10.2,' Mod ',2F10.2,1x,f10.2)
            call report_error('IOSTAT',ierr,'write', outfile,
     .                      1,'writeoutput')
430       continue

      endif

****  Now finish up globk calculations
* MOD TAH 170413: Double output value to give more margin since
*     these values are applied as a constraint in globk.
      if( eq_date(1).eq.0 ) then
          ! No date; make generic string
          eq_date_string = 'YYYY MM DD HR MIN'
      else
          write(eq_date_string, 500) eq_date
 500      format(I4,4(1x,I2.2))
      endif 

      eq_sh = bh/neq * 2
      eq_sv = bv/neq * 2
      if( glong.lt.0 ) glong = glong + 360
      write(300,505) numflts    
 505  format('* DISPLACE Fault summary for ',i3,' faults',/,
     .       '*  F  Long     Lat    Strike    Depth    Dip',
     .       '     Length   Width        u1        u2        u3',/,
     .       '*    (deg)    (deg)   (deg)      (km)   (deg) ',
     .       '   (km)     (km)        (mm)      (mm)      (mm)')
      do m = 1, numflts
         write(300,510,iostat=ierr) m, flong(m), flat(m), fstrike(m), 
     .        depth(m), dip(m), length(m), 
     .        width(m), u1(m), u2(m), u3(m), scale(m)
 510     format('* ',I3,1x,F7.3,1x,F7.3,1x,F7.2,1x,F8.2,1x,F7.2,1x,
     .         F8.2,1x,F8.2,1x,3(F9.1,1x),a)
      end do 
      write(300,'("*")')
      write(300,515) glat, glong, grad/1000.d0, gdepth, 
     .       eq_sh/1000.d0, eq_sv/1000.d0
 515  format('* GLOBK eq parameters:',/,
     .       '* Lat, Long, Radius (km) and Depth (km)',
     .        F10.5,1x,F10.5,1x,F8.2,1x,F8.2,/,
     .       '* Radius Sigmas: Horizontal Vertical (m) ',2F12.3)
*     Re-write the model values

      write(300,520) trim(eq_code), glat, glong, grad/1000.d0, gdepth, 
     .               eq_date_string, eq_sh/1000
 520  format('* eq_def    ',a,1x,F10.5,1x,F10.5,1x,F8.2,1x,F8.2,1x,
     .       a,1x, F10.4)
      write(300,525) trim(eq_code)
 525  format('* eq_rename ',a)
      write(300,530) trim(eq_code), eq_sh/1000, eq_sh/1000, eq_sv/1000
 530  format('* eq_cosei  ',a,6x,' 0.001 0.001 0.001 ',3(F10.4,1x))
      write(300,'("*")')

      close(300)

      return
      end

************************************************************************

CTITLE get_runstring                     Tracadas                 5/4/92

      subroutine get_runstring(filename, dispfn)

*     This procedure decodes the runstring and lists the help file if 
*     the runstring is not correct.

*     variables:rcpar: Function to read the runstring.
*               len_run: length of runstring parameter returned.

      implicit none

      integer rcpar, len_run
      character*256 filename, dispfn

*     Get the name of the input site and setup file.
      len_run = rcpar(1, filename )
      if( len_run.le.0 ) then
        call proper_runstring('displace.hlp', 'displace', 1)
      end if

*     Get the name of the output displacement file.
      len_run = rcpar(2, dispfn )
      if( len_run.le.0 ) then
        call proper_runstring('displace.hlp', 'displace', 1)
      end if

      return
      end

***************************************

 
      subroutine checkinput(flong, flat, fstrike, depth, dip, length, 
     +      width, u1, u2, u3, scale, makegrid, numsta, numflts,
     +      xmin, xinc, xnum, ymin, yinc, ynum, slong, slat, 
     +      sname, infile)

C     See if there are any problems with the input data:

      implicit none
 
      include 'param.displace.h'
c      character*256 infile
      character*5 thescale

      integer
     +         m, n, ierr

      ierr = 0

      thescale = scale(1)
      if (((index(thescale,'mm').gt.0).or.(index(thescale,'dm').gt.0))
     +     .or.((index(thescale,'cm').gt.0).or.
     +     (index(thescale,'m').gt.0))) then
         ierr = 0
      else
         ierr = -1
      endif
      
      do 650, m = 1, numflts
         if ((flong(m).lt.-360).or.(flong(m).gt.360)) then
            ierr = -2
         elseif ((flat(m).lt.-90).or.(flat(m).gt.90)) then
            ierr = -3
         elseif ((dip(m).gt.360).or.(dip(m).lt.-360)) then
            ierr = -4
         elseif (depth(m).lt.(width(m) * sin(dip(m)*PI/180.0))) then
            ierr = -5
         elseif (thescale.ne.scale(m)) then
            ierr = -6
         endif
         
 650  continue
      
      if ((numsta.gt.maxsta).or.(numsta.le.0)) then
         ierr = -7
      endif
      
      if (makegrid.eq.1) then
         do 600 n = 1,numsta
            if (((slong(n).lt.-360).or.(slong(n).gt.360)).or.
     +           (( slat(n).lt.-90 ).or.( slat(n).gt.90 ))) then
               ierr = -8
            end if
 600     continue
      endif      
      
      
      call report_error('IOSTAT',ierr,'input', infile,
     +     1,'checkinput')
      
      return
      end
      
************************************************************************
      
      subroutine writefault(numflts, flong, flat, fstrike, depth, dip,
     +     length, width)
      

      implicit none

      include 'param.displace.h'
      
      integer
     +     m, n, ierr
      
      real
     +     nodelong(4), nodelat(4), nodedep(4),
     +     nodex(4), nodey(4)
      
      open(unit=101,file='fault.out',iostat=ierr)
      write(101,601) 'node long', 'node lat', 'node depth'
      write(101,'(a)') '>'
      
      do 270 m = 1, numflts
         
 601     format(3A15)
         
C     For each fault, 
C     Calculate the 4 nodes of the fault, and their depths, and
C     write the output to a file called fault.out in three
C     columns:  long, lat, and depth in km.
         
         call calcnodes(flong(m), flat(m), fstrike(m), depth(m), 
     +        dip(m), length(m), width(m), 
     +        nodelong, nodelat, nodedep,
     +        nodex, nodey)
         
         do 610 n = 1, 4
            write(101,611) nodelong(n), nodelat(n), nodedep(n)
 611        format (3f15.4)
 610     continue
         write(101,'(a)') '>'
         
 270  continue
      
      close(unit=101)
      
      return
      end
      
************************************************************************
      
      subroutine writetrace(flong, flat, fstrike, dip, depth, length,
     +     width, numflts)

      implicit none
       
      include 'param.displace.h'
      
      integer 
     +     m, n, ierr
      
      real
     +     nodelong(4), nodelat(4), nodedep(4),
     +     nodex(4), nodey(4), 
     +     ftlong(2), ftlat(2)
      
C     Write out the fault trace file
      
      open(unit=105,file='trace.out',iostat=ierr)
      write(105,605) 'fault long', 'fault lat'
      write(105,'(a)') '>'
      do 280 m = 1, numflts
 605     format(2A10)
         
         call calcnodes(flong(m), flat(m), fstrike(m), depth(m), 
     +        dip(m), length(m), width(m), 
     +        nodelong, nodelat, nodedep,
     +        nodex, nodey)
         
         call calctrace(flong(m), flat(m), fstrike(m), dip(m), 
     +        depth(m),ftlong, ftlat, nodex, nodey)
         
         do 650 n = 1,2
            write(105,615) ftlong(n), ftlat(n)
 615        format (2f10.4)
 650     continue
         write(105,'(a)') '>'
         
 280  continue
      close(unit=105)
      
      return
      end
      
************************************************************************
      
      subroutine writevectors(displn, disple, slong, slat, sname, 
     +     numsta, makegrid,infile)
      
C     Write out locations and the displacement vectors

      implicit none
       
      include 'param.displace.h'
      
      integer
     +     n, ierr
      
      open(unit=104,file='vectors.out',iostat=ierr)
      write(104,'(A60)') infile
      write(104,604) 'stn long', 'stn lat','displ e','displ n'
 604  format(4A10)
      if (makegrid.eq.0) then
         do 640 n = 1,numsta
            write(104,614) slong(n), slat(n), disple(n), displn(n),
     +           0,0,0
 640     continue
      else
         do 639 n = 1,numsta
            write(104,614) slong(n), slat(n), disple(n), displn(n),
     +           0,0,0, ' ', sname(n)
 639     continue
         
 614     format (2f10.4, 2f10.4, 3i3, 2a)
      endif
      
      close(unit=104)
      
      return
      end
      
************************************************************************
      
      subroutine writeheights(displu, slong, slat, sname, 
     +     numsta, makegrid,infile)
      
C     Write out locations and the displacement heights

      implicit none
       
      include 'param.displace.h'
      
      integer
     +     n, ierr
      
      open(unit=106,file='heights.out',iostat=ierr)
      write(106,'(A60)') infile
      write(106,604) 'stn long', 'stn lat','dummy ','displ u'
 604  format(4A10)
      if (makegrid.eq.0) then
         do 640 n = 1,numsta
            write(106,614) slong(n), slat(n), 0.0, displu(n),
     +           0,0,0
 640     continue
      else
         do 639 n = 1,numsta
            write(106,614) slong(n), slat(n), 0.0, displu(n),
     +           0,0,0, ' ', sname(n)
 639     continue
         
 614     format (2f10.4, 2f10.4, 3i3, 2a)
      endif
      
      close(unit=106)
      
      return
      end
      
************************************************************************
      subroutine calcnodes(flong, flat, strike, de, di, len, 
     +     wid, nodelong, nodelat, nodedep,
     +     nodex, nodey)
      

      implicit none

      integer numnodes
      parameter (numnodes = 4)
      
      include 'param.displace.h'
      
      real
     +     nodelong(*), nodelat(*), nodedep(*),
     +     deltax, deltay,
     +     whoriz, nodex(*), nodey(*),
     +     di, de, wid, len, strike
      
      
      if (di.lt.1) then
         whoriz = wid
      elseif ((di.gt.89).and.(di.lt.91)) then
         whoriz = 0
      elseif ((di.gt.179).and.(di.lt.181)) then
         whoriz = -wid
      elseif ((di.gt.269).and.(di.lt.271)) then
         whoriz = 0
      elseif (di.gt.359) then
         whoriz = wid
      else
         whoriz =  wid * cos(di * PI / 180) 
      endif
      
      deltay =  whoriz * sin(strike * PI / 180.0)
      deltax = -whoriz * cos(strike * PI / 180.0)
      
      nodex(1) = 0
      nodey(1) = 0
      nodex(2) = len * sin(strike * PI / 180.0)
      nodey(2) = len * cos(strike * PI / 180.0)
      nodex(4) = deltax
      nodey(4) = deltay
      nodex(3) = nodex(2) + nodex(4)
      nodey(3) = nodey(2) + nodey(4)

      nodedep(1) = 0.0 - de
      nodedep(2) = nodedep(1)
      nodedep(3) = nodedep(1) + wid * sin(di * PI / 180.0)
      nodedep(4) = nodedep(3)

      call transcoord(flong, flat, nodelong, nodelat, nodex,  nodey, 
     +             numnodes)

      return
      end

************************************************************************
 
      subroutine calctrace(flong, flat, strike, di, de, 
     +                     ftlong, ftlat, nodex, nodey)

      implicit none
 
      include 'param.displace.h'

      real 
     +     ftlong(*), ftlat(*), 
     +     ftx(2), fty(2), 
     +     nodex(*), nodey(*),
     +     deltax, deltay, whoriz,
     +     di, de, strike
      
      if (di.lt.1) then
         whoriz = 10000 
      elseif ((di.gt.89).and.(di.lt.91)) then
         whoriz = 0
      elseif ((di.gt.179).and.(di.lt.181)) then
         whoriz = 10000
      elseif ((di.gt.269).and.(di.lt.271)) then
         whoriz = 0
      elseif (di.gt.359) then
         whoriz = 10000
      else
         whoriz = de/tan(di * PI / 180.0) 
      endif
      
      deltax = -whoriz * cos(strike * PI / 180.0)
      deltay =  whoriz * sin(strike * PI / 180.0)
      
      ftx(1) = deltax
      fty(1) = deltay
      ftx(2) = nodex(2) + deltax
      fty(2) = nodey(2) + deltay
      
      call transcoord(flong, flat, ftlong, ftlat, ftx, fty, 2)
      return
      end
      
************************************************************************
      
      subroutine adddispl(deltan, deltae, deltau, 
     +     displn,disple,displu, numsta)
      

      implicit none

      include 'param.displace.h'
      
      integer
     +     n
      
      do 350 n = 1, numsta
         displn(n) = displn(n) + deltan(n)
         disple(n) = disple(n) + deltae(n)
         displu(n) = displu(n) + deltau(n)
 350  continue
      
      return
      end
      
************************************************************************
      
      subroutine initdispl(displn, disple, displu, numsta)
      

      implicit none

      include 'param.displace.h'
      
      integer
     +     n
      
      do 370 n = 1, numsta
         displn(n) = 0.0
         disple(n) = 0.0
         displu(n) = 0.0
 370  continue
      return
      end
      
************************************************************************
      subroutine geqmod(slg, slt, de, dn, du, hamp, vamp, edist)

*     Routine to increment the globk solution for the globk
*     earthquake model and to return output information


      implicit none

      include 'param.displace.h' 

* PASSED VARIABLES
      real  slt, slg, de, dn, du 
      real*8 hamp, vamp, edist

* LOCAL VARIABLES
      real*8 geod(3), apart, spos(3)

*     Compute the magnitude of horizontal components
      hamp = sqrt(de**2 + dn**2)
      vamp = abs(du)

*     Get the distance to the eartquake
*     Convert geodetic coordinates to XYZ
      geod(1) = pi/2 - slt*pi/180
      geod(2) = slg*pi/180
      geod(3) = 0
      call GEOD_to_XYZ( geod, spos )
      call eval_dist( spos, eq_pos, edist)

****  Now save the maxiumim distanace with 1 mm displacement
      if( (hamp.gt. 1.0 .or. vamp.gt.1.0 ) .and.
     .     edist.gt.grad ) then
           grad = edist
      end if

****  Now accumulate estimate of 1/r**2 dependence
      apart = (gdepth*1000.d0/edist)**2 
      nume = nume + 1
      neq = neq + apart**2
      bh = bh + hamp*apart
      bv = bv + vamp*apart
 
****  Thats all 
      return
      end
    


