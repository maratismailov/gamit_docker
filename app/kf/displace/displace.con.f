      program displace
C
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
C***********************************************************************
C
C234567890         21234567890         41234567890         61234567890
C        11234567890         31234567890         51234567890         712
C
C     MAIN 

      include 'param.displace.h'
      
      integer m, n
      
C     read the input file and find out whether to generate a grid of 
C     points, or whether specific stations are input.
      call readinput(flong, flat, fstrike, depth, dip, length, width, 
     +     u1, u2, scale, makegrid, numsta, numflts,
     +     xmin, xinc, xnum, ymin, yinc,ynum,slong, slat, 
     +     sname, origlong, origlat)
C     set up an xy grid or transform lat and long into xy coordinates
      if (makegrid.eq.0) then
         call creategrid(xmin, xinc, xnum, ymin, yinc, ynum, x, y)
      elseif (makegrid.eq.1) then
c     set up conical projection using murray's subroutine
         call conetrans(slong,slat,numsta,x,y,origlong,origlat)
cc     set up translation parameters for lambert conformable 
c     conical projection
c         call calctrans(origlong, origlat, maxlon, minlon, maxlat, 
c     +        minlat, parallel1, parallel2, slong(1), slat(1),
c     +        numsta)
cc     do the translation to xy coordinates in km
c         call lamberttrans(origlong, origlat, maxlon, minlon,
c     +        maxlat, minlat, parallel1, parallel2,
c     +        numsta, 0, slong, slat, x, y)
      else 
         write(*,*) 'makegrid must equal 0 or 1'
         stop 'makegrid must equal 0 or 1'
      endif
      
C     Initialize displacement arrays.
      call initdispl(displn, disple, displu, numsta)
C     Find displacement for each fault
      do 180 m = 1, numflts
C     find xy coordinates of fault origin using lambert
C     conformable conical projection
c         call lamberttrans(origlong, origlat, maxlon, minlon,
c     +        maxlat, minlat, parallel1, parallel2,
c     +        one, zero, flong(m), flat(m), fx, fy)
C     find xy coordinates of fault origin using murray's
C     polyconic projection
         write(*,*) 'Longitudes and latitudes:'
         write(*,'(f10.4,f10.4)') flong(m),flat(m)
         write(*,*) 'X and Y coordinates:'
         call poly44(origlong,origlat,flong(m),flat(m),
     .        fx,fy)
C     translate xy coordinates to fault origin
         call translatexy(fx,fy,x,y,tx,ty,numsta)
C     rotate xy coordinates 
         call rotatexy(fstrike(m),tx,ty,rx,ry,numsta)
C     find displacement at all rotated xy points
         call getdispl(rx, ry, depth(m), dip(m), length(m), width(m), 
     +        u1(m), u2(m), deltan, deltae, deltau, numsta)
         
C     remove the rotation of the displacements
         call unrotate(fstrike(m), deltan, deltae, numsta)
         
C     Add in this fault's displacements at each station
         call adddispl(deltan, deltae, deltau, 
     +        displn,disple, displu, numsta)
         
 180  continue
      
C     create an output file for use in gmt plotting routines
      call writeoutput(numflts, flong, flat, fstrike, depth, dip,
     +     length, width, u1, u2, scale, makegrid, numsta,  
     +     x, y, slong, slat, sname, displn, disple, displu,
     +     origlong, origlat)
      
      call writegmt( flong, flat, fstrike, depth, dip, length, 
     +     width, u1, u2, scale, makegrid, numflts, numsta, x, y, 
     +     slong, slat, sname, displn, disple, displu)
      
      stop
      end
      
************************************************************************

      subroutine readinput(flong, flat, fstrike, depth, dip, length, 
     +               width, u1, u2, scale, makegrid, numsta, numflts,
     +               xmin, xinc, xnum, ymin, yinc, ynum, slong, slat, 
     +               sname, origlong, origlat)

      include 'param.displace.h'
      character*256 infile, outfile, line

      integer
     +         m, n, ierr

      call get_runstring(infile, outfile)

C     Open up the input file.  Kill program if there is an error
      open(100, file=infile, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open', infile,
     +                  1,'readinput')

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
         read(line,*, iostat=ierr) depth(m), dip(m), length(m), 
     +        width(m), u1(m), u2(m), scale(m)

         if (ierr.eq.-1) then
            scale(m) = 'mm'
            ierr = 0
         endif
         
 150  continue
      
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
         read(line,*, iostat=ierr) makegrid, numsta
         if ((numsta.gt.MAXSTA).OR.(NUMSTA.LE.0)) then
            write(*,*) 'Displace: number of stations, ',numsta,
     +           ' exceeds maximum allowed, ',MAXSTA
         endif
         do 100  n=1,numsta
            read(100, '(A)', end=100) line
            read(line,*, iostat=ierr) slong(n), slat(n), sname(n)
 100     continue
      end if
      
C     Check for any read errors in the above statements.
      call report_error('IOSTAT',ierr,'read', infile,
     +     1,'readinput')
      
      close(100)
      
      call checkinput(flong, flat, fstrike, depth, dip, length, 
     +     width, u1, u2, scale, makegrid, numsta, numflts,
     +     xmin, xinc, xnum, ymin, yinc, ynum, slong, slat, 
     +     sname,infile)
      
      return
      end

************************************************************************

      subroutine creategrid(xmin, xinc, xnum, ymin, yinc, ynum, x, y)
      
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
      
      subroutine getdispl(rx, ry, de, di, len, wid, us, ud, 
     +     deltan, deltae, deltau, numsta)
      
C     compute displacement for each site using Okada method
      
      include 'param.displace.h'
      integer
     +     n
      
      real
     +     de, di, len, wid, us, ud
      double precision
     +     tempde, tempdi, templen, tempwid, tempus, tempud,
     +     temprx, tempry, tempde, tempdn, tempdu

      tempde = de
      tempdi = di
      templen = len
      tempwid = wid
      tempus = us
      tempud = ud

      do 500 n = 1, numsta
         temprx = rx(n)
         tempry = ry(n)
         call okada(temprx,tempry,tempde,tempdi,templen,tempwid,
     +        tempus,tempud,tempde, tempdn, tempdu)
         deltae(n) = tempde
         deltan(n) = tempdn
         deltau(n) = tempdu
 500  continue

      return
      end
      
************************************************************************
      
      subroutine writegmt( flong, flat, fstrike, depth, dip, length, 
     +     width, u1, u2, scale, makegrid, numflts, numsta, x, y, 
     +     slong, slat, sname, displn, disple, displu)
      
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
     +     numsta, makegrid)
      
      return
      end
************************************************************************

      subroutine writeoutput(numflts, flong, flat, fstrike, depth, dip,
     +     length, width, u1, u2, scale, makegrid, numsta, x, y, 
     +     slong, slat, sname, displn, disple, displu,
     +     origlong, origlat)
      
      include 'param.displace.h'
      character*256 infile, outfile
      
      integer
     +     ierr, n, m
      
      call get_runstring(infile, outfile)
      
C     Open the output file.
      open(300, file=outfile, iostat=ierr)
      call report_error('IOSTAT',ierr,'open', outfile,
     .     1,'writeoutput')
      
C     Write the fault characteristics to the file.
      write(300,440,iostat=ierr) 'flong', 'flat', 'fstrike', 'depth', 
     +     'dip','length',' width', 'u1', 'u2', 'scale'
      do 460 m = 1, numflts
         write(300,450,iostat=ierr) flong(m), flat(m), fstrike(m), 
     +        depth(m), dip(m), length(m), 
     +                             width(m), u1(m), u2(m), scale(m)
460   continue
      write(300,*,iostat=ierr)
440    format(10A8)
450    format(9f8.2, a8)

      if (makegrid.eq.0) then
        call transcoord(origlong, origlat, slong, slat, x, y, numsta)
        write(300, 405, iostat=ierr) 'x','y', 'longitude', 'latitude',
     +          'east','north','up'
        do 420 n = 1, numsta
          write(300, 400, iostat=ierr) x(n), y(n), slong(n), slat(n), 
     +          disple(n), displn(n), displu(n)
405       format(7a10)
400       format(4f10.3,3f10.4)
          call report_error('IOSTAT',ierr,'write', outfile,
     +                      1,'writeoutput')

420     continue

        else 
          write(300, 415, iostat=ierr) 'longitude', 'latitude','name',
     +          'east','north','up'
          do 430 n = 1, numsta
            write(300, 410, iostat=ierr) slong(n), slat(n), sname(n), 
     +          disple(n), displn(n), displu(n)
415         format(6a10)
410         format(2f10.3,a8,3f10.4)
            call report_error('IOSTAT',ierr,'write', outfile,
     +                      1,'writeoutput')
430       continue

        endif

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
     +               width, u1, u2, scale, makegrid, numsta, numflts,
     +               xmin, xinc, xnum, ymin, yinc, ynum, slong, slat, 
     +               sname, infile)

C     See if there are any problems with the input data:

      include 'param.displace.h'
      character*256 infile
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
     +     numsta, makegrid)
      
C     Write out locations and the displacement vectors
      
      include 'param.displace.h'
      
      integer
     +     n, ierr
      
      open(unit=104,file='vectors.out',iostat=ierr)
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
      subroutine calcnodes(flong, flat, strike, de, di, len, 
     +     wid, nodelong, nodelat, nodedep,
     +     nodex, nodey)
      
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
      
      deltay = whoriz * sin(strike * PI / 180.0)
      deltax = whoriz * -cos(strike * PI / 180.0)
      
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
      
      deltax = whoriz * -cos(strike * PI / 180.0)
      deltay = whoriz * sin(strike * PI / 180.0)
      
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

