      Subroutine filopn (versn )

c     Open all the files needed by arc except for the g-file, the scratch integration 
c     files, and the scratch yaw file

      implicit none
      include '../includes/dimpar.h'   
      include '../includes/units.h'
      include '../includes/arc.h'
               
      character*7 stat
      character*16 lowerc       
      character*40 versn
      character*80 line,srpgrdname
      character*256 message

      integer ioerr

      logical fcheck
         
c--------Set the file unit numbers-------------------------------------

      iscrn = 6
      ipole = 32
      iyawtmp = 9
      iyaw =  10 
c     Advanced UCL radiation pressure models require satellite attitude;
c     value reset to 0 below the binary yaw table is not available.
      if( modrad.gt.7 ) then
        iyawtab = 14
      else
        iyawtab = 0 
      endif    
      idebug = 11
      iut = 18
      iug = 29  
      iut1 = 30   
      inut = 31   
      ipole = 32 
      iotide = 33 
c     40+ reserved for scratch t-files, one per SV, set in arc and opened in oscrat 
c     lunar ephemeris  34 set and opened in lib/ephdrd (old code)
c     solar ephemeris  35 set and opened in lib/ephdrd (old code)
c     planetary ephemris 36 set and opened in lib/ephred (new code)
c     svnav.dat  69  set below and in /lib/svnav_read, not kept open 
c     UCL SRP grid files 99  names assigned here, opened in rdsrpgrd
c     UCL SRP fourier series files 5678 - names assigned and opened in srpfmod  


c----------Open the debug print file if requested--------------------
         
      if( idbprt.gt.0 ) then
        dbfname = afname
        dbfname(4:6) = 'bug'
        open(unit=idebug,file=dbfname,status='unknown'
     .      , form='formatted',iostat=ioerr)
        if( ioerr.ne.0 ) call report_stat('FATAL','ARC','arc'
     .         ,'adbfname','Error opening ARC debug print file',ioerr)
      endif  
                    

c----------Open the input files----------------------------------------------

c..... TAI-UT1 table 
      open (unit=iut1,file='ut1.',status='old',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','ARC','filopn','ut1.',
     .  'Error could not open UT1 table',ioerr)
      endif
c
c..... nutation table file  
c     open only if values for Sun, Moon, nutations are from soltab.,luntab., nutabl.
      if( .not.fcheck('nbody') ) then
        open (unit=inut,file='nutabl.',status='old',iostat=ioerr)
        if (ioerr .ne. 0) then
          call report_stat('FATAL','ARC','filopn','nutabl.',
     .    'Error could not open nutation table',ioerr)
        endif 
      endif
c
c..... pole table file
      open(unit=ipole,file='pole.',status='old',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','ARC','filopn','pole.',
     .  'Error could not open pole table',ioerr)
      endif
c
c..... svnav.dat file  -- now opened by svnav_read
c     Check for availability and format the svnav.dat file
c      open(unit=69,file='svnav.dat',status='old',iostat=ioerr)
c      if (ioerr .ne. 0) then
c        call report_stat('FATAL','ARC','filopn','svnav.dat',
c     .  'Error could not open svnav.dat file',ioerr)
c      endif
c      read(69,'(a80)',iostat=ioerr) line
c      read(69,'(a80)',iostat=ioerr) line
c      if (ioerr .ne. 0) then
c        call report_stat('FATAL','ARC','filopn','svnav.dat',
c     .  'Error reading svnav.dat file',ioerr)
c      endif

c......table of ocean tide coefficients  
      if( otidedeg.ne.0 ) then 
        otidfname = 'otides.dat'
        open (unit=iotide,file=otidfname,status='old',iostat=ioerr)   
        if( ioerr.eq.0 ) then
          call report_stat('STATUS','ARC','filopn',otidfname
     .                 , 'Opened ocean tide table ',0 )
        else 
          call report_stat('WARNING ','ARC','filopn',otidfname
     .                  ,'Error could not open ocean tide table',ioerr)
          iotide = 0 
        endif
      else
         iotide = 0 
      endif
               
c.... UCL SRP grid file (and read in the grids)
c      will have 6 grids for each satellite type (xyz bus and xyz panels)
       if (modrad.eq.8) then
        write(message,'(a)')'Opening and reading UCL SRP grid files'
        call report_stat('STATUS','ARC','filopn',' ',message,0)
C      x acc bus IIA i.e. iblock =3
c         srpgrdname= "/home/nejp5/gg/tables/UCLSRPbusX_3IIA_0.grd"
         srpgrdname= "UCLSRPbusX_3IIA_0.grd"
         call rdsrpgrd(srpgrdname,busX3xmin,busX3ymin,busX3cellsize
     .      ,busX3grid,bus3nrows,bus3ncols)
        print*,"busx3grid 1,1",busX3grid(1,1) 
C      y acc bus IIA i.e. iblock = 3
         srpgrdname= "UCLSRPbusY_3IIA_0.grd"
         call rdsrpgrd(srpgrdname,busY3xmin,busY3ymin,busY3cellsize
     .      ,busY3grid,bus3nrows,bus3ncols) 
        print*,"busy3grid 1,1",busY3grid(1,1) 
C      z acc bus IIA i.e. iblock = 3
         srpgrdname= "UCLSRPbusZ_3IIA_0.grd"
         call rdsrpgrd(srpgrdname,busZ3xmin,busZ3ymin,busZ3cellsize
     .      ,busZ3grid,bus3nrows,bus3ncols) 
        print*,"busz3grid 1,1",busZ3grid(1,1) 
C      x acc pnl IIA i.e. iblock = 3
C         srpgrdname= "/home/nejp5/gg/tables/UCLSRPpnlX_3IIA_0.grd"
         srpgrdname= "UCLSRPpnlX_3IIA_0.grd"
         call rdsrpgrd(srpgrdname,pnlX3xmin,pnlX3ymin,pnlX3cellsize
     .      ,pnlX3grid,pnl3nrows,pnl3ncols)
        print*,"pnlx3grid 1,1",pnlX3grid(1,1) 
C      y acc pnl IIA i.e. iblock = 3
         srpgrdname= "UCLSRPpnlY_3IIA_0.grd"
         call rdsrpgrd(srpgrdname,pnlY3xmin,pnlY3ymin,pnlY3cellsize
     .      ,pnlY3grid,pnl3nrows,pnl3ncols) 
        print*,"pnly3grid 1,1",pnlY3grid(1,1) 
C      z acc pnl IIA i.e. iblock = 3
         srpgrdname= "UCLSRPpnlZ_3IIA_0.grd"
         call rdsrpgrd(srpgrdname,pnlZ3xmin,pnlZ3ymin,pnlZ3cellsize
     .      ,pnlZ3grid,pnl3nrows,pnl3ncols) 
        print*,"pnlz3grid 1,1",pnlZ3grid(1,1) 
       endif

c.... input binary yaw file (if required)
             
      ytabfname = ' '
      if( iyawtab.gt.0 ) then  
        ytabfname = 'y'//gfname(2:16)
        ytabfname(6:6) = 't'
c       Check that the arc y-file name is not the same as the ytable file name
C       (yfname is the ascii file produced by ARC and used to get eclipse times.
c        It takes its 6th character from the t-file, so could be anything, but is 
c        usually the year or 'a', 'b', or 'c' for successive orbit iterations.  
c        ytabfile is a binary file of attitude used by MODEL, and potentially now an ARC repeat run.
c        Follow the same naming convention as yfname, but use 't' for 'table' as a reasonably safe 6th character.)
        if( yfname.eq.ytabfname ) then
           call report_stat('FATAL','ARC'
     .   ,'filopn',' ','t is 6th char of y-file and ytable-> same name'
     .   ,0)
        endif
        open( unit=iyawtab,file=lowerc(ytabfname),status='unknown'
     .    , form='formatted',iostat=ioerr)
        if (ioerr .ne. 0) then
           call report_stat('WARNING','ARC','filopn',ytabfname,
     .      'Error opening yawtable file, use nominal yaw',ioerr)
          iyawtab = 0 
        else
          write(iarh,'(a,a16)') ' Input yaw table file: ',ytabfname
        endif
      endif         


c---------Open the output files-----------------------------------

c.... t-file
      iut=18
      stat = 'unknown'
      call topens(tfname,stat,iut,ioerr)
      write(iarh,'(a,a16)')' Output ephemeris file: ',tfname

c.... output ascii yaw file             

      yfname = 'y'//gfname(2:16)
c     first check existence of yaw files
      open( unit=iyaw,file=lowerc(yfname),status='unknown'
     .    , form='formatted',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','ARC','filopn',yfname,
     .  'Error opening output yaw file',ioerr)
      endif
      open( unit=iyawtmp,status='scratch',form='formatted',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','ARC','filopn',' ',
     .  'Error opening temporary yaw file',ioerr)
      endif
      write(iyawtmp,'(2a)') 'PRN  YAW_RAT BIAS YR MO DA  HR MN ECLIPSE'
     .                    , '        BETA'
      write(iyawtmp,'(a)') '     NOMINAL yaw rates'
      write(iarh,'(a,a16)') ' Output yaw file: ',yfname

      return
      end
