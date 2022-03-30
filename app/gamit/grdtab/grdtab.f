      Program GRDTAB

c     Program to read a station list file or interpolate a global grid 
c     files for ocean and atmospheric loading, meterological data,
c     and mapping-function coefficients and to write these into the 
c     u-file for use by MODEL in processing a single session.  
c
c     Written by R. King and P. Tregoning  July 2006
c     Antecedents are utils/octtab.f and an earlier version of utils/grdtab.f, both
c     now obsolete
                                                           
c  **NOTE: Whenever the u-file format is changed, need to update version number (uver)
                               
c     The station-list file is ascii
c     The input grid file is binary, direct access file.   
c
c     A print output file echoing the header lines from the input files, is
c     'grdtab.out'
c
c     The output u-file will be ascii.  Version 1.0 has only ocean-loading tides 
c     and  was written by utils/ottab.f.  Version 2.0 has the format below:
c
c MOD TAH 120502: Changed output of non-tidal loading such that if the file is
c     not filtered (i.e, ' NCEP' rather than 'FNCEP') and the S1/S2 atmospheric
c     tide loading file is specified, then the S1/S2 tide contributions are removed
c     from the non-tidal loading values (so that they are non-tidal loading rather
c     than total loading values).  This should help in interpolation.
c     Changed requires that  Atmospheric (non-tidal) loading code be after  
c     Atmospheric tidal loading code (swapped from original grdtab). 
c MOD TAH 140320: Introduced new center of mass S1/S2 tide model based on analysis
c     of  3-hr loading.u-strasbg.fr/ITRF/CM timeseries IB model.  The rep_S1S2CoM
c     subroutine replaces the CoM model used in the IERS series with the new one
c     developed here.  Use of the new model removes most the diurnal and semidiurnal
c     variations seen in the raw 6-hr time series.  To invoke the new model add a +
c     the name of the atl.grid file (i.e, use  atl.grid+ as name)                 

c * Loading and mapping function values  for dpcnw3.198    
c *  GRDTAB Run on 2005/  2/  9  13:15:38 by rwk 
c * 
c * Values by keyword for each station; epochs are decimal day-of-year
c * Models are 8-characters, followed by two integers giving the number of columns and
c * rows to be read.  The values in the columns are hard-wired according to the 
c * model and number of columns.    
c *
c  Version 2.0 
c  STATION  PVER     PALOS VERDES    40403M001 VLBI  241.5965E   33.7438N    
c  OCEANLOD MODEL CSR4  CE 11  6 
c #           11 constituents   amp (up west south) and phase (deg) 
c #            M2     S2     N2     K2     K1     O1     P1     Q1     Mf     Mm    Ssa
c  OCEANLOD .00597 .00295 .00178 .00090 .01445 .00922 .00475 .00166 .00009 .00004 .00002
c  OCEANLOD .00172 .00018 .00033 .00002 .00374 .00231 .00122 .00041 .00005 .00002 .00000
c  OCEANLOD .00493 .00222 .00110 .00061 .00250 .00163 .00082 .00031 .00003 .00002 .00002
c  OCEANLOD  -28.7  -70.8  -59.6  -76.6   39.3   25.0   38.0   22.0  -80.8 -179.1    0.0
c  OCEANLOD -157.4 -125.3  165.8 -114.6 -131.6 -145.9 -132.8 -149.3  140.2 -169.0 -180.0
c  OCEANLOD   97.8  115.0   80.8  114.9 -174.7  172.9 -175.8  171.2   20.3   33.6 -180.0
c  ATMOSLOD MODEL NCEP CM 
c  ATMOSLOD 199.75     1.00094330    -0.10708474    -6.03734112
c  ATMOSLOD 200.00     0.71372092    -0.48598257    -5.20609474
c  ATMOSLOD 200.25     0.63802755    -0.76032323    -5.21443224
c  ATMOSLOD 200.50     1.39550257    -0.84154856    -5.40882158
c  ATMOSLOD 200.75     1.09849048     0.05170688    -5.66709757
c  ATMTDLOD MODEL NCEP CM  4  1
c #       S1 Ucos  Usin  Ncos  Nsin  Ecos Esin  S1 U cos  U sin  N cos  N sin  E cos E sin  
c  ATMTDLOD 10.12   8.32   2.41   1.12  
c  ATMTDLOD .00212 .00451 .00231 .00113 
c #    
c  ATMOSMET MODEL LOCALPTH 3  2
c  ATMOSMET 199.50 1028.05 25.0  36.0
c  ATMOSMET 200.00 1028.05 25.0  36.0
c #
c  ATMOSMAP MODEL VMF1     4  3
c #                 hydro a     wet a     hzen (mm) wet zen (mm)
c  ATMOSMAP 199.00 0.00122346  0.00053028  2239.1    42.3
c  ATMOSMAP 199.50 0.00122346  0.00053028  2239.1    42.3
c  ATMOSMAP 200.00 0.00122346  0.00053028  2239.1    42.3
c  ENDFILE
c 

c where the first column (# here) non-blank denotes a comment, which may be anywhere
c in the file.  The keyword STATION denotes the beginning of values for each station.
c The keyword MODEL is followed by an 8-character string, and integers indicating the
c number of columns and lines are to be read.  The file must terminate with the 
c string ENDFILE.
                
c There are two modes of running grdtab.  For the GAMIT batch file, the command
c line should read
c
c   grdtab  [d-file]  [yr] [doy] [d-file] [dur] [file1]  [file2] ... [file10]
c      
c     where [d-file] is the name of the d-file containing the list of stations to
c             include and the list of files indicates those that are to be read.
c           [yr] [doy] give the start day of the session
c           [dur] is the duration in days of the sesion.
c     The list can contain from 1 to 10 files, whose local (day-directory) names 
c     are hardwired but may be linked to specific files in another directory.
c       otl.list   ocean tidal loading station files, OSO/IGS format   
c       otl.grid   ocean tidal loading grid file (converted from OSO da files)
c       atml.list  non-tidal atmospheric loading    
c       atml.grid  non-tidal atmospheric loading
c       atl.list   atmospheric tidal loading 
c       atl.grid   atmospheric tidal loading
c       met.list   P, T, H for stations    
c       met.grid   P, T, H from a global circulation model (format not yet defined)
c       map.list   mapping function coefficients  
c       map.grid   mapping function coefficients (IMF, VMF1)
c            
c Also read by the program in this mode is a GAMIT coordinate file, with the name 
c extracted from the d-file. 

c The second mode is meant for stand-alone execution to get OTL components for
c a single station. ATML and met values are not coded for this mode.  Input 
c file is assumed to be only otl.grid, and the command line should read
c 
c   grdtab [lat] [lon] [ht] [site code]     or
c
c   grdtab [x] [y] [z] [site code]
c
c     where [lat] and [lon] are N and E latitude in decimal degrees, and
c              [ht] is ellisoidal height in meters 
c           [x] [y] [z] are Cartesian coordinates in meters
c           [sitecode] is an optional 4-character designation for the site 
c
c This mode outputs the OTL coefficients in u-file format (file named 
c 'ufile.[site code]' and NGS HARPOS format (file named harpos.[site code]'
c If the site code is omitted from the run string, it is given the generic
c name 'SITE'.

c  Both modes produce a documenting file named grdtab.out.


      implicit none     
                                   
      include '../includes/dimpar.h'      
      include '../includes/units.h'
      include '../includes/grdtab.h'   
      include '../includes/model.h'
             

c Values to be written on the u-file for each site ( common /ufcom/ in model.h)   
c
c    Model names to be passed through the c-file to the h-file
c      character*8 otidemod,atmlmod,atidemod,metmod,mapmod  
c    Numbers of components 
c      integer*4 notl,natml,natl,nmet,nmap
c    Numbers of epochs for time-dependent values (atml met map)
c      integer*4 ntatml,ntmet,ntmap 
c    Values and epochs for each model
c      real*4 otides,atml_time,atml_val,atides,met_time,met_val
c     .      ,map_time,map_val     
c    Flag indicating whether model is present
c      logical lotl,latml,latl,lmet,lmap
c
c      common /ufcom/  
c     .   otidemod,atmlmod,atidemod,metmod,mapmod  
c     . , otides(maxotl,6)
c     . , atml_time(maxatml), atml_val(maxatml,3)  
c     . , atides(maxatl,6)    
c     . , met_time(maxmet), met_val(maxmet,3)
c     . , map_time(maxmap), map_val(maxmap,8)
c     . , notl,natl,natml,nmet,nmap    
c     . , ntatml,ntmet,ntmap
c     . , lotl,latml,latl,lmet,lmap   
c   
c            
c   otides(i,6)                 amplitude U N E, phase U N E for i=1,notl tidal components
c   atml_time(i) atml_val(i,3)  doy, U N E for i=1,natml times
c   atides(i,6)                 cos sin for U N E for i=1,natl tidal components
c   met_time(i) met_val(i,3)    doy, U N E for i=1,nmet times
c   map_time(i),map_vel(i,4)    doy, mapping function coefficient for i=1,nmap coefficients


c Parameters from the input list and grid files ( commons of grdtab.h )
                
             
c Version numbers for input files 
c      real*4 otlgver,otllver,atmlgver,atmllver,atllver,atlgver
c     .     , metlver,metgver,maplver,mapgver
c      common/versions/ otlgver,otllver,atmlgver,atmllver,atllver,atlgver
c     .     , metlver,metgver,maplver,mapgver  

        
c Unit numbers - values set in data statement below (avoid block data 
c                because of linking ambiguities with gFortran 
c      integer*4 lud,luc,luotl,luatml,luatl,lumet,lumap,luprnt,luu,luh

c Model id for each list or grid file
c      character*8 otllmod,otlgmod,atmllmod,atmlgmod,atllmod,atlgmod
c     .          , metlmod,metgmod,maplmod,mapgmod
c      common/models/ otllmod,otlgmod,atmllmod,atmlgmod,atllmod,atlgmod
c     .          , metlmod,metgmod,maplmod,mapgmod
 
c Number of components available for each file
c     (For NAO ocean loading grid, this is the number returned by the
c      program, expanded using admittance from the number on the file)
c      integer*4 notll,notlg,natmll,natmlg,natll,natlg,nmetl,nmetg
c     .        , nmapl,nmapg
c      common/filecomp/ notll,notlg,natmll,natmlg,natll,natlg,nmetl,nmetg
c     .        , nmapl,nmapg
        
c Geographic array sizes and byte order for grid files     
c      integer*4 otlglat,otlglon, atmlglat,atmlglon
c     .        , atlglat,atlglon, metglat,metglon
c     .        , mapglat,mapglon    
c      logical otlgswap,atmlgswap,atlgswap,metgswap,mapgswap        
c      common/grids/ otlglat,otlglon, atmlglat,atmlglon
c      .           , atlglat,atlglon, metglat,metglon
c      .           , mapglat,mapglon    
c      .           , otlgswap,atmlgswap,atlgswap,metgswap,mapgswap

c Insititutional source of ocean loading files (OSO or NAO) - determines tides 
c present and order (may be redundant with number of components, but we 
c can't be sure)
c        source of values (OSO or NAO, determines format)         
c      character*3 otllsrc,otlgsrc
c      common/otlsource/otllsrc,otlgsrc
       
c---- Local variables---------------

c Session epoch and duration (days)
      integer*4 syear,sdoy
      real*4 sdur,map_val_save(maxmap)
      character*3 daynum                            

c File names (list, grid, and print file names are hardwired)

      character*10 fname,files(10),dfile,coordfile,ufile
      character*11 harposf

c Flags set true if the file is present
      logical otll,otlg,atmll,atmlg,atll,atlg,metl,metg,mapl,mapg
* MOD TAH 140320: Introduced change in CoM model for S1/S2
      logical atlcom   ! Set true to change CoM model in S1/S2 tide
                       ! Set by adding + to atl.grid file name
      data otll/.false./,otlg/.false./,atmll/.false./,atmlg/.false./
     .    ,atll/.false./,atlg/.false./,metl/.false./,metg/.false./
     .    ,mapl/.false./,mapg/.false./, atlcom /.false./

      character*4 CoMMod   ! either IERS or MIT when atl.grid+ used.
                 

c Station information

c   nstat                number of stations in D-file
c   sitcods(maxsit)      4-character station ids      
c   slat(maxsit)         latitudes of stations (decimal degrees)
c   slon(maxsit)         longitudes of stations(decimal degrees)
c   lstn                 logical true if site information in the station list file
c                        (no need to interpolate from grid)    
      integer*4 nstat   
      real*4 slat(maxsit),slon(maxsit),ell_ht(maxsit)
      character*4 sitcods(maxsit)                           
      logical lstn
       

c Version number for the u-file
      real*4 uver
      data uver/2.0/
                     
c User name, program version, and time of run
      character*16 uname 
      character*40 versn
      integer*4 iyear,imonth,iday,ihr,imn,isec,ihnsec    

c Double precision arguments for geoc_to_geod
      real*8 pi,site_epoch,xsite(3,maxsit),clatd
     .     , glatd,ht,crd1,crd2,crd3
               
      logical eol 
      character*1 latflag
      character*8 arg 
      character*16 eqrnfil
      character*256 message,harposmsg
      integer*4 ioerr,nsessn,nfiles,mode,i,j,k,m
      real*4 alat

c Harmonic identifiers for HARPOS file output (coded for OSO only)
      character*3 harfrq(11)
      data harfrq/'m2 ','s2 ','n2 ','k2 ','k1 ','o1 ','p1 ','q1 '
     .          ,'mf ','mm ','ssa'/
                    
c  Degrees to radians
      real*4 convd
      data convd/0.017453292/

c External functions
      integer*4 nblen,iarg,iclarg
      real*8 decyrs   
      logical fcheck      

c* rwk 131011: These moved to block data at bottom of this file
c      data lud/1/,luc/2/,luprnt/3/,luu/4/,luh/8/
c     .   , luotll/11/,luotlg/12/,luatmll/13/,luatmlg/14/,luatll/15/
c     .   , luatlg/16/,lumetl/17/,lumetg/18/,lumapl/19/,lumapg/20/ 
c*
      data eqrnfil/'eq_rename'/

c MOD TAH 120502: Added variables to support remove atm tidal load from
c     non-filtered non-tidal loading
      real*8 freq(2)   ! S1 and S2 frequencies (from etide.f code) (rad/sec)
      real*8 atmtidune(3)  ! Atm S1/S2 tidal contributions for UNE (system used here)
      real*8 dtatm     ! Time argument for S1/S2 tide (seconds from start of day)
                                                                                         
c--------------------------------------------------------------------------------------

c  Initialization and reading of file headers

      pi = 4.d0*datan(1.d0)  
      freq(1) =  7.27220521664304d-05     ! S1 (rads/sec)
      freq(2) = 14.54441043328608d-05     ! S2 (rads/sec)

      CoMMod = 'IERS'

c  Open the print file and write the program identifiers 

c     exit if a previous step has failed           
      if( fcheck('GAMIT.fatal') )
     .  call report_stat('FATAL','GRDTAB','grdtab',' '
     .                  ,'GAMIT.fatal exists: GRDTAB not executed',0) 
      open (unit=luprnt,file='grdtab.out',form='formatted'
     .     , status='unknown',iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL','GRDTAB','grdtab','grdtab.out',
     .  'Error opening grdtab print file: ',ioerr)
      endif   
      call gversn(versn)
      write(message,'(2a)') 'Program GRDTAB Version '
     .                       ,versn(1:nblen(versn))
      call report_stat('STATUS','GRDTAB','grdtab',' ',message,0) 
      write(luprnt,'(a)') message
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      call getusr(uname)
      write(message,'(a,i4,a,i3,a,i3,2x,i2,a,i2,a,i2,a,a16)')
     .   'GRDTAB Run on ',iyear,'/',imonth,'/',iday,ihr,':',imn,':'
     .     ,isec,' by ',uname         
      harposmsg = message
      call report_stat('STATUS','GRDTAB','grdtab',' ',message,0)
      write(luprnt,'(a)') message

c  Open the eq-rename file if present (not important, but avoids warning in /lib/lread.f)
                    
      eqrnfil = 'eq_rename'
      if( fcheck( eqrnfil ) ) then
c       this variable in common model.h
        iueqrn = 47 
        open(unit=iueqrn,file='eq_rename',STATUS='OLD',iostat=ioerr)
        if ( ioerr .ne. 0 ) then
          call report_stat('WARNING','MODEL','open',eqrnfil
     .                    ,'Error opening rename file: ',ioerr)
        endif   
      else
c       set unit number in model.h to zero to indicate not available
        iueqrn = 0   
      endif
    

c  Read first entry on the command line 

      iarg = iclarg(1,dfile)
      if( iarg.le.0 ) then 
        call report_stat('FATAL','GRDTAB','grdtab',' '
     .                   ,'Missing arguments in command-line',0)     
      else
c       stop if the d-file given on the command line is wrong or does not exist
        if( dfile(1:1).ne."d" ) then
c         assume numeric, Mode 2
          read(dfile,*,iostat=ioerr) crd1
          if( ioerr.ne.0 ) 
     .      call report_stat('FATAL','GRDTAB','grdtab',dfile
     .               ,'First argument not D-file name or numeric:',0)
          mode=2    
        else 
          mode=1
          if( .not.fcheck(dfile) ) 
     .     call report_stat('FATAL','GRDTAB','grdtab',dfile
     .                   ,'D-file not found:',0)   
          open(lud,file=dfile,iostat=ioerr,status='old')   
          if( ioerr.ne.0 ) then
            call report_stat('FATAL','GRDTAB','grdtab',dfile
     .                  ,'Error opening D-file:',ioerr)
          else
            call report_stat('STATUS','GRDTAB','grdtab',dfile
     .                  ,' Opened D-file:',ioerr)
          endif
        endif    
      endif           

c  Mode 1   
      if( mode.eq.1 ) then

c       Read the year, day, and session duration from the command line
           
        iarg = iclarg(2,arg) 
        if( iarg.le.0 ) then 
          call report_stat('FATAL','GRDTAB','grdtab',' '
     .                   ,'Missing second argument in command-line',0)
        else
          read(arg,*,iostat=ioerr) syear
          if( ioerr.ne.0 ) call report_stat('FATAL','GRDTAB','grdtab'
     .       ,' ','Error reading year from command line',ioerr)
        endif    
        iarg = iclarg(3,arg)
        if( iarg.le.0 ) then 
          call report_stat('FATAL','GRDTAB','grdtab',' '
     .                   ,'Missing third argument in command-line',0)
        else
          read(arg,*,iostat=ioerr) sdoy
          if( ioerr.ne.0 ) call report_stat('FATAL','GRDTAB','grdtab'
     .      ,' ','Error reading day-of-year from command line',ioerr)
        endif     
        iarg = iclarg(4,arg)
        if( iarg.le.0 ) then 
          call report_stat('FATAL','GRDTAB','grdtab',' '
     .                   ,'Missing fourth argument in command-line',0)
        else
          read(arg,*,iostat=ioerr) sdur
          if( ioerr.ne.0 ) call report_stat('FATAL','GRDTAB','grdtab'
     .    ,' ','Error reading session duration from command line',ioerr)
        endif     
        write(luprnt,'(/,2a )')  'D-file:                 ',dfile
        write(luprnt,'(a,2i4)' ) 'Date  :                 ',syear,sdoy
        write(luprnt,'(a,f3.1)') 'Session length (days) : ',sdur  
      
c   Mode 2

      elseif( mode.eq.2 ) then
        iarg = iclarg(2,arg)
        if( iarg.le.0 ) then  
          call report_stat('FATAL','GRDTAB','grdtab',' '
     .            ,'Mode 2 missing second argument in command-line',0)
        else
          read(arg,*,iostat=ioerr) crd2
          if( ioerr.ne.0 ) call report_stat('FATAL','GRDTAB','grdtab'
     .   ,' ','Error reading second coordinate from command line',ioerr)
        endif 
        iarg = iclarg(3,arg)
        if( iarg.le.0 ) then  
          call report_stat('FATAL','GRDTAB','grdtab',' '
     .            ,'Mode 2 missing third argument in command-line',0)
        else
          read(arg,*,iostat=ioerr) crd3
          if( ioerr.ne.0 ) call report_stat('FATAL','GRDTAB','grdtab'
     .    ,' ','Error reading third coordinate from command line',ioerr)
        endif 
        iarg = iclarg(4,arg)
        if( iarg.le.0 ) then
          sitcods(1) = 'SITE'
        else
          read(arg,'(a)',iostat=ioerr) sitcods(1)
          if( ioerr.ne.0 ) call report_stat('FATAL','GRDTAB','grdtab'
     .       ,' ','Error site code from command line',ioerr)
        endif 

      else
        call report_stat('FATAL','GRDTAB','grdtab',' '
     .               ,'Mode neither 1 nor 2',0)
      endif

   

c  For Mode 2, may skip the setting of input files, and reading of the 
c  d-file and l-file
            
      if( mode.eq.1 ) then      
        eol = .false.      
        i = 4                
        nfiles = 0
        do while( .not.eol )
          i = i + 1
          iarg = iclarg(i,fname) 
          if( iarg.le.0 ) then
            eol = .true.
          else         
            nfiles = nfiles + 1
            files(nfiles) = fname  
          endif
        enddo             
        do i=1,nfiles
* MOD TAH 140320: See if + adde to atl.grid file name to indicate
*         new CoM model should be used
          if( files(i)(1:nblen(files(i))).eq.'atl.grid+') then
              atlcom = .true.
              files(i) = 'atl.grid'  ! Reset name
              CoMMod = 'MIT'
              atlgmod(8:8) = 'M'     ! Assign last letter to name passed to globk
          endif
          if( files(i)(1:nblen(files(i))).eq.'otl.grid') otlg=.true.
          if( files(i)(1:nblen(files(i))).eq.'otl.list') otll=.true.
          if( files(i)(1:nblen(files(i))).eq.'atml.grid') atmlg=.true.
          if( files(i)(1:nblen(files(i))).eq.'atml.list') atmll=.true.
          if( files(i)(1:nblen(files(i))).eq.'atl.grid') atlg=.true.
          if( files(i)(1:nblen(files(i))).eq.'atl.list') atll=.true.
          if( files(i)(1:nblen(files(i))).eq.'met.grid') metg=.true.
          if( files(i)(1:nblen(files(i))).eq.'met.list') metl=.true.
          if( files(i)(1:nblen(files(i))).eq.'map.grid') mapg=.true.
          if( files(i)(1:nblen(files(i))).eq.'map.list') mapl=.true.
        enddo
        write(luprnt,'(a,10a10)') 'Files : ',(files(i),i=1,nfiles)
                                         

c         Read from the d-file the name of the coordinate file and the i-file 
c         (to create the u-file name);  first skip the solution number and check 
c         the number of sessions
        read( lud,'(1x)')
        read( lud,*) nsessn 
        if( nsessn.ne.1 ) call report_stat('FATAL','GRDTAB','grdtab'
     .    ,' ','GAMIT no longer supports multisessions',0)
        read( lud, '(a)' ) coordfile 
cd        print *,'from d coordfile ',coordfile
cd        print *,'coordfile kfflg ',coordfile,kfflg
        if( .not.fcheck(coordfile))            
     .    call report_stat('FATAL','GRDTAB','grdtab',coordfile
     .     ,'Coordinate file not available',0)   
        call crd_file_type( coordfile,kfflg )
c       first skip the tfile line
        read(lud,'(a)') ufile
        read(lud,'(a)') ufile   
        ufile(1:1) = 'u'
c*      This doesn't always work either since the ifile may (rarely) also have 
c*      6th char = a.  Put in a temporary trap.
        if( ufile(6:6) .eq. 'a' ) 
     .    call report_stat('FATAL','GRDTAB','grdtab',' ',
     .'Temporary logical bug in naming u-file from i-file--see Bob King'
     .    ,0)

c      Read the stations from the d-file and get their coordinates
                                                              
        rewind (lud,iostat=ioerr)
        if( ioerr.ne.0 )   
     .     call report_stat('FATAL','GRDTAB','grdtab',dfile
     .                     ,'Error rewinding d-file:',ioerr)
        call readdf('B',lud,sitcods,daynum,nstat)  
cd        print *,'GRDTAB aft readdf nstat ',nstat  
        do i=1,nstat
          call uppers(sitcods)
        enddo
        open (unit=luc,file=coordfile,form='formatted',status='old'
     .       ,iostat=ioerr)    
        if( ioerr.ne.0 ) then 
             call report_stat('FATAL','GRDTAB','grdtab',coordfile
     .                     ,'Error opening coordinate file:',ioerr)  
        else
           call report_stat('STATUS','GRDTAB','grdtab',coordfile
     .                     ,'Opened coordinate file:',ioerr)  
        endif                           
        iul = luc 
        do i=1,nstat
          site_epoch = decyrs( syear,sdoy,0.d0 )
          call lread( sitcods(i),site_epoch ) 
          do j=1,3
            xsite(j,i) = kpos(j) + kvel(j)*(site_epoch-kepoch0)
          enddo 
          call xyz2sph(xsite(1,i),latr,lonr,radius)  
          clatd = latr*180.d0/pi
c         convert to geodetic latitude
          semi=6378137.d0          
          finv = 298.257223563d0
          call geoc_to_geod( clatd,radius,glatd,ht,semi,finv)
          slat(i) = glatd   
          slon(i) = lonr*180.d0/pi   
          ell_ht(i) = ht
          if( slon(i).lt.0.0 ) slon(i) = slon(i) + 360.
        enddo
        
c   Mode 2
      else  
        otlg = .true.
        ufile = 'ufile.'//sitcods(1)
        harposf = 'harpos.'//sitcods(1) 
        nstat = 1   
c       need both Cartesian and gdetic coordinates for output
        if( dsqrt(crd1**2+crd2**2+crd3**2).lt.6.d6 ) then
c         input is geodetic, compute Cartesian  (geod_to_xyz expects double precision)
          call geod_to_xyz(crd1,crd2,crd3,xsite(1,1)) 
        else
c         input is Cartesian, compute geodetic
          xsite(1,1) = crd1
          xsite(2,1) = crd2
          xsite(3,1) = crd3          
          call xyz_to_geod(xsite(1,1),crd1,crd2,crd3)    
        endif
        slat(1) = crd1
        slon(1) = crd2
        ell_ht(1) = crd3
      endif
                                  

c  Open the input files and read the header information

      if( otll )  call rd_otl_list
      if( otlg )  call rd_otl_grid   
      if( atmll ) call rd_atml_list
      if( atmlg ) call rd_atml_grid    
      if( atll )  call rd_atl_list  
      if( atlg )  call rd_atl_grid       
      if( metl )  call rd_met_list      
      if( metg )  call rd_met_grid 
      if( mapl )  call rd_map_list  
      if( mapg )  call rd_map_grid
                
                                              
c  Open the u-file and write the header lines

      open (unit=luu,file=ufile,form='formatted', status='unknown'
     .     , iostat=ioerr )
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL','GRDTAB','grdtab','grdtab.out'
     .              ,'Error opening grdtab output file: ',ioerr)
      endif  
      write(luu,'(2a)') '* Loading and met values for ',dfile  
      write(message,'(3a,i4,a,i3,a,i3,2x,i2,a,i2,a,i2,a,a16)')
     .   'GRDTAB Version ',versn(1:nblen(versn)),' run on '
     .      ,iyear,'/',imonth,'/',iday,ihr,':',imn,':'
     .     ,isec,' by ',uname 
      write(luu,'(2a)') '* ',message  
      write(luprnt,'(a,/)') message 
      write(luu,'(a,f4.1)') ' Version ',uver
      write(luu,'(a)') '*'
      write(luu,'(a)') '* Values included for each station:'    
                   

c  Optionally, open a HARPOS file and write the header and harmonic lines
         
      if( mode.eq.2 ) then     
        open( unit=luh,file=harposf,form='formatted'
     .      ,status='unknown',iostat=ioerr)
        if(ioerr .ne. 0 ) then
          call report_stat('FATAL','GRDTAB','grdtab','grdtab.out'
     .              ,'Error opening grdtab HARPOS file: ',ioerr)
        endif              
        if( otlg.and.(otlgsrc.ne.'OSO')) 
     .    call report_stat('FATAL','GRDTAB','grdtab','grdtab.out'
     .              ,'HARPOS output coded only for OSO grid',0)
        if( otll.and.(otllsrc.ne.'OSO') ) 
     .   call report_stat('FATAL','GRDTAB','grdtab','grdtab.out'
     .              ,'HARPOS output coded only for OSO list file',0)
        call harpos_head( luh,harposmsg )
      endif
                    
c----------------------------------------------------------------------------------------
c     Report status that we removing tides
      if( (atmll.or.atmlg).and.(atll.or.atlg) .and. 
     .        atmlgmod(1:1).ne.'F' ) then
         message = 'Removing S1/S2 tides from non-filtered load'
         call report_stat('STATUS','GRDTAB','grdtab',' ',message,0)
         write(luprnt,'(a)') trim(message)
      end if


c----------------------------------------------------------------------------------------
c Loop over sites, reading the values from the available files                                      

      do i=1,nstat
                 
c      write the station line of the u-file 
       if( slat(i).ge.0.d0 ) then
         latflag='N'
         alat = slat(i)
       else
         latflag='S'
         alat= -slat(i)
       endif
       write(luu,'(a,a4,1x,f8.4,a1,1x,f8.4,a1)') 
     .     ' STATION ',sitcods(i),slon(i),'E',alat,latflag    
       write(luprnt,'(/,a4,1x,f8.4,a1,1x,f8.4,a1)') 
     .        sitcods(i),slon(i),'E',alat,latflag 
       if( mode.eq.2) write(luh,'(a,/,a,2x,a4,5x,3f14.4,1x,2f9.4,f7.1)')
     .    '#'
     .   ,'S',sitcods(i),(xsite(j,i),j=1,3),slat(i),slon(i),ell_ht(i)


c        Ocean tidal loading 
       if( otll.or.otlg ) then  
        lstn = .false.
        if( otll ) then
          call get_otl_list( slat(i),slon(i),lstn )  
          notl = notll                                                  
          if( lstn )  then
            write(luu,'(a,a8,i3,a3)') 
     .         ' OCEANLOD MODEL ',otllmod,notl,'  6'
            if( mode.eq.2 ) write(luh,'(a,a8)') '# OCEANLOD MODEL ',otllmod
          endif
        endif        
        if( .not.lstn ) then 
          if( otlg ) then      
            call get_otl_grid( slat(i),slon(i) ) 
            notl = notlg
            write(luu,'(a,a8,i3,a3)') ' OCEANLOD MODEL ',otlgmod,notl
     .            ,'  6'  
           if(mode.eq.2) write(luh,'(a,a8)') '# OCEANLOD MODEL ',otlgmod
          else                                 
            write(message,'(a,f9.5,2x,f9.5,a)') 
     .         'No site match within 10km for (',slat(i),slon(i)
     .         ,' ) in otl.list and otl.grid missing' 
            call report_stat('FATAL','GRDTAB','grdtab',' ',message,0)
          endif
        endif 
        if( notl.eq.11 ) then
          write(luu,'(2a)')  '#          '
     .       ,' 11 constituents   amp (up west south) and phase (deg)'
          write(luu,'(2a)')  '#             M2     S2     N2     K2'
     .        ,'     K1     O1     P1     Q1     Mf     Mm    Ssa'
c         for now, omit comment for NAO 54-component model
        endif   
        do j=1,3      
c        convert from m to mm for new u-file   
         write(luu,'(a,54f7.2)') ' OCEANLOD',(otides(m,j)*1.e3,m=1,notl)
        enddo
        do j=4,6
         write(luu,'(a,54f7.1)') ' OCEANLOD',(otides(m,j),m=1,notl)
        enddo   
       endif  
       if( mode.eq.2 ) then 
         do m=1,11
          write(luh,'(a,2x,a3,7x,a4,4x,2(2x,3f9.5))') 
     .      'D',harfrq(m),sitcods(i)
c            HARPOS is cos U E N, sin U E N (U W S for OSO and u-file)
     .     ,   otides(m,1)*cos(convd*otides(m,4))
     .     ,  -otides(m,2)*cos(convd*otides(m,5))
     .     ,  -otides(m,3)*cos(convd*otides(m,6))
     .     ,   otides(m,1)*sin(convd*otides(m,4))
     .     ,  -otides(m,2)*sin(convd*otides(m,5))
     .     ,  -otides(m,3)*sin(convd*otides(m,6))
         enddo
       endif 

c-------------------------------------------------------------------------
c      Atmospheric tidal loading 
       if( atll.or.atlg ) then
        lstn = .false.
        if( atll ) then
          call get_atl_list( sitcods(i),slat(i),slon(i),lstn )
          natl = natll                         
          if( lstn ) write(luu,'(a,a8,a,i3)') 
     .          ' ATMTDLOD MODEL ',atllmod,'  6',natl
        endif
        if( .not.lstn ) then
          if( atlg ) then
            call get_atl_grid( slat(i),slon(i),sdoy ) 
            natl = natlg
* MOD TAH 140320: If alternative CoM model used, replace last letter of
*           tide S1/S2 model
            if( atlcom ) atlgmod(8:8) = 'M'     ! Assign last letter to name passed to globk
            write(luu,'(a,a8,a,i3,1x,a4)') ' ATMTDLOD MODEL ',atlgmod,
     .        '  6',natl, CoMMod 
          else
            call report_stat('FATAL','GRDTAB','grdtab',' '
     .      ,'No site match in atl.list and atl.grid missing',0)
          endif
        endif
* MOD TAH 140320: Modify results to better represent the S1/S2 center
*       of mass variations.  Remove the IERS CoM and replace with values
*       from Strasberg 3-hr series (loading.u-strasbg.fr/ITRF/CM)
        if( atlcom ) call rep_S1S2CoM(atides, slat(i),slon(i) ) 
 
        write(luu,'(2a)')  '# S1 S2 '
     .             , '  U cos  U sin  N cos  N sin  E cos  E sin (mm)' 
        do k=1,2 
c PT060825: output in U N E rather than N E U as stored
          write(luu,'(a,6f7.2)') 
     .      ' ATMTDLOD',(atides(k,j),j=1,6)
        enddo           
       endif

c MOD TAH Moved non-tidal code after tidal code so that atides values
c     are available.

c-------------------------------------------------------------------------
c      Atmospheric (non-tidal) loading 
       if( atmll.or.atmlg ) then
        lstn = .false.
        if( atmll ) then  
          call get_atml_list( syear,sdoy,sdur,sitcods(i),slat(i),slon(i)
     .                      , lstn)
          natml = natmll 
          if( lstn ) write(luu,'(a,a8,2i3)') 
     .        ' ATMOSLOD MODEL ',atmllmod,natml,ntatml
        endif
        if( .not.lstn ) then   
          if( atmlg ) then
            call get_atml_grid( syear,sdoy,sdur,slat(i),slon(i) )
            natml = natmlg
            write(luu,'(a,a8,i3,i5)')
     .        ' ATMOSLOD MODEL ',atmlgmod,natml,ntatml
          else
            call report_stat('FATAL','GRDTAB','grdtab',' '
     .      ,'No site match in atml.list and atml.grid missing',0)
          endif
        endif
        write(luu,'(a)') '#          doy     U      N      E  (mm)'
        do k=1,ntatml 
c PT060825, explicitly order the output to be U N E
c MOD TAH 120502: See if we should remove tide from value before output.
c         We need atm S1/S2 tide be specified and non-smoothed non-tidal
c         loading.
          if( (atll.or.atlg) .and. atmlgmod(1:1).ne.'F' ) then
*             Tide being used and not filtered so remove tide from values
*             Care with order here (atml_val is NEU, atides is UNE)
              dtatm = mod(atml_time(k),1.0)*86400.0
              do m = 1,3
                 atmtidune(m) = 0.0d0
                 do j=1,2         
                     atmtidune(m) =   atmtidune(m)
     .                   + atides(j,(m-1)*2+1)*dcos(freq(j)*dtatm)
     .                   + atides(j,(m-1)*2+2)*dsin(freq(j)*dtatm)
                 enddo
              enddo                        
*             Now remove from values
              atml_val(k,3) =  atml_val(k,3) - atmtidune(1)  ! U
              atml_val(k,1) =  atml_val(k,1) - atmtidune(2)  ! N
              atml_val(k,2) =  atml_val(k,2) - atmtidune(3)  ! E
          endif
          write(luu,'(a,4f7.2)')  
     .      ' ATMOSLOD',atml_time(k),atml_val(k,3),(atml_val(k,j),j=1,2)
        enddo
       endif              


c-------------------------------------------------------------------------
c      Met data for delay 
       if( metl.or.metg ) then
        lstn = .false.
        if( metl ) then                     
          call get_met_list( syear,sdoy,sdur,sitcods(i),slat(i),slon(i)
     .                     , lstn )
          if( lstn ) write(luu,'(a,a8,i3,i2)') 
     .        ' ATMOSMET MODEL ',metlmod,nmet,' 3'     
        endif
        if( .not.lstn ) then
          if( atlg ) then
            call get_met_grid( syear,sdoy,sdur,slat(i),slon(i) )
            nmet = nmetg
            write(luu,'(a,a8,2i3)')
     .        ' ATMOSMET MODEL ',metgmod,nmet,ntmet 
          else
            call report_stat('FATAL','GRDTAB','grdtab',' '
     .      ,'No site match in met.list and met.grid missing',0)
          endif
        endif
        write(luu,'(a)') '#          doy  P (hPa) T (C)  RH(%)'
        do k=1,ntmet
          write(luu,'(a,f7.2,3f8.1)')
     .       ' ATMOSMET',met_time(k),(met_val(k,j),j=1,3)
        enddo
       endif 
 
c        Mapping Function coefficients
       if( mapl.or.mapg ) then
         lstn = .false.
         if( mapl ) then                     
           call get_map_list( syear,sdoy,sdur,sitcods(i),slat(i),slon(i)
     .                      , lstn )  
           if( lstn ) then
c PT060825: for the moment, only the VMF1 has a list file. We will be 
c           outputting 5 variables (in addition to the time stamp)  
c RWK060829: Which 5? I've coded for all eight but PT can change. 
c RWK070404: I've added a 9th value (orthometric hgt) for current VMF1 files
             if( maplver.eq.0.5 ) then
               write(luu,'(a,a8,2i3,a)') ' ATMOSMAP MODEL '
     .                 ,maplmod,nmapl,ntmap,' AH AW ZH ZW '
               write(luu,'(4a)') '#          doy ','  A hydro    A wet'
     .           ,'       dzen (mm)  wzen (mm) '
             elseif( maplver.eq.1.0 ) then
               write(luu,'(a,a8,2i3,a)') ' ATMOSMAP MODEL '
     .                 ,maplmod,nmapl,ntmap,' AH AW ZH ZW TM PR TP WP'
               write(luu,'(4a)') '#          doy ','  A hydro    A wet'
     .           ,'       dzen (mm)  wzen (mm)  Tmean (C)  Pres (hPa) '
     .           ,' Temp (C)  WV Pres (hPa) '   
             elseif( maplver.eq.1.1 ) then
               write(luu,'(a,a8,2i3,a)') ' ATMOSMAP MODEL '
     .                ,maplmod,nmapl,ntmap,' AH AW ZH ZW TM PR TP WP HT'
               write(luu,'(4a)') '#          doy ','  A hydro    A wet'
     .           ,'       dzen (mm)  wzen (mm)  Tmean (C)  Pres (hPa) '
     .           ,' Temp (C)  WV Pres (hPa)  Hgt (m)'   
             endif
           endif
           mapmod=maplmod 
           nmap = nmapl
         endif
         if( .not.lstn ) then  
           if( mapg ) then                   
c            trap the case where the first day is beyond the validity of the grid
* MOD TAH 140320: Check that data falls in range.  Leap years will over estimate
*            a little but should OK to verfiy end of file.
!            if( ((sdoy-.01).gt.mapg_end(2) )) then 
             if( syear+sdoy/365.25.lt. 
     .                       mapg_start(1)+mapg_start(2)/365.25 .or.
     .           syear+sdoy/365.25.gt. 
     .                           mapg_end(1)+mapg_end(2)/365.25 ) then
               write(message,'(a,2(f5.0,1x,f8.2,1x))') 
     .              'Requested day beyond grid: Range '
     .              ,mapg_start, mapg_end
               call report_stat('FATAL','GRDTAB','grdtab',' ',message,0)
             endif
             call get_map_grid( syear,sdoy,sdur,slat(i),slon(i)
     .          ,ell_ht(i) )
             write(luu,'(a,a8,2i3,1x,a)')
     .        ' ATMOSMAP MODEL ',mapgmod,nmapg,ntmap,'AH AW ZH'
             mapmod=mapgmod   
             nmap = nmapg   
           else
             call report_stat('FATAL','GRDTAB','grdtab',' '
     .       ,'No site match in map.list and map.grid missing',0)
           endif
           if( mapmod(1:3).eq.'VMF'.and.ntmap.gt.0 ) then
              write(luu,'(2a)') '#           doy  '
     .          ,'   A hydro     A wet    dzen (mm)'
            elseif( mapmod(1:3).eq.'IMF' ) then
              write(luu,'(2a)') '#           doy  '
     .                                ,'   a      b        c'
            endif  
         endif
         do k=1,ntmap  
c          The format assumes AH AW are the first two, others can be 
c          combinations of zenith delays, pressure, and temperatures   
c          if we have tried to read values beyond the end of the
c          valid grid, substitute the previous values for 12 hrs
c          and issue a warning; after that, stop the program,  This
c          works except for when the first value is beyond the grid 
c          end, but this case trapped above and the program stopped
* MOD TAH: Add year to calculation.  Probably note needed since tested above
*          at least for start time
           if( (syear*365+map_time(k)-(mapg_end(1)*365+mapg_end(2))).gt.
     .                                                     0.6 ) then
             write(message,'(2a,f7.2)') 'Requested time > 12 hrs after '
     .      ,'last map grid epoch',mapg_end(2)    
             call report_stat('FATAL','GRDTAB','grdtab',' ',message, 0)
* MOD TAH: Same mod as above
           elseif( (syear*365+map_time(k)-
     .                 (mapg_end(1)*365+mapg_end(2))).gt. 0.2 ) then
             write(message,'(a,f7.2,a)') 
     .         'Requested time 6-12 hrs after last map grid epoch'
     .          ,mapg_end(2),' substitute last valid value'
             call report_stat('WARNING','GRDTAB','grdtab',' '
     .          ,message, 0)           
             do j=1,nmap               
               map_val(k,j) = map_val_save(j)
             enddo
           endif
           write(luu,'(a,f7.2,2f11.8,3f11.1,3f11.2,f11.1)')
     .       ' ATMOSMAP',map_time(k),(map_val(k,j),j=1,nmap)     
           do j=1,nmap
             map_val_save(j) = map_val(k,j)
           enddo  
         enddo
c      endif on map files
       endif  
             
c     end loop on stations
      enddo

c    C'est tout

      write(luu,'(a)') ' ENDFILE'
      if(mode.eq.2) then
         write(luh,'(a)') 'HARPOS Format version of 2002.12.12'
         write(message,'(4a)') 'Normal stop in GRDTAB - created '
     .        , ufile,' and ',harposf    
      else
        write(message,'(2a)')'Normal stop in GRDTAB - created ',ufile
      endif
      call report_stat('STATUS','GRDTAB','grdtab',' ',message,0)
      stop
      end
     
c***************************************************************************************
  
C     BLOCK DATA for GRDTAB commons

      block data grdtab_bd

      include  '../includes/grdtab.h'

 
      data lud/1/,luc/2/,luprnt/3/,luu/4/,luh/8/
     .   , luotll/11/,luotlg/12/,luatmll/13/,luatmlg/14/,luatll/15/
     .   , luatlg/16/,lumetl/17/,lumetg/18/,lumapl/19/,lumapg/20/ 
                
      end


