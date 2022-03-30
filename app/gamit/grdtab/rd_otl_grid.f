      Subroutine rd_otl_grid
                                  
c     R. King 25 July 2004, modified for version 2 10 March 2018

c     Read the header of an ocean tidal loading direct-access grid file, 
c     either a Scherneck-style 'daf' file created by daf2da, or a version
c     2 file created by toc2da from acsii toc files written from a netCDF
c     grid via MatLab (M. Floyd).  

      implicit none     
                                   
      include '../includes/dimpar.h'
      include '../includes/model.h'
      include '../includes/grdtab.h'  

c Input unit numbers from grdtab.h                                        
c   luprnt  print file  
c   luotlg  OTL grid file


c Output to grdtab.h      
c   otlgver   r*4      version number of grid file
c   otlgswap  logical  true if the grid file needs to be byte-swapped
c   otlgmod   c*8      model name (last 2 char are cntr-of-figure CE or cntr-ofmass CM)
c   otlgsrc   c*3      institutional source of tides (OSO or NAO
c   otlglat   i*4      number of latitude values on grid
c   otlglon   i*4      number of longitude values on grid 
c   notlg     i*4      number of tidal components available by interpolating
c                      the grid.  For OSO, =11; for NAO, 21 are read from the
c                      grid but admittance is used to compute 33 minor 
c                      constituents, giving a total of 54.                   
c   otlgns    c*1      order of latitudes: 'N' for N-to-S, 'S' for S-toN
c   otlgunits r*4      units if i*2 values on version 2 grids 

c  Local variables                               

      integer*4 grecl
c   ngtid   i*4  number of values per lat/lon on grid 
c         OSO ver 1 complex amp/phs:  4 x #waves = 4 x 11 = 44   
c         NAO ver 1 complex amp/phs:  4 x #waves = 4 x 21 = 84
c         OSO ver 2 primary sin/cos UNE: 6 x 11 =  66
c         OSO ver 2 all sin/cos UNE    : 6 x 36 = 216
c   grecl i*4 record length of direct-access file  
c         OSO ver 1: real*4 x ngtid = 4 x 44 = 176 bytes
c         NAO ver 1: real*4 x ngtid = 4 x 84 = 336 bytes
c         OSO ver 2: integer*2 x ngtid = 2 x 66 = 132 bytes or 2 x 216 = 432 bytes
      integer*4 ngtid
      integer*4 ioerr,i
     
      logical debug/.false./      
                                                                        

c  Set byte-swapped status to false (will be checked below)
      otlgswap = .false.

c  Open the file with dummy recl to determine the type of grid 

      open(luotlg,file='otl.grid',status='old',access='direct'
     .        ,form='unformatted',recl=432,iostat=ioerr)
      if( ioerr.eq.0 ) then
        call report_stat('STATUS','GRDTAB','rd_otl_grid'
     .   ,'otl.grid','Open the OTL grid file',0)
      else
        call report_stat('FATAL','GRDTAB','rd_otl_grid'
     .   ,'otl.grid','Error opening the OTL grid file',ioerr)
      endif    
      read(luotlg,rec=1,iostat=ioerr) ngtid
      if(debug) print *,'RD_OTL_GRID checking first record ngtid=',ngtid
c     Use the header to determine source and hence the version and record length
c      Version 1  (ngtid otlglon otlglat otlgmod otlgns)
c         44 2880 1441 FES2004 N 
c      Version 2  (otlgmod otlgunits otlns otlglon otlglat ngtid waves)
c         FES2012 1.e-5 n 2880 1441 11 M2  S2 N2 ...       
c     If the first variable is an expected integer, assume version 1
      if( ioerr.ne.0 ) then
c       otlgver = 2.0 
      else 
c       check # of values for OSO and NAO respectively
        if( ngtid.eq.44.or.ngtid.eq.84 ) then 
          otlgver = 1.0          
        else
c         but may be byte-swapped 
          call swap_bytes(4,ngtid,1)
          if(debug) print *,'after byte-swap ngtid=',ngtid 
          if( ngtid.eq.44.or.ngtid.eq.84 ) then 
            otlgver = 1.0  
            otlgswap = .true. 
          else
            otlgver = 2.0 
          endif 
        endif   
        if( debug) print *,'2 otlgver ',otlgver
        if( otlgver.eq.1.0 ) then
          if(ngtid.eq.44) then
            otlgsrc = 'OSO'
            notlg = 11
            grecl = 176       
          elseif(ngtid.eq.84) then
            otlgsrc = 'NAO'
            notlg = 54
            grecl = 336
          endif                            
          read(luotlg,rec=1,iostat=ioerr) 
     .            ngtid,otlglon,otlglat,otlgmod,otlgns 
          if(otlgswap) then
            call swap_bytes(4,ngtid,1)
            call swap_bytes(4,otlglon,1)
            call swap_bytes(4,otlglat,1)  
            if(debug) print *,'rec 1 bytes swapped '
          endif
          if(debug)  print *
     .     ,'version 1  ngtid,otlglon,otlglat,otlgmod otlgns '
     .                , ngtid,otlglon,otlglat,otlgmod,otlgns
c         set defaults for old grids missing model name
          if( otlgmod.eq.'       ' ) then
            if( ngtid.eq.44 ) otlgmod = 'CSR4   E'
            if( ngtid.eq.84 ) otlgmod = 'NAO99b E'
          endif  
c         Need to translate some names to match the OSO names in otlcmc.dat
c         So far, all grids are CE 
          if( otlgmod(8:8).ne.'M' ) otlgmod(8:8)='E'
          if( otlgmod(1:5).eq.'FES04') otlgmod(1:7) = 'FES2004'
c         if order not on header, set this from the model name 
c         (all grids are N->S now except the original CSR3/4 and NAO99)
          if( otlgns.eq.' ' ) then
            if( otlgmod(1:3).eq.'CSR' .or. otlgmod(1:3).eq.'NAO' ) then
              otlgns = 'S'
            elseif (otlgmod(1:3).eq.'TPX'.or.otlgmod(1:3).eq.'FES'. or.
     .              otlgmod(1:3).eq.'GOT' ) then
              otlgns = 'N'
            else
              call report_stat('FATAL','GRDTAB','rd_otl_grid','otl.grid'
     .          ,'Cannot determine latitude order for o-tide grid',0) 
            endif
          endif 
        endif
      endif        

      if( otlgver.eq.2.0 ) then 
        if( debug) print *,'at ver 2'
c       1st variable not an integer, but rather the character model name 
        otlgver = 2.0 
        read(luotlg,rec=1,iostat=ioerr) otlgmod,otlgns,otlgunits
     .            ,otlglon,otlglat,notlg,(otlwaves(i),i=1,notlg)
        if(ioerr.ne.0) call report_stat('FATAL','GRDTAB','rd_otl_grid',
     .  'otl.grid','Error reading first record of version 2 file',ioerr)
        if(debug)  print *, 'version 2 '
     .        ,'otlgmod otlgns otlgunits otlglon otlglat notlg waves '
     .        , otlgmod,otlgns,otlgunits,otlglon,otlglat,notlg
     .        , (otlwaves(i),i=1,notlg)
c       record length = #waves x #values x i*2    
        grecl = notlg*6*2 
      endif 

      close(luotlg)
   
c  Reopen the file with the correct length for the data records
        
      if(debug) print *,'reopening ver ',otlgver,'with recl ',grecl                           
      open(luotlg,file='otl.grid',status='old',access='direct'
     .        ,form='unformatted',recl=grecl,iostat=ioerr)
      return
      end







