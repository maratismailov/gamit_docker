      Subroutine rd_atl_grid

c     Read the header of an ANU-style or LU-style grid file for atmospheric tidal loading 

c     Written by R. King  10 August 2006 for ANU, LU added 31 January 2013
c     using P. Tregoning code from  old utils/grdtab.f.  
                            
      implicit none   

      include '../includes/grdtab.h'
                         

c Version 1.0 (ANU file) has constant coefficients, 3 header records, and
c    a record length of 48
c Version 2.0 (LU file ) has monthly coefficients, 1 header record, and
c    a record length of 288

c Input unit numbers from grdtab.h                                        
c   luprnt  print file  
c   luatlg  ATL grid file


c Output to grdtab.h      
c   atlgver   r*4      version number of grid file
c   atlgswap  logical  true if the grid file needs to be byte-swapped
c   atlgmod   c*8      model name 
c   atlglat   i*4      number of latitude values on grid
c   atlglon   i*4      number of longitude values on grid    
c   natlg     i*4      number of tidal components available by interpolating
c                      the grid.  Currently 2 (S1 and S2)
c   atlg_int  r*4      interval on time-dependent grid (~30 days for LU montly averages)
                                             
C Local variables                        

c  grecl  i*4  record length of direct-access file 
      integer*4 grecl    
c  ngtid   i*4  number of values per lat/lon on grid 
c                = 6 (UNE sin cos) x 2 components = 12
      integer*4 ngtid
c  ngval   i*4  number of values per record (same as ngtid for ANU, 144 for LU)
      integer*4 ngval 
c  nhead  i*4  number of header records (3 for ANU, 1 for LU)
      integer*2 nhead
      integer*4 ioerr            

c  integer version number from file
      integer*2 nversion

      character*3 ref_frame,source
      character*5 press_src,tide   
      character*8 atidemod 
                       
c  Open the grid file provisionally to determine whether it's an ANU-
c  style file with constant coefficients (version 1.0) or an LU-style
c  file with monthly coefficients (version 2.0), and whether byte-swapping
c  is needed for a BIG ENDIAN / little endian number storage scheme
                    
      atlgswap = .false.
      grecl = 2
      open(unit=luatlg,file='atl.grid',form='unformatted'
     .        ,access='direct',status='old',recl=grecl,iostat=ioerr)
      if(ioerr.ne.0)then
          call report_stat('FATAL','GRDTAB','rd_atl_grid','atl.grid'
     .           ,'Error opening atm tidal loading grid file:',ioerr)
      else
          call report_stat('STATUS','GRDTAB','rd_atl_grid','atl.grid'
     .                  ,'Opened atm tidal loading grid file:',ioerr)
      endif
      read(luatlg,rec=1,iostat=ioerr)  nversion
cd      print *,'nversion ',nversion
c     see if needs byte-swapping
      if(nversion.gt.100.or.nversion.lt.0) then
c       If the grid size does not seem be correct, swap the bytes and see if OK.   
        write(*,'(a,i4,a)') '**nversion =',nversion,' swapping**'
        call swap_bytes(2,nversion,1) 
cd        print *,'swapped nversion ',nversion
c       Now see if OK size.
        if(nversion.gt.100.or.nversion.lt.0) then 
          write(*,'(a,i4,a)') 'Invalid version (',nversion
     .            ,') of ATML grid file '
c         Still a problem so kill.   
          stop
        else
c         File seems to be swapped.  Set the swapped status and continue.
c         (We now need to swap every record read)  
          write(*,'(a)') 'Byte-swapping input grid file '
          atlgswap = .true.
        endif 
      endif     
c     close the file and reopen with the correct record length
      close(unit=luatlg)


c  Now read all the header records
                 
      if( nversion.eq.10 ) then
       grecl = 48
       open(unit=luatlg,file='atl.grid',form='unformatted'
     .        ,access='direct',status='old',recl=grecl,iostat=ioerr)
c         Record 1     
        read(luatlg,rec=1,iostat=ioerr)  
     .      nversion,press_src,ref_frame,source,tide
c            I*2        C*5       C*3      C*3   C*5 (these not used?)
        if( atlgswap ) call swap_bytes(2,nversion,1)
c         Record 2 is blank - left for future expansion 
c         Record 3  -  array limits
        read(luatlg,rec=3,iostat=ioerr) ngtid,atlglon,atlglat
        if( atlgswap )then
          call swap_bytes(4,ngtid,1)
          call swap_bytes(4,atlglon,1)
          call swap_bytes(4,atlglat,1)
        endif
        natlg = 2                     
        nhead = 3
        atlgver = nversion
        atlgver = atlgver/10.           
        write(luprnt,'(/,a,f5.1,a,3i6)') 'Opened ANU atl.grid version '
     .     ,atlgver,'   # val lat lon ',ngtid,atlglat,atlglon
        write(luprnt,'(a)') 'constant S1 and S2 tides only'
        write(luprnt,'(a,1x,a5,1x,a3)') 
     .   'Pressure source  frame ',press_src,ref_frame
        atlgmod(1:5) = press_src
        atlgmod(6:8) = ref_frame

      elseif( nversion.eq.20 ) then
c       it's an LU file with monthly values
        grecl = 288
        natlg = 2
        open(unit=luatlg,file='atl.grid',form='unformatted'
     .        ,access='direct',status='old',recl=grecl,iostat=ioerr)
c        Record 1      
c         version  nhead  model frame  ngval  nglat nglon 
c           I*2     I*2    C*8   C*3    I*4    I*4   I*4
        read(luatlg,rec=1,iostat=ioerr) nversion,nhead,atlgmod,ref_frame
     .      ,ngval,atlglat,atlglon
cd        print *,'swapped ',atlgswap 
cd        print *,'header ',nversion,nhead,atlgmod,ref_frame
cd     .    ,ngval,atlglat,atlglon
        if(atlgswap)then                   
          call swap_bytes(4,ngval,1)
          call swap_bytes(4,atlglat,1)
          call swap_bytes(4,atlglon,1)
        endif              
        atlgver = float(nversion)/10.
        write(luprnt,'(/,a,f5.1,/,a,3i5,2x,a)') 
     .     'Opened LU atl.grid version ',atlgver   
     .    , ' # val lat lon ',ngval,atlglat,atlglon
     .    , 'monthly coefficients for S1 and S2 tides '
      endif

      return
      end



