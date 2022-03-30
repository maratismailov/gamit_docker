c     Program to convert a set of gridded ascii atmospheric tidal loading files 
c     produced by Tonie VanDam (Luxembourg U) to a single direct-access file to be 
c     used by GAMIT program grdtab, which in turn interpolates the files to a 
c     particular latitude and longitude.  

c     Written by R. King  23 January 2013 from older OTL program daf2da

c     The LU global model comes as 12 files, one for each month of the year
c     and gridded at 1-deg intervals.  A typical list is 
c
c       s1_s2_def_apr_cm.dat  s1_s2_def_jan_cm.dat  s1_s2_def_may_cm.dat
c       s1_s2_def_aug_cm.dat  s1_s2_def_jul_cm.dat  s1_s2_def_nov_cm.dat
c       s1_s2_def_dec_cm.dat  s1_s2_def_jun_cm.dat  s1_s2_def_oct_cm.dat
c       s1_s2_def_feb_cm.dat  s1_s2_def_mar_cm.dat  s1_s2_def_sep_cm.dat
c
c     There are no headers on the files, only a list of numbers in the format
c     
c   0.000   90.000  0.0774  0.1079 -0.1076  0.0190  0.0529  0.6618  0.0286  0.0608 -0.6999 -0.0382 -0.0778  0.0280
c   0.000   89.000  0.0809  0.1007 -0.1143  0.0193  0.0545  0.6629  0.0273  0.0606 -0.7001 -0.0378 -0.0775  0.0275
c     ...
c
c     where the first two values are long and lat, and the remainder are the coefficients 
c     of CosS1, SinS1, CosS2, SinS2 for dR, dN, dE.
c
c     This program will read the 12 files and write a binary, direct-access file
c     of the form used for the OTL, ATML, and VMF1 grids, comprising three header
c     records:
c             
c  Header record
c     version  nhead  model frame  ngval  nglat nglon 
c       I*2     I*2    C*8   C*3   I*4    I*4   I*4
c
c    and data records comprising integer*2 values (units 1.e-4 mm) representing, for
c    each grid point (lat,lon) the 12 coefficients grouped by months, making the size 
c    of each record 2x12x12=288 bytes.  The grid order for the ascii input files is 
c    lon,lat, with latitude increasing more rapidly, but the grid file is written
c    with longitude increasing more rapidly to be consistent with the ATML,VMF1, and
c    OTL grid files.
c
c    To run atl_asc2bin, create an atl_list file which will be a command
c    line argument for atl_asc2bin.  The first line has the name of the 
c    output grid file, the model name (a8), and the frame (cm/ce/cf)
c                                       
c                           
c        RP130131_grid.atl  RP130131 CM   
c    col:  1-20              22-29   31-33   
c
c    The 12 monthly files should be in time-order in the list.
      
      implicit none

c   Files
c     Input:  lulist= 1  List of GAMIT .atl file to write and .dat files to read 
c             ludat=  2  VanDam .dat files, opened one at a time
c     Output: luatl=  3  GAMIT .atl file (direct access), 
c             lprnt=  4  atl_asc2bin.out, print record (not yet used)
                      

c        maxlon - number of longitude entries (nominally 0-360 at 1 deg inclusive = 361)
c        maxlat - number of latitude entries (nominally 90 to -90 at 1 deg = 181)
c        maxtid - number of tidal coefficients (nominally 2 cos + 2 sin)*3 = 12 )
c        maxdat - number of dat files or months of coefficients (12)
      integer*4 maxtid,maxdat,maxlon,maxlat,maxgval
      parameter(maxtid=12,maxdat=12,maxlon=361,maxlat=181)
      parameter(maxgval=maxtid*maxdat)

c       header values
      integer*2 version,nhead
      integer*4 nglat,nglon,ngval
      character*8 atidemod
      character*3 ref_frame 

c       values from one line of an input file
      real*4 values(12)
c       grid for output direct-access file (file is grid-size + 1 record)
      common grid(maxgval,maxlon,maxlat)
      integer*2 grid

      logical fcheck,eof
                          
      integer*4 lulist,ludat,luatl,lprnt
      integer*4 nblen,mchkey,ioerr,iday,imonth,iyear,ihr,imn,isec,ihnsec
     .        , iarg,iclarg,idat,ndat,ntid
     .        , irec,i,j,k
      real*4 interval

      character*1 ftyp,nsflag
      character*4 vartyp,headfmt
      character*16 uname 
      character*20 listfile,atlfile,datfiles(maxdat),filename
      character*33 buf33
      character*40 versn   
      character*261 pad261
      character*256 message


      data lulist/1/,ludat/2/,luatl/3/,lprnt/4/

c     Set and write the program identifiers and open the print file


c     exit if a previous step has failed           
      call gversn(versn)
      write(message,'(2a)') 'Program ATL_ASC2BIN Version '
     .                       ,versn(1:nblen(versn))
      call report_stat('STATUS','ATL_ASC2BIN','atl_asc2bin',' '
     .                 ,message,0)
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      call getusr(uname)
      write(message,'(a,i4,a,i3,a,i3,2x,i2,a,i2,a,i2,a,a16)')
     .  ' ATL_ASC2BIN Run on ',iyear,'/',imonth,'/',iday,ihr,':',imn,':'
     .     ,isec,' by ',uname
      call report_stat('STATUS','ATL_ASC2BIN','atl_asc2bin',' '
     .                 ,message,0)
      open (unit=lprnt,file='atl_asc2bin.out',form='formatted'
     .     , status='unknown',iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin'
     .    ,'atl_asc2bin.out','Error opening daf2da print file: ',ioerr)
      endif


c  Get the name of the list file from the command-line argument 
                      
      iarg = iclarg(1,listfile)
      if  (iarg .le. 0) then    
        call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',' '
     .                   ,'Missing list file in command-line',0)
      else
c       stop if the list-file given on the command line does not exist
        if( .not.fcheck(listfile) ) 
     .   call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',listfile
     .                   ,'List-file not found:',0)
      endif                      
      open(lulist,file=listfile,iostat=ioerr,status='old')   
      if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',listfile
     .                     ,'Error opening list-file:',ioerr)

       
c  Read the name of the direct-access output file (a20) and the GAMIT model ID 
      
      buf33 = ' '  
      atlfile = ' '   
      atidemod = ' '
      ref_frame = '   '                                      
      read(lulist,'(a)',iostat=ioerr) buf33
      if( ioerr.gt.0 ) 
     .   call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',' ' 
     .          ,'Error reading 1st line of list file',ioerr)
      read(buf33(1:20),'(a20)',iostat=ioerr) atlfile
      if( ioerr.gt.0 ) 
     .    call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',listfile
     .         ,'Error reading output file name from list file',ioerr)  
      read(buf33(23:30),'(a8)',iostat=ioerr) atidemod
      if( ioerr.gt.0 ) 
     .   call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',atidemod
     .         ,'Error reading model name list file',ioerr)  
      read(buf33(31:33),'(a3)',iostat=ioerr) ref_frame
      if( ioerr.gt.0 ) 
     .   call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',ref_frame
     .        ,'Error reading reference frame from list file',ioerr)  
        
        
c  Open the output grid file 

      open(luatl,file=atlfile,status='unknown',access='direct'
     .    ,form='unformatted',recl=288,iostat=ioerr)
      if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',atlfile
     .                 ,'Error opening output direct-access file',ioerr)
       write(lprnt,'(/,2a)') 'Output direct access file: ',atlfile


c  Initialize the output array
                
      do k=1,maxlat
        do j=1,maxlon
          do i=1,maxgval
            grid(i,j,k) = 0
          enddo
        enddo
      enddo                             


c  Initialize the array limits and header values to the nominals
      ntid = maxtid 
      ndat = maxdat
      ngval = ntid*ndat
      nglat = maxlat
      nglon = maxlon
      nhead = 1
      version = 20
                       
c  Read the names of the input files and see if they all exist
                                    
      eof = .false.
      ndat = 0
      do while (.not.eof)             
        read(lulist,'(a)',iostat=ioerr) filename    
        if( ioerr.gt.0 ) 
     .     call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',listfile
     .             ,'Error reading dat file name from list file',ioerr)
        if (ioerr.eq.-1 ) eof =.true.   
        if( .not.eof ) then
          ndat = ndat + 1 
          if( ndat.gt.maxdat ) then
            write(message,'(a,i2,a,i2,a)') 'Number of dat files (',ndat
     .          ,') > maxdat (',maxdat,')'
            call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',' '
     .                      ,message,0)    
          endif
          datfiles(ndat) = filename
          if( .not.fcheck(datfiles(ndat)) ) 
     .    call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin'
     .      ,datfiles(ndat),'Input dat file not found:',ioerr)    
        endif
      enddo
                          

c  Open each dat file and move the values to the direct-access file
                                               
c      write(lprnt,'(/,a)') 'Input file           Tide  #   nlon  nlat '
c      write(lprnt,'(a)')   '------------------------------------------'
      do idat = 1,ndat
        open(ludat,file=datfiles(idat),status='old',iostat=ioerr)
        call report_stat('STATUS','ATL_ASC2BIN','atl_asc2bin'
     .                  ,datfiles(idat),'Opened input grid file',ioerr)    
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin'
     .       ,datfiles(idat),'Error opening input file',ioerr)
        do k=1,nglon
          do j=1,nglat
            read(ludat,'(18x,12f8.0)',iostat=ioerr) (values(i),i=1,12)
            do i=1,ntid
              grid(ntid*(idat-1)+i,k,j) = nint(values(i)*1.e4)  
cd              if( ((j.eq.6.and.k.eq.4).or.(j.eq.1.and.k.eq.2)).and.
cd     .             (i.eq.1.or.i.eq.2) )  then 
cd               print *,'lon lat idat i index grid'
cd     .           , k,j,idat,i,ntid*(idat-1)+i
cd     .           ,grid(ntid*(idat-1)+i,k,j)    
cd              endif  
            enddo                                
cd            print *,'grid(1-12,1,1) ',(grid(i,1,1),i=1,12)
cd            print *,'grid(13-24,1,1) ',(grid(i,1,1),i=13,24)
cd            if( nglat.gt.1 ) stop
cd              if( idat.le.2 ) print *,'DEBUG k j  grid(13,4,6) '
cd     .              ,k,j,grid(13,4,6)
          enddo               
        enddo             
cd        print *,'DEBUG grid(13,4,6) ',grid(13,4,6)
        close(ludat)
c     end of loop on files 
      enddo                                     
cd      print *,'DEBUG grid(13,4,6) ',grid(13,4,6)
cd      print *,'DEBUG  lon 0 lat 90 Jan'
cd      write(*,*) (grid(i,1,1),i=1,12)
cd      print *,'DEBUG lon 3 lat 85 Feb '
cd      write(*,*) (grid(i,4,6),i=13,24)    
cd      print *,'DEBUG lon 200 lat -85 Feb '
cd      write(*,*) (grid(i,201,176),i=13,24)    


      call report_stat('STATUS','ATL_ASC2BIN','atl_asc2bin',' '
     .       ,'Filled grid array',0)   
     
    
   
c  Write the header line of the output file, defining its size
             
c     version  nhead  model  frame  ngval  nglat nglon 
c       I*2     I*2    C*8    C*3    I*4    I*4   I*4       27 bytes
c
      
      pad261 = ' ' 
      write(luatl,rec=1,iostat=ioerr) version,nhead,atidemod,ref_frame
     .                               ,ngval,nglat,nglon,pad261
      if( ioerr.ne.0 ) then
        call report_stat('FATAL','ATL_ASC2BIN','atl_asc2bin',atlfile
     .     ,'Error writing 1st record to direct-access file',ioerr)
      else                            
        write(message,'(a,2i2,a8,a3,3i4)') 
     .    'Record 1 of grid file: '
     .       ,version,nhead,atidemod,ref_frame,ngval,nglat,nglon
        call report_stat('STATUS','ATL_ASC2BIN','atl_asc2bin',' '
     .         ,message,0)
cx        write(lprnt,'(/,a)') message 
      endif

c  Write the array, starting at record 2
c     note that the lat/lon order is reversed from Tonie's files in order to be
c     consistent with our other grid files; i.e., latitude is the more slowly
c     changing value.  
     
      irec = 1      
      do k=1,nglat                         
        do j=1,nglon
          irec = irec + 1
          write(luatl,rec=irec) (grid(i,j,k),i=1,ngval)
cd          write(*,*) 'irec  grid(1-12,1,1) ',irec,(grid(i,1,1),i=1,12)
        enddo 
      enddo
      write(message,'(a,i7,a)') 'Output file written with ',irec
     .                       ,' records'
      call report_stat('STATUS','ATL_ASC2BIN','atl_asc2bin',atlfile
     .       ,message,0)
                         
      stop
      end

      
                                             









