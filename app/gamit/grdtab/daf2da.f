c     Program to convert a set of gridded ascii ocean tide files produced by
c     Hans-Georg Scherneck (Onsala) to a single direct-access file to be used
c     by GAMIT program grdtab, which in turn interpolates the files
c     to a particular latitude and longitude.  

c     Written by R. King  24 April 2000

c     The Scherneck global models comes as 11 files, one for each tidal
c     constituent, and each gridded at 1-deg, 0.5-deg, or 0.25-deg in 
c     latitude and longitude.   A typical list is
c
c       csr4-ptanP1.daf   csr4-radiP1.daf      
c       csr4-ptanK1.daf   csr4-radiK1.daf       
c       csr4-ptanK2.daf   csr4-radiK2.daf      
c       csr4-ptanM2.daf   csr4-radiM2.daf       
c       csr4-ptanN2.daf   csr4-radiN2.daf      
c       csr4-ptanO1.daf   csr4-radiO1.daf                
c       csr4-ptanQ1.daf   csr4-radiQ1.daf
c       csr4-ptanS2.daf   csr4-radiS2.daf    
c       eqlt-ptanSsa.daf  eqlt-radiSsa.daf    
c       schw-ptanMf.daf   schw-radiMf.daf
c       schw-ptanMm.daf   schw-radiMm.daf
c
c     where the 'radi' group have the (scalar) radial displacements and 
c     the 'ptan' group have the horizontal (tangential) potential from 
c     which the north and east components can be derived (in grdtab) by
c     finite differencing.

c    Each ascii input file has a one-line header giving the number of
c    lon and lat values, an 'ftyp' variable ('M'), the representation 
c    of the entries (CMPX for complex), the constituent, and an integer
c    ('1' or '2').  In the files Scherneck sent in March and July of 2000, 
c    the format of the headers is
c
c  720  360 M CMPX M2          2   
c
c    with the '720' in columns 3-5, but in the files received in March 2005, 
c    the format is
c
c        720       361 M CMPX M2          1       
c
c    with the '720' in columns 8-10.  I'll distinguish the two formats by
c    the location of the 'CMPX' entry.     
c                             

c    To run daf2da, create a daf_list file which will be a command
c    line argument for daf2da.  The first line has the name of the 
c    output grid file, the model name, and a flag indicating whether
c    the input daf files are ordered S->N (older CSR3, CSR4 and NAO99)
c    or N->S (new CSR4, GOT00, TPX7, FES99, and FES2004).  The format is
c                                                          
c        FES2004_grid.otl    FES2004   N
c    col:  1-20              22-29     31 
c
c    The output file is direct access with record length 176 bytes.
c    Each record contains 4 x 11 complex numbers representing the radial
c    and tangential components for 11 tidal constituents, in the order
c    M2  S2  N2  K2  K1  O1  P1  Q1  Mf  Mm  Ssa    

c    Note: 2000-format files have the long-period tides designated with
c    a lower-case second character; 2005-format files used uppercase for
c    both.  Change the data list to use uppercase and convert the file header
c    constituent to uppercase before comparing.
          
c    Before starting, the input .daf files are opened and ordered 
c    so that they match the order of the constituents in the output file.
      
      implicit none

c   Files
c     Input:  lulist= 1  List of GAMIT .otl file to write and .daf files to read 
c             ludaf=  2  Scherneck .daf files, opened one at a time
c     Output: luotl=  3  GAMIT .otl file (direct access), 
c             lprnt=  4  daf_to_da.out, print record
                      
      logical fcheck,eof
                          
      integer*4 lulist,ludaf,luotl,lprnt
      integer*4 nblen,mchkey,iyear,imonth,iday,ihr,imn,isec,ihnsec,ioerr
     .        , iarg,iclarg,idaf,ndaf,nlon,nlat,nz,ntid,ix,icomp,itide
     .        , ipointr,ilat,ilon,iline,nlines,irec,ii,i,j,k
      integer*4 maxtid,maxdaf,maxz,maxlon,maxlat  
c**  rwk 060830: maxlon maxlat need to be 2880 1441 for FES2004 but too large 
c                for g77; use the Intel version
      parameter(maxtid=11,maxlon=1440,maxlat=721)
      parameter(ntid=11)
      parameter(maxdaf=2*maxtid,maxz=4*maxtid)
      parameter(nz=4*ntid)
c      maxtid = number of tidal constituents 
c      maxdaf = number of input files = 2 * maxtid
c      maxz = complex radi (2) * complex ptan (2) * tangential * 11 constituents = 44
c      maxlon = 360 degrees / 0.5-degree interval
c      maxlat = 180 degrees / 0.5-degree interval
  
c       values from one line of an input file
      real*4 values(10)
c       grid for output direct-access file (file is grid-size + 1 record
      common grid(maxz,maxlon,maxlat)
      real*4 grid
                                         
      character*1 ftyp,nsflag
      character*2 atide,tides(maxtid)                   
      character*4 vartyp,headfmt
      character*8 otidemod
      character*16 uname 
      character*20 listfile,otlfile,daffiles(maxdaf),filename
      character*32 buf32
      character*40 versn   
      character*164 pad164
      character*256 message


      data lulist/1/,ludaf/2/,luotl/3/,lprnt/4/
      data tides/'M2','S2','N2','K2','K1','O1','P1','Q1','MF','MM','SS'/
                                               

c     Set and write the program identifiers and open the print file

c     exit if a previous step has failed           
      call gversn(versn)
      write(message,'(2a)') 'Program DAF2DA Version '
     .                       ,versn(1:nblen(versn))
      call report_stat('STATUS','DAF2DA','daf2da',' ',message,0)
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      call getusr(uname)
      write(message,'(a,i4,a,i3,a,i3,2x,i2,a,i2,a,i2,a,a16)')
     .   ' DAF2DA Run on ',iyear,'/',imonth,'/',iday,ihr,':',imn,':'
     .     ,isec,' by ',uname
      call report_stat('STATUS','DAF2DA','daf2da',' ',message,0)
      open (unit=lprnt,file='daf2da.out',form='formatted'
     .     , status='unknown',iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL','DAF2DA','daf2da','daf2da.out',
     .  'Error opening daf2da print file: ',ioerr)
      endif
      write(lprnt,'(a,/)') message

c  Get the name of the list file from the command-line argument 
                      
      iarg = iclarg(1,listfile)
      if  (iarg .le. 0) then    
        call report_stat('FATAL','DAF2DA','daf2da',' '
     .                   ,'Missing list file in command-line',0)
      else
c       stop if the list-file given on the command line does not exist
        if( .not.fcheck(listfile) ) 
     .   call report_stat('FATAL','DAF2DA','daf2da',listfile
     .                   ,'List-file not found:',0)
      endif                      
      open(lulist,file=listfile,iostat=ioerr,status='old')   
      if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','DAF2DA','daf2da',listfile
     .                     ,'Error opening list-file:',ioerr)

       
c  Read the name of the direct-access output file (a20), the GAMIT model ID (a8),
c  and the latitude order (a1, 'N' or 'S')
      
      buf32 = ' '  
      otlfile = ' '   
      otidemod = ' '                                      
      nsflag = ' '
      read(lulist,'(a)',iostat=ioerr) buf32
      if( ioerr.gt.0 ) 
     .   call report_stat('FATAL','DAF2DA','daf2da',' ' 
     .          ,'Error reading 1st line of list file',ioerr)
      if (nblen(buf32).le.20 ) then
c       file name only
        read(buf32,'(a20)',iostat=ioerr) otlfile
        if( ioerr.gt.0 ) 
     .     call report_stat('FATAL','DAF2DA','daf2da',listfile
     .          ,'Error reading output file name from list file',ioerr)  
        call report_stat('WARNING','DAF2DA','daf2da',listfile
     .    ,'Model name and latitude order missing from command file',0)  
      elseif( nblen(buf32).le.24 ) then
c       file name plus model name  
        read(buf32,'(a20,1x,a8)',iostat=ioerr) otlfile,otidemod  
c       make some assumptions about the latitude order
        if( otidemod(1:3).eq.'CSR'.or.otidemod(1:3).eq.'NAO' ) then
          nsflag = 'S'
        else
          nsflag = 'N'
        endif
        if( ioerr.gt.0 ) 
     .     call report_stat('FATAL','DAF2DA','daf2da',listfile
     . ,'Error reading output file name and model from list file',ioerr)  
        write(message,'(a,a1)') 
     .      'Latitude order missing from command file, assume ',nsflag
        call report_stat('WARNING','DAF2DA','daf2da',listfile
     .      ,message,0)
      else
c       N-S flag also present
        read(buf32,'(a20,1x,a8,1x,a1)',iostat=ioerr) 
     .      otlfile,otidemod,nsflag  
        if( ioerr.gt.0 ) 
     .     call report_stat('FATAL','DAF2DA','daf2da',listfile
     . ,'Error reading output file name and model from list file',ioerr)   
      endif

        
c  Open the output grid file 

      open(luotl,file=otlfile,status='new',access='direct'
     .    ,form='unformatted',recl=176,iostat=ioerr)
      if( ioerr.ne.0 ) 
     .     call report_stat('FATAL''DAF2DA','daf2da',otlfile
     .                 ,'Error opening output direct-access file',ioerr)
       write(lprnt,'(/,2a)') 'Output direct access file: ',otlfile


c  Initialize the output array
          
      do k=1,maxlat
        do j=1,maxlon
           do i=1,maxz
             grid(i,j,k) = 0.
           enddo
        enddo
      enddo                             

                       
c  Read the names of the input files and see if they all exist
                                    
      eof = .false.
      ndaf = 0
      do while (.not.eof)             
        read(lulist,'(a)',iostat=ioerr) filename    
        if( ioerr.gt.0 ) 
     .     call report_stat('FATAL','DAF2DA','daf2da',listfile
     .             ,'Error reading daf file name from list file',ioerr)
        if (ioerr.eq.-1 ) eof =.true.   
        if( .not.eof ) then
          ndaf = ndaf + 1 
          if( ndaf.gt.maxdaf ) then
            write(message,'(a,i2,a,i2,a)') 'Number of daf files (',ndaf
     .          ,') > maxdaf (',maxdaf,')'
            call report_stat('FATAL','DAF2DA','daf2da',' '
     .                      ,message,0)    
          endif
          daffiles(ndaf) = filename
          if( .not.fcheck(daffiles(ndaf)) ) 
     .    call report_stat('FATAL','DAF2DA','daf2da',daffiles(ndaf)
     .                    ,'Input daf file not found:',ioerr)    
        endif
      enddo
                          

c  Open each daf file and move the values to the direct-access file
                                               
      write(lprnt,'(/,a)') 'Input file           Tide  #   nlon  nlat '
      write(lprnt,'(a)')   '------------------------------------------'
      do idaf = 1,ndaf
        open(ludaf,file=daffiles(idaf),status='old',iostat=ioerr)
        call report_stat('STATUS','DAF2DA','daf2da',daffiles(idaf)
     .                 ,'Opened input grid file',ioerr)    
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','DAF2DA','daf2da',daffiles(idaf)
     .                 ,'Error opening input daf file',ioerr)
c       get grid spacing and tide type from first line of file
c       there are (at least two possible formats--see comments at top of this file)
        read(ludaf,'(a30)',iostat=ioerr) buf32
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL''DAF2DA','daf2da',daffiles(idaf)
     .                 ,'Error reading header of input file',ioerr)    
        if( buf32(14:17).eq.'CMPX' ) then
          headfmt = '2000'
        elseif( buf32(24:27).eq.'CMPX' ) then
          headfmt = '2005'
        else
          call report_stat('FATAL','DAF2DA','daf2da',daffiles(idaf)
     .            ,'Cannot recognize header format of input file',ioerr)
        endif  
        rewind( ludaf )
        if( headfmt.eq.'2000') then
          read(ludaf,'(2i5,1x,a1,1x,a4,1x,a2)',iostat=ioerr)  
     .         nlon,nlat,ftyp,vartyp,atide   
          if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','DAF2DA','daf2da',daffiles(idaf)
     . ,'Error reading values from fmt-2000 header of input file',ioerr) 
        elseif ( headfmt.eq.'2005' ) then
          read(ludaf,'(2i10,1x,a1,1x,a4,1x,a2)',iostat=ioerr)  
     .         nlon,nlat,ftyp,vartyp,atide     
          if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','DAF2DA','daf2da',daffiles(idaf)
     . ,'Error reading values from fmt-2005 header of input file',ioerr) 
        endif       
        if( nlon.gt.maxlon ) then 
          write(message,'(a,i4,a,i4,a)') 'Number of lon values(',nlon
     .         ,') > maxlon (',maxlon,')' 
          call report_stat('FATAL','DAF2DA','daf2da',' ',message,0)
        elseif (nlat.gt.maxlat) then
          write(message,'(a,i3,a,i3,a)') 'Number of lat values(',nlat
     .        ,') > maxlat (',maxlat,')' 
          call report_stat('FATAL','DAF2DA','daf2da',' ',message,0)    
        endif
c       component-type not in file, must get from file name
        ix = nblen(daffiles(idaf))
        if( mchkey(daffiles(idaf),'radi',ix,4).gt.0 ) then
          icomp = 1
        elseif( mchkey(daffiles(idaf),'ptan',ix,4).gt.0 ) then
          icomp = 2
        else 
          icomp = 0
          call report_stat('FATAL','DAF2DA','daf2da',daffiles(idaf)
     .       ,'Component type in file name neither radi nor ptan',0)
        endif      
c       compute the constituent index
        do i=1,11
          call uppers(atide)
          if(atide.eq.tides(i)) then 
             itide = i
          endif
        enddo       
c       ipointr points to the slot in the tide index (1-44) of the grid array
c       immediately prior to where the component read from daf file should be written
        ipointr = (itide-1)*4
        write(lprnt,'(a20,2x,a2,1x,i3,2i6)') 
     .        daffiles(idaf),atide,itide,nlon,nlat   
        write(message,'(3a)') 'Writing ',daffiles(idaf),'values to grid'
        call report_stat('STATUS','DAF2DA','daf2da',' ',message,0)
c       now read through the values and put them in the .otl file array 
c         # lines in file is nlon x nlat x 2(complex) / 10 values per line
        nlines = 2*nlon*nlat/10                     
        ilat = 0
        ilon = 0
        do iline=1,nlines
c         every nlines/nlat (144) lines, increment latitude index and reset longitude
          if( mod(iline-1,nlines/nlat).eq.0 ) then
            ilat = ilat + 1
            ilon = 0
          endif
          read(ludaf,'(10e12.4)',iostat=ioerr) (values(i),i=1,10)
          do i=1,10
c           if odd, increment longitude 
            if( mod(i,2).ne.0 ) ilon = ilon +1    
c           radial are first two, tangential second two of each four elements
            if     ( icomp.eq.1 .and. mod(i,2).ne.0 ) then 
              ii = ipointr + 1
            elseif ( icomp.eq.1 .and. mod(i,2).eq.0 ) then 
              ii = ipointr + 2
            elseif ( icomp.eq.2 .and. mod(i,2).ne.0 ) then 
              ii = ipointr + 3
            elseif ( icomp.eq.2 .and. mod(i,2).eq.0 ) then 
              ii = ipointr + 4
            else
              ii = 0
              call report_stat('FATAL','DAF2DA','daf2da',' '
     .                        ,'Something wrong in logic',0)
            endif   
            grid(ii,ilon,ilat) = values(i)
          enddo
c        if( iline.gt.150 ) stop
c       end of loop on lines
        enddo 
c     end of loop on files
      enddo     
      call report_stat('STATUS','DAF2DA','daf2da',' '
     .       ,'Filled grid array',0)
    
   
c  Write the header line of the output file, defining its size
                   
      pad164 = ' ' 
      pad164(1:8) = otidemod  
      pad164(9:9) = nsflag
      write(luotl,rec=1,iostat=ioerr) nz,nlon,nlat,pad164  
      if( ioerr.ne.0 ) then
        call report_stat('FATAL','DAF2DA','daf2da',otlfile
     .     ,'Error writing 1st record to direct-access file',ioerr)
      else                            
        write(message,'(a,3i4,a)') 
     .    'Record 1 of grid file: ',nz,nlon,nlat,pad164(1:nblen(pad164))
        call report_stat('STATUS','DAF2DA','daf2da',' ',message,0)
        write(lprnt,'(/,a)') message 
      endif

c  Write the array, starting at record 2
                    
      irec = 1    
      do k=1,nlat 
        do j=1,nlon 
          irec = irec + 1
          write(luotl,rec=irec) (grid(i,j,k),i=1,nz)
        enddo 
      enddo
      write(lprnt,'(/,a,i7,a)') 'Output file written with ',irec
     .                       ,' records'
      call report_stat('STATUS','DAF2DA','daf2da',otlfile
     .       ,'Wrote direct-access tide grid file',0)
                         
      stop
      end

      
                                             









