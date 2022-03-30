c     Program to convert a set of gridded ascii ocean tide files produced by
c     Hans-Georg Scherneck (Onsala) to a single direct-access file to be used
c     by GAMIT program utils/octtab, which in turn interpolates the files
c     to a particular latitude and longitude.  

c     Written by R. King  24 April 2000
c     Modified by K. Matsumoto 10 Feb 2003
 
c     The NAO.99b global models comes as 21 files, one for each tidal
c     constituent, and each gridded at 0.5-degrees in latitude and
c     longitude.   A typical list is
c
c       nao-ptanM2.daf   nao-radiM2.daf       
c       nao-ptanS2.daf   nao-radiS2.daf    
c       nao-ptanN2.daf   nao-radiN2.daf      
c       nao-ptanK2.daf   nao-radiK2.daf      
c       nao-ptan2N2.daf  nao-radi2N2.daf      
c       nao-ptanMu2.daf  nao-radiMu2.daf       
c       nao-ptanNu2.daf  nao-radiNu2.daf      
c       nao-ptanL2.daf   nao-radiL2.daf    
c       nao-ptanT2.daf   nao-radiT2.daf
c       nao-ptanK1.daf   nao-radiK1.daf       
c       nao-ptanO1.daf   nao-radiO1.daf                
c       nao-ptanP1.daf   nao-radiP1.daf      
c       nao-ptanQ1.daf   nao-radiQ1.daf
c       nao-ptanJ1.daf   nao-radiJ1.daf       
c       nao-ptanM1.daf   nao-radiM1.daf      
c       nao-ptanOO1.daf  nao-radiOO1.daf                
c       nao-ptanMtm.daf  nao-radiMtm.daf                
c       nao-ptanMf.daf   nao-radiMf.daf                
c       nao-ptanMm.daf   nao-radiMm.daf                
c       nao-ptanSsa.daf  nao-radiSsa.daf                
c       nao-ptanSa.daf  nao-radiSa.daf                
c
c     where the 'radi' group have the (scalar) radial displacements and 
c     the 'ptan' group have the horizontal (tangential) potential from 
c     which the north and east components can be derived (in octtab) by
c     finite differencing.

c    The output file is direct access with record length 336 bytes.
c    Each record constains 4 x 21 complex numbers representing the radial
c    and tangential components for 21 tidal constituents, in the order
c    M2  S2  K1  O1  N2  P1  K2  Q1
c    M1  J1  OO1 2N2 Mu2 Nu2 L2  T2
c    Mtm Mf  Mm  Ssa Sa    

c    Before starting, the input .daf files are opened and ordered 
c    so that they match the order of the constituents in the output file.
      
      implicit none

c   Files
c     Input:  lulist= 1  List of GAMIT .oct file to write and .daf files to read 
c             ludaf=  2  Scherneck .daf files, opened one at a time
c     Output: luoct=  3  GAMIT .oct file (direct access), 
c             lprnt=  4  nao2da.out, print record
                      
      logical fcheck,eof
                          
      integer*4 lulist,ludaf,luoct,lprnt
      integer*4 nblen,mchkey,iyear,imonth,iday,ihr,imn,isec,ihnsec,ioerr
     .        , iarg,iclarg,idaf,ndaf,nlon,nlat,ix,icomp,itide,ii
     .        , ipointr,ilat,ilon,iline,nlines,irec,i,j,k
      integer*4 maxtid,maxdaf,maxz,maxlon,maxlat
      parameter(maxtid=21,maxlon=720,maxlat=360)
      parameter(maxdaf=2*maxtid,maxz=4*maxtid)
c      maxtid = number of tidal constituents 
c      maxdaf = number of input files = 2 * maxtid
c      maxz = complex radi (2) * complex ptan (2) * tangential * 21 constituents = 84
c      maxlon = 360 degrees / 0.5-degree interval
c      maxlat = 180 degrees / 0.5-degree interval
  
c       values from one line of an input file
      real*4 values(10)
c       grid for output direct-access file (file is grid-size + 1 record
      real*4 grid(maxz,maxlon,maxlat)
                                         
      character*1 ftyp
      character*3 atide,tides(maxtid)                   
      character*4 vartyp
      character*6 module
      character*16 uname 
      character*20 listfile,octfile,daffiles(maxdaf),filename
      character*40 versn   
      character*164 pad164
      character*256 message


      data lulist/1/,ludaf/2/,luoct/3/,lprnt/4/
      data tides/'M2 ','S2 ','K1 ','O1 ','N2 ','P1 ','K2 ','Q1 ',
     .           'M1 ','J1 ','OO1','2N2','Mu2','Nu2','L2 ','T2 ',
     .           'Mtm','Mf ','Mm ','Ssa','Sa '/
                                               

c     Set and write the program identifiers and open the print file

      module = 'NAO2DA'
c     exit if a previous step has failed           
      call gversn(versn)
      write(message,'(2a)') 'Program NAO2DA Version '
     .                       ,versn(1:nblen(versn))
      call report_stat('STATUS',module,'utils/nao2da',' ',message,0)
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      call getusr(uname)
      write(message,'(a,i4,a,i3,a,i3,2x,i2,a,i2,a,i2,a,a16)')
     .   ' NAO2DA Run on ',iyear,'/',imonth,'/',iday,ihr,':',imn,':'
     .     ,isec,' by ',uname
      call report_stat('STATUS',module,'utils/nao2da',' ',message,0)
      open (unit=lprnt,file='nao2da.out',form='formatted'
     .     , status='unknown',iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL',module,'utils/nao2da','nao2da.out',
     .  'Error opening nao2da print file: ',ioerr)
      endif
      write(lprnt,'(a,/)') message

c  Get the name of the list file from the command-line argument 
                      
      iarg = iclarg(1,listfile)
      if  (iarg .le. 0) then    
        call report_stat('FATAL',module,'utils/nao2da',' '
     .                   ,'Missing list file in command-line',0)
      else
c       stop if the list-file given on the command line does not exist
        if( .not.fcheck(listfile) ) 
     .   call report_stat('FATAL',module,'utils/nao2da',listfile
     .                   ,'List-file not found:',0)
      endif                      
      open(lulist,file=listfile,iostat=ioerr,status='old')   
      if( ioerr.ne.0 ) 
     .     call report_stat('FATAL',module,'utils/nao2da',listfile
     .                     ,'Error opening list-file:',ioerr)

       
c  Read the name of the direct-access output file and open it

      read(lulist,'(a)',iostat=ioerr) octfile  
      if( ioerr.gt.0 ) 
     .   call report_stat('FATAL',module,'utils/nao2da',listfile
     .          ,'Error reading output file name from list file',ioerr)
      open(luoct,file=octfile,status='new',access='direct'
     .    ,form='unformatted',recl=336,iostat=ioerr)
      if( ioerr.ne.0 ) 
     .     call report_stat('FATAL',module,'utils/nao2da',octfile
     .                 ,'Error opening output direct-access file',ioerr)
       write(lprnt,'(/,2a)') 'Output direct access file: ',octfile


c  Write the header line of the output file, defining its size
                   
      pad164 = ' ' 
      write(luoct,rec=1,iostat=ioerr) maxz,maxlon,maxlat,pad164  
      if( ioerr.ne.0 ) 
     .     call report_stat('FATAL',module,'utils/nao2da',octfile
     .       ,'Error writing 1st record to direct-access file',ioerr)
      write(lprnt,'(/,a,3i4)') 'First record of direct-access file: '
     .        ,maxz,maxlon,maxlat


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
     .     call report_stat('FATAL',module,'utils/nao2da',listfile
     .             ,'Error reading daf file name from list file',ioerr)
        if (ioerr.eq.-1 ) eof =.true.   
        if( .not.eof ) then
          ndaf = ndaf + 1
          daffiles(ndaf) = filename
          if( .not.fcheck(daffiles(ndaf)) ) 
     .    call report_stat('FATAL',module,'utils/nao2da',daffiles(ndaf)
     .                    ,'Input daf file not found:',ioerr)    
        endif
      enddo
                          

c  Open each daf file and move the values to the direct-access file
                                               
      write(lprnt,'(/,a)') 'Input file           Tide  #   nlon  nlat '
      write(lprnt,'(a)')   '------------------------------------------'
      do idaf = 1,ndaf
        open(ludaf,file=daffiles(idaf),status='old',iostat=ioerr)
        call report_stat('STATUS',module,'utils/nao2da',daffiles(idaf)
     .                 ,'Opened input grid file',ioerr)    
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL',module,'utils/nao2da',daffiles(idaf)
     .                 ,'Error opening input daf file',ioerr)
c       get grid spacing and tide type from first line of file
        read(ludaf,'(2i5,1x,a1,1x,a4,1x,a3)',iostat=ioerr)  
     .         nlon,nlat,ftyp,vartyp,atide     
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL',module,'utils/nao2da',daffiles(idaf)
     .                 ,'Error reading header of input file',ioerr)    
c       component-type not in file, must get from file name
        ix = nblen(daffiles(idaf))
        if( mchkey(daffiles(idaf),'radi',ix,4).gt.0 ) then
          icomp = 1
        elseif( mchkey(daffiles(idaf),'ptan',ix,4).gt.0 ) then
          icomp = 2
        else 
          icomp = 0
          call report_stat('FATAL',module,'utils/nao2da',daffiles(idaf)
     .       ,'Component type in file name neither radi nor ptan',0)
        endif      
c       compute the constituent index
        do i=1,maxtid
          if(atide.eq.tides(i)) then 
             itide = i
          endif
        enddo       
c       ipointr points to the slot in the tide index (1-44) of the grid array
c       immediately prior to where the component read from daf file should be written
        ipointr = (itide-1)*4
        write(lprnt,'(a20,2x,a3,1x,i3,2i6)') 
     .        daffiles(idaf),atide,itide,nlon,nlat   
        write(message,'(3a)') 'Writing ',daffiles(idaf),'values to grid'
        call report_stat('STATUS',module,'utils/nao2da',' ',message,0)
c       now read through the values and put them in the .oct file array 
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
              call report_stat('FATAL',module,'utils/nao2da',' '
     .                        ,'Something wrong in logic',0)
            endif   
            grid(ii,ilon,ilat) = values(i)
          enddo
c        if( iline.gt.150 ) stop
c       end of loop on lines
        enddo 
c     end of loop on files
      enddo     
      call report_stat('STATUS',module,'utils/nao2da',' '
     .       ,'Filled grid array',0)
    

c  Write the array, starting at record 2
                    
      irec = 1    
      do k=1,maxlat 
        do j=1,maxlon 
          irec = irec + 1
          write(luoct,rec=irec) (grid(i,j,k),i=1,maxz)
        enddo 
      enddo
      write(lprnt,'(/,a,i7)') 'Output file written with ',irec
     .                       ,' records'
      call report_stat('STATUS',module,'utils/nao2da',octfile
     .       ,'Wrote direct-access tide grid file',0)
                         
      stop
      end

      
                                             









