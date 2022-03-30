c     Program to convert a set of gridded ascii ocean tide files produced by 
c     a MatLab dump of a netCDF file (Mike Floyd) to a single direct-access file 
c     to be used  by GAMIT program grdtab, which in turn interpolates the files
c     to a particular latitude and longitude.  

c     Written by R. King  15 February 2018 from earlier program daf2da.f
                 
c     The original netCDF file for 2012 is toc_fes2012_grid_d359.nc. From
c     this Mike Floyd produced 216 ascii files, a sin and cos file (2) for
c     Up, North, and East (3) for each of 36 tidal constituents (waves) 
c     (14 primaryplus + or - side-lobes for some constituents).  Each ascii
c     file from the MatLab dump has a tab-delineated list by lon,lat of the
c     displacement in meters, e.g. for toc_fes2012_grid_d359-M2_cos_Up.dspl

c   0	    90	-0.00239999988116
c   0.25	90	-0.00239999988116
c   0.5	90	-0.00239999988116
c ....
c  359.75	90	-0.00239999988116
c  0   	89.75	-0.00251999986358
c  0.25	89.75	-0.00249999994412
c  0.5 	89.75	-0.00249999994412
                                                       
c Command-line:

c  toc2da [input-listfile] [output-file] [model-name] [N-S-order] [grid-size] [tides]
c
c   where 
c     input-listfile : name of fiie containing list of all the toc_ ascii files to be read
c     output-file    : name of the GAMIT binary direct-access file (e.g. otl_FES2014b.grid)
c     model-name     : 8-character name, e.g. 'FES2014b'
c     N-S order      : 'N' for north-to-sout, 'S' for south-to-north in the grid order
c     grid-size      : lon/lat divisions in degrees, e.g. 0.25           
c      Optional
c     tides          : 'all' if a complete set (36 for FES2012) otherwise only 11 primary 

c    The output file is direct access with record length 432 bytes. There is a header 
c    record comprising the model name (a8), a flag (1x,a1) for the latitude order 
c    ('N' or 'S'), the units of the displacement values (e6.0), the number of longitude, 
c    bands (i5), the  number of latitude bands (i5), the number of  tidal constituents 
c    (i4), and a list of the constituents (36(1x,a3));  e.g
c
c    FES2012 N 1.e-5 1440  721  36 M2  M2- S2  N2  N2- ....
c
c    The remaining records each have 6 x nv integer*2 values of the displacements, 
c    ordered sin U  cos U   sin E  cos E  sin N  cos N 
                             
      implicit none

c   Files
c     Input:  lulist= 1  List of GAMIT .otl file to write and toc_ files to read 
c             lutoc = 2  MatLab 'toc' files, opened one at a time
c     Output: luotl = 3  GAMIT .otl file (direct access), 
c             luprnt= 4  toc2da.out, print record
                      
      integer*4 lulist,lutoc,luotl,luprnt,reclen
      integer*4 iyear,imonth,iday,ihr,imn,isec,ihnsec,ioerr
     .        , iarg,ntoc,itoc,nlon,nlat,nval,nwave
     .        , icomp,iwave,islot,ilat,ilon,iline,nlines,icol,indx,irec
     .        , ii,i,j,k
      integer*4 maxwav,maxtoc,maxval,maxlon,maxlat  
c**  rwk 060830: maxlon maxlat need to be 2880 1441 for FES2004 but too large 
c                for g77; use the Intel version
      parameter(maxwav=36,maxlon=1440,maxlat=721)
      parameter(maxtoc=6*maxwav,maxval=4*maxwav)
c      maxwav = number of tidal constituents 
c      maxtoc = number of input files = 6 * maxwav
c      maxval = sin/cos (2) * U/E/N (3) * 36 constituents = 216 
c               sin/cos/(2) * U/E/N (3) * 11 constituents =  66 
c      maxlon = 360 degrees / 0.25-degree interval  = 1440
c      maxlat = 180 degrees / 0.25-degree interval  + 1   = 721
  
c       values from one line of an input file
      real*4 lon,lat,disp,gridsize 

c       units of the i*2 values on the output grid file (hardwired to 1.e-5)
      real*4 units

c       grid for output direct-access file (file is grid-size + 1 record
c       448 Mb for 36 constituents at 0.25-deg intervals
      integer*2 grid(maxval,maxlon,maxlat)

      character*1 ftyp,nsflag,wflag,comp
      character*2 atide       
      character*3 wave,waves(maxwav),trig
      character*4 vartyp,headfmt
      character*5 asize 
      character*8 otidemod
      character*16 uname 
      character*20 listfile,otlfile
      character*32 buf32
      character*40 versn,filename,tocfiles(maxtoc)   
      character*164 pad164
      character*256 message
cc      character*432 header 

      logical fcheck,eof 
      logical debug/.false./

c     functions       
      integer iclarg,mchkey,nblen,index

      data lulist/1/,lutoc/2/,luotl/3/,luprnt/4/
      data units/1.e-5/
      data waves/
     . 'M2 ','S2 ','N2 ','K2 ','K1 ','O1 ','P1 ','Q1 ','MF ','MM ','SSA'
     .,'J1 ','L2 ','M4 ','M6 ','M8 ','N4 ','S1 ','S4 ','T2'
     .,'2N2','J1+','K1-','K1+','K2+','LA2','M2-','MF+','MSF'
     .,'MTM','MU2','N2-','O1-','Q1-','SSA','   '/


c Set and write the program identifiers and open the print file

      call gversn(versn)
      write(message,'(2a)') 'Program TOC2DA Version '
     .                       ,versn(1:nblen(versn))
      call report_stat('STATUS','TOC2DA','toc2da',' ',message,0)
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      call getusr(uname)
      write(message,'(a,i4,a,i3,a,i3,2x,i2,a,i2,a,i2,a,a16)')
     .   ' TOC2DA Run on ',iyear,'/',imonth,'/',iday,ihr,':',imn,':'
     .     ,isec,' by ',uname
      call report_stat('STATUS','TOC2DA','ml2da',' ',message,0)
      open (unit=luprnt,file='toc2da.out',form='formatted'
     .     , status='unknown',iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL','TOC2DA','toc2da','toc2da.out',
     .  'Error opening toc2da print file: ',ioerr)
      endif
      write(luprnt,'(a,/)') message

c  Read the command line

      iarg = iclarg(1,listfile)
      if  (iarg .le. 0) then    
        call report_stat('FATAL','TOC2DA','toc2da',' '
     .                   ,'Missing list file in command-line',0)
      else
c       stop if the list-file given on the command line does not exist
        if( .not.fcheck(listfile) ) 
     .   call report_stat('FATAL','TOC2DA','toc2da',listfile
     .                   ,'List-file not found:',0)
      endif                      
      iarg = iclarg(2,otlfile)
      iarg = iclarg(3,otidemod)
      iarg = iclarg(4,nsflag)
      iarg = iclarg(5,asize)
      read(asize,'(f4.0)',iostat=ioerr) gridsize
      if(ioerr.ne.0 ) call report_stat('FATAL','TOC2DA','toc2da'
     .   ,' ','Error reading grid size from command-line',ioerr)
      iarg = iclarg(6,wflag)
      if( iarg.eq.0 ) then
         call report_stat('STATUS','TOC2DA','toc2da',' '
     .      ,'Use by default only 11 constituents',0) 
      elseif( wflag.eq.'all') then
         call report_stat('STATUS','TOC2DA','toc2da',' '
     .      ,'Use all (usually 36) constituents',0)       
      endif

c  Open the list file and the output grid file

      open(lulist,file=listfile,iostat=ioerr,status='old')   
      if( ioerr.ne.0 ) then
        call report_stat('FATAL','TOC2DA','toc2da',listfile
     .                   ,'Error opening list-file:',ioerr)  
      else
        call report_stat('STATUS','TOC2DA','toc2da',listfile
     .                     ,'Opened list-file',ioerr)
      endif
      if( wflag.eq.'all' ) then
        nwave = 36  
        nval = 216 
        reclen = 432
      else 
        nwave = 11   
        nval = 66
        reclen = 132       
      endif
      open(luotl,file=otlfile,status='unknown',access='direct'
     .    ,form='unformatted',recl=reclen,iostat=ioerr)
      if( ioerr.ne.0 ) then
        call report_stat('FATAL','TOC2DA','toc2da',otlfile
     .             ,'Error opening output direct-access file',ioerr)
      else
        call report_stat('STATUS','TOC2DA','toc2da',otlfile
     .                 ,'Opened output direct-access file',ioerr)
        write(luprnt,'(/,2a)') 'Output direct access file: ',otlfile
      endif

c  Set the number of lon/lat entries from the grid size
     
      nlon = int(360./gridsize)
      if( nlon.gt.maxlon ) then
         write(message,'(a,i4,a)') 'Number of lon values (',nlon
     .     ,') > maxlon (1440)'
         call report_stat('FATAL','TOC2DA','toc2da', ' '
     .                   ,  message,0)           
      endif
      nlat = int(180/gridsize) + 1
      if( nlat.gt.maxlat ) then
        write(message,'(a,i4,a)') 'Number of lat values (',nlat
     .     ,') > maxlat (721)'
         call report_stat('FATAL','TOC2DA','toc2da',' '
     .                   ,message,0)           
      endif

c  Write the parameters to the print file

      write(luprnt,'(9a,f5.1,3(a,i5))' ) 
     .  'List file: ',listfile,'   Grid file: ',otlfile
     . ,'Model: ',otidemod,'  NS flag: ',nsflag,'  Grid size: ',gridsize
     . ,'# lon : ',nlon,'  #lat :',nlat ,' # tides :',nwave

c  Initialize the output array
          
      do k=1,maxlat
        do j=1,maxlon
           do i=1,maxval
             grid(i,j,k) = 0
           enddo
        enddo
      enddo                             

                       
c  Read the names of the input files and see if they all exist
                                    
      eof = .false.
      ntoc = 0
      do while (.not.eof)             
        read(lulist,'(a)',iostat=ioerr) filename    
        if( ioerr.eq.-1 ) then
          eof = .true.
          if( ntoc.ne.216 ) then
            write(message,'(a,i4,a)') 'toc files (',ntoc
     .          ,') in listfile ne 216'
            call report_stat('WARNING','TOC2DA','toc2da',listfile
     .          ,message,0)    
          endif                   
        elseif( ioerr.ne.0 ) then 
           call report_stat('FATAL','TOC2DA','toc2da',listfile
     .           ,'Error reading toc file name from list file',ioerr)
        else
          ntoc = ntoc + 1 
          if( ntoc.gt.maxtoc ) then
            write(message,'(a,i2,a,i2,a)') 'Number of toc files (',ntoc
     .          ,') > maxtoc (',maxtoc,')'
            call report_stat('FATAL','TOC2DA','toc2da',' '
     .                      ,message,0)    
          endif
          tocfiles(ntoc) = filename 
          print *,'DEBUG ntoc tocfiles ',ntoc,tocfiles(ntoc)
          if( .not.fcheck(tocfiles(ntoc)) ) 
     .    call report_stat('FATAL','TOC2DA','toc2da',tocfiles(ntoc)
     .                    ,'Input toc file not found:',ioerr)    
        endif
      enddo
                          

c  Open each toc file and move the values to the direct-access file
                                               
      write(luprnt,'(/,2a)') 'Input file                             '
     .      ,' Wave  #   nlon  nlat '
      write(luprnt,'(2a)')   '---------------------------------------'
     .                      ,'---------------------------------------'
      do itoc = 1,ntoc
        open(lutoc,file=tocfiles(itoc),status='old',iostat=ioerr)
        call report_stat('STATUS','TOC2DA','toc2da',tocfiles(itoc)
     .                 ,'Opened input grid file',ioerr)    
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','TOC2DA','toc2da',tocfiles(itoc)
     .                 ,'Error opening input toc file',ioerr)
  
c      Get the constituent, component, and whether sin or cos from the file name
                                                     

c        Assume the file names have the form 
c            toc_fes2012_grid_d359-M2_cos_Up.dspl   or 
c            toc_fes2012_grid_d359-MF+_sin_East.dspl  
c                                
        if( mchkey(tocfiles(itoc),'sin',40,3).gt.0 ) trig= 'sin' 
        if( mchkey(tocfiles(itoc),'cos',40,3).gt.0 ) trig= 'cos' 

        if( mchkey(tocfiles(itoc),'Up',40,2).gt.0 )    comp= 'U' 
        if( mchkey(tocfiles(itoc),'East',40,4).gt.0 )  comp= 'E' 
        if( mchkey(tocfiles(itoc),'North',40,5).gt.0 ) comp= 'N' 
         
        icol = index(tocfiles(itoc),'-')
        read(tocfiles(itoc)(icol+1:icol+3),'(a3)',iostat=ioerr) wave
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL''TOC2DA','toc2da',tocfiles(itoc)
     .                 ,'Error reading wave from toc file name',ioerr)    
c       some wave names are only 2-characters long
        if(wave(3:3).eq.'_') wave(3:3) = ' '
c       now find the slot number for the wave according to the pre-determined order
        do i=1,nwave
          if(wave.eq.waves(i)) islot = i
        enddo                     
c       and the component index within the wave (i.e. sinU=1, cosU=2, sinE=3,....)
        if( trig.eq.'sin' ) then
          if( comp.eq.'U' ) icomp = 1
          if( comp.eq.'E' ) icomp = 3
          if( comp.eq.'N' ) icomp = 5
        elseif( trig.eq.'cos' ) then  
          if( comp.eq.'U' ) icomp = 2
          if( comp.eq.'E' ) icomp = 4
          if( comp.eq.'N' ) icomp = 6
        endif 
        
        write(luprnt,'(a40,2x,a2,1x,i3,2i6)') 
     .        tocfiles(itoc),wave,islot,nlon,nlat   
        write(message,'(3a)') 'Writing ',tocfiles(itoc),'values to grid'
        call report_stat('STATUS','TOC2DA','toc2da',' ',message,0)
c       now read through the values and put them in the .otl file array 
c       NEW CODE
c       # toc files have one value per line: lon  lat val 
        nlines = nlon*nlat
c              = 1,038,240 lines for 0.25 deg grid 
        ilon = 1
        ilat = 1
        do iline = 1,nlines
c         There should never be an EOF if the logic is correct
          read(lutoc,*,iostat=ioerr) lon,lat,disp
          if(debug.and.iline.eq.1 ) print *,'iline lon lat disp '
     .          ,iline,lon,lat,disp
          if(ioerr.ne.0 ) then
            write(message,'(a,i7)') 'Error reading line ',iline
            call report_stat('FATAL','TOC2DA','toc2da',tocfiles(itoc)
     .                        ,message,ioerr) 
          endif 
c         The order of the displacements for each wave is ( sin U  cos U sin E cos E sin N cos N )
c         and the waves follow the order in the data statement.  For 11 waves, there will be
c         66 values; for 36 waves, 216 values
          indx = 6*(islot-1)+1 + (icomp-1)
          if(debug) print *,'wave trig comp islot icomp indx '
     .                     , wave,trig,comp,islot,icomp,indx
          if(debug) print *,'lon lat ilon ilat ',lon,lat,ilon,ilat 
          grid(indx,ilon,ilat) = int(disp/units)
          if ( debug.and.iline.le.3 ) then
             print *,'iline,grid indx ilon ilat disp units float int '
     .       ,iline,indx,ilon,ilat,disp,units,disp/units,int(disp/units)
cd          else
cd          stop 
          endif 
          ilon = ilon + 1 
          if( (mod(ilon,1440)).eq.1 ) then
            ilat = ilat + 1
            ilon = 1
          endif 
c        if( iline.gt.150 ) stop
c       end of loop on lines
        enddo 
c     end of loop on files
      enddo     
      call report_stat('STATUS','TOC2DA','toc2da',' '
     .       ,'Filled grid array',0)
    
   
c  Write the header line of the output file

c      header = ' '                                
c      header(1:8) = otidemod
c      header(10:10) = nsflag
c      header(12:16) = '1.e-5'
c      write(header(18:21),'(i4)') nlon
c      write(header(23:26),'(i4)') nlat  
c      write(header(27:30),'(i3)') nwave
c      if( nwave.eq.11 ) then
c        write(header(32:75),'(11(a3,1x))') (waves(i),i=1,11)
c        write(luotl,rec=1,iostat=ioerr) header(1:132)       
c        write(luprnt,'(/,a/,a)') 
c     .         'Header of OTL grid file:',header(1:132)
c      elseif( nwave.eq.36) then
c        write(header(32:143),'(36(a3,1x))') (waves(i),i=1,36)
c        write(luotl,rec=1,iostat=ioerr) header(1:432)
c        write(luprnt,'(/,a/,a)') 
c     .      'Header of OTL grid file:',header(1:432)
c      endif    
       write(luotl,rec=1,iostat=ioerr) otidemod,nsflag,units,nlon,nlat
     .      ,nwave,(waves(i),i=1,nwave)
       write(luprnt,'(/,a,a8,1x,a1,1x,e6.1,1x,3i4,36(1x,a2))')
     .      'Header of OTL grid file:'
     .     , otidemod,nsflag,units,nlon,nlat,nwave
     .     , (waves,(i),i=1,nwave)

c     FES2012  N 1.e-5 1440  721  36 M2  M2- S2  N2  N2- ....
      if(ioerr.ne.0) call report_stat('FATAL','TOC2DA','toc2da',otlfile
     .     ,'Error writing 1st record to direct-access file',ioerr)

c  Write the array, starting at record 2
                    
      irec = 1    
      do k=1,nlat 
        do j=1,nlon 
          irec = irec + 1
          write(luotl,rec=irec) (grid(i,j,k),i=1,nval)
          if( irec.ge.1.and.irec.le.5 )
     .       write(*,'(a,4i5,2x,66i5)') 'irec j k nval val '
     .                        ,  irec,j,k,nval,(grid(i,j,k),i=1,nval)
        enddo 
      enddo
      write(luprnt,'(/,a,i7,a)') 'Output file written with ',irec
     .                       ,' records'
      call report_stat('STATUS','TOC2DA','toc2da',otlfile
     .       ,'Wrote direct-access tide grid file',0)
                         
      stop
      end

      
                                             









