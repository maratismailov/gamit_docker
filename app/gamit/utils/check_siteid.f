c     Program to compare 4-character IDs and coordinates between two lists and
c     identify different IDs for the same site and conflicting IDs for different sites
c     

c     R. King 3 November 2008; last modified 13 August 2010

c      Supports:  glist
c                 sopac colon-delmited (csv) file from Paul Jamason 
c                 norcal space-delimited file from http://csrc.ucsd.edu/projects/pgm/norcal2004.html
c                 sinex (NRCan)  
c                 magnet space-delimited from Geoff Blewitt      
c                 cart  id followed by free-format cartesian (INGV)
c                 -------
c                 idf   RWK-defined list with source identifier (for merging)  
c                       ID in columns 1-4, source in col 6-9, then lat, lon (free format)


c       The program is command-line driven with the form

c           check_siteid [file1] [ftype1] [src1] [file2] [ftype2] [src2] [mfile] [tol]

c        Required:

c          file1   Input file1 with site IDs and coordinates  
c          ftype1  Type of file (idf, glist, sopac, norcal, sinex, magnet, usgs)

c         Optional:

c           src1    3-character source of file1 (default 'MIT'; 
c                     blank to get from input idf file 
c           file2   Input file 2 with site IDs and coordinates; if omitted, 
c                     sites within file1 with will checked        
c           ftype2  Type of file (idf, glist, sopac, noral, sinex, magnet, usgs)
c           src2    3-characater source of file2 (default 'MIT'; 
c                    blank to get from input idf file)  
c           mfile   Output file merging the two input files                  
c           tol     Tolerance for match (default 10m) 

c     Examples:    check_siteid  pnw.glist  glist MIT
c                  check_siteid  pnw.glist  glist MIT  scec.glist glist SCE
c                  check_siteid  merged.idf idf  ' '   usgs.txt usgs USG


      implicit none
       
      integer*4 maxsit,num1,num2,numm,ifile1,ifile2,ifilem
     .        , iout,i,j,ioerr,iarg
      parameter( maxsit=15000) 
                 
      character*1 ans 
      character*3 source1(maxsit),source2(maxsit),sourcem(maxsit)
     .            ,src1,src2
      character*8 id1(maxsit),id2(maxsit),idm(maxsit) 
      character*6 ftype1,ftype2   
      character*13 span1(maxsit),span2(maxsit)
      character*24 file1,file2,filem,filemalph,filemlon     
      character*32 arg 
 
      real*8 lat1(maxsit),lat2(maxsit),lon1(maxsit),lon2(maxsit)
     .     , latm(maxsit),lonm(maxsit) 
     .     , convd,dtol,lattol,lontol
           

c* Functions
      integer*4 iclarg,nblen
           
      data convd/0.017453292519943296d0/


c*  Set units and initialize file names
      ifile1 = 1
      ifile2 = 2
      iout   = 3 
      ifilem = 4
      file1  = ' '    
      file2 =  ' ' 
      src1 = 'MIT'  
      src2 = 'MIT' 
      do i=1,maxsit
        id1(i) = ' '
        id2(i) = ' ' 
        idm(i) = ' ' 
        source1(i) = '   '
        source2(i) = '   '
      enddo
      filem = ' ' 
      filemalph = ' ' 
      filemlon = ' '             
                        
c*  Get the run-string arguments
 
      iarg = iclarg(1,file1)
      if( iarg.eq.0 ) then
        write(*,'(a)')
        write(*,'(2a)') 'check_siteid [file1] [ftype1] [src1] [file2]'
     .    ,' [ftype2] [src2] [mfile] [tol]'
        write(*,'(a)') 
        write(*,'(a)') 'Required: '
        write(*,'(a)') ' file1  Input file1 w/ site IDs and coordinates'   
        write(*,'(a)') ' ftype1 (see below)'
        write(*,'(/,a)') 'Optional:'
        write(*,'(2a)') '    srcid1  3-char source of file1 (default'
     .   ,' MIT; blank or omit to get from input idf file)'
        write(*,'(2a)') '    file2   Input file2 w/  IDs and coord; '
     .    ,'if omitted, sites within file1 within will checked'      
        write(*,'(a)') '    ftype2  (see below)'
        write(*,'(2a)') '    srcid2  3-char source of file2 (default'
     .    ,' MIT; blank or omit to get from input idf file)'
        write(*,'(2a)') '    mfile   output file merging the'
     .    ,' inputs (idf format)'         
        write(*,'(a)') '    tol     distance (m) tolerance default 10 m'
        write(*,'(/,a)') 'File types supported: '
        write(*,'(a)') ' idf glist sopac noral sinex magnet usgs'
        stop
      endif
      iarg = iclarg(2,ftype1) 
      if( iarg.eq.0.or.ftype1(1:1).eq.' ') then
        write(*,'(a)') 'File type required, stop '
        stop 
      endif 
      iarg = iclarg(3,arg)   
      if( iarg.ne.0 ) then
        src1 = arg(1:3)  
        do i=1,maxsit
         source1(i) = src1 
        enddo
      endif
      iarg = iclarg(4,file2)           
      iarg = iclarg(5,ftype2)  
      if( file2(1:1).ne.' '.and.iarg.eq.0 ) then
        write(*,'(a)') 'File type required for file2, stop'
        stop
      endif
      iarg = iclarg(6,arg) 
      if( iarg.ne.0 ) then
         src2 = arg
         do i=1,maxsit
           source2(i) = src2
         enddo
      endif        
      iarg = iclarg(7,filem)
      iarg = iclarg(8,arg)
      if( iarg.ne.0 ) then
        read(arg,*) dtol
      else
        dtol = 10.d0
      endif
      lattol = (dtol/6378137.d0)/convd
c     longitude tolerance gotten from this once latitude is known

             
* Read the site ids and coordinates into storage
           
      open(unit=ifile1,file=file1,status='old',iostat=ioerr) 
      if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening input file ',file1,ioerr   
        stop
      else
        write(*,'(a,a)') ' Opened input file ',file1 
      endif
      if( ftype1.eq.'idf   ' ) then
        call read_idf(ifile1,id1,lat1,lon1,source1,maxsit,num1)     
      elseif( ftype1.eq.'glist ' ) then
        call read_glist(ifile1,id1,lat1,lon1,span1,maxsit,num1)
      elseif( ftype1.eq.'cart  ') then 
        call read_cart(ifile1,id1,lat1,lon1,maxsit,num1)
      elseif( ftype1.eq.'sopac ') then 
        call read_sopac(ifile1,id1,lat1,lon1,maxsit,num1)  
      elseif( ftype1.eq.'norcal' ) then
        call read_norcal(ifile1,id1,lat1,lon1,maxsit,num1)  
      elseif( ftype1.eq.'sinex ' ) then
        call read_sinex(ifile1,id1,lat1,lon1,maxsit,num1)
      elseif( ftype1.eq.'magnet' ) then
        call read_magnet(ifile1,id1,lat1,lon1,maxsit,num1)
      elseif( ftype1.eq.'usgs  ' ) then
        call read_usgs(ifile1,id1,lat1,lon1,maxsit,num1)
      else
        write(*,*) 'File type 1 not valid ',ftype1                                              
        stop
      endif                        
      if( source1(1)(1:1).eq.' ' ) then
         do i=1,num1
           source1(i) = src1
         enddo
      endif  
      write(*,'(a,i5,a,a3)') ' #sites = ',num1,'  source=',src1    
      if( file2(1:1).eq.' ' ) then
         num2 = num1
         do i=1,num1
           id2(i) = id1(i)
           lat2(i) = lat1(i)
           lon2(i) = lon1(i)   
           source2(i) = source1(i)
           span2(i) = span1(i)
         enddo   
      else         
        open(unit=ifile2,file=file2,status='old',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening input file ',file2,ioerr   
          stop
        else  
          write(*,'(a,a)') ' Opened input file ',file2 
        endif
        if( ftype2.eq.'idf  ' ) then
          call read_idf(ifile2,id2,lat2,lon2,source2,maxsit,num2)
        elseif( ftype2.eq.'glist ' ) then
          call read_glist(ifile2,id2,lat2,lon2,span2,maxsit,num2)  
        elseif( ftype2.eq.'cart  ') then
          call read_cart(ifile2,id2,lat2,lon2,maxsit,num2)
        elseif( ftype2.eq.'sopac ') then 
          call read_sopac(ifile2,id2,lat2,lon2,maxsit,num2)
        elseif( ftype2.eq.'norcal' ) then
          call read_norcal(ifile2,id2,lat2,lon2,maxsit,num2)
        elseif( ftype2.eq.'sinex ' ) then
          call read_sinex(ifile2,id2,lat2,lon2,maxsit,num2)   
        elseif( ftype2.eq.'magnet' ) then
          call read_magnet(ifile2,id2,lat2,lon2,maxsit,num2)   
        elseif( ftype2.eq.'usgs  ' ) then
          call read_usgs(ifile2,id2,lat2,lon2,maxsit,num2)
        else
          write(*,*) 'File type 2 not valid ',ftype2
          stop
        endif     
        if( source2(1)(1:1).eq.' ' ) then
           do i=1,num2
             source2(i) = src2
           enddo
        endif   
       write(*,'(a,i5,a,a3)') ' #sites = ',num2,'  source=',src2    
      endif

cd      print *,'DEBUG file 1'
cd      do i=1,num1
cd         print *,i,id1(i),lat1(i),lon1(i),span1(i)
cd      enddo                 
cd      print *,'DEBUG file 2'
cd      do i=1,num2
cd         print *,i,id2(i),lat2(i),lon2(i),span2(i)
cd      enddo

                                                             
* Open the output report file 
                                    
      open(unit=iout,file='check_siteid.out',status='unknown'
     .     ,iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(2a)') 'Error opening output file ','check_sited.out'
     .    ,ioerr
        stop  
      endif    
      if( file2(1:1).eq.' ' ) src2 = '   '
      write(iout,'(9a)') 
     .  ' Files: ',file1(1:nblen(file1)),' (',src1,')  '
     .            ,file2(1:nblen(file2)),' (',src2,')'
      dtol = convd*lattol*6378137.d0
      write(iout,'(/,a,f7.2,a)') ' Distance tolerance is ',dtol,' m'
                 

* Replace any _GPS extents with blanks for output

      do i=1,num1
       if( id1(i)(5:8).eq.'_GPS' ) id1(i)(5:8)='    '
      enddo
      do i=2,num2
       if( id2(i)(5:8).eq.'_GPS' ) id2(i)(5:8)='    '
      enddo


* Check for duplicate IDs for different sites

      write(iout,'(/,a)') '* Duplicate IDs'
c      print *,'DEBUG num1 num2 ',num1,num2
      do i=1,num1
        do j=1,num2 
c         print *,'id1 id2 ',id1(i),id2(j)                     
c          print *,id1(i),lat1(i),lon1(i),id2(j),lat2(j),lon2(j)
          if( id1(i)(1:4).eq.id2(j)(1:4) ) then   
             lontol = lattol*dcos(lat1(i)*convd) 
             if( dabs(lat1(i)-lat2(j)).gt.lattol  .or.
     .         dabs(lon1(i)-lon2(j)).gt.lontol ) then 
                 write(iout,'(/,2x,a8,2x,a3,2f11.6,2x,a13)')
     .              id1(i),source1(i),lon1(i),lat1(i),span1(i)
                 write(iout,'(2x,a8,2x,a3,2f11.6,2x,a13)')
     .              id2(j),source2(j),lon2(j),lat2(j),span2(j)
             endif
           endif 
c          if( num2.gt.20 ) stop 20
         enddo                    
c         if( num1.gt.20 ) stop 20
      enddo        

* Check for duplicate coordinates

      write(iout,'(//,a)') '* Duplicate coordinates '
      do i=1,num1         
         lontol = lattol*dcos(lat1(i)*convd)
         do j=1,num2     
           if( dabs(lat1(i)-lat2(j)).lt.lattol  .and.
     .         dabs(lon1(i)-lon2(j)).lt.lontol ) then
             if( id1(i)(1:4).ne.id2(j)(1:4) ) then    
c                print *,'dupc ',id1(i),id2(j)
c     .             ,lat1(i),lat2(j),lon1(i),lon2(j),span1(i),span2(j)
                write(iout,'(/,2x,a8,2x,a3,2f11.6,2x,a13)')
     .              id1(i),source1(i),lon1(i),lat1(i),span1(i)
                 write(iout,'(2x,a8,2x,a3,2f11.6,2x,a13)')
     .              id2(j),source2(j),lon2(j),lat2(j),span2(j)
             endif 
           endif
         enddo
      enddo        
         
* Check for common sites
                       
      write(iout,'(//,a)') '* Common sites '
      do i=1,num1         
         lontol = lattol*dcos(lat1(i)*convd)
         do j=1,num2     
           if( dabs(lat1(i)-lat2(j)).lt.lattol  .and.
     .         dabs(lon1(i)-lon2(j)).lt.lontol ) then 
               write(iout,'(/,2x,a8,2x,a3,2f11.6,2x,a13)')
     .              id1(i),source1(i),lon1(i),lat1(i),span1(i)
               write(iout,'(2x,a8,2x,a3,2f11.6,2x,a13)')
     .              id2(j),source2(j),lon2(j),lat2(j),span2(j)
           endif
         enddo
      enddo   
             
* Write out a merged file if requested
           
      if( filem(1:1).ne.' ') then
        numm = num1 + num2
c        print *,'DEBUG num1 num2 numm ',num1,num2,numm
c        print *,'DEBUG2 source1(1,1841) ',source1(1),source1(1841)
        do i=1,num1
          idm(i) = id1(i)
          latm(i) = lat1(i)
          lonm(i) = lon1(i)       
          if( source1(i)(1:1).ne.' ') then
            sourcem(i) = source1(i)
          else
            sourcem(i) = src1
          endif
        enddo
        do i=1,num2 
          j = num1 + i   
          if( j.gt.maxsit ) then
            write(*,'(a,i6)') 'Site count > maxsit ',maxsit
            stop
          endif
          idm(j) = id2(i)
          latm(j) = lat2(i)
          lonm(j) = lon2(i)       
          if( source2(i)(1:1).ne.' ') then
            sourcem(j) = source2(i)
          else
            sourcem(j) = src2
          endif
        enddo  
c        print *,'calling write_merge maxsit numm ',maxsit,numm   
c        print *,'filem ',filem
c        print *,'1st file 1st,last ',idm(1),idm(num1)
c        print *,'2nd file 1st,last ',idm(num1+1),idm(num1+num2) 
c        print *,'DEBUG3 source1(1,1841) ',source1(1),source1(1841)  
c        print *,'DEBUG3 sourcem(1,1841) ',sourcem(1),sourcem(1841)
        call write_merge(ifilem,filem,idm,latm,lonm,sourcem,maxsit,numm) 
c        print *,'DEBUG4 sourcem(1,1841) ',sourcem(1),sourcem(1841)
        endif      

      write(*,*) 'Normal end '
      stop
      end                

c-------------------------------------------------------------------------------

      Subroutine read_idf(ifile,site,lat,lon,source,maxsit,numsit)     
         
      implicit none
                 
      integer*4 ifile,maxsit,numsit,ioerr,is,lastcol,nblen
     
      real*8 lat(maxsit),lon(maxsit)

      logical eof
                                          
      character*3 alat,alon,source(maxsit),src
      character*8 site(maxsit)
      character*256 line
                           
*     Format is simple: ID in columns 1-4, source in col 10-12 then lon, lat (free format)
                                 
      eof = .false. 
      is = 0 
      do while (.not.eof ) 
        read(ifile,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then 
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error reading line of IDF file ',ioerr
        else       
          is = is + 1
*         read the site id and source
          read(line(1:12),'(a8,1x,a3)',iostat=ioerr) site(is),src
          call uppers(site(is)) 
          call uppers(src)     
          if( source(is).eq.'   ' ) source(is) = src    
*         read lon and lat
          read(line(13:38),*,iostat=ioerr) lon(is),lat(is)
          if( lon(is).lt.0 ) lon(is) = lon(is) + 360.d0
        endif
      enddo
      numsit = is
      return
      end

c-------------------------------------------------------------------------------

      Subroutine read_glist(ifile,site,lat,lon,span,maxsit,numsit)     
         
      implicit none
                 
      integer*4 ifile,maxsit,numsit,ioerr,is,i
                                    
      real*4 yr1,yr2
      real*8 lat(maxsit),lon(maxsit)

      logical eof,found
           
      character*8 site(maxsit) 
      character*13 span(maxsit)
      character*256 line
       
*     Find the first line of the single-column coordinates section   

      eof = .false.    
      found = .false.
      do while (.not.eof.and..not.found ) 
         read(ifile,'(a)',iostat=ioerr) line 
         if( ioerr.ne.0 ) then
            write(*,*) 'Error looking for SUMMARY ', ioerr 
            if( ioerr.eq.-1 ) eof = .true.
         endif 
         if( line(3:9).eq.'SUMMARY' ) found = .true. 
      enddo        
c     skip the next two lines
      read(ifile,'(1x)')
      read(ifile,'(1x)') 

*     Read the values into storage 

      is = 0  
      eof = .false.
      do while (.not.eof)            
        line = ' ' 
        read(ifile,'(a)',iostat=ioerr) line  
        if( ioerr.eq.-1 .or. line(3:5).eq.'  ' ) then
           eof = .true.
        elseif ( ioerr.ne.0 ) then 
           write(*,*) 'Error reading line, ioerr: ',ioerr
        else
          is = is + 1   
          if( is.le.maxsit) then   
c            Sometime between Nov 2004 and Sep 2005, the glist format was
c            changed to add a third decimal place to the start/stop times 
c**             read(line,'(1x,2f9.4,14x,2f9.3,8x,a8)',iostat=ioerr) 
cc** rwk 161027 fixed for current format           
            read(line,'(1x,2f9.4,15x,2f9.0,8x,a8)',iostat=ioerr) 
     .             lon(is),lat(is),yr1,yr2,site(is)    
c               if( ioerr.ne.0 .or. stop(is).lt.1985.) then
c                   write(*,*) 'Error decoding line, ioerr: ',ioerr 
c                   write(*,*) 'Old-format glist file? , stop'
c                   stop
c               endif     
                write(span(is),'(f6.1,a,f6.1)') yr1,' ',yr2
          else
            write(*,*) '# sites > maxsit, STOP'
            stop 
          endif
        endif
      enddo  
      numsit = is   
c      print *,'DEBUG read glist file numsit ',numsit       
      
      return
      end                              

c---------------------------------------------------------------------------

      Subroutine read_cart(ifile,site,lat,lon,maxsit,numsit)     
         
      implicit none
                 
      integer*4 ifile,maxsit,numsit,ioerr,is,lastcol,nblen
     
      real*8 lat(maxsit),lon(maxsit),geod_pos(3),x(3),rot(3,3),pi

      data pi/3.1415926535897932D0/

      logical eof         

      character*3 alat,alon
      character*8 site(maxsit)
      character*256 line
                           
*     Format is simple: ID in columns 1-4, x y z follows in free format
*     Skip lines with blank column 1
                                 
      eof = .false. 
      is = 0 
      do while (.not.eof ) 
        read(ifile,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then 
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error reading line of CART file ',ioerr
        elseif( line(1:1).eq.' ') then
          continue                              
        else       
          is = is + 1
*         read the site id
          read(line(1:4),'(a4,1x,a3)',iostat=ioerr) site(is)(1:4)
          call uppers(site(is)) 
*         read x y z
          read(line(6:42),*,iostat=ioerr) x 
          if(ioerr.ne.0) then
            write(*,'(a,i4)') 'Error reading XYZ coordinates',ioerr
            write(*,'(a)') line
            stop
          endif
          call XYZ_to_GEOD(rot,x,geod_pos) 
          lat(is) = 90.d0 - geod_pos(1)*180.d0/pi 
          lon(is) = geod_pos(2)*180.d0/pi   
          if( lon(is).lt.0.d0) lon(is) = lon(is) + 360.d0
        endif
      enddo
      numsit = is
      return
      end

c--------------------------------------------------------------------------------

      Subroutine read_sopac(ifile,site,lat,lon,maxsit,numsit)     
         
      implicit none
                 
      integer*4 ifile,maxsit,numsit,ioerr,is,ix,field,i
     
      real*8 lat(maxsit),lon(maxsit)

      logical eof
                                          
      character*3 alat,alon
      character*8 site(maxsit)
      character*256 line
                           
      alat = 'lat'
      alon = 'lon'

*     Skip the first line
      read(ifile,'(a)') line 
                                 
*     Routine assumes lat,lon,site are in fields 2, 11, and 12, colon-delimited
     
      eof = .false. 
      is = 0 
      do while (.not.eof ) 
        read(ifile,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then 
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error reading line of SOPAC file ',ioerr
        else  
*         find the first semi-colon and read the site id
          ix = index(line,';')   
          is = is + 1              
          read(line(ix+1:ix+4),'(a4)',iostat=ioerr) site(is)(1:4)
          call uppers(site(is))  
          site(is)(5:8) = '    '
          if( ioerr.ne.0 ) then
            write(*,'(a,i4)') 'Error reading site id ',ioerr
            stop
          endif  
          field = 1
*         find the 11th field (lat)
          do while (field.lt.11) 
            if( line(ix:ix).eq.';' ) then
              field = field + 1      
            endif
            ix = ix + 1 
          enddo 
*         Lat and lon can have as few as 7 ( .123456) or as many as 11 (-123.456789)
*         characters, presenting a problem for formatted reads (unformatted doesn't work
*         because of the semi-colon); hence the complicate read in sb read_latlon.         
          call read_latlon( line,site(is),'lat',ix,lat(is) )  
*         find the 12th field (lon)
          do while(field.lt.12)
            if( line(ix:ix).eq.';' ) field = field + 1
            ix = ix + 1   
          enddo          
          call read_latlon( line,site(is),'lon',ix,lon(is) )
          if( ix.gt.600 ) stop
          if( lon(is).lt.0 ) lon(is) = lon(is) + 360.d0
        endif
      enddo
      numsit = is
                                              
      return
      end                                

c-------------------------------------------------------------------------

      Subroutine read_magnet(ifile,site,lat,lon,maxsit,numsit)     
         
      implicit none
                 
      integer*4 ifile,maxsit,numsit,ioerr,is,lastcol,nblen
     
      real*8 lat(maxsit),lon(maxsit)

      logical eof
                                          
      character*8 site(maxsit)
      character*256 line
                           
*     Format is simple: ID in columns 1-4, followed by lat, long, and height

                                 
      eof = .false. 
      is = 0 
      do while (.not.eof ) 
        read(ifile,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then 
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error reading line of MAGNET file ',ioerr
        else       
          is = is + 1
*         read the site id 
          read(line(1:4),'(a4)',iostat=ioerr) site(is)    
          call uppers(site(is))                       
          site(is)(5:8) = '    '
*         read lat and lon
          read(line(5:34),*,iostat=ioerr) lat(is),lon(is)
          if( lon(is).lt.0 ) lon(is) = lon(is) + 360.d0
        endif
      enddo
      numsit = is
      return
      end
          

c-------------------------------------------------------------------------

      Subroutine read_usgs(ifile,site,lat,lon,maxsit,numsit)     
         
      implicit none
                 
      integer*4 ifile,maxsit,numsit,ioerr,is,lastcol,nblen
     
      real*8 x(3),lat(maxsit),lon(maxsit),rot(3,3),geod_pos(3),pi

      logical eof
                                          
      character*8 site(maxsit)
      character*256 line    
        
      data pi/3.1415926535897932D0/
                           
*     Format begins in column 2:
* 0001 1985 01 01 00:00:00.00 9999999.00   -2582012.8647  -4256382.1787   3973935.9940    
                             
      eof = .false. 
      is = 0 
      do while (.not.eof ) 
        read(ifile,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then 
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error reading line of USGS file ',ioerr
        else       
          is = is + 1
*         read the site id 
          read(line(2:5),'(a4)',iostat=ioerr) site(is)(1:4)
          call uppers(site(is))                            
          site(is)(5:8) = '    '
*         read xyz
          read(line(40:85),*,iostat=ioerr) x
*         convert to lat, long
          call XYZ_to_GEOD(rot,x,geod_pos)  
          lat(is) = 90.d0 - geod_pos(1)*180.d0/pi 
          lon(is) = geod_pos(2)*180.d0/pi   
          if( lon(is).lt.0.d0 )  lon(is) = lon(is) + 360.d0
        endif
      enddo
      numsit = is
      return
      end

c-----------------------------------------------------------------------

      Subroutine read_norcal(ifile,site,lat,lon,maxsit,numsit)     
         
      implicit none
                 
      integer*4 ifile,maxsit,numsit,ioerr,is,lastcol,nblen
     
      real*8 lat(maxsit),lon(maxsit)

      logical eof
                                          
      character*3 alat,alon
      character*8 site(maxsit)
      character*256 line
                           
      alat = 'lat'
      alon = 'lon'
               
*     Each line has an variable number of entries, but the site ID is 
*     always in column 7 and the last two entries are always latitude
*     and longitude in format (f8.3,4x,f8.3) beginning in the last column
*     minus 20.
                                 
      eof = .false. 
      is = 0 
      do while (.not.eof ) 
        read(ifile,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then 
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error reading line of SOPAC file ',ioerr
        else       
          is = is + 1
*         read the site id from the second field
          read(line(7:10),'(a4)',iostat=ioerr) site(is)(1:4)
          call uppers(site(is))
          site(is)(5:8) = '    '
*         find the last non-blank column
          lastcol = nblen(line)
*         read the latitude and longitude 
          read(line(lastcol-19:lastcol),'(f8.3,4x,f8.3)',iostat=ioerr)
     .                 lat(is),lon(is)
          if( lon(is).lt.0 ) lon(is) = lon(is) + 360.d0
        endif
      enddo
      numsit = is
      return
      end

      subroutine uppers(string)
      character*(*) string
      integer*4 i,j,ichar

c     return an upper case string

      j = ichar('A')-ichar('a')
      do 10 i=1,len(string)
         if (lge(string(i:i),'a') .and. lle(string(i:i),'z')) then
            string(i:i)=char(ichar(string(i:i))+j)
         endif
  10  continue

      return
      end


c--------------------------------------------------------------------
 
      Subroutine read_sinex(usnx,site,lat,lon,maxsit,nss)     
 
*     Routine to get station coordinates from a SINEX file
*     rwk  May 2009 from subroutine upd_from_sinex in kf/htoglb/mstinf.f
                                                              
*       iusnx is the unit number, set in the main program
*       lat,lon in decimal degrees
*       nss is number of sites read                                                       

      implicit none
                              

      integer*4 usnx,maxsit,ioerr,nss
     .       , ilondeg,ilonmin,ilatdeg,ilatmin,i
     
      real*8 lat(maxsit),lon(maxsit),seclon,seclat

      logical eof,block_end,finished,found
                                          
      character*8 site(maxsit)
      character*256 line
                           
      finished = .false. 
      do while( .not.finished )
               
        read(usnx,'(a)', iostat=ioerr) line  
        if( ioerr.eq.-1 ) then
          finished = .true.
        elseif( ioerr.ne.0 ) then 
          write(*,*) 'Error reading SINEX line',ioerr
        elseif( line(1:8).eq. '+SITE/ID' ) then
cd          print *,'Found +SITE/ID'
          block_end = .false.  
          found = .false.    
          nss = 0
          do while (. not.block_end .and. .not.found ) 
            read(usnx,'(a)', iostat=ioerr) line
cd             print *,'SINEX line ',line
            if( ioerr.ne.0 ) then
              write(*,*) 'Error reading SITE ID line ',ioerr
            elseif( line(1:8).eq.'-SITE/ID' ) then 
cd              print *,'found -SITE/ID'
              block_end = .true. 
              finished = .true.
            elseif( line(1:1).eq.' ' ) then
              nss = nss + 1
              read(line,'(1x,a4,38x,2(i4,1x,i2,f5.1))')  site(nss)(1:4)
     .             ,ilondeg,ilonmin,seclon,ilatdeg,ilatmin,seclat    
              site(nss)(5:8) = '    '
              if( ioerr.ne.0 )  then
                 write(*,*)'Error decoding SINEX site coord line ',ioerr
              else  
                lon(nss)=ilondeg+ilonmin/60.d0+seclon/3600.d0
                lat(nss)=ilatdeg+ilatmin/60.d0+seclat/3600.d0
              endif
            endif
          enddo
        endif
      enddo                        
cd      print *,'nss ',nss
cd      if( nss.gt.0.and.nss.lt.maxsit ) then
cd         do i=1,nss
cd           print *,site(i),lat(i),lon(i)
cd         enddo
cd      endif
      return
      end      


c-------------------------------------------------------------------------------

      Subroutine write_merge(ifile,file,id,lat,lon,source,maxsit,numsit)

      implicit none

      integer*4 ifile,ifilea,ifilel,maxsit,numsit,ioerr,i  

      real*8 lat(maxsit),lon(maxsit) 
                
      character*3 source(maxsit)        
      character*8 id(maxsit)
      character*24 file,filealph,filelon

*     function
      integer*4 nblen
               
c      print *,'In write_merge maxsit numsit ',maxsit,numsit
*     Write the unsorted file
                             
cd      print *,'file ifile ',file,ifile                                          
      open(unit=ifile,file=file,status='unknown',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening output merged file '
     .         ,file,ioerr  
         stop  
      else
        write(*,'(2a)') 'Opened merge file: ',file
      endif    
cd      print *,'In write_merge numsit 1st,last ',numsit,id(1),id(numsit)
      do i = 1,numsit
        write(ifile,'(a8,1x,a3,2f13.7)') id(i),source(i),lon(i),lat(i)
      enddo

*     Sort alphabetically and write the file again
                                                  
      ifilea = ifile + 10   
cd      print *,'ifilea ',ifilea
      filealph = file(1:nblen(file))//'.sortalph'
      ifilea = ifile + 1
      open(unit=ifilea,file=filealph,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening output merged file '
     .           ,filealph,ioerr  
        stop  
      else
        write(*,'(2a)') 'Opened merge file: ',filealph
      endif
      call sorta(maxsit,numsit,id,source,lat,lon)
cd      print *,'From sorta maxsit numsit id source lat lon '
cd     .       , maxsit,numsit,id(1),source(1),lat(1),lon(1)
      do i = 1,numsit
        write(ifilea,'(a8,1x,a3,2f13.7)') id(i),source(i),lon(i),lat(i)
      enddo    

*     Sort by longitude and write the file again
    
      ifilel = ifile + 11
      filelon = file(1:nblen(file))//'.sortlon'
      open(unit=ifilel,file=filelon,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening output merged file '
     .           ,filelon,ioerr  
        stop  
      else
        write(*,'(2a)') 'Opened merge file: ',filelon
      endif  
      call sortn( maxsit,numsit,id,source,lat,lon ) 
      do i = 1,numsit
       write(ifilel,'(a8,1x,a3,2f13.7)') id(i),source(i),lon(i),lat(i)
      enddo                                  
    
      return
      end
           
   
c------------------------------------------------------------------------------------
                 
      integer function nblen (string)
c
c     given a character string, nblen returns the length of the string
c     to the last non-blank character, presuming the string is left-
c     justified, i.e. if string = '   xs  j   ', nblen = 8.
c
c     called non-library routines: none
c     language: standard fortran 77
c
      integer ls,i
      character*(*) string
      character*1 blank,null
      data blank /' '/
c
      null = char(0)
      nblen = 0
      ls = len(string)
      if (ls .eq. 0) return
      do 1 i = ls, 1, -1
         if (string(i:i) .ne. blank .and. string(i:i) .ne. null) go to 2
    1    continue
      return
    2 nblen = i
      return
      end

c--------------------------------------------------------------------------------------

      subroutine read_latlon( line,site,acrd,ix,coord )
      implicit none
                              
      character*(*) line,site,acrd                 
      integer*4 ix,ioerr
      real*8 coord

*     Lat and lon can have as few as 7 ( .123456) or as many as 11 (-123.456789)
*     characters, presenting a problem for formatted reads (unformatted doesn't work
*     because of the semi-colon).  So I've added some laborious logic to figure out
*     the number of characters well enough to read 5 decimal places (1 m).
*     this works because lat/lon are at least 9 characters line with 6 d.p 

*     Try first reading 6 characters; if the value is larger, read more
      read(line(ix:ix+5),'(f6.0)',iostat=ioerr) coord   
      if( ioerr.ne.0 ) then
         write(*,'(a,a4,1x,a3,a)') 'Error reading ',site,acrd,' f6.0'
         stop
      endif
      if( coord.lt.0.d0 ) then 
        if( coord.le.-100.d0 ) then 
          read(line(ix:ix+10),'(f11.0)',iostat=ioerr) coord  
          if( ioerr.ne.0 ) then
           write(*,'(a,a4,1x,a3,a)') 'Error reading ',site,acrd,' f11.0'
           stop 
          endif
        elseif( coord.le.-10.d0 ) then
          read(line(ix:ix+9),'(f10.0)',iostat=ioerr) coord  
          if( ioerr.ne.0 ) then
           write(*,'(a,a4,1x,a3,a)') 'Error reading ',site,acrd,' f10.0'
            stop 
          endif
        elseif( coord.le.-1.d0 ) then
          read(line(ix:ix+8),'(f9.0)',iostat=ioerr) coord  
          if( ioerr.ne.0 ) then
            write(*,'(a,a4,1x,a3,a)') 'Error reading ',site,acrd,' f9.0'
            stop 
          endif
        else
          continue
        endif
      else           
        if( coord.ge.100.d0 ) then
          read(line(ix:ix+9),'(f10.0)',iostat=ioerr) coord   
          if( ioerr.ne.0 ) then
           write(*,'(a,a4,1x,a3,a)') 'Error reading ',site,acrd,' f10.0'
           stop    
          endif
        elseif( coord.ge.10.d0 ) then
          read(line(ix:ix+8),'(f9.0)',iostat=ioerr) coord    
          if( ioerr.ne.0 ) then
            write(*,'(a,a4,1x,a3,a)') 'Error reading ',site,acrd,' f9.0'
            stop 
          endif
        elseif( coord.gt.1.d0 ) then
          read(line(ix:ix+7),'(f8.0)',iostat=ioerr) coord 
          if( ioerr.ne.0 ) then
            write(*,'(a,a4,1x,a3,a)') 'Error reading ',site,acrd,' f8.0'
            stop      
          endif
        else
          continue
        endif
      endif 

      return
      end


c----------------------------------------------------------------------------------

      Subroutine sorta( maxsit,n,id,src,lat,lon )
      
      implicit none
 
      integer*4 maxsit,n,i,j 

      real*8 lat(maxsit),lon(maxsit),latbuf1,latbuf2,lonbuf1,lonbuf2
                    
      
      character*8 id(maxsit),idbuf1,idbuf2  
      character*3 src(maxsit),srcbuf1,srcbuf2

c     Sort alphabetically 
cd      print *,'In sorta maxsit n id src lat lon '
cd     .   ,maxsit,n,id(1),src(1),lat(1),lon(1)
               
      do i = 1,n-1 
        do j = 1,n-i
           idbuf1 = id(j)
           idbuf2 = id(j+1)  
           srcbuf1 = src(j)
           srcbuf2 = src(j+1)
           latbuf1 = lat(j)
           latbuf2 = lat(j+1)
           lonbuf1 = lon(j)
           lonbuf2=  lon(j+1) 
           if( lle(idbuf1,idbuf2) ) then  
             id(j) = idbuf1
             id(j+1) = idbuf2  
             src(j) = srcbuf1
             src(j+1) = srcbuf2
             lat(j) = latbuf1
             lat(j+1) = latbuf2 
             lon(j) = lonbuf1
             lon(j+1) = lonbuf2 
           else    
             id(j) = idbuf2
             id(j+1) = idbuf1  
             src(j) = srcbuf2
             src(j+1) = srcbuf1
             lat(j) = latbuf2
             lat(j+1) = latbuf1 
             lon(j) = lonbuf2
             lon(j+1) = lonbuf1     
           endif 
         enddo
      enddo   
      return
      end


      Subroutine sortn( maxsit,n,id,src,lat,lon )
      
      implicit none
 
      integer*4 maxsit,n,i,j 

      real*8 lat(maxsit),lon(maxsit),latbuf1,latbuf2,lonbuf1,lonbuf2
                 
      character*3 src(maxsit),srcbuf1,srcbuf2
      character*8 id(maxsit),idbuf1,idbuf2
                 
c      print *,'In sortn maxsit n id src lat lon '
c     .    ,maxsit,n,id(1),src(1),lat(1),lon(1)


c     Sort by longitude
               
      do i = 1,n-1 
        do j = 1,n-i
           idbuf1 = id(j)
           idbuf2 = id(j+1)  
           srcbuf1 = src(j)
           srcbuf2 = src(j+1)
           latbuf1 = lat(j)
           latbuf2 = lat(j+1)
           lonbuf1 = lon(j)
           lonbuf2=  lon(j+1) 
           if( lonbuf1.le.lonbuf2 ) then  
             id(j) = idbuf1
             id(j+1) = idbuf2  
             src(j) = srcbuf1
             src(j+1) = srcbuf2
             lat(j) = latbuf1
             lat(j+1) = latbuf2 
             lon(j) = lonbuf1
             lon(j+1) = lonbuf2 
           else    
             id(j) = idbuf2
             id(j+1) = idbuf1  
             src(j) = srcbuf2
             src(j+1) = srcbuf1
             lat(j) = latbuf2
             lat(j+1) = latbuf1 
             lon(j) = lonbuf2
             lon(j+1) = lonbuf1     
           endif 
         enddo
      enddo   

      return
      end

