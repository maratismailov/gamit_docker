      Subroutine READU 

                                                         
c  R. W. King 19 July 2006
c  Read the u-file to get values for empirical models by matching
c  station code and keywrds:
c
c    OCEANLOD  ocean tidal loading 
c    ATMOSLOD  atmopsheric (non-tidal) loading
c    ATMTDLOD  atmospheric tidal loading
c    ATMOSMET  metereological values for atmospheric signal delay 
c    ATMOSMAP  coefficients for site- and time-dependent mapping functions


c   Output written to common /ufcom/ in model.h :
c
c     lotl          logical  true if values of ocean tidal loading read
c     otidemod      char*8   ocean tidal loading model used, e.g. FES2004
c     notl          int*4    number of components of ocean tidal loading
c     otides(i,6)   real*4   ocean tide loading values for i=notl components, 
c                            U W S amplitude,  U W S phase
c
c     latml         logical  true if values of atmospheric loading read
c     atmlmod       char*8   source/character of loading 
c     natml         integer  number of values read
c     atml_time(i)  real*4   times (doy) of values
c     atml_val(i,3) real*4   U N E loading values    
c
c     latl          logical  true if values of atmospheric tidal loading read
c     atidemod      char*8   source of tidal loading
c     natl          integer  number of components (hardwired to 2; S1 and S2)
c     atides(2,6)   real*4   loading values for S1 and S2 components
c                            N cos sin  E cos sin  U cos sin 
c               
c     lmet          logical  true if met values read (from RINEX if available, 
c                            otherwise from u-file if available)
c     metmod        char*8   string indicating source of pressure, temp, humidity
c     nmet          integer  number of values read
c     met_time(i)   real*4   i=1,nmet, times (doy) of values
c     met_val(i,3)  real*4   i=1,nmet pressure (hPa) temp (C) relative humidity (%) 
c
c     lmap          logical  true if mapping function coefficients read
c     mapmod        char*8   mapping function model
c     nmap          integer  number of values read
c     map_time(i)   real*4   i=1,nmap, times (doy) of values
c     map_val(i,4)  real*4   i=1,nmap, values of hydrostatic 'a', wet 'a' coefficients
c                            and hydrostatic zenith delay (m) and wet zenith delay (m)

       
      implicit none

      include '../includes/dimpar.h' 
      include '../includes/model.h' 

      logical eoh,eof,found_src,found_site,newblock,end_site
      integer*4 ioerr,nline,i,j

      real*4 uvers 

      real*8 tmptide(3)
                       
      character*4 upperc
      character*8 keywrd,keywrd1
      character*80 message
      character*709 line


c Initialize all u-file values
                            
      lotl = .false.
      latml = .false.
      latl = .false.
      lmet = .false.
      lmap = .false.
      otidemod = '        '
      atmlmod = '        '
      atidemod = '        '
      metmod = ' '
      mapmod = ' ' 
      notl = 0
      natml = 0
      ntatml = 0
      natl = 0
      nmet = 0  
      ntmet  = 0
      nmap = 0  
      ntmap = 0
      do j=1,6
        do i=1,maxotl
          otides(i,j) = 0.
        enddo
      enddo 
      do j=1,3
        do i=1,maxatml
          atml_val(i,j) = 0.
        enddo
      enddo
      do j=1,6
        do i=1,maxatl
          atides(i,j) = 0.
        enddo
      enddo
      do j=1,3
        do i=1,maxmet
          met_val(i,j) = 0.
        enddo
      enddo
      do j=1,8
        do i=1,maxmap
          map_val(i,j) = 0.
        enddo
      enddo


c Read the version number from the header record    
c    Old-style 
  
      eoh =.false.
      do while (.not.eoh)
        read(iuu,'(a)',iostat=ioerr) line 
c       if first non-comment line is not ' Version', assume Version 1.0 (OTL only)
        if( line(1:1).ne.' ') then
          continue
        else
          if( line(2:8).eq.'Version' ) then
            read(line(10:13),'(f4.1)',iostat=ioerr)  uvers 
            if(ioerr.ne.0)  call report_stat('FATAL','MODEL','readu',' '
     .                ,'Error reading u-file version number',ioerr)
          else
            uvers = 1.0
          endif
          eoh = .true.   
        endif
      enddo                  
      rewind(iuu)

c If Version 1.0 file, read the OTL in the old way  (code from setup)
      
c      print *,'READU uvers ',uvers
      if( uvers.eq.1.0 ) then  
        lotl = .true.
        otidemod= '        '
c       read the first six lines of the header to determine which model we have    
        nline = 0                              
        found_src = .false.
        do while (nline.lt.6 .and. .not.found_src )
          read(iuu,'(a)') line  
          nline = nline + 1
          if( index(line,'Scherneck').ne.0 ) then  
            notl = 11  
            found_src = .true.
          elseif ( index(line,'Matsumoto').ne.0 ) then
            notl = 54     
            found_src = .true.
          endif
        enddo
        if( .not.found_src ) then
            call report_stat('WARNING','MODEL','readu'
     .            ,' ','Ocean tide source not known, assuming OSO',0)
            notl = 11
        endif                
        rewind(iuu)
        found_site = .false.
        do while (.not.found_site) 
          read(iuu,'(a)',iostat=ioerr) line 
          if(ioerr.eq.-1 ) then
            call report_stat('FATAL','MODEL','readu'
     .         ,' ','EOF before finding site on old-style u-file',ioerr)
          elseif( ioerr.ne.0 ) then 
            call report_stat('FATAL','MODEL','readu'
     .     ,' ','Error reading ocean tides from old-style u-file',ioerr)
          endif
          if( upperc(line(2:5)).eq.upperc(sitecd)) then
            found_site = .true.  
c           old-style u-file will be missing model name but have 'lon' in cols 58-60
            if( line(58:60).eq.'lon' ) then   
c             model will probably be CSR3SCHW but mark 'OSO' to be safe
              otidemod = 'OSO     '
            else
              otidemod = line(58:65)
            endif
c           old amplitude units are meters, no conversion
            do j=1,6  
              if ( notl.eq.11 ) then                    
                read(iuu,'(11(1x,f6.0))',iostat=ioerr)
     .                 (otides(i,j),i=1,11)  
              elseif ( notl.eq.54 ) then
                read(iuu,'((10f7.0))',iostat=ioerr)
     .               (otides(i,j),i=1,54)     
              endif
c            print *,'OTIDES otidemod ',otidemod 
c            print *,(otides(i,j),i=1,54)
             if(ioerr.ne.0) call report_stat('FATAL','MODEL','setup'
     .           ,' ','Error reading ocean tide table',ioerr)
            enddo
          endif
        enddo 
        return
      endif
          

c  For Version > 1.0, read in the new way, checking keywrds for values present
                       
c    First look for the right site

      found_site = .false.
      eof = .false.
      do while (.not.found_site .and. .not.eof ) 
 
        read(iuu,'(a)',iostat=ioerr) line 
c        print *,'line 1 ',line   
        if( ioerr.ne.0 ) then    
          call report_stat('FATAL','MODEL','readu',' '
     .           ,'Error reading line of u-file',ioerr)  
        elseif( line(1:1).ne.' ') then
          continue
        elseif( line(2:8).eq.'ENDFILE' ) then
          eof = .true.    
        else 
          if( line(2:8).eq.'STATION'.and.
     .          upperc(line(10:13)).eq.upperc(sitecd) ) then    
            found_site = .true.     
          endif
        endif
      enddo  
      
c     Now read the block for each model type

c       The assumption in this loop is that within each model block,
c       all of the available records will be read, so that the 
c       next non-comment line should be a new model, a new station,
c       or an ENDFILE. 
                      
      end_site = .false.
      do while ( .not.end_site )
        read(iuu,'(a)',iostat=ioerr) line 
c        print *,'line 2 ',line
        if( ioerr.ne.0 ) then 
          call report_stat('FATAL','MODEL','readu',' '
     .             ,'Error reading MODEL line of u-file',ioerr) 
        elseif( line(1:1).ne.' ') then
          keywrd = ' '   
c          print *,'comment '
        elseif( line(2:8).eq.'STATION'.or.line(2:8).eq.'ENDFILE ') then
          end_site = .true.
        else
          keywrd = line(2:9)  
          if( line(11:15).ne.'MODEL' )  
     .      call report_stat('FATAL','MODEL','readu',' '
     .             ,'Keyword in first block line not MODEL',0)  
c          print *,'DEBUG end_site keywrd ',end_site,keywrd
        endif   

c       Read values for the block

        if( keywrd.eq.'OCEANLOD' .and. .not.end_site ) then  
c         read the model line, some comments, and then 6 data lines  
          read(line,'(16x,a8,i3)',iostat=ioerr)  otidemod,notl    
c          print *,'otidemod notl ',otidemod,notl
          if( ioerr.ne.0 ) 
     .       call report_stat('FATAL','MODEL','readu',' '
     .          ,'Error reading OCEANLOD model from u-file',ioerr)
c         read six non-comment lines of values    
          nline = 1
          do while ( nline.le.6 )
            read(iuu,'(a)',iostat=ioerr) line  
c            print *,'line 3 ',line
            if( ioerr.ne.0 ) then
               call report_stat('FATAL','MODEL','readu',' '
     .                ,'Error reading OCLEANLOD line of u-file',ioerr)
            elseif( line(1:1).ne.' ' ) then
c rwk 180531: temporarily keep the constituent line as a comment but read in the tokens
c             fix this when we start using FES2014
               if( line(15:16).eq.'M2') then
                 read(line,'(9x,54(4x,a3))',iostat=ioerr) 
     .              (otlwaves(i),i=1,notl)
                 if( ioerr.ne.0 )  
     .            call report_stat('FATAL','MODEL','readu',' '
     .                ,'Error reading OCEANLOD waves from u-file',ioerr)
               endif 
            else 
              read(line,'(1x,a8,100f7.0)',iostat=ioerr) 
     .             keywrd1,(otides(i,nline),i=1,notl) 
              if( ioerr.ne.0 ) then   
                   call report_stat('FATAL','MODEL','readu',' '
     .           ,'Error reading OCLEANLOD values ',ioerr)   
              elseif( keywrd1.ne.keywrd ) then   
                call report_stat('FATAL','MODEL','readu',' '
     .         ,'Unexpected keywrd reading OCEANLOD values from u-file'
     .         ,ioerr)
              else         
c                print *,'ok nline ',nline
                lotl = .true.
                if(nline.le.3) then
c                  convert amplitudes from mm to m
                   do i=1,notl
                      otides(i,nline) = otides(i,nline)/1000.  
                   enddo  
                endif
                nline = nline + 1   
              endif
            endif
          enddo   
                                                           
        elseif( keywrd.eq.'ATMOSLOD' .and. .not.end_site) then
c         read the model line and then an arbitrary number of epochs
          read(line,'(16x,a8,2i3)',iostat=ioerr) atmlmod,natml,ntatml  
c          print *,'atmlmod natml,ntatml ',atmlmod,natml,ntatml
          if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','readu',' '
     .        ,'Error reading MODEL line for ATMOSLOD in u-file',ioerr)
          ntatml = 0     
          newblock = .false.
          do while( .not.newblock )   
            read (iuu,'(a)',iostat=ioerr) line   
c            print *,'line 3 ',line
            if( ioerr.ne.0 ) then  
                call report_stat('FATAL','MODEL','readu',' '
     .           ,'Error reading ATMOSLOD MODEL from u-file',ioerr)    
            elseif( line(1:1).ne.' ') then
              continue
            elseif( line(2:9).ne.keywrd ) then
c              print *,'keywrd wrong ',keywrd
              newblock = .true.   
              backspace(iuu)    
c              print *,'backspaced ntatml ',ntatml
            else
              ntatml = ntatml + 1
              if( ntatml.gt.maxatml ) then
                write(message,'(a,i5,a)') 
     .          '# atm loading values on u-file (',ntatml ,') > maxatml'
                call report_stat('FATAL','MODEL','readu',' ',message,0)
              endif          
              latml = .true.
              read(line,'(9x,4f7.0)',iostat=ioerr ) 
     .             atml_time(ntatml),(atml_val(ntatml,j),j=1,3) 
c              print *,'latml ntatml time val '
c     .       ,latml,ntatml,atml_time(ntatml),(atml_val(ntatml,j),j=1,3)
              if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','readu'
     .           ,' ','Error reading ATMOSLOD values from u-file',ioerr)
            endif
          enddo
        
        elseif( keywrd.eq.'ATMTDLOD' .and. .not.end_site) then
c         read the model line and then always 1 data line   
          read(line,'(16x,a8,1x,i2)',iostat=ioerr) atidemod,natl  
c          print *,'atidemod natl ',atidemod,natml
          if( ioerr.ne.0)  call report_stat('FATAL','MODEL','readu',' '
     .       ,'Error reading ATMTDLOD MODEL from u-file',ioerr)
c         assume only S1 and S2 for now
          nline = 1  
          do while (nline.le.2)
            read(iuu,'(a)',iostat=ioerr) line      
c            print *,'line 3 ',line
            if( ioerr.ne.0 ) then  
              call report_stat('FATAL','MODEL','readu',' '
     .          ,'Error reading ATMTDLOD values line from u-file',ioerr)
            elseif( line(1:1).ne.' ' ) then
              continue   
            elseif (line(2:9).ne.keywrd ) then  
              call report_stat('FATAL','MODEL','readu',' '
     .          ,'Unexpected end of ATMTDLOD block in u-file',0)
            else 
              latl = .true.
              read(line,'(9x,12f7.0)',iostat=ioerr) 
     .             (atides(nline,j),j=1,6) 
              if (ioerr.ne.0) call report_stat('FATAL','MODEL','readu'
     .         ,' ','Error reading ATMTDLOD values from u-file,ioerr',0) 
c              print *,'i atides ',i,(atides(nline,j),j=1,6)   
c PT100827: model/etide is expecting NEU but the U-file has UNE. Need to
c             swap these around. This bug seems to have been introduced when
c             the u-file was revamped in 2006.
              tmptide(1) = atides(nline,1)
              tmptide(2) = atides(nline,2)
              do j=1,4
                atides(nline,j) = atides(nline,j+2)
              enddo
              atides(nline,5) = tmptide(1)
              atides(nline,6) = tmptide(2)
              nline = nline + 1
            endif
          enddo   

        elseif( keywrd.eq.'ATMOSMET' .and. .not.end_site) then
c       read the model line and then an arbitrary number of epochs
          read(line,'(16x,a8,2i3)',iostat=ioerr) metmod,nmet,ntmet
          if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','readu',' '
     .        ,'Error reading MODEL line for ATMOSMET in u-file',ioerr) 
          ntmet = 0   
          newblock = .false.
          do while( .not.newblock )   
            read (iuu,'(a)',iostat=ioerr) line    
c            print *,'line 3 ',line
            if( ioerr.ne.0 ) then  
                call report_stat('FATAL','MODEL','readu',' '
     .           ,'Error reading ATMOSMET MODEL from u-file',ioerr)    
            elseif( line(1:1).ne.' ') then
              continue
            elseif( line(2:9).ne.keywrd ) then
              newblock = .true.     
              backspace(iuu)  
              if( nmet.ne.3 ) call report_stat('FATAL','MODEL','readu'
     .         ,' ','# values on ATMOSMET u-file line must be 3',ioerr)
            else
              ntmet = ntmet + 1
              if( ntmet.gt.maxmet ) then
                write(message,'(a,i5,a)') 
     .           '# met values on u-file (',ntmet ,') > maxmet'
                call report_stat('FATAL','MODEL','readu',' ',message,0)
              endif        
              lmet = .true.   
c              print *,' ntmet nmet ',ntmet,nmet
              read(line,'(9x,4f8.0)',iostat=ioerr ) 
     .             met_time(ntmet),(met_val(ntmet,j),j=1,nmet)      
              if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','readu'
     .           ,' ','Error reading ATMOSMET values from u-file',ioerr)
            endif
          enddo
              
        elseif( keywrd.eq.'ATMOSMAP' .and. .not.end_site) then
c       read the model line and then an arbitrary number of epochs
          read(line,'(16x,a8,2i3,20(1x,a2))',iostat=ioerr) 
     .                mapmod,nmap,ntmap,(map_name(i),i=1,nmap)      
          if( metmod(1:4).eq.'    ' ) metmod = mapmod 
c          print *,'READU nmap map_name ',nmap,map_name
          if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','readu',' '
     .        ,'Error reading MODEL line for ATMOSMET in u-file',ioerr)
          ntmap = 0   
          newblock = .false.
          do while( .not.newblock )   
            read (iuu,'(a)',iostat=ioerr) line     
c            print *,'line 3 ',line
            if( ioerr.ne.0 ) then  
                call report_stat('FATAL','MODEL','readu',' '
     .           ,'Error reading ATMOSMAP MODEL from u-file',ioerr)    
            elseif( line(1:1).ne.' ') then
              continue
            elseif( line(2:9).ne.keywrd ) then
              newblock = .true. 
              backspace(iuu)                      
c              print *,'wrong keyword, backspacing '
            else
              ntmap = ntmap + 1
              if( ntmap.gt.maxmap ) then
                write(message,'(a,i5,a)') 
     .      '# mapping function  values on u-file (',ntmap ,') > maxmap'
                call report_stat('FATAL','MODEL','readu',' ',message,0)
              endif        
              lmap = .true.
c              print *,'lmap, reading values ',lmap
              read(line,'(9x,f7.0,9f11.0)',iostat=ioerr ) 
     .             map_time(ntmap),(map_val(ntmap,j),j=1,nmap)       
              if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','readu'
     .           ,' ','Error reading ATMOSMAP values from u-file',ioerr)
            endif
          enddo

        else
          continue 
c         endif for selecting block (MODEL line) 
        endif 
                    
c     enddo on site block
      enddo
            
      return
      end

