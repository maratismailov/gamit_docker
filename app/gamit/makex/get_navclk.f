      Subroutine GET_NAVCLK( gnss_sel,iprn,iprns,nepc,jd,t
     .                     , clock,use,debug )

c     Read a RINEX navigation file and store all of the clock values
c     R. King December 2015

      implicit none 

      include '../includes/dimpar.h'
      include '../includes/makex.h'
         
      integer*4 iprn,iprns(maxsat),nepc,jd(maxepc)
     .        , iprnx,iyear,month,iday,idn,ihr,min,isec
     ,        , itcorr,iwkcorr,iutc,iwkn,ioerr,i

      real*4 version  
      real*8 t(maxepc),clock(maxepc)
     .     , xeaf0,xeaf1,xeaf2,sec,a0,a1,utcoff    
     .     , TauN,GammaN,aGf0,aGf1,frame_time,sow 

      character*1 gnss_sel, gnss_type,gnss
c       gnss_sel : requested GNSS
c       gnss_type: type of nav file (may be 'M' for mixed)
c       gnss     : gode for PRN read, each record
      character*4 corrtyp
      character*5 corrsrc
      character*80 line    
      character*256 message 

      logical eoh,eof,use(maxepc),debug

c     Function
      integer*4 julday,idoy

      if(debug) print *,'GET_NAVCLK gnss_sel iprn ',gnss_sel,iprn

      eof = .false.   
      eoh = .false.      
      nepc = 0    

c  Read the first line to get the RINEX version number and GNSS 
  
      read (unav,'(a80)',iostat=ioerr) line
      rewind (unav)                                                   
      if( index(line,'RINEX') .eq. 0)  call report_stat('FATAL'
     .    ,'MAKEJ','makex/get_navclk',' ','Input nav file not RINEX ',0)
c       prior to v2.10, the version number could be an integer in column 6,
c       so to avoid an error reading the integer under f9.2, read this
c       value in free format.
       read(line(1:9),*,iostat=ioerr) version
       if(version.lt.2.0.or.version.gt.3.9) call report_stat('WARNING'
     .   ,'MAKEJ','makex/get_navclk',' '
     .   ,'Unsupported RINEX nav version',0)
      if(version.lt.3.0) then
        if( index(line,'NAVIGATION').gt.0 .or.
     .      index(line,'GPS').gt.0 ) then
          gnss_type = 'G'
        elseif( index(line,'GLONASS').gt.0 ) then
          gnss_type = 'R'
        endif 
        if(debug) print *,'V2 gnss_type ',gnss_type
      else
c       assume version 3
        read(line,'(40x,a1)',iostat=ioerr) gnss_type
        if(debug) print *,'V3 gnss_type ',gnss_type
c           with version 3, may have G R etc., or M for mixed 
      endif

c  Read the rest of the header to get the time correction (skip other records)

      eoh = .false.
      do while (.not.eoh )
        read(unav,'(a)',iostat=ioerr) line 
        if( ioerr.ne.0 ) then
          call report_stat('FATAL','MAKEJ','makex/get_navclk',' '
     .         ,'Error reading RINEX nav header',ioerr)
        else
          if(line(61:73).eq.'END OF HEADER') then
            eoh = .true.
          elseif(line(61:69).eq.'DELTA-UTC' ) then
c           this label for version 2 only
            read(line,'(3x,2d19.2,2i9)',iostat=ioerr) 
     .           a0,a1,itcorr,iwkcorr
          elseif(line(61:76).eq.'TIME SYSTEM CORR' ) then
            read(line,'(a4,1x,d17.10,d16.9,i7,i5,1x,a5,1x,i2)'
     .        ,iostat=ioerr) corrtyp,a0,a1,itcorr,iwkcorr,corrsrc,iutc
          endif   
        endif
      enddo

c  Read the clock values into storage

      do while (.not.eof )
        read(unav,'(a)',iostat=ioerr) line
        if( ioerr.eq.-1 ) then 
           eof = .true.
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','MAKEJ','makex/get_navclk'
     .        ,' ','Error reading nav-file PRN line',ioerr)
        else       
c         Record 1 PRN/EPOCH/CLK
         if( version.lt.3.0 ) then
c          need to use the header info to determine the GNSS: mixed not allowed
           gnss = gnss_type
           if(gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then 
             read(line,'(i2,5i3,f5.1,3d19.12)' ,iostat=ioerr)
     .               iprnx,iyear,month,iday,ihr,min,sec
     .             , xeaf0,xeaf1,xeaf2    
c            skip through the ephemeris records
             do i=1,7
               read(unav,'(1x)') 
               if( ioerr.ne.0 ) call report_stat('FATAL','MAKEJ'
     .            ,'makex/get_navclk',' '
     .            ,'Error skipping ephemeris records for non-Glonass'
     .            ,ioerr)
             enddo
           elseif(gnss.eq.'R' ) then 
             read(line,'(i2,5i3,f5.1,3d19.12)' ,iostat=ioerr)
     .            iprnx,iyear,month,iday,ihr,min,sec
     .          , TauN,GammaN,frame_time 
             xeaf0 = TauN  
             xeaf1 = GammaN
cd            print *,'V2 read prn date xeaf0 xeaf1 '
cd     .             , nprn,month,iday,ihr,min,sec,xeaf0,xeaf1
c            skip through the ephemeris records
             do i=1,3
               read(unav,'(1x)') 
               if( ioerr.ne.0 ) call report_stat('FATAL','MAKEJ'
     .            ,'makex/get_navclk',' '
     .            ,'Error skipping ephemeris records for Glonass',ioerr)
             enddo
           endif
           call fix_y2k(iyear)
         else
c          assume version 3 - can get the GNSS from the record itself
           read(line,'(a1)') gnss     
cd           print *,'V3 BUFF80 gnss',line,gnss 
           if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .         gnss.eq.'J'.or.gnss.eq.'I' ) then         
             read(line,'(a1,i2.2,1x,i4,5(1x,i2.2),3d19.12)'
     .         ,iostat=ioerr) gnss,iprnx,iyear,month,iday,ihr,min,isec
     .                     , xeaf0,xeaf1,xeaf2        
c            skip through the ephemeris records
             do i=1,7
               read(unav,'(1x)') 
               if( ioerr.ne.0 ) call report_stat('FATAL','MAKEJ'
     .            ,'makex/get_navclk',' '
     .            ,'Error skipping ephemeris records for non-Glonass'
     .            ,ioerr)
             enddo
           elseif( gnss.eq.'R' ) then             
             read(line,'(a1,i2.2,1x,i4,5(1x,i2.2),3d19.12)'
     .          ,iostat=ioerr) gnss,iprnx,iyear,month,iday,ihr,min,isec
     .                     , TauN,GammaN,frame_time   
             xeaf0 =  TauN  
             xeaf1 =  GammaN
cd              print *,'V3 read prn date xeaf0 xeaf1 '
cd     .             , gnss,nprn,month,iday,ihr,min,sec,xeaf0,xeaf1    
c            skip through the ephemeris records
             do i=1,3
               read(unav,'(a)') line
               if( ioerr.ne.0 ) call report_stat('FATAL','MAKEJ'
     .            ,'makex/get_navclk',' '
     .            ,'Error skipping ephemeris records for Glonass',ioerr)
             enddo
           elseif( gnss.eq.'S' ) then
             read(line,'(a1,i2.2,1x,i4,5(1x,i2.2),3d19.12)'
     .          ,iostat=ioerr) gnss,iprnx,iyear,month,iday,ihr,min,isec
     .                       , aGf0,aGf1,frame_time   
             xeaf0 = aGf0 
             xeaf1 = aGf1 
c            skip through the ephemeris records
             do i=1,3
               read(unav,'(1x)') 
               if( ioerr.ne.0 ) call report_stat('FATAL','MAKEJ'
     .            ,'makex/get_navclk',' '
     .            ,'Error skipping ephemeris records for SBAS',ioerr)
             enddo
                                       
           else
             write(message,'(a,a1,a)') 'GNSS id ',gnss,' not recognized'
             call report_stat('FATAL','MAKEJ','makex/get_navclk',' '
     .                      , message,0)
           endif  
           if( ioerr.ne.0 ) call report_stat('FATAL','MAKEJ'
     .      ,'makex/get_navclk',' ','Error decoding 1st nav line',ioerr)
           sec = dfloat(isec)
           endif                  
         endif
         if( gnss.eq.gnss_sel.and.iprnx.eq.iprns(iprn) ) then 
           nepc = nepc + 1               
           if( nepc.gt.maxepc ) then
             write(message,'(a,i4)') '# epochs > maxepc ',maxepc
             call report_stat('FATAL','MAKEJ','makex/get_navclk',' '
     .                      , message,0)
           endif   
c          Get the time tags in PEP JD,sec-of-day
c          (for Glonass doy/hr/min/sec are UTC, so need to convert the time
c          tags used for the polynomial fit.     
           if( gnss.eq.'R' ) then
             idn = idoy(iyear,month,iday)
             call timcon(-2,iwkn,sow,iyear,idn,ihr,min,sec,utcoff) 
           endif 
           jd(nepc) = julday(month,iday,iyear)
           t(nepc) = ihr*3600.d0 + min*60.d0 + sec     
c          use this utcoff to correc the Glonass clock offset from GPST: No don't do 
c          this since not correcting the PRs
c*           if( gnss.eq.'R') then
c*             xeaf0 = xeaf0 - utcoff 
c*           endif
           clock(nepc) = xeaf0 
           use(nepc) = .true.
         endif
      enddo     
      if(debug) then
        print *,'iprn nepc ',iprn,nepc
        do i=1,nepc
          write(*,'(i7,f10.2,1x,l1,d15.4)') jd(i),t(i),use(i),clock(i) 
        enddo
      endif
      return
                
      end

