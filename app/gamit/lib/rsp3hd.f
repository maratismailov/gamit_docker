      Subroutine rsp3hd( iungs,gnss_sel
     .                 , iyear,imon,iday,ihr,imin,sec
     .                 , delt,mjd,fmjd,nepoch 
     .                 , numsp3sv,numsat,itsat, accsat 
     .                 , pvsigb,clksigb,otlmod )
c
c	P. Fang  March 1993, adopted from rsp1hd (R.W. King March 1988)    
c  R. King  December 2002: Mods to accommodate sp3-c format          
c  R. King  October 2014: Allow selecting any or all satellites
c  R. King  November 2015: Restrict selection to one GNSS 
c  
c      Read the header records for NGS Standard Product #3 orbit format
c      and save the PRN #s for the GNSS to be used 

c      Reference:  B.W. Remondi, "Extending the National Geodetic Survey
c                  Standard GPS Orbit Formats", NOAA Tech. Rep. NOS 133 NGS 46
c                  Rockville, MD, November 1989. 
c                  http://igscb.jpl.nasa.gov/newformats.html
c                  ftp://igscb.jpl.nasa.gov/pub/data/format/sp3c.txt.gz
                        
c    Input
c        iungs     Unit number of SP3 file   
c        gnss_sel  GNSS code to select from SP3 file (G R E C J I) 
c                                                     
c    Output 
c        spver    'a' 'b' 'c' or 'd' for SP3
c        iyear,imon,iday, ihr, imin, sec : Start date from header
c        mjd      Modified Julian Day of start
c        fmjd     Fractional day (from midnight)
c        nepoch   Number of epochs                      
c        numsp3sv Number of SVs on original sp3 file 
c        numsat   Number of PRNs matching the requested GNSS
c        itsat    PRN numbers SVs matching GNSS (dimension maxsat)  
c        accsat   Accuracy of SV from header (m)
c        pvsigb   Base for position/velocity
c        clksigb  Base for clock
c        otlmod   c*8 name of ocean tide model; matched with listing in otlcmc.dat
c                 to apply corrections to get from the CE frame of sp3 to CM frame of the t-file;
c                 if missing from header, assume no model (CM=CE frame)
c
c    Other quantities from the SP3 header
c        orbtyp   Orbit Type (F=fitted, E=extrapolated or predicted
c                            B=Broacast)
c        org      Agency source of ephemeris
c        crdsys   Coordinate System of ephemeris
c                    WGS72
c                    WGS84
c                    IER85 = Earth-fixed 1985 (IERS)
c                    ITR90
c                    ITR91  

      implicit none

      include '../includes/dimpar.h'

      logical eoh
                     
      character*1 pvflag,spver,gdum,otlflg,atlflg,gnss_sel
     .          , sp3gnss(maxsp3sv)
      character*2 linesym,ftype
      character*3 orbtyp,time_sys,clkmod,orbmod
      character*4 org
      character*5 crdsys
      character*8 otlmod,atlmod 
      character*10 pcvmod
c rwk 160408: line changed from c*69 to c*80 to accomodate sp3-d 
      character*80 line 
      character*80 prog_name
      character*256 message
      integer*4 iungs,mjd,jd,julday,iyear,iday,imon,ihr,imin,nepoch
     .        , numsp3sv,isp3sv(maxsp3sv),idumy(maxsp3sv)
     .        , numsat,itsat(maxsat),iaccsat(maxsp3sv),igpswk
     .        , ioerr,len,i,is
      real*8 fmjd,fmjdt,delt,ts,sec,accsat(maxsat),pvsigb,clksigb

c     function
      integer*4 rcpar 

      logical done,debug/.false./

c  initialization

      otlmod = ' '

c  Read the first and second records of the SP3 file header

      read(iungs,'(a)',iostat=ioerr) line
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .                 ,' ','Error reading 1st line of sp3 file',ioerr)  
      read(line,'(a2,a1,i4,4i3,1x,f11.8,3x,i5,7x,a5,1x,a3,1x,a4)'
     .    ,iostat=ioerr) linesym,pvflag,iyear,imon,iday,ihr,imin,sec
     .                 , nepoch,crdsys,orbtyp,org   
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .              ,' ','Error decoding 1st line of sp3 file',ioerr)  
      if( linesym.ne.'# '.and.linesym.ne.'#a'.and.linesym.ne.'#b'.and.
     .    linesym.ne.'#c'.and.linesym.ne.'#d'  ) then
        call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .        ,' ','Invalid version type on line 1 of sp3 file',0)
      else
        spver = linesym(2:2)
      endif
      call check_y2k(iyear) 
      if( debug ) print *,'RSP3HD linesym pvflag iyear imon iday nepoch'
     .                   ,        linesym,pvflag,iyear,imon,iday,nepoch
      read(iungs,'(3x,i4,17x,f14.8,1x,i5,1x,f15.13)',iostat=ioerr)
     .            igpswk,delt,mjd,fmjd                                 
      if( debug ) print *,'  igpswk delt mjd fmjd ',igpswk,delt,mjd,fmjd
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .                 ,' ','Error reading 2d line of sp3 file',ioerr)  

c  Set the start time

      ts = (dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec)/86400.d0
      jd = julday(imon,iday,iyear)
      fmjdt = dble(jd) + ts - 2400001.d0 
c     replace value off header with more precise computed value
      fmjd = ts

c  Read the SV numbers and their accuracy codes
            
c     lines 3-7 (sp3-c) or 3-11 (sp3-d): PRNSs 
c     rwk 160408 : it looks like the existing code works for any number of SVs
      read(iungs,'(a60)',iostat=ioerr) line
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .             ,' ','Error reading line 3 of sp3 file',ioerr)  
      read(line,'(4x,i2)') numsp3sv    
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .             ,' ','Error decoding line 3 of sp3 file',ioerr)  
      if( numsp3sv.gt.maxsp3sv ) then
        write(message,'(a,i2,a)') 
     .     'Number of satellites on SP3 file (',numsp3sv,') >',maxsp3sv
        call report_stat('FATAL',prog_name,'lib/rsp3hd',' '
     .                  ,message,0) 
      endif 
      backspace(iungs)          
      if( spver.eq.' ' .or. spver.eq.'a' .or. spver.eq.'b' .or.
     .    spver.eq.'c' ) then 
c       this code works when the #SVs allowed is fixed at 85
        if( numsp3sv.gt.85 ) call report_stat('FATAL',prog_name
     .      ,'lib/rsp3hd','SP3 version = a, b or c but # SVs > 85',0)
        read(iungs,'(4x,i2,3x,17(a1,i2),4(/,9x,17(a1,i2)))'
     .              ,iostat=ioerr)
     .     numsp3sv,(sp3gnss(i),isp3sv(i),i=1,numsp3sv),
     .               (gdum,isp3sv(i),i=numsp3sv+1,85)
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .            ,' ','Error reading lines 3-7 of sp3[abc] file',ioerr)
c       lines 8-12  accuracy codes  
        read(iungs,'(9x,17i3,4(/,9x,17i3))',iostat=ioerr)
     .    (iaccsat(i),i=1,numsp3sv),(idumy(i),i=numsp3sv+1,85)  
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .           ,' ','Error reading lines 8-12 of sp3[abc] file',ioerr)
      elseif( spver.eq.'d' ) then                           
c       the # SVs is now 999 (unlimited), so new logic is required
c       read the 1st SV line (sp3 line 3) for the number of SVs, then
c       use the column one symbol to read ther rest
        read(iungs,'(a)',iostat=ioerr) line
        if( ioerr.ne.0)  call report_stat('FATAL', prog_name
     .           ,'lib/rsp3hd',' ','Error reading line 3 ',ioerr)
        read(line(1:6),'(a2,1x,i3)') linesym,numsp3sv
        if( ioerr.ne.0.or.linesym.ne.'+ ') call report_stat('FATAL'
     .     ,prog_name,'lib/rsp3hd',' ','Error reading line 3 ',ioerr)
        read(line,'(9x,17(a1,i2))',iostat=ioerr) 
     .     (sp3gnss(i),isp3sv(i),i=1,17)
        if(ioerr.ne.0) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .     ,' ','Error decoding line 3 ',ioerr)
        is = 18              
        done = .false.
        do while (.not.done ) 
          read(iungs,'(a)',iostat=ioerr) line
          if(ioerr.ne.0)  call report_stat('FATAL',prog_name
     .       ,'lib/rsp3hd',' ','Error reading PRN line ',ioerr)   
          if(line(1:2).eq.'+ ') then 
            read(line,'(9x,17(a1,i2))',iostat=ioerr)
     .         (sp3gnss(i),isp3sv(i),i=is,is+16)
            if(ioerr.ne.0) call report_stat('FATAL',prog_name
     .          ,'lib/rsp3hd',' ','Error decoding PRN line',ioerr)
            is = is + 17
          else
            done = .true.
          endif
        enddo
c       next line should be the start of accuracy codes (++)
        read(line,'(a2,7x,17i3)',iostat=ioerr) linesym
     .        , (iaccsat(i),i=1,17)
        if(ioerr.ne.0) call report_stat('FATAL',prog_name
     .       ,'lib/rsp3hd',' ','Error decoding 1st accuracy line',ioerr)
        is = 18
        done = .false.
        do while (.not.done ) 
          read(iungs,'(a)',iostat=ioerr) line
          if(ioerr.ne.0)  call report_stat('FATAL',prog_name
     .        ,'lib/rsp3hd',' ','Error reading accuracy line ',ioerr)   
          if(line(1:2).eq.'++') then 
            read(line,'(9x,17(a1,i2))',iostat=ioerr) 
     .         (iaccsat(i),i=is,is+16)
            if(ioerr.ne.0) call report_stat('FATAL',prog_name
     .          ,'lib/rsp3hd',' ','Error decoding accuracy line',ioerr)
            is = is + 17
          else
           done = .true.
          endif
        enddo
      endif
         
c Select the GNSS SVs to be used

      numsat = 0 
      do i=1,numsp3sv    
        if( debug ) print *,'RSP3HD i sp3gnss ',i,sp3gnss(i)
        if(sp3gnss(i).eq.gnss_sel .or. 
     .     (gnss_sel.eq.'G' .and. sp3gnss(i).eq.' ') ) then
           numsat = numsat + 1   
          if( numsat.gt.maxsat ) then
            write(message,'(a,i2,a,i2,a)') 
     .        'Number of satellites selected from SP3 file ('
     .           ,numsat,') > maxsat (',maxsat,')'
            call report_stat('FATAL',prog_name,'lib/rsp3hd',' '
     .                      ,message,0)
          endif 
          itsat(numsat) = isp3sv(i)
c         convert the accuracy codes from 2**n mm to floating point m
          accsat(numsat) = 1.d-3*2.d0**iaccsat(i)
        endif
      enddo    

cd    print *,'numsp3sv isp3sv '
cd     .     ,numsp3sv,(sp3gnss(i),isp3sv(i),i=1,numsp3sv) 
cd    print *,'numssat itsat ',numsat',numsp3sv,numsat


c Read the file and time types and exponent base for time-dependent position and velocity sigmas    
                                                    
c     line 13 (c) or undetermined (d)
      read(iungs,'(a2,1x,a2,1x,a3)',iostat=ioerr) linesym,ftype,time_sys 
c     these will be dummy but readable for sp3-a
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .               ,' ','Error reading time-type line sp3 file',ioerr)
c    line 14 (c) or 16 (d)
      read(iungs,'(1x)',iostat=ioerr) 
c    line 15 (c)  or 17 (d)
      read(iungs,'(3x,f10.7,1x,f12.9)',iostat=ioerr) pvsigb,clksigb
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .   ,' ','Error reading position and velocity exponent base',ioerr)

c Rewind the file and look for models in the comment line

      rewind(iungs,iostat=ioerr)
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .     ,' ','Error rewinding sp3 file after reading header',ioerr)
      eoh = .false.
      pcvmod = ' '
      otlmod = ' '
      atlmod = ' '
      clkmod = ' '
      orbmod = ' '
      do while (.not.eoh ) 
        read(iungs,'(a)',iostat=ioerr) line   
        if( ioerr.ne.0 ) then
          call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .                ,' ','Error reading header for models line',ioerr)
        elseif( line(1:1).eq.'*') then 
c         end of comments
          eoh = .true.
        elseif( line(14:15).eq.'TL' ) then 
c         read comment line for models   
c         this group for temporary proposed format (Aug - 21 Nov 2006) 
          pcvmod = line(8:12) 
          otlmod = line(17:24) 
          otlflg = line(28:28)
          atlmod = line(30:37)
          atlflg = line(41:41) 
          clkmod = line(47:49)
          orbmod = line(55:57) 
c         this group for erroneous MIT files 29 Nov 2006 - 14 Apr 2011
        elseif( line(22:24).eq.'ALL') then
          pcvmod = line(8:17) 
          otlmod = line(26:33) 
c          otlflg = line(37:37) ?
          atlmod = line(39:46)  
          otlflg = line(55:55)
          atlflg = line(56:56) 
          orbmod = line(59:61)  
          clkmod = line(67:69)
c         this group for correct files according to IGS Mail 5490 22 Nov 2006
        elseif( line(22:24).eq.'AL:') then
          pcvmod = line(8:17) 
          otlmod = line(25:32) 
          atlmod = line(34:41)
          otlflg = line(43:43)
          atlflg = line(44:44)
          orbmod = line(50:52)
          clkmod = line(58:60)
        endif                                                  
      enddo

c  File should now be positioned at the first epoch; backspace one
  
      backspace(iungs,iostat=ioerr) 
      if( ioerr.ne.0 )  call report_stat('FATAL',prog_name,'lib/rsp3hd'
     .   ,' ','Error backspacing after reading comments',ioerr)

c  Echo the SP-3 file headers

      if( debug ) then 
        write(6,30) spver,orbtyp,crdsys,org,igpswk
     1          , iyear,imon,iday,ihr,imin,sec,mjd,fmjd,nepoch,delt
   30 format(//,1x,'Header records from NGS Standard Product file:'
     1      , // 
     a      ,1x,'SP version: ',a1,//
     2      ,1x,'Orbit Type: ',a3,'   Coordinate System : ',a5
     3      , '  Organization : ',a4,'   GPS Week : ',i4,/
     4      ,1x,'Start epoch (GPST): ',i4,4i3,1x,f10.7,/
     5      ,1x,'            MJD   : ',i5,1x,f15.14,/
     6      ,1x,'Number epochs     : ',i6,'  Interval :',f7.2,' sec')
        write(*,'(a,i3)')   ' Total number of SVs = ',numsp3sv                                                                      
        write(*,'(a,i3)')   ' Num SVs selected =',numsat
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
        write(*,'(a,50i3)') ' PRN #s           =',(itsat(i),i=1,numsat)
        write(*,'(a,50i3)') ' Accuracy (2**n mm)='
     .                          ,(iaccsat(i),i=1,numsat)
        if( pcvmod(1:1).eq.' ') then
          write(6,'(a)') 'No comment line for models'
        else
        write(6,32) pcvmod,otlmod,otlflg,atlmod,atlflg,clkmod,orbmod
   32     format(1x,'PCV model : ',a5,/
     .          ,1x,'OTL model : ',a8,'  applied = ',a1,/
     .          ,1x,'ATL model : ',a8,'  applied = ',a1,/
     .          ,1x,'CLK corrections : ',a3,/
     .          ,1x,'ORB corrections : ',a3)
        endif
      endif

      return
      end



