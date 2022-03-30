      Subroutine j_from_sp3 (afmt,gnss_sel,nprns,iprns,debug,iprndb )

c     Write a J-file by fitting an offset and rate (optionally an acceleration)
c     to values from the SV3 file; patterned after j_from_c.f, Feigl 1990.

c     King May/July 2015

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'

   
c     Unit number for SP3 file (opened in makej)
      integer*4 usv3
c     --unit number for J-file is usvclk in makex.h  
c     --unit number for nav file is unav in makex.h
                       
c        J-file format, gnss requested, and gnss codes for all SVs
       character*60 afmt
       character*1 gnss_sel

c        Quantities from the SP3 file
            
      character*1 spver,asvid
      character*2 linesym
      character*8 otlmod 
      integer*4 yr,mon,day,hr,min,doy,mjds,nepchs,numsp3sv,nprns
     .        , id,iprns(maxsat),clksigexpon,nepc,jd(maxepc,maxsat)
      real*8  sec,delts,frcts,clk,accsat(maxsat),pvsigb,clksigb
     .       ,clock(maxepc,maxsat),clkoff(maxseg),rate(maxseg)
     .       ,t(maxepc,maxsat),err(maxepc,maxsat)
      logical use(maxepc,maxsat)
                            
c        Quantities from the navigation file
      real*8 xeaf0,xeaf1,xeaf2
                               
c       Quantities from the clock fit
      integer*4 jdstart,npseg(maxseg),jd0(maxseg,maxsat),iwkn,nseg
      real*8    t0(maxseg,maxsat),sow,rms(maxseg)
      logical goodseg(maxseg)
                    
c        Local
      integer*4 iep,iprn,iprnx,isatcnt,sod,ioerr,i,j
      real*8 utcoff
      logical goodest(maxsat),batch,eoh,eof
      character*3 source
      character*80 line
      character*256 message 

c     function
      integer*4 julday,idoy

c     satellite to debug
      integer*4 iprndb
      logical  debug       

c Initialization
                   
      do i=1,maxsat
        iprns(i)=0
      enddo
        
c Read the SP3 file header to get times and satellites
          
      batch = .true. 
c     determine the file format using the first two characters of the first line
c       SP1    ' # '
c       SP3-a  '# '  or '# a'
c       SP3-b  '#b'
c       SP3-c  '#c'
c       SP3-d  '#d'
      read(usp3,'(a2)') linesym
      if ( linesym.eq." #") then  
        spver = '1'
      elseif ( linesym.eq.'# '.or.linesym.eq.'#a' ) then
        spver = 'a'
      elseif ( linesym.eq.'#b' ) then
        spver = 'b'
      elseif ( linesym.eq.'#c' ) then
        spver = 'c'
      elseif ( linesym.eq.'#d' ) then
        spver = 'd'
      else
        call report_stat('FATAL','MAKEJ','makex/j_from_sp3',' '
     .                  ,'Unrecognized SP orbit format',0 )
      endif
      rewind(usp3)
      call rsp3hd( usp3,gnss_sel
     .           , yr,mon,day,hr,min,sec
     .           , delts,mjds,frcts,nepchs 
     .           , numsp3sv,nprns,iprns,accsat
     .           , pvsigb, clksigb, otlmod )
      if( debug ) then  
        print *,'J_FROM_SP3 aft rsp3hd maxsp3sv numsp3sv nprns '
     .                        ,        maxsp3sv,numsp3sv,nprns   
        print *,'pvsigb,clksigb ',pvsigb,clksigb 
        print *,'nprns maxsat ',nprns,maxsat
        print *,'iprns ',(iprns(i),i=1,nprns)
      endif
      if( spver.ne.'c'.and.spver.ne.'d' ) 
     .   call report_stat('FATAL','MAKEJ','makex/j_from_sp3',' '
     .                   ,' SP3 file for clocks must be c or d',0)
                                                
c Convert start time to JD, sec-of-day (for debug only)
       if( debug ) print *,'read header ymd ',yr,mon,day 
       jdstart = julday(mon,day,yr)
c      SP3 is GPS time, convert to (GPS) JD, second-of-day
       doy = idoy(yr,mon,day)
       call timcon(-4,iwkn,sow,yr,doy,hr,min,sec,utcoff)  
       if( debug ) print *,'J_FROM_SP3 jdstart doy0 ',jdstart,doy 
      
            
c       Read the clock values into storage
         
      if(nepchs.gt.maxepc ) then 
        write(message,'(a,i4)') '# epochs on sp3 file > maxepc ',maxepc
        call report_stat('FATAL','MAKEJ','makex/j_from_sp3'
     .                        ,' ',message,0)
      endif
      do i = 1,nepchs
        read(usp3,'(3x,i4,4(1x,i2),1x,f10.7)',iostat=ioerr) 
     .      yr,mon,day,hr,min,sec
        if(ioerr.ne.0 ) call report_stat('FATAL','MAKEJ'
     .       ,'makex/j_from_sp3',' '
     .       ,'Error reading time tag from sp3 file ',ioerr)  
        if( debug) print *,'epoch ymd hm ',i,yr,mon,day,hr,min
        isatcnt = 0     
        
        do j = 1,numsp3sv 
          read(usp3,'(a)',iostat=ioerr) line
c         assume that sp3-c has the full 80 columns, though some can be blank
          read(line,'(1x,a1,i2,43x,f13.0,10x,i3)',iostat=ioerr) 
     .       asvid,id,clk,clksigexpon
          if(ioerr.ne.0 ) call report_stat('FATAL','MAKEJ'
     .            ,'makex/j_from_sp3'
     .            ,' ','Error reading clock value from sp3 file ',ioerr)
          if( asvid.eq.gnss_sel ) then
            isatcnt = isatcnt + 1        
            if(isatcnt.gt.maxsat ) then 
              write(message,'(a,i3)') 
     .          '# selected SVs on sp3 file > maxsat ',maxsat
              call report_stat('FATAL','MAKEJ','makex/j_from_sp3'
     .                              ,' ',message,0)
            endif
            jd(i,isatcnt) = julday(mon,day,yr)
            t(i,isatcnt) =  3600.d0*hr+60.d0*min+sec      
            if( id.ne.iprns(isatcnt)) then 
              write(message,'(a,i2,a,i2,a,i5)') 'Unexpected PRN (',id
     .          ,') at slot ',isatcnt,' for epoch ',i
              call report_stat('FATAL','MAKEJ','makex/j_from_sp3'
     .                        ,' ',message,0)
            endif
            clock(i,isatcnt) = clk      
            if(debug) print *,'  jd t i prn clk ',jd(i,isatcnt)
     .            ,t(i,isatcnt),isatcnt,iprns(isatcnt),clock(i,isatcnt)
c           if a value missing (999999.999999), flag it to not be used
            if(clock(i,isatcnt).gt.1.d5 ) then
               use(i,isatcnt) = .false. 
            else                         
              use(i,isatcnt) = .true.
              clock(i,isatcnt) = clock(i,isatcnt)*1.d-6
              err(i,isatcnt) = (clksigb*1.d-06)**clksigexpon
            endif  
            if(debug.and.iprns(isatcnt).eq.iprndb ) 
     .        print *,'prn clk clksig use ',iprns(isatcnt)
     .           ,clock(i,isatcnt),err(i,isatcnt),use(i,isatcnt) 
          endif
        enddo
      enddo
      
c Estimate an offset and rate for each satellite from the SP3 values
c allowing for jumps, and therefore multiple segments/j-file-entries
c If no good estimate available from the SP3 values, use the navigation file

      do  iprn=1,nprns
                                                 
        if(debug.and.iprns(iprn).eq.iprndb) print *
     .     ,'Fitting SP3 for PRN ',iprns(iprn)
         call clkfit( nepchs,jd(1,iprn),t(1,iprn)
     .              , clock(1,iprn),use(1,iprn)
     .              , nseg, npseg,goodseg, jd0(1,iprn), t0(1,iprn)
     .              , clkoff, rate, rms, iprns(iprn),debug, iprndb )
        if(debug.and.iprns(iprn).eq.iprndb ) 
     .     print *,'nseg goodseg jd0 t0 clkoff rate '
     .         ,nseg,goodseg(1),jd0(1,iprn),t0(1,iprn),clkoff(1),rate(1)
c       see if at least one good segment     
        goodest(iprn) = .false.
        do i=1,nseg
          if(goodseg(i)) goodest(iprn) = .true.
        enddo 
        if( goodest(iprn) ) then
c         write the values for all good segments
          source = 'SP3 '
          do i=1,nseg
            if(goodseg(i)) then    
              call writej( usvclk,jd0(i,iprn),t0(i,iprn),gnss_sel
     .                   , iprns(iprn),clkoff(i),rate(i),0.d0,rms(i)
     .                   , npseg(i),source,afmt )
            endif
          enddo  
        else
c         read the navigation file for this PRN    
          write(message,'(a,i3,a)') 'No good SP3 segment for PRN ',iprn
     .          ,', use nav-file'
          call report_stat('WARNING','MAKEJ','makex/j_from_sp3',' '
     .                    , message,0)
          rewind(unav) 
          call get_navclk( gnss_sel,iprn,iprns,nepc,jd(1,iprn),t(1,iprn)
     .                   , clock(1,iprn),use(1,iprn),debug ) 
          if( debug ) print *,'From GET_NAVCLK gnss iprn nepc 1st-last '
     .              , gnss_sel,iprn,nepc,clock(1,iprn),clock(nepc,iprn)
cd          if(debug.and.iprn.eq.iprndb) print *,'Fitting BRDC for PRN '
cd    .         ,iprns(iprn)
          call clkfit( nepc, jd(1,iprn), t(1,iprn)
     .               , clock(1,iprn), use(1,iprn)
     .               , nseg, npseg,goodseg, jd0(1,iprn), t0(1,iprn)
     .               , clkoff, rate, rms, iprns(iprn),debug, iprndb ) 
          if(debug) then 
            print *,'Aft CLKFIT nseg  ',nseg
            do i=1,nseg
              print *,'seg goodseg npseg jd0 t0 '
     .               ,nseg,goodseg(i),npseg(i),jd0(i,iprn),t0(i,iprn)
            enddo
          endif 
c         see if at least one good segment     
          goodest(iprn) = .false.
          do i=1,nseg
            if(goodseg(i)) goodest(iprn) = .true.
          enddo   
          if (goodest(iprn) ) then    
            write(message,'(a,i2,a)') 'Nav-file used for PRN '
     .              ,iprns(iprn),' clock '
            call report_stat('WARNING','MAKEJ','makex/j_from_sp3'
     .             ,' ',message,0)             
c           write the values for all good segments
            source = 'NAV'
            do i=1,nseg
              if(goodseg(i)) then    
                if(debug) print *
     .              ,'Calling WRITEJ nseg iprn prn jd0 t0 npseg'
     .             , i,iprn,iprns(iprn),jd0(i,iprn),t0(i,iprn),npseg(i)
                call writej( usvclk,jd0(i,iprn),t0(i,iprn)
     .                     , gnss_sel,iprns(iprn),clkoff(i),rate(i),0.d0
     .                     , rms(i),npseg(i),source,afmt)
              endif
            enddo
          else 
* MOD TAH 201003: Modified to make missing clock as warning and to write
*           zero for the clock entries.  (Normally if no clock something
*           is bad with satellite and it is deleted else where).
            write(message,'(a,i3)') 
     .       'No good clock fit from sp3 or nav file for PRN ',iprn
            call report_stat('WARNIG','MAKEJ','makex/j_from_sp3',' '
     .          ,message,0)
            call writej( usvclk,jdstart,0.0d0
     .                     , gnss_sel,iprns(iprn),0.d0,0.d0,0.d0
     .                     , 99.999d-9, 0 ,source,afmt)

          endif
        endif     
c       endif on goodest for SP3 file
      enddo 
c------ end of loop on SVs


c Report a a summary

      write(message,'(a,i3,a)') 'J-file written for ',nprns
     .     ,' SVs from SP3 file'
      call report_stat('STATUS','MAKEJ','makex/j_from_sp3',' '
     .   ,message,0)

      return
      end

c-------------------------------------------------------------------
      Subroutine writej(usvclk,jd,t,gnss,iprn,a0,a1,a2,rms,nepc,source
     .                 ,afmt)
                      
      implicit none 

c     Write a record to the J-file 

c     yr/doy/hms in UTC; wkno/sow in GPST

      integer*4 usvclk,jd,iprn,yr,doy,hr,min,itflag,iwkn,nepc

      real*8 t,a0,a1,a2,sec,utcoff,sow,rms                        

      character*1 gnss
      character*3 source                            
      character*60 afmt                   
      character*75 afmt2
                                   
c     get GPST in yr/doy/h/m/s for timcon
      call dayjul (jd,yr,doy)
      call ds2hms (yr,doy,t,hr,min,sec) 
c     convert to GPS week/sow
      call timcon(-4,iwkn,sow,yr,doy,hr,min,sec,utcoff)
c     convert to UTC yr/doy/h/m/s
      call timcon(1,iwkn,sow,yr,doy,hr,min,sec,utcoff)       
      afmt2 = afmt(1:59)//',f9.3,i5,2x,a3)'
cd      print *,'afmt ',afmt
cd      print *,'afmt2',afmt2 
c      write(usvclk,'(i4,1x,i4,2i3,1x,f10.7,2x,i4,1x,f14.7,2x
c     .              ,i2.2,2x,3d16.8,f12.3,i5,2x,a3)') yr,doy,hr,min,sec
c     .            , iwkn,sow,iprn,a0,a1,a2,rms*1.d9,nepc,source
       write(usvclk,afmt2) yr,doy,hr,min,sec,iwkn,sow,gnss,iprn
     .               , a0,a1,a2,rms*1.d9,nepc,source
      return
      end


