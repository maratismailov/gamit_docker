      subroutine reade( lu,icall,gnss
     .                , iflag,trans_sow,nprn,iewkn,bephem,bclock,subfr1)

C     Read the Broadcast Ephemeris from a file on unit lu (assumed open).
c     On first call (icall=0), the header is read and checked against
c     the requested GNSS (input with icall=0). On  subsequent calls, 
c     the next record (epoch/PRN) is returned, with gnss an output 
c     variable, and the calling program decides whether to use the record
c     (i.e., a useful epoch and the correct GNSS). 
c
c     Now accepts only RINEX formats (no longer FICA [GAMIT e-file']):

c     Input 
c       lu      : logical unit number (assumed open)
c       icall   : Called with file rewound (read headers)
c                    or not (read data records)
c       gnss    : Input variable if icall=0: requested GNSS
c
c     Output
c       iflag     : = -1 end of broadcast input file 
c                   =  0 record ok
c                   =  1 bad record  
c       trans_sow : Transmission time seconds of GPS week 
c       gnss     : System (G R E C J I S ) from the record 
c       nprn     : Satellite ID number (PRN)
c       iewkn    : Week number at reference epoch
c       bephem   : Broadcast ephemeris elements  
c       bclock   : Broadcast clock parameters (block(1) is sow at ref epoch)
c       subfr1   : Subframe 1 parameters
          
  
c     CLOCK PARAMETERS bclock(i)   
c      GPS
c      1   sow
c      2   xeaf0
c      3   xeaf1
c      4   xeaf2
c      GLONASS
c      1   sow
c      2   xeaf0 (-TauN)
c      3   xeaf1 (GammaN)
c      4   xfrmtim 


C     EPHEMERIS PARAMETERS bephem(i)
c      GPS
C      1   xetoe
C      2   xem0
C      3   xedn
C      4   xeart     square root of semi-major axis, in m**1/2
C      5   xeecc
C      6   xei0
C      7   xeidt
C      8   xeom0
C      9   xeomd
C     10   xew
C     11   xecuc
C     12   xecus
C     13   xecrc
C     14   xecrs
C     15   xecic
C     16   xecis
c     GLONASS
c      1   xpos (km, convert to meters when assigned)
c      2   ypos 
c      3   zpos 
c      4   xvel  (km/s, convert to m/sc)
c      5   yvel
c      6   zvel
c      7   health    (0 ok)
c      8   freqno
c      9   aoe (age of ephemeris in days) 
c
c
c     SUBFRAME 1 PARAMETERS subfr(i)
c      1     L2 code flag cflgl2
c      2     SV accuracy
c      3     SV health
c      4     Age of Data Clock
c      5     L2 P data flag
c      6     Group delay differential
c      7     Age of Data (ephemeris) 
c      8     SOW transmission time

c     TRANSMISSION TIME trans_sow in RINEX2 line 7 derived from 
c       Z-count in Hand Over Word (HOW)

c MOD TAH 200605: Added feature to check year is valid in GLONASS
c     entries.

      implicit none
       
      include '../includes/makex.h'

      character*1 gnss,gnss_type
      character*20 buff20
      character*80 buff80,prog_name
      character*256 message

C     Number of elements in array
      integer numf,numi,numc
C     Number of the block
      integer iblkid

      integer*4 ioerr
      real*8 bephem(16),bclock(6),subfr1(8)

      integer*4 it,lu,iday,idoy,jdoy,nprn,imonth,ihr
     .        , min,iflag,iyear,iewkn,idumb,icall,i
     .        ,indx,lerr,len,rcpar

c     RINEX NAVIGATION FILE VARIABLES
      real*4 version 
      real*8
     .         sow,sec,utcoff
     .       , xeaf0,xeaf1,xeaf2,aode,xecrs,xedn,xem0,xecuc,xeecc
     .       , xecus,xeart,xetoe,xecic,xeom0,xecis,xei0,xecrc,xew,xeomd
     .       , xeidt,cflgl2,weekno,pflgl2,svaccr,svhlth,tgd,aodc
     .       , trans_sow,spare1,spare2,spare3,aoe
c     Glonass
      real*8 xpos,xvel,xacc,health,ypos,yvel,yacc,freqno
     .     , zpos,zvel,zacc,age,TauN,GammaN,aGf0,aGf1,frame_time                                           
c     RINEX 3
      integer*4 isec 

c     Time correction correction 
      character*4 corrtyp
      character*5 corrsrc
      integer*4 itcorr,iwkcorr,iutc
      real*8    a0,a1,wkcorr
 
c     local variables
      integer inqerr
      character*16 fname
      logical eoh 
      logical debug/.false./

* MOD TAH 200605: Add save year
      integer*4 save_year   ! Year read in first record. Initialized
                            ! to zero so can be set.

      data save_year / 0 / 

      save version,gnss_type, save_year

c Get the name of the calling program and the nav file for error messsages

      len = rcpar(0,prog_name)         
      if( prog_name(1:2).eq.'BC' .or. prog_name(1:2).eq.'bc' )
     .    prog_name = 'ORBITS'  
      inquire ( unit=lu,name=fname,iostat=inqerr )

c Set the EOF flag  

      iflag = 0 
                 

c Read the header
    
      if (icall.eq. 0) then   
c       Read the first line to get the RINEX version number and GNSS 
        read (lu,'(a80)',iostat=ioerr) buff80
        rewind (lu)                                                   
        if( index(buff80,'RINEX') .eq. 0)  call report_stat('FATAL'
     .      ,prog_name,'lib/reade',fname,'Input file not RINEX ',0) 
c       prior to v2.10, the version number could be an integer in column 6,
c       so to avoid an error reading the integer under f9.2, read this
c       value in free format.
        read(buff80(1:9),*,iostat=ioerr) version
        if(version.lt.2.0.or.version.gt.3.9) call report_stat('WARNING'
     .       ,prog_name,'lib/reade',fname,'Unsupported RINEX version',0)
        if(version.lt.3.0) then
          if( index(buff80,'NAVIGATION').gt.0 .or.
     .        index(buff80,'GPS').gt.0 ) then
             gnss_type = 'G'
          elseif( index(buff80,'GLONASS').gt.0 ) then
            gnss_type = 'R'
          endif 
          if(debug) print *,'version < 3 gnss_type ',gnss_type
        else
c         assume version 3
          read(buff80,'(40x,a1)',iostat=ioerr) gnss_type
          if(debug) print *,'gnss_type ',gnss_type
c         with version 3, may have G R etc., or M for mixed 
        endif

c       Look for the time correction and the end-of header (skip other records)
        eoh = .false.
        do while (.not.eoh )
          read(lu,'(a)',iostat=ioerr) buff80 
          if( ioerr.ne.0 ) then
            call report_stat('FATAL',prog_name,'lib/reade',fname
     .       ,'Error reading RINEX nav header',ioerr)
          else
            if(buff80(61:73).eq.'END OF HEADER') then
              eoh = .true.
            elseif(buff80(61:69).eq.'DELTA-UTC' ) then
c             this label for version 2 only
              read(buff80,'(3x,2d19.2,2i9)',iostat=ioerr) 
     .           a0,a1,itcorr,wkcorr
            elseif(buff80(61:76).eq.'TIME SYSTEM CORR' ) then
              read(buff80,'(a4,1x,d17.10,d16.9,i7,i5,1x,a5,1x,i2)'
     .          ,iostat=ioerr) corrtyp,a0,a1,itcorr,iwkcorr,corrsrc,iutc
            endif
          endif   
        enddo
           
c     endif for initial header read
      endif     

c  Now read the ephemeris values
                                
      buff80 = ' '
      read (lu,'(a80)',iostat=ioerr) buff80
      if( ioerr.eq.-1.or.buff80(1:2).eq.' ' ) then
        iflag = -1                
        return        
      elseif( ioerr.ne.0 ) then
         call report_stat('WARNING',prog_name,'lib/reade',' '
     .     ,'Error reading record of navigation file',ioerr)
      else                         
c       Record 1 PRN/EPOCH/CLK
        if( version.lt.3.0 ) then
c         need to use the header info to determine the GNSS: mixed not allowed
          gnss = gnss_type       
          if(gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .       gnss.eq.'J'.or.gnss.eq.'I' ) then 
            read(buff80,'(i2,5i3,f5.1,3d19.12)' ,iostat=ioerr)
     .              nprn,iyear,imonth,iday,ihr,min,sec
     .            , xeaf0,xeaf1,xeaf2   
            elseif(gnss.eq.'R' ) then 
            read(buff80,'(i2,5i3,f5.1,3d19.12)' ,iostat=ioerr)
     .              nprn,iyear,imonth,iday,ihr,min,sec
     .            , TauN,GammaN,frame_time 
            xeaf0 = TauN  
            xeaf1 = GammaN
          endif
          call fix_y2k(iyear)
          if(debug)  print *,'V2 read prn date xeaf0 xeaf1 '
     .             , nprn,imonth,iday,ihr,min,sec,xeaf0,xeaf1
        else
c         assume version 3 - can get the GNSS from the record itself
          read(buff80,'(a1)') gnss     
          if(debug)  print *,'V3 BUFF80 gnss',buff80,gnss 
          if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .        gnss.eq.'J'.or.gnss.eq.'I' ) then         
             read(buff80,'(a1,i2.2,1x,i4,5(1x,i2.2),3d19.12)'
     .        ,iostat=ioerr) gnss,nprn,iyear,imonth,iday,ihr,min,isec
     .                     , xeaf0,xeaf1,xeaf2  
          if(debug) print *
     .       ,'READE gnss nprn iyear imonth iday ihr min isec '
     .       ,       gnss,nprn,iyear,imonth,iday,ihr,min,isec
          elseif( gnss.eq.'R' ) then             
            read(buff80,'(a1,i2.2,1x,i4,5(1x,i2.2),3d19.12)'
     .        ,iostat=ioerr) gnss,nprn,iyear,imonth,iday,ihr,min,isec
     .                     , TauN,GammaN,frame_time 

* MOD TAH 200605: Check that the year makes sense.  Some brdc files have
*     2042 years which fatals in leap.sec.  Assumption here is that first
*     record is not bad so that we know what year it is.
            if( save_year.eq.0 ) then
               save_year = iyear
            elseif ( abs(iyear-save_year).gt.1 ) then  ! Allow possible year
                                                       ! boundary
                write(message,210) iyear, save_year
 210            format('BAD BRDC year: ',I5,' Reset to ',i4)
                call report_stat('WARNING',prog_name,'lib/reade',fname,
     .                  message, 0)
                iyear = save_year
            endif
  
* MOD TAH 180316: There is no sign change of the clock.  Both Rinex 2
*     and Rinex 3 nav-files have the same sign for clock. The time
*     tag is UTC.  Clock rates seem to ~1e-11 and 14 second offset will
*     change clock error by < 1ns.
C           xeaf0 = -TauN  
            xeaf0 =  TauN  
            xeaf1 =  GammaN
          if(debug)    print *,'V3 read prn date xeaf0 xeaf1 '
     .             , gnss,nprn,imonth,iday,ihr,min,sec,xeaf0,xeaf1 
          elseif( gnss.eq.'S' ) then
            read(buff80,'(a1,i2.2,1x,i4,5(1x,i2.2),3d19.12)'
     .        ,iostat=ioerr) gnss,nprn,iyear,imonth,iday,ihr,min,isec
     .                     , aGf0,aGf1,frame_time   
            xeaf0 = aGf0 
            xeaf1 = aGf1
          else
            write(message,'(a,a1,a)') 'GNSS id ',gnss,' not recognized'
            call report_stat('FATAL',prog_name,'lib/reade',fname
     .                      , message,0)
          endif  
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/reade',fname,'Error decoding 1st nav line',ioerr)
          sec = dfloat(isec)
        endif                  
     
cRWK 160218/160219: Flag the one instance so far of a bogus entry
        if(iyear.lt.1980.or.iyear.gt.2050) then 
          iflag = 1
          write(message,'(a,a1,i2,1x,i4,4i3)') 'Bad year on nav-file: '
     .                 ,gnss,nprn,iyear,imonth,iday,ihr,min
            call report_stat('WARNING',prog_name,'lib/reade',fname
     .                      , message,0)
        endif

c       Record 2 ORBIT-1
        if( version.lt.3.0 ) then
          if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .        gnss.eq.'J'.or.gnss.eq.'I' ) then 
             read(lu,'(3x,4d19.12)',iostat=ioerr) aode,xecrs,xedn,xem0
          elseif(gnss.eq.'R'.or.gnss.eq.'S' ) then 
            read(lu,'(3x,4d19.12)',iostat=ioerr) xpos,xvel,xacc,health
          endif
        else
c         V 3 
          if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .        gnss.eq.'J'.or.gnss.eq.'I' ) then 
             read(lu,'(4x,4d19.12)',iostat=ioerr) aode,xecrs,xedn,xem0
          elseif(gnss.eq.'R'.or.gnss.eq.'S' ) then
             read(lu,'(4x,4d19.12)',iostat=ioerr) xpos,xvel,xacc,health
          endif
        endif
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/reade',fname,'Error decoding 2nd nav line',ioerr)


c       Record 3 ORBIT-2
        if( version.lt.3.0 ) then
          if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .        gnss.eq.'J'.or.gnss.eq.'I' ) then 
            read(lu,'(3x,4d19.12)',iostat=ioerr) xecuc,xeecc,xecus,xeart
          elseif(gnss.eq.'R'.or.gnss.eq.'S' ) then 
            read(lu,'(3x,4d19.12)',iostat=ioerr) ypos,yvel,yacc,freqno

          endif
        else
c         V 3 
          if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .        gnss.eq.'J'.or.gnss.eq.'I' ) then 
            read(lu,'(4x,4d19.12)',iostat=ioerr) xecuc,xeecc,xecus,xeart
          elseif(gnss.eq.'R'.or.gnss.eq.'S' ) then
             read(lu,'(4x,4d19.12)',iostat=ioerr) ypos,yvel,yacc,freqno
          endif
        endif
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/reade',fname,'Error decoding 3rd nav line',ioerr)

c       Record 4 ORBIT-3
        if( version.lt.3.0 ) then
          if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .        gnss.eq.'J'.or.gnss.eq.'I' ) then 
            read(lu,'(3x,4d19.12)',iostat=ioerr) xetoe,xecic,xeom0,xecis
          elseif(gnss.eq.'R'.or.gnss.eq.'S' ) then 
            read(lu,'(3x,4d19.12)',iostat=ioerr) zpos,zvel,zacc,aoe
          endif
        else
c         V 3 
          if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .        gnss.eq.'J'.or.gnss.eq.'I' ) then 
            read(lu,'(4x,4d19.12)',iostat=ioerr) xetoe,xecic,xeom0,xecis
          elseif(gnss.eq.'R'.or.gnss.eq.'S' ) then
             read(lu,'(4x,4d19.12)',iostat=ioerr) zpos,zvel,zacc,age
          endif
        endif
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/reade',fname,'Error decoding 4th nav line',ioerr)

        if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .      gnss.eq.'J'.or.gnss.eq.'I' ) then
c         Glonass and SBAS have only 4 records

c         Record 5 ORBIT-4
          if( version.lt.3.0 ) then
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then 
              read(lu,'(3x,4d19.12)',iostat=ioerr) xei0,xecrc,xew,xeomd
            endif
          else
c           V 3 
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then 
               read(lu,'(4x,4d19.12)',iostat=ioerr) xei0,xecrc,xew,xeomd
            endif
          endif
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/reade',fname,'Error decoding 5th nav line',ioerr)

c         Record 6 ORBIT-5
          if( version.lt.3.0 ) then
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then 
              read(lu,'(3x,4d19.12)',iostat=ioerr) xeidt,cflgl2, weekno
            endif
          else
c           V 3 
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then 
              read(lu,'(4x,4d19.12)',iostat=ioerr) xeidt,cflgl2, weekno
     .                                            , pflgl2
            endif
          endif
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/reade',fname,'Error decoding 6th nav line',ioerr)

c         Record 7 ORBIT-6
          if( version.lt.3.0 ) then
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then 
             read(lu,'(3x,4d19.12)',iostat=ioerr) svaccr,svhlth,tgd,aodc
            endif
          else
c           V 3 
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then 
             read(lu,'(4x,4d19.12)',iostat=ioerr) svaccr,svhlth,tgd,aodc
            endif
          endif
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/reade',fname,'Error decoding 7th nav line',ioerr)

c         Record 8 ORBIT-7
          if( version.lt.3.0 ) then
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     ,          gnss.eq.'J'.or.gnss.eq.'I' ) then 
              read(lu,'(3x,4d19.12)',iostat=ioerr) trans_sow
     .                                            ,spare1,spare2,spare3
            endif
          else
c           V 3 
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then 
              read(lu,'(4x,4d19.12)',iostat=ioerr) trans_sow
     .                                            ,spare1,spare2,spare3 
            endif
          endif
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/reade',fname,'Error decoding 8th nav line',ioerr)
c       endif on GPS,Beidou,Galileo extra records
        endif

c     endif on reading data records
      endif

c  Convert calender day, hr, min, sec to GPS week, sec-of-week
c     (For Glonass, this will include a time-system conversion from
c      UTC to GPST) 
c RWK 160218: We've found at least one instance of a bad year on a nav-file;
c     for this case, avoid a fatal in julday or taiutc by skipping the time 
c     conversion since the values will not be used by the calling programs 
c     (makexp, makex, bctot)
      if( iflag.eq.0 ) then 
        if( gnss.eq.'R') then
          it = -2
        else 
          it = -4  
        endif
        jdoy = idoy(iyear,imonth,iday)   
        call timcon (it,iewkn,sow,iyear,jdoy,ihr,min,sec,utcoff)     
        if(debug) 
     .     print *,'READE timcon it iyear jdoy ihr min sec iewkn sow '
     .            ,              it,iyear,jdoy,ihr,min,sec,iewkn,sow
      endif 

c  Assign the clock values
      bclock(1)   =  sow
      bclock(2)   =  xeaf0  
      bclock(3)   =  xeaf1
      if( gnss.eq.'R'.or.gnss.eq.'S' ) then
        bclock(4) = frame_time
      else
        bclock(4)   =  xeaf2
      endif                          
      bclock(5) = 0.d0
      bclock(6) = 0.d0                         
      if(debug) print *,'bclock ',bclock 

c  Assign the ephemeris values

      if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .        gnss.eq.'J'.or.gnss.eq.'I' ) then 
        subfr1(7 )  =  aode
C       units are meters
        bephem(14)  =  xecrs
        bephem( 3)  =  xedn
        bephem( 2)  =  xem0
        bephem(11)  =  xecuc
        bephem( 5)  =  xeecc
        bephem(12)  =  xecus
        bephem( 4)  =  xeart
        bephem( 1)  =  xetoe
        bephem(15)  =  xecic
        bephem( 8)  =  xeom0
        bephem(16)  =  xecis
        bephem( 6)  =  xei0
        bephem(13)  =  xecrc
        bephem(10)  =  xew
        bephem( 9)  =  xeomd
        bephem( 7)  =  xeidt
        subfr1(1)   =  cflgl2
        iewkn       =  weekno
c** RWK temporary code to fix incorrect week number for IRNSS
        if( gnss.eq.'I'.and.weekno.lt.1024) iewkn = weekno + 1024
c** RWK Beidou seems to start at week 1356 
        if( gnss.eq.'C'.and.weekno.lt.1356 ) iewkn = weekno + 1356
        subfr1(5)   =  pflgl2
        subfr1(2)   =  svaccr
        subfr1(3)   =  svhlth
        subfr1(6)   =  tgd
        subfr1(4)   =  aodc
      elseif( gnss.eq.'R'.or.gnss.eq.'S' ) then
        bephem(1) = xpos*1.d3
        bephem(2) = ypos*1.d3
        bephem(3) = zpos*1.d3
        bephem(4) = xvel*1.d3
        bephem(5) = yvel*1.d3
        bephem(6) = zvel*1.d3 
        bephem(7) = health
        bephem(8) = freqno
        bephem(9) = aoe 
        do i=10,16
         bephem(i) = 0.d0
        enddo
      endif
      if(debug) print *
     .  ,'READE trans_sow nprn iewkn bephem bclock subfr1 '
     .        , trans_sow,nprn,iewkn,bephem,bclock,subfr1
      return
      end 

