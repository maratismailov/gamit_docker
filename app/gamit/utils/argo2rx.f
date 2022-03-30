      Program argo2rx

c     translate NGS ARGO (CIGNET) dat format to RINEX format

c     Kurt Feigl MIT/IPG 1991-1992
c     M.Burc Oral 5/12/1992
c     M.Burc Oral 5/20/1992   v. 8.2
c     R. King    12/12/1997   v. 9.1 to read station.info
                           
c     Normally run by com/sh_argo2rx (see)
c     Fortran input stream:
c       argo-filename
c       Y             (or N) to indicate RINEX file overwrite (or not)
c       apr-filename  (optional GLOBK apr file for header coordinates)
c     Will also attempt to read station.info for header info
* MOD TAH 200203: Added AntDAZ to list of values from station.info


      implicit none
       
c     --dimensions (maxsat)
      include '../includes/dimpar.h'  
c     --dimensions (maxchn,maxobt,,maxlin)
      include '../includes/makex.h'
             
      real*8 sec,sow,data(6,maxsat)
      logical first(maxsat)
c     first phase values of each satellite
      real*8 phase1(maxsat), phase2(maxsat)

c     line buffers and error flags
      character*256 message,line
      character*80 buff80
      character*4 upperc  
      integer ioerr
c     ioerr is generic, ioerr1 for ARGO epoch line, ioerr2 for ARGO data line

c     overwrite flags
      character*1 overwrite
      logical rinex_ok

c     ARGO file items
      integer iyr,iyr2,imo,idy,imn,prn(maxsat)
      character*3 rcvrsw
      character*20 stnnam20
      real*4 swver

c     RINEX defined items
      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
c     comment
      character*60 rxcom(maxlin)
c     mark name
      character*60 rxmrk
c     observer
      character*20 rxobs
c     agency
      character*40 rxagy
c     receiver serial number, type and SW version
c      character*20 rcvrsn - model.h
      character*20 rctype
c      character*20 rcvers - in model.h
c     antenna serial number and type
c      character*20 antsn - in model.h
      character*20 anttyp
c     aproximate coordinates
      real*8 apx,apy,apz
c     antenna offsets
      real*8 anth,ante,antn
c     wavelength factors
      integer nwave1,nwave2
c     observation types
      integer nobtyp
      character*2 rxobtp(maxobt)
c     data interval in seconds
      real*8 rxint 
c     number of satellites 
      integer*4 numsv
c     --satellite numbers for header
      integer*4 totsv(maxsat)
c     --number of obs of each type for each satellite for header
      integer*4 numobs(maxobt,maxsat)
          
c     station.info and anthi.dat items             
      character*4 blank4,prject,orbit
      character*5 hgtcod,radome
      character*6 rcvcod,antcod   
      character*16 stnnam16,amodel
      character*20 rcvers,rcvrsn,antsn
      integer*4 icall,jsessn,jstart(5),jstop(5)
      real*8 raw_up,raw_north,raw_east,offarp(3),dhpab
      real*8 antdaz  ! Antenna aligment from True N (deg).

c     data start time
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0
c     data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec1
      character*4 stncod,aprstn
                
      character*1 pcncod
      character*2 ayr0
      character*16 buff16

c     external function
      integer*4 iarray
      integer*4 idflag,idoy,jdoy,jsod,nsat,nerr,nprn,isessn,span
     .        , irxcom,issi0,issi1,issi2,isnr1,isnr2,nweek,nblen
     .        , ift,ihr,jmo,jdy,kmo,kdy,jyr,kyr,i,j    
c    debug
c     . ,k

      logical fcheck,apr_found,found_prn,old_stinf,warnings


c     file names and units
      character*80 argofile,rinexfile,aprfile
      integer uargo,ustnfo,uapr
c     urinex and uanthi declared in ../includes/makex.h
      parameter (uargo=1)
      parameter (ustnfo=3)
      parameter (uapr=8)  
      urinex = 2
      uanthi = 4
           
      do i=1,maxsat
         totsv(i) = 0
         phase1(i) = 0.d0
         phase2(i) = 0.d0
      enddo                             
      found_prn = .false.
 


***** Remove old versions of the status, warning, and error files

      call report_stat('CLEAR','ARGO2RX',' ',' ', ' ',0)
      call report_stat('CLEAR','LIB',' ',' ', ' ',0)
 
           
****** Read the input stream

c     ARGO file name
      argofile = ' '
      read(*,'(a)') argofile
      stncod = argofile(1:4)
      overwrite = 'n'
      read(*,'(a)') overwrite 
      call lowers(overwrite)
      aprfile = ' '
      read(*,'(a)') aprfile 
           
******* Open the ARGO file

      open(unit=uargo,file=argofile,status='old'
     .    ,form='formatted',iostat=ioerr) 
      if (ioerr.ne.0) then 
         call report_stat('FATAL','ARGO2RX','argo2rx',argofile
     .      ,'Cannot open input ARGO file',0) 
      else
        call report_stat('STATUS','ARGO2RX','argo2rx',argofile
     .                  ,'Opening: ',0)
      endif

       
******* Read the ARGO file header to get station name and date
       
      read(uargo,'(a80)',iostat=ioerr) buff80 
      if( ioerr.ne.0 ) call report_stat('FATAL','ARGO2RX','argo2rx',' '
     .  ,  'Error reading ARGO header line ',ioerr)
      read(buff80,'(a20,5i5)',iostat=ioerr)
     .     stnnam20,nweek,jdoy,jyr,jmo,jdy
      if( ioerr.ne.0 ) call report_stat('FATAL','ARGO2RX','argo2rx',' '
     .  ,  'Error decoding ARGO header line ',ioerr)
      jsod = 0  
      write(message,'(a,2i4,a,i4,a,i2,a,i2,a)') 'ARGO header epoch : '
     .     ,jyr,jdoy,' (',jyr,'-',jmo,'-',jdy,')' 
      call report_stat('STATUS','ARGO2RX','argo2rx',' ',message,0)

      
****** Check to make sure there is no frequency offset not accounted for by this program     

**    Week 573: 1st week of Jan 1991   arbitrary
      if ( nweek .lt. 460 ) then                                    
         call report_stat('FATAL','ARGO2RX','argo2rx',' '
     .          ,'Frequency plan used before week 460--run ARGO2FIC',0)
      else if (nweek .gt. 460 .and. nweek .lt. 573 ) then  
         call report_stat('WARNING','ARGO2RX','argo2rx',' '
     .     ,'Frequency plan MAY be used until week 573--run ARGO2FIC',0)
      endif
       
******  Get the RINEX header information 

c     read the data epoch line to get starting epoch for station.info
      read(uargo,'(i4,4i3,f11.7,f16.7,i4,11i3)',iostat=ioerr)
     .  irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,sow,nprn,
     .  (prn(i),i=1,nprn) 
      if( ioerr.ne.0 ) call report_stat('FATAL','ARGO2RX','argo2rx',' '
     .  ,   'Error reading ARGO epoch line for start epoch',ioerr)
      call fix_y2k(irxyr0)
      backspace (uargo)     
c     ARGO headers typically have generic (non-RINEX/IGS, non-GAMIT) receiver names:
c        TI4100, MINIMAC, ROGUE, ASHTECH, TRIMBLE
c        Use these only if station.info is not available and set the first three
c        characters of the GAMIT code for uniform identification of observable types 
      rctype = buff80(51:57)//'             '  
      rcvcod = rctype(1:3)//'   '
      if( rcvcod.eq.'TRI   ' ) rcvcod = 'TRM   '      
      rcvrsn = ' '
      rcvers = ' ' 
      anttyp = ' ' 
      antsn = ' '
c     GAMIT variable
      swver = 0.0  
c     sampling interval not known from header (and not necessary for RINEX)
      rxint = 0.
      anth = 0.
      ante = 0.
      antn = 0.  
      antdaz = 0.d0 
      rxpgm = 'ARGO2RX v 9.56'
      rxagy = 'CIGNET'
      rxobs = ' ' 
      rxcom(1) = 'NGS ARGO .dat file translated into RINEX'
      irxcom = 1                                            
      rxver = 2.10
      call getusr(buff16)
      rxusr = buff16//'    '
      call getdat(kyr,kmo,kdy)
      write (rxdat,'(i4.4,"/",i2.2,"/",i2.2)') kyr,kmo,kdy

c     Read the apr file if available

      if( aprfile(1:1).eq.' ' ) then
        call report_stat('WARNING','ARGO2RX','argo2rx',' '
     .,'No apr file input--approx coordinates omitted from RINEX header'
     .         , 0)       
        apx = 0.d0
        apy = 0.d0
        apz = 0.d0
      else
        open(unit=uapr,file=aprfile,status='OLD',iostat=ioerr) 
        if (ioerr .ne. 0)      
     .     call report_stat('FATAL','ARGO2RX','argo2rx',aprfile
     .                     ,'Error opening apr file',ioerr) 
        apr_found = .false.
        do while (ioerr.eq.0)
          read(uapr,'(a)', iostat=ioerr) line
          if( ioerr.eq.0 .and. nblen(line).gt.0 ) then
c             see if first character is non-blank
              if( line(1:1).eq.' ' ) then 
                 read(line,'(1x,a4)') aprstn
                 if( upperc(aprstn).eq.upperc(stncod) ) then
                   read(line,'(9x,a)') buff80
                   read(buff80,*) apx,apy,apz 
                   apr_found = .true.
                 endif
              endif
          endif
        enddo
        if( .not.apr_found ) 
     .   call report_stat('WARNING','ARGO2RX','argo2rx',stncod
     .,'Station id not found on apr file',0 )
      endif 
 

c     Read station.info if available

      if( fcheck('station.info') ) then
        open(unit=ustnfo,file='station.info',status='OLD',iostat=ioerr)
        if (ioerr .ne. 0) then     
           call report_stat('FATAL','ARGO2RX','argo2rx','station.info'
     .                     ,'Error opening  station.info file',ioerr)
        endif 
* MOD TAH 200203: Added AntDAZ to list of values from station.info
        call rstnfo( ustnfo
     .         , stncod, jyr, jdoy, jsod, span
     .         , stnnam16, raw_up, raw_north, raw_east, antdaz
     .         , rcvcod, antcod
     .         , hgtcod, radome, swver, rcvers, rcvrsn, antsn
     .         , jstart, jstop )
        rxmrk = ' '
        rxmrk(1:16) = stnnam16   
        rxcom(2) = 'Header information from station.info'
        irxcom = 2 
        if (fcheck('rcvant.dat')) then 
c        get the antenna full name 
         call read_rcvant(1,1,antcod,anttyp,radome,rcvcod,rctype,pcncod)
c        get the receiver full name
         call read_rcvant(1,2,antcod,anttyp,radome,rcvcod,rctype,pcncod)  
        else
           call report_stat('WARNING','ARGO2RX','argo2rx',' '
     .  ,'File rcvant.dat not linked, cannot get full rcvr/ant names',0)
        endif
c       if height not DHPAB (RINEX ARP), convert
        if( hgtcod.ne.'DHPAB' ) then 
          if( fcheck('hi.dat') ) then 
            open(unit=uanthi,file='hi.dat',status='OLD'
     .         ,iostat=ioerr)
            if (ioerr .ne. 0) then
              call report_stat('FATAL','ARGO2RX','argo2rx','antmod.dat'
     .                     ,'Error opening hi.dat',ioerr)
            endif       
            warnings = .true.
            call hisub(uanthi,raw_up,raw_north,raw_east,antcod
     .      , hgtcod,stncod,iyr,jdoy,isessn,offarp,warnings )
            anth = dhpab
            antn = offarp(2)
            ante = offarp(3)   
          else
            call report_stat('FATAL','ARGO2RX','argo2rx',' '
     . ,'File hi.dat must be linked to get RINEX antenna offsets',0)
          endif     
        else
          anth = raw_up
          antn = raw_north
          ante = raw_east
        endif
      else
         call report_stat('WARNING','ARGO2RX','argo2rx',' '
     .    ,'No station.info, using program defaults for RINEX header',0) 
        rxcom(2) = 'Most header values dummied'
        irxcom = 2                                            
      endif  


****** Set observable information 
                                               
      if ( rcvcod(1:3).eq.'MIN') then
         rcvrsw = 'MIN'
         nwave1 = 1
         nwave2 = 2
         nobtyp = 3
         rxobtp(1) = 'L1'
         rxobtp(2) = 'C1'
         rxobtp(3) = 'L2'
      elseif (rcvcod(1:3).eq.'ROG'.or.rcvcod(1:3).eq.'TRB') then
         rcvrsw = 'ROG'
         nwave1 = 1
         nwave2 = 1
         nobtyp = 4
         rxobtp(1) = 'L1'
         rxobtp(2) = 'P1'
         rxobtp(3) = 'L2'
         rxobtp(4) = 'P2'
      elseif (rcvcod(1:3).eq.'ASH') then  
         rcvrsw = 'ASH'
         nwave1 = 1
         rxobtp(1) = 'L1'
         rxobtp(3) = 'L2'
         if( swver.eq.0. ) then
           call report_stat('WARNING','ARGO2RX','argo2rx',' ',
     .  'Ashtech firmware version unknown, assuming P-code receiver',0)
            nobtyp = 4
            rxobtp(2) = 'P1'
            rxobtp(4) = 'P2'
            nwave2 = 1
         elseif( swver.lt.2.0 ) then 
            nobtyp = 3
            rxobtp(2) = 'C1'
            nwave2 = 2
         else 
            nobtyp = 4
            rxobtp(2) = 'P1' 
            rxobtp(4) = 'P2'
            nwave2 = 1
         endif 
      elseif (rcvcod(1:3).eq.'TRM') then    
         rcvrsw = 'TRM'
         nwave1 = 1     
         rxobtp(1) = 'L1' 
         rxobtp(3) = 'L2'
         if( swver.ge.5.50 ) then 
           nobtyp = 4 
           rxobtp(2) = 'P1' 
           rxobtp(4) = 'P2'
           nwave2 = 1 
         else 
           nobtyp = 3
           rxobtp(2) = 'C1'
           nwave2 = 2 
           if (swver.eq.0 ) 
     .         call report_stat('WARNING','ARGO2RX','argo2rx',' '
     .       ,'No firmware version, assuming codeless Trimble (SST)',0)
         endif 
      elseif ( rcvcod(1:3).eq.'TI4' ) then  
c        rcvrsw used only for mapamp call, in which GESAR and CORE are both 'GES'
         rcvrsw = 'GES'
         nobtyp = 4
         nwave1 = 1
         nwave2 = 1
         nobtyp = 4
         rxobtp(1) = 'L1'
         rxobtp(2) = 'P1'
         rxobtp(3) = 'L2'
         rxobtp(4) = 'P2'
      else     
         write(message,'(a,1x,a3,1x,a3)') 'Receiver code (',rcvcod(1:3)
     .        ,') unknonwn: cannot get wavelength factors'               
         call report_stat('FATAL','ARGO2RX','argo2rx',' ',message,0)
      endif                   


******  Come here at each change of day to open a new RINEX output file

c     read the data epoch line to get starting epoch
  10  read(uargo,'(i4,4i3,f11.7,f16.7,i4,11i3)',iostat=ioerr)
     .  irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,sow,nprn,
     .  (prn(i),i=1,nprn) 
      call fix_y2k(irxyr0)
      if( ioerr.ne.0 ) call report_stat('FATAL','ARGO2RX','argo2rx',' '
     .  ,   'Error reading ARGO epoch line for start epoch',ioerr)
      backspace (uargo)     
      jdoy = idoy(irxyr0,irxmo0,irxdy0)
      write(message,'(a,2i4,1x,2i3,f4.0,a,i4,a,i2,a,i2,a)') 
     .     'First epoch of day: ',irxyr0,jdoy,irxhr0,irxmn0,rxsec0
     .     ,' (',irxyr0,'-',irxmo0,'-',irxdy0,')'
      call report_stat('STATUS','ARGO2RX','argo2rx',' ',message,0) 
c     save alphameric year to use with epoch-line validity test 
      write(ayr0,'(i2)') mod(irxyr0,100)
c     set flag to avoid overwriting an existing RINEX file
      rinex_ok = .false.
      isessn = 1   
      do while (.not.rinex_ok)
         write (rinexfile,'(a4,i3.3,i1,".",i2.2,"o")')  
     .                    stncod,jdoy,isessn,mod(irxyr0,100)
         call lowers(rinexfile)  
 
c        if the file exists and overwrite has not been specified, increment the session number   
         if( overwrite.eq.'y') then
            rinex_ok = .true.
         else  
            if( fcheck(rinexfile) ) then
               isessn=isessn+1
            else
               rinex_ok = .true.
            endif
         endif
      enddo      
      open(unit =urinex,file=rinexfile,status='unknown',form='formatted'
     .    ,iostat = ioerr)
      if (ioerr.ne.0) then 
            call report_stat('FATAL','ARGO2RX','argo2rx',rinexfile
     .         ,'Cannot open output RINEX file',0) 
      else   
            call report_stat('STATUS','ARGO2RX','argo2rx',rinexfile
     .                       ,'Opening: ',0)
      endif

      jdy = irxdy0

c*     guess at last epoch
c      irxyr1 = irxyr0
c      irxmo1 = irxmo0
c      irxdy1 = irxdy0
c      irxhr1 = 23
c      irxmn1 = 59
c*     No, set month = 0 to key no valid entry --rwk 970321
       irxmo1 = 0  
c      set span for station.info check to 1 day arbitrarily--may generate warning
       span = 0

c      set header number of satellites = 0 to omit entry
       numsv = 0

c     write the RINEX header    
      call wrxhed (urinex,
     .   rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .   rcvrsn,rctype,rcvers,antsn,anttyp,apx,apy,apz,
     .   anth,ante,antn,nwave1,nwave2,nobtyp,rxobtp,rxint,
     .   irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .   irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1,
     .   numsv,totsv,numobs)  
 

***** every time a new day's file opened, reset first = true
      do ift=1,maxsat
          first(ift) = .true.
      enddo                   


****** start reading the file *********
c     read until end of file
      nerr = 0

c     **** come back here for each new epoch line (except with new day go to 10) ****
  20  continue
c     read time tag line
      read(uargo,'(a80)',end=1000,iostat=ioerr) buff80    
      if (ioerr.ne.0 .or. buff80(3:4).ne.ayr0 )  then 
         write(message,'(a,1x,a80)') 'Bad epoch line: ',buff80 
         call report_stat('WARNING','ARGO2RX','argo2rx',' ',message,0)
         nerr = nerr + 1 
         if( nerr.gt.100 ) goto 500
         goto 20
      endif
      read(buff80,'(i4,4i3,f11.7,f16.7,i4,11i3)',iostat=ioerr)
     .      iyr,imo,idy,ihr,imn,sec,sow,nprn,(prn(i),i=1,nprn)  
      call fix_y2k(iyr)
      if (ioerr .ne. 0) then 
        write(message,'(a,1x,a80)') 'Error decoding epoch line: ',buff80
        call report_stat('WARNING','ARGO2RX','argo2rx',' ',message,0) 
        nerr = nerr + 1 
        if( nerr.gt.100 ) goto 500
        goto 20
      endif

c     if day change, open a new RINEX file
      if (idy .ne. jdy) then
         jdy = idy
         backspace (uargo)
         close (urinex)
         goto 10
      endif

c     set flag to OK and write the RINEX epoch line  
      idflag = 0
      iyr2 = mod(iyr,100)
      write(urinex,'(5i3,f11.7,i3,12(i3))',iostat=ioerr)
     .      iyr2,imo,idy,ihr,imn,sec,idflag,nprn,(prn(i),i=1,nprn)    
         
c     add the PRN to the SV list if not there already
      do i = 1,nprn
         if( prn(i).ne.0 ) then
           found_prn = .false. 
           do  j = 1, numsv 
             if( prn(i).eq.totsv(j) ) found_prn = .true.
           enddo 
c          if this is new SV; if so add it to the list
           if( .not.found_prn ) then
             numsv = numsv + 1
             totsv(numsv) = prn(i) 
c            print *,'prn numsv totsv ',prn(i),numsv,(totsv(k),k=1,numsv)
           endif  
c           numobs(1,iarray(prn(i),totsv,numsv)) = 
c     .            numobs(1,iarray(prn(i),totsv,numsv)) + 1
c           idum = iarray(prn(i),totsv,numsv)
c          print *,'i prn numsv idum numobs',i,prn(i),numsv,idum
c     .          , numobs(1,iarray(prn(i),totsv,numsv))
         endif
      enddo 

c     read and write the data lines
c     one line per sat
      nsat = nprn
      do j = 1,nprn
        read(uargo,'(a80)',end=1000,iostat=ioerr) buff80 
        if( ioerr.ne.0 ) then 
          write(message,'(a,1x,a80)') 'Bad data line: ',buff80 
          call report_stat('WARNING','ARGO2RX','argo2rx',' ',message,0)
          nerr = nerr + 1 
          if( nerr.gt.100 ) goto 500
          goto 20
        endif
        read(buff80,'(4f17.5,2f5.1)',iostat=ioerr) (data(i,j),i=1,6)   
        if (ioerr .ne. 0) then 
          write(message,'(a,1x,a80)') 'Error decoding epoch line: '
     .        ,buff80
          call report_stat('WARNING','ARGO2RX','argo2rx',' ',message,0) 
          nerr = nerr + 1    
          if( nerr.gt.100 ) goto 500
          do i=1,6
             data(i,j) = 0.0d0
          enddo
          goto 20
        endif
      enddo
      issi0 = 0
c     remove the large phase.. use the first L1/L2 phase
c     ARGO format is always  data1 : L1     data2 : L2
      do j=1,nprn
        if( first(iarray(prn(j),totsv,numsv)) ) then  
c           phase1(prn(j)) = dint( data(1,j) )
c           phase2(prn(j)) = dint( data(3,j) )
          phase1(iarray(prn(j),totsv,numsv)) = dint( data(1,j) )  
          phase2(iarray(prn(j),totsv,numsv)) = dint( data(3,j) )  
          first(iarray(prn(j),totsv,numsv)) = .false.
c     debug
c      write(6,*)prn(j),data(1,j),phase1(iarray(prn(j),totsv,numsv))
c     .                ,data(2,j),phase2(iarray(prn(j),totsv,numsv))
        endif
c        data(1,j) = data(1,j) - phase1(prn(j))
c        data(3,j) = data(3,j) - phase2(prn(j))  
        data(1,j) = data(1,j) - phase1(iarray(prn(j),totsv,numsv))
        data(3,j) = data(3,j) - phase2(iarray(prn(j),totsv,numsv))
      enddo 
c     ARGO does not have loss-of-lock indicator
      issi0 = 0
c     write the RINEX data lines
      do j =1,nprn
c        map SNR into signal strength
         call mapamp (1,rcvrsw,swver,1,data(5,j),isnr1,issi1)
         call mapamp (1,rcvrsw,swver,2,data(6,j),isnr2,issi2) 
cdebug   write(6,*) '::',data(5,j),isnr1,issi1,'==',data(6,j),isnr2,issi2
         if (nobtyp .eq. 4) then
           write (urinex,31)
     .        data(1,j),issi1,
     .        data(2,j),' ',
     .        data(3,j),issi2,
     .        data(4,j),' '
   31      format (f14.3,1x,i1,f14.3,1x,1a,f14.3,1x,i1,f14.3,1x,1a)
        else if (nobtyp .eq. 3) then
             write (urinex,32)
     .           data(1,j),issi1,
     .           data(2,j),' ',
     .           data(3,j),issi2
   32        format (f14.3,1x,i1,f14.3,1x,1a,f14.3,1x,i1)
        endif
      enddo
      goto 20   
c     go read another epoch line

      
 500  call report_stat('FATAL','ARGO2RX','argo2rx',' '
     .    , 'Too many errors in ARGO file',0 )

1000  call report_stat('STATUS','ARGO2RX','argo2rx',' '
     .    ,'EOF on ARGO file--normal termination',0)

      stop
      END

