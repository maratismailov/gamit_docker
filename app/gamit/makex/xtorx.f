      program XTORX

c     Translate GAMIT X-files into RINEX format
c
c     Written by R. King from Kurt Feigl program CIG2RX  18 June 1990
c     and completed by Chris Vigny 28 June 1991
c     Mods by Burc Oral Feb and May 1992
c          Set bit 2 of the lli indicator to denote deleted, outlier, or
c            unweighted observations (GAMIT error flags ignone, igchop,
c            igoutl, and igunwt.
c     Mods by Kurt Feigl Jan 93:
c          Modify call to new xhdred subroutine
c          Call ../makex/settyp to determine important nwave factors.
c          Disable pseudorange warning if makex version cannot be found.
c          Dummy receiver type and antenna type if we cannot figure it out.
c          Handle Sercel receivers and antennas.
c          Catch typo bug in ../makex/settyp.f
c     Mods by Burch Oral and Tom Herring Jul 92 merged by Bob King with
c          Jan 93 Feigl changes:
c          Dimension jchan with maxsat, not maxchn.
c          Put a bias flag at the first epoch if file 'put_bias_rinex' exists
c          Display the allowed receiver types and/or software, and allow
c            manual input
c          Display a warning that bit 2 of the lli indicator is being set
c            (not RINEX standard);  add low amplitude (iglamp), low elevation
c            (igloel) and too-few points (ig2few) to conditions for setting bit 2.
c     Mods Tom Herring 930127:
c          (1) Fixed problem with the rinex header data types (L1,L2,P1,P2,C1)
c              not matching the data as written out.
c          (2) Removed setting Bit 2 in the rinex flag (not official use of this
c              bit, and replaced with not writting marginal data.
c     Mods Bob King 930215:
c          Calculate the RINEX antenna height either by reading station.info
c          or, if that file doesn't exist, from the antenna type and values on
c          the X-file.  If neither exists, zero is written.
c          Clean-up the statement numbers and comments.
c     Mods Peng Fang 930407:
c          Add Turbo-Rogue
c     Mods Peng Fang 930422
c          Add extra variable to ORDSIT calls to match new LIB routine
c     Mods Yehuda Bock 950211
c          Fix time conversion for GPST X-files
c          Hard wire X-file name day number for RINEX file name
c     Mods Yehuda Bock 950213
c          Fix time conversion for UTC X-files near day boundary
c          Add LEICA 200
c     Mods S. McClusky 951020
c          New call to hisub. Must open the antmod.dat file now to get correct
c          ARP heights in the rinex files.
c     Mods R. King 970111
c          Incorporate call to read_rcvant to get full names if not on X-file
c          Add comments and restructure some.
c     Mods R. King 970303
c          Set RINEX version number to 2.  
c     Mods R. King 970321
c          Moved from /utils to /makex and use makex.h file variables
c          Subtract the intial integer to conform to RINEX standards
c     Mods R. King 980720
c          Replace calls to gamit/lib sb with comlib/ function pickfn.
c          Fix bug in double conversion of UTC to GPST
c          Print warning and header comment about non-standard pseudoranges
c            with MAKEX ver < 7.19 ONLY for MiniMacs. 
c     Mods R. King 030320
c          Use doy from x-file name rather than start time to name RINEX file.
c     Mods R King  050618
c          Add radome to read_rcvant (but don't do anything with it)
c     Mods R King 050928
c          Add option to read table using hisub2 (eventually remove hisub)
c     Mods R King 100908
c          Add rcvers, rcvnum, antsn to rstnfo calling arguments
c     Mods R King 140828
c          Add satnam to xhdred calling arguments and RINEX output
c MOD TAH 200203: Added AntDAZ token for antenna  Alignment from True N

      implicit none
                     
c     --maximum dimensions for satellites in session
       include '../includes/dimpar.h'  

c     --file names and unit numbers
      include '../includes/makex.h'

c     --error flag codes
      include '../includes/errflg.h' 

c     --unit numbers not contained in makex.h
      integer*4 ustnfo,usestbl,udfile
       
      integer ioerr
     .      , iyr,iwkn0,jdoy0,im0,iwkn1,jdoyr1,jdoyr0,jdoyrx,jdoyx
     .      , jdobs0,id0,julday
     .      , jchan(maxchn),msat,prn(maxsat),issi0,issi1,issi2
     .      , idum,icount,iepoch,mepoch,iwkn,idoy,jdoy,istart,istop
     .      , kyr,kmo,kdy,itflag,imiss,len,span
     .      , istat,nstat,imapamp,ioc_ver,ioc_trans,nbatch,iclarg
     .      , iwarning,lli,isod,ii,i,j

      real*8 sec,sow0,sow,sow1,data(maxobt,maxsat),delta,utcoff
     .      , tobs0,dum8

c     not used:
      real*8 phs10(maxsat),phs20(maxsat)

* MOD TAH 930127: New variables to be used with the modified code.
*     iss_rx(maxobt) - Rinex definition signal strength by each
*                      data type.
*     obs_rx(maxobt) - Rinex observables.  The ordering here is
*                      same as the X-file order
*     msat_rx - Number of satellites observed at an epoch after
*               removing the data not used in the SOLVE solution
*     prn_rx(maxchn) - PRN numbers of the msat_rx satellites left
*               after checking for good data.
*     iremoved       - Number of data removed because of lgood
*                      not true.

      integer*4 iss_rx(maxobt), msat_rx, prn_rx(maxchn),
     .          iremoved
      real*8    obs_rx(maxobt)
*

      character*1  latflag,lonflag,pcncod
      character*3  rcvrsw,buff3,rxobtyp(maxdat)
      character*16 fname,xfname,rfname,wildcard,pickfn
      character*16 sitnam,amodel,satnam(maxsat)
      character*16 uname
      character*109 version
      character*80 buff80
      character*256 message

      logical first(maxsat), fcheck,sestbl_reqd,old_stinf

      character*1 upperc,fileid,xchar
      character*3 daynum
      character*4 snames(maxsit)

c     string ('yes' or 'no ') read from trivial file 'put_bias_rinex' to
c     indicate whether a bias flag is inserted at the first point
      character*3 put_bias_first

c     X-file defined items

      integer ischan(maxsat),ier(maxsat),lambda(maxsat,maxobt)
     .      , dattyp(maxobt),nchan,iy,im,id,ndat,mtime
     .      , ihr,min,inter,isnr(maxobt,maxsat),nepoch
     .      , latd,lond,latm,lonm,ntext,icol,ircint,isessn
     .      , numobsx(maxsat)
      character*1 gnss
      real*4 swver


c     items from station.info

      character*4 sitcod
      character*5 hgtcod,radome
      character*6 rcvcod,rxxcod,antcod,antcodx
      character*16 stanam
      integer*4 icall,istarts(5),istops(5)   
      real*4 swvers
      real*8 anthgt,offstn,offste,dhpab
      real*8 antdaz  ! Antenna aligment from True N (deg).


      real*4 MAKEX_ver
      real*8 offarp(3),seclat,seclon,height
      character*80 text(maxtxt)

      logical lbias,lgood,stnfo_ok,warnings
c     not used:  logica lmarg


c     RINEX defined items

      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
c     comment
      integer irxcom
      character*60 rxcom(maxtxt)
c     mark name
      character*60 rxmrk
c     observer
      character*20 rxobs
c     agency
      character*40 rxagy
c     receiver serial number,type and SW version
      character*20 rcvnum
      character*20 rctype   
      character*20 rcvers
c     antenna serial number and type
      character*20 antnum
      character*20 anttyp
c     aproximate coordinates
      real*8 apx,apy,apz
c     antenna offsets
      real*8 anth,ante,antn
c     wavelength factors
      integer nwave1,nwave2
c     observation types
      integer nobtyp
c     number of obs of each type for each satellite
      integer*4 numobs(maxobt,maxsat)
c     data interval in seconds
      real*8 rxint
c     data start time
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0
c     data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,iyr2
      real*8 rxsec1
      character*4 stcode
c     data flag
      integer idflag


c------------------------------------------------------------------------ 
           
c     Remove old versions of the status, warning, and error files

      call report_stat('CLEAR','XTORX',' ',' ', ' ',0)
      call report_stat('CLEAR','LIB',' ',' ', ' ',0)

c     Initialization

c     these stored in ../includes/makex.h
      uinfor = 8
      uscren = 6
      usceno = 10
      urinex = 11
      uxfile = 20
      uanthi = 32        
c     these local  (though makex.h/usited could be used for ustnfo)
      ustnfo  = 13
      usestbl = 15  
      udfile = 16        
      

c     Identify XTORX, user, and agency
          
c     open the print file  
      open (unit=uinfor,file='xtorx.out',form='formatted'
     .     , status='unknown',iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL','XTORX','xtorx','xtorx.out',
     .  'Error opening XTORX print file: ',ioerr)
      endif
      write(uinfor,'(/,a,/)') 'Program XTORX'
c     get the utilities version number and write it to the screen or log file
      call mversn(version)      
c     form is char*40: '9.37 of 97/01/13 07:50:00 (SunOS)'
      write(uinfor,'(a1)') ' '
c     begin the GAMIT.status file  
      write(message,'(3a)')  'Started XTORX (Utilities v. ',version,')'
      call report_stat('STATUS','XTORX','xtorx',' ',message,0) 
      call report_stat('STATUS','XTORX','xtorx',' '
     .  ,'Summary output written to file xtorx.out',0)

c     Fill the RINEX header lines that are common to all files

c    -END OF HEADER makes this version 2.  Optional records not written.
      rxver = 2.10
c     identify XTORX version in 20 characters
      write(rxpgm,'(a6,a5,a9)') 'XTORX ',version(1:5),version(9:17) 
c     XTORX user name, agency, and run date
      call getusr(uname) 
      rxusr = '    '//uname
      if( fcheck('sestbl.')) then
         open (unit=usestbl,file='sestbl.',status='old',iostat=ioerr)
         if( ioerr.ne.0) call report_stat('FATAL','XTORX','xtorx'
     .                 ,'sestbl.','Error opening file ',ioerr)
         sestbl_reqd = .true.
         call rdsest( 17, 'Processing agency',3,buff3,usestbl
     .              , sestbl_reqd, ioerr )
         rxagy = buff3 
         if( ioerr.ne.0 ) call report_stat('WARNING','XTORX','xtorx',' '
     .      ,'Agency missing from sestbl., set blank in RINEX header',0)
         rxagy = ' '
      else
         call report_stat('WARNING','XTORX','xtorx',' '
     .     ,'No sestbl. available, set agency blank in RINEX header',0)
          rxagy = ' '
      endif
      call getdat(kyr,kmo,kdy)   
      write (rxdat,'(i4.4,"/",i2.2,"/",i2.2)') kyr,kmo,kdy
c     put the GAMIT message into the first two comment lines of all RINEX headers
      irxcom = 2
      rxcom(1) = 'GAMIT (clean) X-files translated into RINEX'
      rxcom(2) = '  '  


c     Warn the user of limitations of the output RINEX file
  
      call report_stat('WARNING','XTORX','xtorx',' '
     . ,'Some items in RINEX may be missing or changed from original',0)
      call report_stat('WARNING','XTORX','xtorx',' '
     .  ,'--Observer not yet passed from X-file, set blank',0)  
      call report_stat('WARNING','XTORX','xtorx',' '
     .  ,'--Firmware version from GAMIT codes--see manual',0)  
      call report_stat('WARNING','XTORX','xtorx',' '
     .,'--Computation of RINEX-std ant hts not checked for all types',0)

         
c     Original (Burc Oral) code for loss-of-lock indicators and the adding of
c     bias flags at the first epoch (file 'put_bias_rinex') has been removed 
c     but is saved at MIT in /data13/rwk/old_gamit_active. --rwk 970113 
                 
c     first epoch flagged with loss-of-lock indicator (RINEX standard?)
      put_bias_first = 'yes'


c     Read the command line, if no argument is given then
c     Instruct to read a D-file with all sites included, or one X-file at a time
c

c     get input file name from the command line argument
      ii=iclarg(1,fname)

      if(ii.le.0) then
         write(uscren,101)
 101     format(/,1X,'Enter a X-file name, or',/
     1           ,1x,'      a D-file name to read a group of X-files :')
         wildcard = '[xXdD]*.*' 
         len = 9
         fname = pickfn( wildcard,len )
         fname = fname(1:len)
      endif
      fileid = upperc(fname(1:1))
      if (fileid.eq.'X') then
         nstat=1
         xfname = fname
C        Find the period in the filename
         icol = index(xfname,'.')
         snames(1) = xfname(icol-5:icol-2)
         daynum    = xfname(icol+1:icol+3)
         xchar     = xfname(icol-1:icol-1)
      else if (fileid.eq.'D') then
         open (unit=udfile,file=fname,status='old',iostat=ioerr)
         if( ioerr.ne.0 ) call report_stat('FATAL','XTORX','xtorx'
     .      ,fname,'Error opening D-file ',ioerr)
         call readdf( nbatch,udfile,snames,daynum,nstat )
         ii=iclarg(2,xchar)
         if(ii.le.0) then
            write(uscren,102)
 102        format(/,
     1      ' Enter series id (6th character) for input X-File names') 
c           the following is a kluge to avoid calling ftell (not available
c           in g77) to reposition the stdin buffer after calling pickfn.
c           --Peng Fang / Bob King  990525
            read(5,103) xchar,xchar
 103        format(a1)
         endif
      else        
         call report_stat('FATAL','XTORX','xtorx',fname
     .           ,'Input file neither X nor D',0)
      endif

c     Begin loop over files
      do 200 istat=1,nstat

c        construct the X-file name
         if( fileid.eq.'D' ) then
            xfname= 'X'//snames(istat)//xchar//'.'//daynum
         endif
         call lowers(xfname)

c        open the X-file (input)
         open(unit=uxfile,
     1        file=xfname,
     2        status='old',
     3        form='formatted',
     4        iostat=ioerr)
         if (ioerr.ne.0) then 
            write(message,'(a,a16,a)') 
     .          'Cannot open file ',xfname,' --continue'
            call report_stat('WARNING','XTORX','xtorx',' ',message,0) 
            write(uinfor,'(a,//)') message
            goto 200
         else     
           write(message,'(a,a16)') 'Opened file ',xfname 
           call report_stat('STATUS','XTORX','xtorx',' ',message,0)
         endif

c------------------------------------------------------------------------


c **** Determine the RINEX header information
                        
c        Re-initialize the comment line counter 
c        (The first two lines have the GAMIT comments, common to all files)
         irxcom=2

c        Scan the X-file to get the actual start and stop times and the
c        number of observations for each SV.
         write(uinfor,'(a)')
     .     ' Scanning the X-file for start and stop times'
         irxyr0 = 0       
cd         print *,'DEBUG XTORX calling scanx '
         call scanx ( istart,istop
     1              , irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0
     2              , irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1
     3              , nchan,numobsx )
c         write(*,*) irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0
c         write(*,*) irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1
         if( irxyr0.eq.0 ) then 
         write(message,'(a)') 'No valid data on X-file: go to next file'
           call report_stat('WARNING','XTORX','xtorx',' ',message,0)
           goto 200
         endif
         rewind uxfile                       

c        Scan the X-file to look for old fashion amplitude maping
         call passheader(uxfile)
         imapamp=0
         do 110 idum=1,1000
c           read the epoch, # sats, and observation time
            read ( uxfile,
     1             fmt = '(/,2I4,I5,I4,2I3,F11.7)',
     2             end=111,iostat=ioerr)
     3             iepoch,msat,iyr,jdoy,ihr,min,sec
            if( ioerr.ne.0 ) call report_stat('FATAL','XTORX','xtorx'
     .           ,' ','Error reading epoch line ',ioerr) 
c
            if( jdoy.gt.0 ) call fix_y2k(iyr)
c           read the data records
            do 110 j=1,msat
               read( uxfile,
     1               fmt = '(10X,2I2,2(1X,D22.15,1X,I3),2(2X,D22.15))',
     2               end=111)
     3               ier(j),jchan(j),(data(i,j),isnr(i,j),i=1,2)
     4               ,(data(i,j),i=3,ndat)   
               if( ioerr.ne.0 ) call report_stat('FATAL','XTORX','xtorx'
     .           ,' ','Error reading data line ',ioerr) 
                        
               if( isnr(1,j).gt.10.or.isnr(2,j).gt.10 ) then
                  imapamp=1
                  rxcom(irxcom+1)=
     1            'WARNING: Old fashion signal amplitudes on X-file'
                  rxcom(irxcom+2)=
     1            '         translation done by mapamp'
                  irxcom=irxcom+2
                  goto 111
               endif
 110     continue

 111     rewind uxfile  
c        patch to avoid iyr = 0 after a no-obs record - fang/rwk 940818
         if( iyr.eq.0 ) iyr = irxyr0

c        check MAKEX version number,in order to detect MiniMac pseudo-ranges
c        not meeting RINEX specification.
c        The MAKEX key is not always there, though.
         ioc_ver=0
         ioc_trans=0
         iwarning=0
  115    continue
           read(uxfile,'(a80)') buff80
           ioc_ver=index(buff80,'MAKEX v. ')
           ioc_trans=index(buff80,'X-File written from X-File')
           if (ioc_ver.ne.0) read(buff80,fmt='(9x,f5.2)') MAKEX_ver
           if (ioc_trans.ne.0) iwarning=1
         if (buff80.ne.'END') goto 115
         rewind uxfile  
c
c        Read the X-file header

         call xhdred ( uxfile,uinfor,uscren
     1               ,nepoch,inter,ircint,isessn
     2               , mtime,iy,im,id,ihr,min,sec
     3               , nchan,ischan,satnam
     4               , ndat,dattyp,rxobtyp,lambda
     5               , offarp,sitnam,rcvrsw,swver,antcodx
     6               , rctype,rcvnum,anttyp,antnum
     7               , latflag,latd,latm,seclat
     8               , lonflag,lond,lonm,seclon,height
     9               , ntext,text,gnss )       
                
         if( rcvrsw.eq.'MIN' .and. MAKEX_ver.lt.7.19.and.MAKEX_ver.gt.0
     .     .and. iwarning.eq.0 ) then
           rxcom(irxcom+1)=
     1'WARNING: Pseudoranges DO NOT meet RINEX definition !!!!!!!!!'
           rxcom(irxcom+2)=
     1'         Satellite clock offset has been removed by receiver'
           rxcom(irxcom+3)=
     1'         Phase data conform to RINEX specification.'
           irxcom=irxcom+3
         endif
                

c        Get the GPST week and second-of-week for start
 
         jdoy0= idoy(iy,im,id)
c        Case of X-file UTC time
         if (mtime.eq.1) then
          itflag = -2
c        Convert to GPS time
          call timcon (itflag,iwkn0,sow0,iy,jdoy0,ihr,min,sec,utcoff)
c        Convert to GPS calendar day
          itflag = 4
          call timcon (itflag,iwkn0,sow0,iy,jdoy0,ihr,min,sec,utcoff)
c        Case of X-file GPS time
         elseif (mtime.eq.2) then
          itflag = -4
          call timcon (itflag,iwkn0,sow0,iy,jdoy0,ihr,min,sec,utcoff)
         endif        
         isod = ihr*3600 + min*60 + dint(sec)  
c        jd, sod used when epoch missing   
         call monday(jdoy0,im0,id0,iy)
         jdobs0 = julday( im0,id0,iy )
         tobs0 = isod

c        Get the sampling interval and span

         rxint = dfloat(ircint)
         if( ircint.eq.0 .or. ircint.lt.inter ) rxint = dfloat(inter)
         span = nepoch * inter
                 
       
c        Read station.info (primary for antenna offsets; check for rcvr and ant types)
                            
         stnfo_ok = .false.
         if( fcheck('station.info') ) then 
            open( unit=ustnfo,file='station.info',status='old'
     .           , iostat=ioerr)        
           write(uinfor,'(//,a,//)')
     .            '**Reading station.info for antenna height'   
           sitcod = snames(istat) 
* MOD TAH 200203: Added AntDAZ to list of values from station.info
           call rstnfo(ustnfo,sitcod,iyr,jdoy0,isod
     .         , span, stanam,  anthgt, offstn, offste, antdaz
     .         , rcvcod, antcod
     .         , hgtcod, radome, swvers, rcvers, rcvnum, antnum
     .         , istarts, istops ) 
           stnfo_ok = .true.
         endif 

c        Receiver name -- If non-blank from X-file, assume passed
c            from original RINEX and copy to the new one.  If blank, create
c            from the firmware version and/or station.info

         if( rctype(1:2).eq.'  ' ) then  
    
c          Get the receiver code from the X-file firmware codes
           if (rcvrsw.eq.'COR' ) then
              rxxcod = 'TI4100'   
           else if (rcvrsw.eq.'GES' ) then
              rxxcod = 'TI4100'
           else if (rcvrsw.eq.'TRM') then
              rxxcod = 'TRMSST'
              do i=1,maxobt
                if(dattyp(i).eq.4) rxxcod = 'TRMSSE'
              enddo   
c             **Note: Since MAKEX renames C2 to P2, we cannot distinguish
c                     an SSE (or SSi) from a serial P-code SST (SLD?) 
c                     Is there a way based on firmware version?
           else if (rcvrsw.eq.'ASH') then
              rxxcod = 'ASHTEC'
           else if (rcvrsw.eq.'ROG' ) then  
              rxxcod ='ROGSNR'
c             cannot distinguish the various Rogue and Mini-Rogue receivers  
           else if (rcvrsw.eq.'TRB' ) then
              rxxcod = 'TRBROG'
           else if (rcvrsw.eq.'MIN' ) then
              rxxcod = 'MIN6AT'
           else if (rcvrsw.eq.'MAC' ) then
              rxxcod = 'MAC_II'
           else if (rcvrsw.eq.'SRT' ) then
              rxxcod = 'SRTR5S'
           else if (rcvrsw.eq.'SRN' ) then
              rxxcod = 'SRNR52'
c          No Leica receivers yet coded in MAKEX
           endif
c          compare with station,info if available
           if( stnfo_ok ) then
             if( rxxcod.ne.rcvcod ) then
               write(message,'(a,a6,a,a6,a)')
     .                   'Rcvr type inferred from X-file (',rxxcod
     .                  ,') differs from station.info (',rcvcod
     .                  ,') -- use the latter' 
               call report_stat('WARNING','XTORX','xtorx',' ',message,0)
              endif
              if( abs(swver-swvers).gt.0.0001 ) then
                write(message,'(a,f5.2,a,f5.2,a)') 
     .                   'Rcvr firmware from X-file (',swver
     .                  ,') differs from station.info (',swvers
     .                  ,') -- use the latter'     
               call report_stat('WARNING','XTORX','xtorx',' ',message,0)
               swver = swvers         
              endif           
           else
              rcvcod = rxxcod
           endif     
c          get the full receiver name from the rcvant.dat table   
           call read_rcvant(1,2,antcod,anttyp,radome,rcvcod,rctype
     .                      ,pcncod )
c          end if on blank receiver name from X-file
         endif
             
c       Firmware version -- GAMIT, not original version passed from X-file

c          write the (real) swver into the RINEX (character) variable  
           rcvers = ' '         
           write(rcvers(7:11),'(f5.2)') swver  
c          TI is the only case meriting a firmware name
           if( rcvcod.eq.'TI4100') then
              if( rcvrsw.eq.'COR' ) rcvers(1:4) = 'CORE'
              if( rcvrsw.eq.'GES' ) rcvers(1:5) = 'GESAR'
           endif 


c        Determine wavelength factors and observable types
c        --always the same for x-files except for C2/P2

         nobtyp = ndat	   
         do i=1,nobtyp
	   if( dattyp(i).eq.1 ) then
	     rxobtyp(i) = 'L1'
	   elseif( dattyp(i).eq.2) then
	     rxobtyp(i) = 'L2'
	   elseif( dattyp(i).eq.3 ) then
	     rxobtyp(i) = 'P1'
	   elseif( dattyp(i).eq.4) then
	      rxobtyp(i) = 'P2'
	   elseif( dattyp(i).eq.5) then
	      rxobtyp(i) = 'C1'
	   endif
	 enddo         
         nwave1 = iabs(lambda(1,1))
	 nwave2 = iabs(lambda(1,2))

c        set the number of valid observations for all observable types to be the same
c        (X-file makes no distinction)
         do i = 1,nchan
           do j = 1,nobtyp
             numobs(j,i) = numobsx(i)
           enddo
         enddo
              

c        Antenna name  - We must get this one right, so check and use the 
c                        station.info or X-file 6-char codes rather than 
c                        copying the 20-char description
          

         if( stnfo_ok ) then 
             if( antcodx.ne.antcod )  
     .          write(message,'(a,a6,a,a6,a)')
     .                   'Antenna type inferred from X-file (',antcodx
     .                  ,') differs from station.info (',antcod
     .                  ,') -- use the latter'   
               call report_stat('WARNING','XTORX','xtorx',' ',message,0)
         else
            antcod = antcodx
            write(*,'(/,a)') 
     .         'No station.info, antenna code taken from X-file'
         endif
c        get the full antenna name from table rcvant.dat
         call read_rcvant(1,1,antcod,anttyp,radome,rcvcod,rctype,pcncod)

          
c        Antenna offsets
          
         if( antcod.ne.'      ' ) then  
            anth = offarp(1)
            antn = offarp(2)
            ante = offarp(3) 
         else
            write(*,'(/,2a,//)') '**No antenna code on X-file and no'
     .                     , '  station.info file; set antenna ht = 0.'
            anth = 0.
            anth = 0.
            ante = 0.
         endif
         
       
c        Monument name and coordinates
   
         rxmrk = sitnam
         call geodms_xyz ( latd,latm,seclat,latflag,lond,lonm,seclon
     1                    , lonflag,height,apx,apy,apz )

          
c        RINEX program headers   
         
c        program name (XTORX), user, agency, date, RINEX version, and 
c        first two comment lines are set at beginning of run
c        observer info not yet passed from in X-file
         rxobs = ' '


c        Construct the RINEX file name
         jdoyrx = idoy(irxyr0,irxmo0,irxdy0)  
         read(daynum,'(i3)') jdoyx
         if( jdoyrx.ne.jdoyx ) then  
           write(message,'(a,i3,a,i3,a)') 
     .       'DOY from x-file name (',jdoyx
     .       ,') differs from actual start day ('
     .       ,jdoyrx,'), use x-file DOY for RINEX file name'
            call report_stat('WARNING','XTORX','xtorx',xfname
     .                      ,message,0)   
            jdoyrx = jdoyx
         endif
         stcode = xfname(2:5)  
         write (rfname,139) stcode,jdoyx,0,mod(irxyr0,100)
 139     format (a4,i3,i1,'.',i2.2,'o')
         if ( rfname(5:5) .eq. ' ' ) rfname(5:5) = "0"
         if ( rfname(6:6) .eq. ' ' ) rfname(6:6) = "0"
         call lowers(rfname)
c        open the RINEX file (output)  
c         --don't allow overwriting of RINEX files--too dangerous 
         if( fcheck(rfname) ) 
     .       call report_stat('FATAL','XTORX','xtorx',rfname
     .                      ,'RINEX file exists',ioerr)
         open(unit=urinex,
     1        file=rfname,
c*   2        status='unknown',
c*        Don't allow overwriting of RINEX files--too dangerous
     2        status='new',
     3        form='formatted',
     4        iostat=ioerr)
         if (ioerr.ne.0) then    
            write (uinfor,*) 'Error opening file : ', rfname
            call report_stat('FATAL','XTORX','xtorx',rfname
     .                      ,'Error opening file',ioerr)
         else   
            call report_stat('STATUS','XTORX','xtorx',rfname
     .                      ,'Opened file',0)
            write (uinfor,*)  'Opened file: ', rfname
         endif


c        Write the RINEX header

         call wrxhed (urinex,
     1      rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     2      rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz,
     3      anth,ante,antn,nwave1,nwave2,nobtyp,rxobtyp,rxint,
     4      irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     5      irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1,
     6      nchan,ischan,numobs )

c-----------------------------------------------------------------------


c *****  Loop over all data records of the X-file

         icount = 0
         imiss  = 0
         iremoved = 0

         do  i=1,nchan
           first(i) = .true.
         enddo

         do 170 idum=1,nepoch

c          read the epoch, # sats, and observation time
           read ( uxfile,
     1           fmt = '(/,2I4,I5,I4,2I3,F11.7)')
     2           iepoch,msat,iyr,jdoy,ihr,min,sec
           if( jdoy.gt.0 ) call fix_y2k(iyr)
                        
c          initialize variables
           do j = 1, nchan
              ier(j) = 1
              do i = 1, ndat
                 isnr(i,j) = 0
                 data(i,j) = 0.D0
              enddo
           enddo

c          read the data records
           do j=1,msat
              read( uxfile,
     1              fmt = '(10X,2I2,2(1X,D22.15,1X,I3),2(2X,D22.15))')
     2              ier(j),jchan(j),(data(i,j),isnr(i,j),i=1,2)
     3              ,(data(i,j),i=3,ndat)
              prn(j) = ischan(jchan(j))

c             put a bias flag at the channel when it first appears
c             pick the FIRST GOOD OBSERVATION..........
              if( first(jchan(j)) .and.  lgood(ier(j)) ) then
                if(put_bias_first .eq. 'yes' ) ier(j) = 10
c               code to remove the initial integer phase from all observations
c               at a given channel--implemented by rwk 970321
                 phs10(jchan(j))= dint( data(1,j) )
                 phs20(jchan(j))= dint( data(2,j) )
                first(jchan(j)) = .false.
              endif
              data(1,j) = data(1,j) - phs10(jchan(j))
              data(2,j) = data(2,j) - phs20(jchan(j))
           enddo


c          write the RINEX record if the epoch has good data
* MOD TAH 930127: Set the number of used satellites to zero
*          here incase their are no satellites at all and
*          the checking code below will not be executed.
           msat_rx = 0

           if( msat.gt.0 ) then

              icount= icount + 1
              if( iyr.gt.0 ) then
c             convert UTC to GPS time
                if(mtime.eq.1) itflag = -2
c             convert GPS to GPS time
                if(mtime.eq.2) itflag = -4
c               UTC yr,doy,hms to GPST wk,sow
                call timcon(itflag,iwkn,sow,iyr,jdoy,ihr,min,sec,utcoff)
              else
c               compute the observation time if it's not there
                delta = (iepoch-1)*dble(inter)
                call secsum( iwkn0,sow0,delta,iwkn,sow )
              endif
              itflag = +4
c             GPST wk,sow to GPST yr,doy,hms
              call timcon( itflag,iwkn,sow,iyr,jdoy
     1                   , irxhr1,irxmn1,rxsec1,utcoff)
              call monday( jdoy,irxmo1,irxdy1,irxyr1 )

* MOD TAH 920127: Scan all of the current xfile data at this epoch
*             and find out which PRN's have bad data.  NOTE:
*             msat_rx set to zero before check on msat so that
*             it will be zero if msat=0
              do j = 1, msat
                if( lgood(ier(j)) ) then
*                   This is a good observation so save in
*                   prn_rx list
                    msat_rx = msat_rx + 1
                    prn_rx(msat_rx) = prn(j)
                else
                    iremoved = iremoved + 1
                end if
              end do
*                         ! end if on msat > 0, Now check to
*                         ! see if msat_rx > 0
           end if

*          See if we still have data
           if( msat_rx.gt.0 ) then

c             set flag to OK
              idflag = 0

*             Change here to write out the number of good data
*             only  
              iyr2 = mod(irxyr1,100)   
              gnss = satnam(i)(1:1)
              write(urinex,145,iostat=ioerr)
     1              iyr2,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1
     2            , idflag,msat_rx,(gnss,prn_rx(i),i=1,msat_rx)
 145          format (5i3,f11.7,i3,14(a1,i2))


              do 150 j=1,msat

c                translate GAMIT flags into RINEX flags
                 lli=0
                 issi0 = 0
                 issi1 = isnr(1,j)
                 issi2 = isnr(2,j)

c----------------
* MOD TAH 930127 / RWK 930215:  Code to propagate GAMIT bias flags using
*     bit 2 of the RINEX loss-of-lock indicator has been removed for con-
*     formity with RINEX standards.  XTORX will simply not write out any
*     data that are not flagged as 'good'.  The bit-2 scheme worked out by
*     Burc Oral is preserved in the following comments:

c                if(lgood(ier(j))) then
c     dummy end if
c                if( 1 .eq. 1) then
c ** GAMIT bias flags ier      RINEX loss-of-lock indicator lli
c          (iggood =  0)    :  0, " "   OK
c          (igbias = 10)    :  1
c          (ignone =  1)    :  4   [set bit 2 ]
c          (igloel =  4)    :  4   [set bit 2 ]
c          (igchop =  2)    :  4   [set bit 2 ]
c          (igoutl =  6)    :  4   [set bit 2 ]
c          (igunwt = -1)    :  4   [set bit 2 ]
c          (iglamp =  3)    :  4   [set bit 2 ]
c          (ig2few =  5)    :  4   [set bit 2 ]

c                if (lbias(ier(j))) lli=1
c                if (lmarg(ier(j))) lli=4
c                if (ier(j) .eq. ignone ) lli=4
c                if (ier(j) .eq. igchop ) lli=4
c--------------

*                If we are mapping xfile amplitudes to rinex amplitudes
*                then call mapamp
                 if( imapamp.eq.1 ) then
                   call mapamp (2,rcvrsw,swver,1,dum8,isnr(1,j),issi1)
                   call mapamp (2,rcvrsw,swver,2,dum8,isnr(2,j),issi2)
c                  check that mapamp did the job right
                   if(issi1.lt.0.or.issi1.ge.10) issi1=1
                   if(issi2.lt.0.or.issi2.ge.10) issi2=1
                 endif

* MOD TAH 930127:
* The code below is not needed.  Both bad and marginal data
* will not be written to the rinex file.
C               else if (lmarg(ier(j))) then
C                issi1=1
C                issi2=1
C               endif

*               Check to see if we need to set the loss of lock flag
*               (X-file bias flag)

                 if( lbias(ier(j))) lli=1

c                if the phase observable is too large to fit, subtract a large no.
c                if the pseudorange is too large, set it to zero
                 do i = 3,4
c                   if the pseudorange is too large, set it to zero
                    if( dabs(data(i,j)).gt. 999999999.999d0 )
     .                          data(i,j)=0.d0
                 enddo

* MOD TAH 930127: Write out the X-file data to the rinex files.
*                Only write the line if the data is considered good.
                 if( lgood(ier(j)) ) then

*                    Copy the X-file data to the arrays to be written
*                    set the signal strength values.
                     do i = 1, nobtyp
                        obs_rx(i) = data(i,j)
                        iss_rx(i) = issi0

*                       check the rinex data types to see if phase.  If
*                       it is then the sign needs to be changed.
                        if( rxobtyp(i).eq.'L1' ) then
                           obs_rx(i) = -obs_rx(i)
                           iss_rx(i) = issi1
                        end if
                        if( rxobtyp(i).eq.'L2' ) then
                           obs_rx(i) = -obs_rx(i)
                           iss_rx(i) = issi2
                        end if
                     end do

*                    Now write the line out to the rinex file
                     write(urinex,146) (obs_rx(i),lli,iss_rx(i),
     .                               i = 1, nobtyp)
 146                 format (4(f14.3,i1.0,i1))
                 end if

 150         continue
c            End loop over satellites

          else

             imiss = imiss + 1

          endif

 170     continue
c        End loop over epochs


c        Write out the epochs and times written

         if(mtime.eq.1) itflag = -2
         if(mtime.eq.2) itflag = -4
         jdoyr0= idoy(irxyr0,irxmo0,irxdy0)
         call timcon ( itflag,iwkn0,sow0,irxyr0,jdoyr0
     1               , irxhr0,irxmn0,rxsec0,utcoff )
         jdoyr1= idoy(irxyr1,irxmo1,irxdy1)
         call timcon ( itflag,iwkn1,sow1,irxyr1,jdoyr1
     2               , irxhr1,irxmn1,rxsec1,utcoff )
         mepoch= istop-imiss

* MOD TAH 930127: The following check changed (it was wrong).
*        istop should be number of epoch with data and should
*        equal icount (accumulated as the xfile is read if
*        msat>0 at the epoch.
* SOMETHING STILL WRONG IN BOOK-KEEPING BUT SEEMS TO UNIMPORTANT
* ALL DATA NEEDED SEEMS TO GET THROUGH AND THE START AND STOP
* TIMES IN THE RINEX FILES ARE CORRECT.
C        if( icount.ne.mepoch ) then
C        if( icount.ne.istop  ) then
C          write(uscren,114) icount,istop
C114       format(1x,'Counted epochs (',i5,')  .ne. computed epochs (',
C    1                i5,')',/,'---something wrong')
C        endif
         write(uinfor,185) icount,istart,istop
     1                 , irxyr0,jdoyr0,irxhr0,irxmn0,rxsec0
     2                 , irxyr1,jdoyr1,irxhr1,irxmn1,rxsec1
 185     format(/,1x,'Copied ',i4,' epochs (',i4,'-',i4,')  Start: '
     1           ,3i4,i3,f7.3,/
     2           ,34x,'Stop: ',3i4,i3,f7.3)
         if( imiss.gt.0 ) write(uinfor,186) imiss
         if( iremoved.gt.0 ) write(uinfor,187) iremoved
 186     format(/,34x,i4,' empty X-file epochs skipped')
 187     format(33x,i5,' removed data due to marginal data')

c     end of loop over X-files
 200  continue  
      call report_stat('STATUS','XTORX','xtorx',' '
     .                ,'Normal end of XTORX',0)
      stop
      end

c======================================================================

      subroutine passheader(uxfile)
c     Pass the X-file header to position the pointer where starts the actual data
c     Written by C. Vigny 02 July 1991
      character*6 word
      integer*4 uxfile,ioerr             
      logical eof,found
                        
      eof =.false.
      found = .false.
      do while (.not.eof. and. .not.found)
        read(uxfile,'(a6)',iostat=ioerr) word
        if( ioerr.eq.-1 ) then
           eof = .true.
        elseif( ioerr.ne.0 ) then 
            call report_stat('FATAL','XTORX','xtorx/passheader'
     .        ,' ','Read error looking for EPOCH #',ioerr)
        else
          if( word.eq.' EPOCH' ) found = .true.
        endif
      enddo  
        if( .not.found ) call report_stat('FATAL','XTORX',
     .     'xtorx/passheader',' ','Cannot find EPOCH line',ioerr)
      return
      end

c=========================================================================

