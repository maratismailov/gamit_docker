       Subroutine OPEN( trash ,scratch_dir )

c     Open all of input and output files.   R. King 8 April 1987.
c     Update for IGRF11 Remove dipole variable EJP/Andrzej A. 21 Sept 2010

      implicit none    

      include '../includes/dimpar.h'
      include '../includes/units.h'
      include '../includes/model.h' 

      LOGICAL fcheck,stnfo_exist,sesfo_exist
c**     LOGICAL lask
C
      CHARACTER*1 CFILE,XFILE,sfile,TRASH
      CHARACTER*5 NONE  
      CHARACTER*16 TMPNAM,UNAME,eqrnfil
      character*1  upper1            
      CHARACTER*16 NUTFIL,UT1FIL,POLFIL,datfil,hifil,pcvfil
      character*16 solfil,lunfil,dcbfil
      CHARACTER*22 BUF22,DATE
      CHARACTER*40 VERSN
      CHARACTER*80 which_line,scratch_dir
      character*128 tmp_cfile_name,cvalue
      character*256 message,inline 

       
c     function
      character*4 lowerc
                                                  
      real*8 value

      INTEGER*4 IOERR,DAYOYR
     .        , IYEAR,IHR,IMN,M,IHNSEC,ISEC,IDAY,IMONTH,nblen
     .        , trimlen,indx,getpid,getuid,pid,uid 

      data nutfil/'nutabl.         '/
      data ut1fil/'ut1.            '/
      data polfil/'pole.           '/
      data solfil/'soltab.         '/
      data lunfil/'luntab.         '/
      data hifil/'hi.dat          '/
      data datfil/'gdetic.dat     '/
      data pcvfil/'antmod.dat      '/    
      data dcbfil/'dcb.dat      '/    
      data eqrnfil/'eq_rename'/
      data none /'none'/
* MOD TAH 060317: Changed cfile to C from c because test below is casefolded up.
* MOD RWK 120727: Changed back for consistency, changing the test below
      data cfile,xfile,sfile /'c','x','s'/
c     the following may be necessary if the routine is converted to uppercase
c      call uppers(cfile)
c      call uppers(xfile)
c      call uppers(sfile)
c      call uppers(earthf)

C
C         Unit Numbers of Files:  10  Receiver-clock (I-File)    iui
C                                 11  Ephemeris (T-) File        iut
C                                 12  Input X, C, or S File      iobs
C                                 13  Output C-File              iuc
C                                 20  Nutation Table             inut
C                                 21  UT1 Table                  iut1
C                                 22  Pole Position Table        ipole
C                                 23  Input Met Table (W-File)   iuw
C                                 24  Output Met Table (Z-File)  iuz
C                                 25  Geodetic Datum (D-File)    iud
C                                 26  SV Clock File (J-File)     iuj
c                                 27  station.info               istnfo
C                                 28  Tides and loading (U-file) iuu 
c                                 29  session.info               isesfo
c                                 30  Archive (print) File       iprnt
C                                 31  Site coordinates(L-file)   iul
c                                 32  Ant PCV model (antmod.dat) ipcv
C                                 33  luntab.                    ilun
c                                 34  soltab.                    isun
c                                 36  n-body ephemerides         ibody 
c                                 37  SV attitude (y-file)       iuy   
c                                 38  simulation displacement file - may be opened in simred   
c                                 41  GPT2 grid file             igpt2 (opened in setup)
c                                 42  Diff. code biases (dcb.dat) idcb   
c                                 43  autcln.cmd                  (opened in setup)  
c                                 44  Antenna mechanical offsets (hi.dat) ihi
c                                 45  Met rinex file             (opened in metmod)   
c                                 46  Ocean tide center-of-mass correction (opened in /lib/otlcmc) 
c                                 47  Coordinate rename file (eq_rename)  
c                                 48  IONEX (F-) file             iuf 
c                                 49  Site-specific antenna model iepcv
c
c
c     ITERM=  5  - assigned in MODEL
c     ISCRN=  6  - assigned in MODEL
      iui  = 10
      IUT  = 11
      IOBS = 12
      IUC  = 13
      INUT = 20
      IUT1 = 21
      IPOLE= 22
      IUW  = 23
      IUZ  = 24
      IUD  = 25
      IUJ  = 26
      istnfo  = 27  
      iuu  = 28
      isesfo= 29 
      iprnt= 30
      iul  = 31
      ipcv = 32
      ilun = 33
      isun = 34  
      ibody = 35 
      iuy = 37     
      iudcb = 42 
      ihi = 44 
      iueqrn = 47                                                      
      iuf = 48
      iepcv = 49 
      
c  Open the output archive (print) file

      OPEN (UNIT=IPRNT,FILE=pfiln,FORM='FORMATTED',status='unknown'
     .     ,iostat=ioerr)
      if ( ioerr .ne. 0 ) then
        call report_stat('FATAL','MODEL','open',pfiln
     .                  , 'Error opening printfile: ',ioerr)
      endif  
      call uppers(pfiln)
c     find the period in the file name
      M=INDEX(pfiln,'.')
      SITECD=pfiln(M-5:M-1)
      READ (pfiln(M+1:M+3),'(I3)') DAYOYR
      call lowers(pfiln)

c  Write the site name and MODEL version number to the status, warning, and p-files

      call mversn(versn)
      write(message,'(a,a4,a,a40)') 'Site ',sitecd
     .                             ,': Started MODEL version ',versn
      call report_stat('STATUS','MODEL','open',' ',message,0)
      call report_stat('WARNING','MODEL','open',' ',message,0)
      write(iprnt,'(a,/,a)')
     .     '----------------------------------------------------'
     .    ,' ** MODEL VERSION AND INPUT FILES SUMMARY **'
      write(iprnt,'(//,10x,2a,//)') 
     .     'Program MODEL Version ',versn(1:nblen(versn))
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      call getusr(uname)
      write(iprnt,'(a,i4,a,i3,a,i3,2x,i2,a,i2,a,i2,a,a16)')
     .  ' MODEL Run on ',iyear,'/',imonth,'/',iday,ihr,':'
     .    ,imn,':',isec,' by ',uname
      write(buf22,'(i4,a,i2,a,i3,1x,i2,a,i2,a,i2)')
     .    iyear,'/',imonth,'/',iday,ihr,':',imn,':',isec
      read(buf22,'(a20)') date
      ntext = ntext + 1
      text(ntext) = 'Model '//versn
      ntext = ntext + 2
      text(ntext)=' Run on '//date//' by '//uname
    

c  Open the Receiver-clock (I-) File for the Experiment
c  temporarily allow for S-file instead
          
      if( ifiln.ne.'NONE' ) then 
        OPEN (UNIT=iui,file=ifiln,STATUS='OLD',iostat=ioerr)
        if ( ioerr .ne. 0 ) then
          call report_stat('FATAL','MODEL','open',ifiln
     .                    ,'Error opening ifile: ',ioerr)
        endif
        if( ifiln(1:1).eq.'s' ) then
          write(iscrn,'(/,a)') 
     .      'Warning: S-file rather than I-file input'
          write(iprnt,'(/,a,a16)') ' Session (S-) File :   ',ifiln
        else    
          write(iprnt,'(/,a,a16)') ' Receiver-clock (I-) File : ',ifiln
        endif
      else
        iui = 0
      endif

c  Open the coordinate file 

      OPEN (UNIT=iul,file=lfiln,STATUS='OLD',iostat=ioerr)
      if ( ioerr .ne. 0 ) then
        call report_stat('FATAL','MODEL','open',lfiln
     .                  ,'Error opening lfile: ',ioerr)
      endif  
      write(iprnt,'(a,a16)') ' Site Coordinate File     : ',lfiln
c     see if l-file or apr file and set the flag in model.h
      call crd_file_type( lfiln, kfflg )
        
c  Open the coordinate rename file 

      if( fcheck( eqrnfil ) ) then
        open(unit=iueqrn,file=eqrnfil,STATUS='OLD',iostat=ioerr)
        if ( ioerr .ne. 0 ) then
          call report_stat('FATAL','MODEL','open',eqrnfil
     .                    ,'Error opening rename lfile: ',ioerr)
        endif   
        write(message,'(a,a16)' ) ' Site rename File         : ',eqrnfil
        call report_stat('STATUS','MODEL','open', ' ',message,0)
        write(iprnt,'(a,a16)') ' Eq-rename file           : ',eqrnfil
      else
c       set unit number in model.h to zero to indicate not available
        iueqrn = 0   
      endif


c     Open the Input Observation File
                                      
      if (fcheck(obfiln)) then
         IF (obfiln(1:1).eq. cfile) then
            CALL COPENS (obfiln,'OLD',IOBS,IOERR)
         else if (obfiln(1:1).eq.xfile .or. obfiln(1:1).eq.sfile ) then
            OPEN (UNIT=IOBS,FILE=obfiln
     .           , FORM='FORMATTED',iostat=ioerr,STATUS='OLD')
         ELSE
            WRITE(ISCRN,'(A)') 'Input observation file of unknown type'
         ENDIF
         if (ioerr .eq. 0) then
          write(message,'(a,a16)') ' Input Observation File   : ',obfiln
          call report_stat('STATUS','MODEL','open', ' ',message,0)
          write(iprnt,'(a)') message
         else
           call report_stat('FATAL','MODEL','open',obfiln
     .                ,'Error opening input observation file: ',ioerr)
         endif
      else
         call report_stat('FATAL','MODEL','open',obfiln
     .                   ,'Could not find input file: ',ioerr)
         write (iscrn,*) 'Could not find inputfile: ',obfiln
      endif

c  Open the Output C-file 
                
c     Create the scratch c-file name (use uid and pid to make it unique)
      if ( scratch_dir .ne. 'NONE' ) then 
        pid = getpid()
        uid = getuid()
        write(tmp_cfile_name,'(a,a1,a,a1,i6.6,i5.5)')
     .     scratch_dir(1:trimlen(scratch_dir)),'/',
     .     cfiln(1:trimlen(cfiln)),'.', pid,uid
      else
        tmp_cfile_name = cfiln
      endif    
      call copens (tmp_cfile_name(1:trimlen(tmp_cfile_name)),
     .             'unknown',iuc,ioerr)
      if (ioerr.eq.0 ) then    
         write(message,'(a,a)' ) ' Output C-file            : '
     .      ,tmp_cfile_name(1:trimlen(tmp_cfile_name))
         call report_stat('STATUS','MODEL','open', ' ',message,0)
         write(iprnt,'(a)') message
      else
         call report_stat('FATAL','MODEL','open',tmp_cfile_name,
     .   'Error opening output C-file: ',ioerr)
      endif
C Debug EJP Dec 2007
C      Print*,'MODEL/open: output c-file name', cfiln, 
C     . tmp_cfile_name
             

c  Open the Ephemeris (T-) File

      call topens (tfiln,'old',iut,ioerr)
      write(message,'(a,a16)' ) ' Ephemeris (T-) File      : ',tfiln
      call report_stat('STATUS','MODEL','open', ' ',message,0)
      write(iprnt,'(a)') message


c  Open the IONEX files
                            
      if( ionfiln(1:4).eq.'NONE' .or. ionfiln(1:2).eq.'  ' .or.
     .    ionfiln(1:2).eq.'I ') then
c       set unit number to zero to indicate no ion file
        iuf = 0      
        ionsrc = 'NONE'
      elseif( ionfiln(1:1).eq.'f') then
        open(unit=iuf,file=ionfiln,status='old',iostat=ioerr)
cd        Print*,'Model\open iuf',iuf,'ioerr',ioerr
        if( ioerr.ne.0 ) then
          call report_stat('FATAL','MODEL','open',ionfiln
     .                    , 'Error opening IONEX file:',ioerr)
         else
          write(message,'(a,a16)') ' IONEX File              : ',ionfiln
          write(iprnt,'(a)') message
        endif 
      endif            

c  Open the UT1 table file

      OPEN (UNIT=IUT1,FILE=UT1FIL,STATUS='OLD',iostat=ioerr)
      if (ioerr .ne. 0) then
         call report_stat('FATAL','MODEL','open',ut1fil,
     .      'error opening UT1 table:',ioerr)
      endif


c  Open the Nutation table file
                           
c     open only if it is needed (no 'nbody' for lunar-solar ephemerides)
      if( .not.fcheck('nbody') ) then
        OPEN (UNIT=INUT,FILE=NUTFIL,STATUS='OLD',iostat=ioerr)
        if (ioerr .ne. 0) then
           call report_stat('FATAL','MODEL','open',nutfil,
     ,     'error opening nutation table:',ioerr)
        endif 
      endif


c  Open the Polar Motion file

      OPEN (UNIT=IPOLE,FILE=POLFIL,STATUS='OLD',iostat=ioerr)
      if (ioerr .ne. 0) then
         call report_stat('FATAL','MODEL','open',polfil,
     .     'error opening nutation polar motion:',ioerr)
      endif
                        

c  Open the Ocean tide table file
                       
c     since we don't know yet whether ocean loading, atm loading, u-file met will
c     be used, don't require this file, but set the unit # to zero if it is missing
c*      tmpnam = lfiln 
c*  add this temporarily to avoid l-file naming change  --rwk 001220
      tmpnam = ifiln   
c*  this doesn't always work either since the ifile may also have 6th char = a
c*  put in a temporary trap, but move this code later in the MODEL eventually
      if( tmpnam(6:6) .eq. 'a' ) 
     .    call report_stat('FATAL','MODEL','open',' ',
     .'Temporary logical bug in naming u-file from i-file--see Bob King'
     .  ,0)                                                                  
c*  rwk 130213: For simulations, i-file is normally NONE, so use the s-file name
      if( upper1(obfiln(1:1)).eq.upper1('S') ) tmpnam = obfiln
      tmpnam(1:1) = 'u'
      call lowers(tmpnam) 
      if(fcheck(tmpnam)) then
        open (unit=iuu,file=tmpnam,status='old',iostat=ioerr)
        if( ioerr.eq.0 ) then 
         write(message,'(a,a16)' ) ' Loading/Met (U-) File    : ',tmpnam
         call report_stat('STATUS','MODEL','open', ' ',message,0)
         write(iprnt,'(a)') message
        else
          call report_stat('FATAL','MODEL','open',tmpnam
     .               ,'error opening ocean loading input file:',ioerr)  
        endif
      else
         iuu = 0
      endif  

c        Open the antenna offset file (hi.dat)
      OPEN (UNIT=ihi,FILE=hifil,STATUS='OLD',iostat=ioerr)
      if (ioerr .ne. 0) then   
         call report_stat('FATAL','MODEL','open',hifil,
     .     'error opening hi.dat file:',ioerr)
      endif    
 
             
c  Open the antenna phase centre model file

      OPEN (UNIT=ipcv,FILE=pcvfil,STATUS='OLD',iostat=ioerr)
      if (ioerr .ne. 0) then
         call report_stat('FATAL','MODEL','open',pcvfil,
     .   'error opening antmod.dat file:',ioerr)
      endif 

                                
c  Open the empirical (site-specific) antenna model file
      if( upper1(epcvflg).eq.'Y' ) then
        epcvfiln = 'antmod.'//lowerc(sitecd)
        OPEN (UNIT=iepcv,FILE=epcvfiln,STATUS='OLD',iostat=ioerr)
        if (ioerr .ne. 0) then
           call report_stat('FATAL','MODEL','open',epcvfiln,
     .     'error opening site-specific ANTEX file:',ioerr)
        endif                  
      else
        iepcv = 0 
      endif
          

c  Open the pseudorange correction file

      OPEN (UNIT=iudcb,FILE=dcbfil,STATUS='OLD',iostat=ioerr)
      if (ioerr .ne. 0) then
         call report_stat('WARNING','MODEL','open',dcbfil,
     .   'error opening C/P code bias file:',ioerr)
      endif
                       

c  Open the Solar table file
                             
      if( .not.fcheck('nbody') ) then
        OPEN (UNIT=isun,FILE=SOLFIL,STATUS='OLD',iostat=ioerr)
        if (ioerr .ne. 0) then
           call report_stat('FATAL','MODEL','open',solfil,
     ,     'error opening solar ephemeris file:',ioerr)
        endif
      endif

c  Open the Lunar table file
           
      if( .not.fcheck('nbody') ) then
        OPEN (UNIT=ILUN,FILE=LUNFIL,STATUS='OLD',iostat=ioerr)
        if (ioerr .ne. 0) then
           call report_stat('FATAL','MODEL','open',lunfil,
     1     'error opening lunar ephemeris file:',ioerr)
        endif
      endif 


C  Open the Datum file

      OPEN (UNIT=IUD,FILE=datfil,STATUS='OLD',iostat=ioerr)
      if (ioerr .ne. 0) then
         call report_stat('FATAL','MODEL','open',datfil,
     1   'error opening datum file:',ioerr)
      endif


c  Open station.info

      inquire( file='station.info', exist=stnfo_exist )
      if( stnfo_exist ) then
        OPEN (UNIT=istnfo,FILE='station.info',STATUS='OLD',iostat=ioerr)
         if (ioerr .ne. 0) then
           call report_stat('FATAL','MODEL','open','station.info',
     .     'Error opening station.info file:',ioerr)
        endif
      else     
         write(message,'(2a)')
     .     'No station.info files antenna offsets taken from'
     .   , ' the X-file:'
         call report_stat('WARNING','MODEL','open','station.info'
     .                   ,message,ioerr)
c        set unit to zero as a flag for SETUP not to read (this is easy but a little dangerous)
         istnfo = 0
      endif
  
 
c Open session.info for simulations

      if (obfiln(1:1).eq.sfile ) then
       inquire( file='session.info', exist=sesfo_exist )
       if( sesfo_exist ) then
        OPEN (UNIT=isesfo,FILE='session.info',STATUS='OLD',iostat=ioerr)
        if (ioerr .ne. 0)  call report_stat('FATAL','MODEL','open'
     .    ,'session.info','Error opening session.info file:',ioerr)
        else 
         call report_stat('FATAL','MODEL','open','session.info.',
     .   'No session.info for simulation',ioerr)
        endif
      endif  
               
             
c   Open the input met (weather) file (now required to be a RINEX m-file)
         
      if( metfiln(1:1) .ne. ' ' ) then
        open(unit=iuw,file=metfiln,status='old',iostat=ioerr)
        if (ioerr .ne. 0) then
          call report_stat('FATAL','MODEL','open',metfiln
     .                     ,'Error opening input met file:',ioerr)
        else
         write(message,'(a,a16)') ' RINEX Met File           : ',metfiln
         call report_stat('STATUS','MODEL','open', ' ',message,0)
         write(iprnt,'(a)') message 
        endif  
      else
c       line is blank, set iuw=0 to indicate no met file available
        iuw = 0 
      endif 


c  Open the output (Z-) met print file
                  
      if( zfiln(1:3).ne.'   '.and.zfiln(1:3).ne.'N  ')  then  
c       test on 'N   ' is to allow use of old-style batch files
        open(unit=iuz,file=zfiln,status='unknown',iostat=ioerr)
        if (ioerr .ne. 0) then
          call report_stat('FATAL','MODEL','open',zfiln,
     .     'Error opening met print (z-) file:',ioerr)   
        else       
         write(message,'(a,a16)' ) ' Output Met (Z-) File     : ',zfiln
         call report_stat('STATUS','MODEL','open', ' ',message,0)
         write(iprnt,'(a)') message
        endif   
      else
c       use unit number 0 to tell atmdel.f there is no output file
        iuz = 0
      endif

                                 
c  Open the SV Clock (J-)file

      call lowers(jfiln)
      call lowers(none)
      IF(index(jfiln,none) .gt. 0) THEN
         IUJ=0
      ELSE
         OPEN(UNIT=IUJ,FILE=jfiln,STATUS='OLD',iostat=ioerr)
         if ( ioerr .ne. 0 ) then
            call report_stat('FATAL','MODEL','open',jfiln,
     1      'Error opening J-file:',ioerr)
         endif
         write(iprnt,'(a,a16)') ' SV clock (J-) File       : ',jfiln
         write(iprnt,'(/)')
      ENDIF                       

      RETURN
      END
