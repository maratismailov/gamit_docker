      Subroutine MDMAKE( irun, bfil2, tfile, ifile, lfile, jfile
     ,                 , yfile, ffile, xfile, cfile1
     .                 , klock, lsite, idatum, ietide, year, lexp
     .                 , delmod, etidemod, ionsrc, magfield, gnss
     .                 , lsess, sfile, ierr_sestbl, ierr_sittbl
     .                 , atmlflag, scratch_dir, metsource )

C     subroutine to make MODEL batch file
C
C     S.Shimada              original at NRCDP
C     S.Shimada    01/26/90  modified at IGPP
C     S.Shimada    09/21/90  modified at IGPP
C     Y. Bock      12/90     modified at IGPP
C     MIT/Scripps mods from 12/90 : see FVERSN 
c     P. Tregoning/R.King  060815  mods for met models and data 
C
C     arguments
C        IRUN   : run mode
C                 =0 : for quick solution
C                 =1 : for full  solution
C        BFIL2  : secondary batch file name        input
C        TFILE  : T-file name                      input
C        IFILE  : I-file name                      input
C        LFILE  : L-file name                      input
C        JFILE  : SV clock polynomial file name    input
c        YFILE  : SV yaw file name                 input
c        FFILE  : F-file name (IONEX)              input
C        XFILE  : X-(C-)file name                  input
C        CFILE1 : new C-file name                  input
C        KLOCK  : clock model type                 input
C        LSITE  : unit number of site table        input
C        IDATUM : geodetic datum                   input
C        year   : 4-digit integer year             input
C        LEXP   : experiment type                  input
C        GNSS   : GNSS system (G R C E I)          input
C        IERR_SESTBL, IERR_SITTBL : return codes for sestbl. and sittbl. errors    output
c        atmlflag : char*1 (Y or blank) to indicate whether atmospheric loading
c                 should be applied                            
c        metsource : sestbl. line giving met-source hierarchy (c*20)
C
C     sample MODEL secondary batch file :
C        S                   Static Mode
C        pmath0.034          Print file
C        iscal0.034          Station clock polynomial (I-) file
C        lscal0.034          Coordinates (L) file
C        xmath0.034          Input X or C-file, or simulation control (S-file)
C        cmath0.034 /tmp     Output C-file / Scratch directory (blank if day directory)
C        N                   Delete input C-file?
C        tscal0.034          T-file
C        fscal0.034 IGRF12   Ionosphere source  (formerly 'I' for Inertial frame; 'E' not supported
C        wmath0.034          W-file (blank if none)
C        zmath0.034          Z-file (blank if none)
C        jscal0.034          Satellite clock polynomial (J-) file   
c        0  15 7 IERS92 Y N  Datum / Tides applied / SP EOP / E-tide model / Atm load / Hydrol load 
C        N     NONE NONE     Use site-specific antenna model (Y/N) / antenna model / SV antenna model
C        3                   Clock model / Yaw file
C        RNX ufile GPT 50    Met options (source hierarchy + humidity) or P T H
C        SAAS SAAS NMFH NMFW Met models (dryzen wetzen drymap wetmap)

c PT060612: mods to include in the model batch file the name of a rinex
c           file containing met temp, press, humidity. If requested, the name
c           will replace the values in the second-last line of the batch file.

      implicit none

      REAL   * 4  WETH(3)

      LOGICAL  reqd, in_sestbl, kbit, sittbl_met, last

      INTEGER* 4  year,lsite,lexp,irun,nout,idatum,lfntop
     .         ,  ill,ierr_sestbl,ierr_sittbl,i0,i1,nblen,klock
     .         ,  lsess,ietide,isptide,nyr2,ioerr,indx,imet,ival,i
      integer*4 jerr  ! IOSTAT error reading Earth Rotation sestbl. entry

      CHARACTER* 1  XORC, delmod, atmlflag, hydrflg, epcvmod, gnss
      character* 2  aerp,cyr,buf2 
      CHARACTER* 4  DZEN,WZEN,DMAP,WMAP,STNNAM,antmod,svantmod,buf4
     .           ,  ionsrc
      character* 4  metmod(4), sitcod    
      character* 5  cval
      character* 6  etidemod
c**   probably temporary:
      character* 6  magfield
      character*12  rnxmetfile
      character*16  bfil2,tfile,ifile,jfile,xfile,pfile,cfile1,yfile
     .             ,zfile,LFILE,sfile,xfile1,ffile
      character*20  metsource
      character*34  rnxmetline
      character*71  pthline     
      character*51  metcomment
      CHARACTER*80  line    
      character*80 scratch_dir
      character*256 message
      character*1  kinflg

c      External functions
     
      logical fcheck 
                      
c--------------------------------------------------------------------

c  Open the MODEL batch file and write a line of the primary batch file
                   

      OPEN( 21, FILE=BFIL2, STATUS='UNKNOWN' )
      WRITE( 17, '(A,A16)' )  'model  < ', BFIL2
                        
c  Write first line: GNSS code replaces obsolete SKD (static, kinematic, or dynamic)
      if( gnss.eq.' ' ) gnss = 'G'
      write( 21, '(a1,19x,a)') gnss,'GNSS code'

c  Get the names of the primary data files and write the next five lines

C     is it C-, or X-file ?
C     check first character of input file name
      I0 = LFNTOP( XFILE )
      I1 = NBLEN( XFILE )
      XORC = XFILE(I0:I0)  
      STNNAM = XFILE(I0+1:I0+4)
C     set up default P-file name
      IF( IRUN .EQ. 0 )  PFILE = 'p'//XFILE(I0+1:I1)
      IF( IRUN .EQ. 1 )  PFILE = 'p'//CFILE1(2:NBLEN(CFILE1))
c     if simulation, replace the x-file name with the s-file name  
      if( sfile(1:1).ne.' ') then
        if( sfile(2:5).eq.'site' ) then
           xfile1 = xfile 
           xfile1(1:1) = 's'
         else
           xfile1  = sfile
         endif
      else
         xfile1 = xfile
      endif 
      WRITE( 21, 100 )  PFILE, IFILE, LFILE, xfile1
 100  FORMAT ( A16, 4X, 'Print file', /,
     .         A16, 4X, 'Station clock polynomial (I-) file', /,
     .         A16, 4X, 'Coordinates (L) file',/,
     .         A16, 4X, 'Input X, C, or S file')

      if( scratch_dir(1:1).eq.' ') then
        write(21,'(a16,4x,a)') cfile1,'Output C-file'
      elseif( scratch_dir(2:4).eq.'tmp'.and.nblen(scratch_dir).eq.4 ) 
     .    then
        write(21,'(a10,1x,a4,5x,a)') cfile1(1:10)
     .    , scratch_dir(1:nblen(scratch_dir))
     .    , 'Output C-file / Scratch directory'
      else
        write(21,'(a10,1x,a,2x,a)') cfile1(1:10)
     .       ,scratch_dir(1:nblen(scratch_dir))
     .    ,'  Output C-file / Scratch directory'
      endif

C Write the c-file delete-input question line

      IF( XORC .EQ. 'x') then
         write( 21, '(A1,19X,A)' ) 'N', 'Delete input C-file?'
      ELSEIF( XORC .EQ. 'c') then
         if( delmod.eq.'Y' ) then
             write( 21, '(A1,19X,A)' ) 'Y', 'Delete input C-file?'
         else
             write( 21, '(A1,19X,A)' ) 'N', 'Delete input C-file?'
         endif
      endif
                 
c Write the t-file line

      WRITE( 21, '(A16,4X,A)' )  TFILE, 'T-file'

c Determine the source of ionospheric corrections (2nd and 3rd order only, for now)   
c     (NONE, blank, or 'I' [old-style batch file] means no corrections)
                                  
      if( ionsrc.eq.'GMAP'.and.ffile(1:1).eq.'f') then           
c**       temporarily allow (undocumented choice of dipole field (blank default is IGRF)
c**         write( 21, '(a10,10x,a)' )  ffile,'Ionosphere source'
          write( 21, '(a10,1x,a6,3x,a)' )  ffile,magfield
     .        ,'Ionosphere source'
      else
         write( 21, '(a)' ) 'NONE                Ionosphere source'
      endif
                                    

c Determine the source of met values
  
c        Read the sestbl for preferred sources    
      rnxmetline = '                    RINEX met file'  
      metcomment = 'Met options (source hierarchy + humidity) or P T H'
      pthline    = metsource//metcomment    
c       Tokens are hierrarchical, with GPT or STP followed by a humidity 
c       value; e.g. 'RNX ufile GPT 50'.  Look for RNX first for rnxmetline,
c       then copy other tokens  in order to the batch file (changing 'ufile'
c       to 'UFL'.  If RINEX met is specified but not as first choice, 
c       this scheme will not work.  GPT or STP must  appear in the line to 
c       get a graceful exit. 
      if( pthline(1:3).eq.'RNX' ) then
c       create the met-file name
        nyr2 = mod(year,100)
        write(buf2,'(i2.2)') nyr2
        read(buf2,'(a2)') cyr
        rnxmetfile=cfile1(2:5)//cfile1(8:10)//'0.'//cyr//'m'   
        if( fcheck(rnxmetfile) ) then
          rnxmetline(1:12) = rnxmetfile  
        else                       
          sitcod = rnxmetfile(1:4)
          call uppers(sitcod)
          write(message,'(a,a4)') 
     .       'No RINEX met file available for ',sitcod
          call report_stat('WARNING','FIXDRV','mdmake',' ',message,0)
        endif
      endif
      do i=1,4
       metmod(i) = ' '
      enddo
      indx = 1           
      imet = 0 
      last = .false.     
      do while (.not.last)                             
        call read_line(metsource,indx,'CH',ioerr,ival,cval)   
        if( ioerr.ne.0 ) then
           call report_stat('WARNING','FIXDRV','mdmake'
     .     ,' ','Error decoding  Met Obs Source entry in sestb',ioerr)   
           call report_stat('FATAL','FIXDRV','mdmake'
     .     ,' ','Must include GPT or STP as last value',ioerr)   
c**        elseif( cval(1:3).eq.'RNX') then 
c**          continue    *** rwk 060826: keep RNX even though redundant
        else               
          if( cval(1:5).eq.'UFILE' ) cval(1:5) = 'UFL  '
          imet = imet + 1
          metmod(imet)(1:nblen(cval)) = cval(1:nblen(cval))
          if( metmod(imet).eq.'GPT '.or.metmod(imet).eq.'STP ') 
     .        last = .true.   
        endif
      enddo
c     this read should be the humidity value--store it also as a charcter string
      call read_line(metsource,indx,'CH',ioerr,ival,cval)
      imet = imet + 1
      metmod(imet)(1:nblen(cval)) = cval(1:nblen(cval))
c     now write the model tokens (except RNX) into a string written later to 
c     line 16 of the batch file
      write(pthline(1:20),'(5a4)') (metmod(i),i=1,imet)
c     If values present in the sittbl., they override the sestbl commands
      reqd = .false.
      call rdsitt(stnnam(1:4),10,'MET. VALUE',nout,line,lsite, reqd,ill)
      if( ill.ne.0 ) ierr_sittbl = ill
      if( line(1:1).ne.' ') then
        READ( LINE(1:NOUT),*,iostat=ioerr )  WETH  
        if( ioerr.ne. 0 ) then 
          write(message,'(2a)') 
     .      'Error reading P, T,H from sittbl.; line=',line
          call report_stat('FATAL','FIXDRV','mdmake',' ',message,ioerr)
        else
         call report_stat('WARNING','FIXDRV','mdmake',' '
     .     ,'Sittbl entries for met data will override sestbl models',0)
         if( weth(3).gt.99.9 ) then
          write(pthline(1:17),'(f7.2,1x,f4.1,1x,f4.0)') (weth(i),i=1,3)
         else
          write(pthline(1:17),'(f7.2,1x,f4.1,1x,f4.1)') (weth(i),i=1,3)
         endif
        endif
      endif

                      
c  Write the input met (RINEX) file line

      write(21,'(a)') rnxmetline


c  Write the ouput met (Z-file) file line

      reqd = .false.       
      call rdsest( 10,'Output Met',4,buf4,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      if( buf4(1:1).eq.'Y' ) then 
        zfile = xfile
        zfile(1:1) = 'z'
        write( 21, '(a16,4x,a)' )  zfile,'Z-file' 
      else
         write(21,'(a1,19x,a)' )  ' ', 'Z-file'
      endif
        

c  Write the satellite clock (J-) file line  

       CALL LOWERS(JFILE)
       IF( INDEX(JFILE,'jnone') .EQ. 0 )  THEN
          WRITE( 21, '(A16,4X,A)' )
     .        JFILE, 'Satellite clock polynomial (J-) file'
       ELSE
         WRITE( 21, '(A5,15X,A)' )
     .       'jnone', 'Satellite clock polynomial (J-) file'
       ENDIF

               
c  Write the line for geodetic datum, tide model, short-period earth-rotation,
c  and atmospheric and hydrological loading      

c     datum and tide model read in bmake
* MOD TAH 200505: default of IERS2010 + UT1LIB signal.  After
*     repro3 change to Desai and Sibios ( 3+16+64 )  
      isptide = 11 + 16
      call rdsest( 14,'Earth Rotation',2,aerp,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      if( aerp(1:1).ne.' ' ) then 
* MOD TAH 200504: Updated processing of options to make sure a 
*        unqiue model is selected and setting the default of
*        IERS2010 + UT1LIB signal.
         read(aerp,*,iostat=jerr) isptide
         if ((kbit(isptide,1)).or.(kbit(isptide,2))) then
* MOD TAH 200504: Added Desai and Sibious model (bit 7).
*           Adopt the latest model and set all lower bits to zero.
            if( kbit(isptide,7)) then   !Desai and Sibous
                call sbit(isptide,3,0)
                call sbit(isptide,4,0)
                call sbit(isptide,6,0)
            elseif( kbit(isptide,6) ) then
                call sbit(isptide,3,0)
                call sbit(isptide,4,0)
            elseif( kbit(isptide,4) ) then
                call sbit(isptide,3,0)
            else   ! Old Ray model; force to IERS2010
                call report_stat('WARNING','ORBITS','trot',' ',
     .             'old diurnal EOP model--override with IERS10 model',
     .             0)
* MOD TAH 200504: Use sbit to set IERS model as default.  This might
*                 change later.
                  call sbit(isptide,4,1)
            endif
         else
            isptide = 0   ! No model to apply.
         endif
      else  
*        No string returned.  Set default
         isptide = 11 + 16  ! IERS2010 + UT LIBRATION terms
      endif

* 
      if( atmlflag.eq.' ' ) atmlflag = 'N'
c     hydrological loading not yet supported
      hydrflg = 'N'
* MOD TAH 200219: Changed format to allow for larger ietide value
*     when IERS20 mean pole used I2->I3.  Change read_batch in model
*     to read with free fomtat to handle old and new formats.
      write( 21,'(i1,1x,i3,1x,i2,1x,a6,1x,a1,1x,a1,1x,2a)')
     .     idatum,ietide,isptide,etidemod,atmlflag,hydrflg
     .    ,'Datum / Tides applied / SP EOP / E-tide model'
     .    ,' / Atm load / Hydrol load'

         
c Write the line for receiving and transmitting antenna phase-center models
                      
c     The first argument was formerly the elevation cutoff angle, but since this
c     is no longer used in MODEL, we've supplanted it by a Y/N question of whether
c     there will be a site-specific ANTEX file read.
c     set default
      epcvmod = 'N'
      antmod = 'AZEL' 
      svantmod = 'ELEV'
      reqd = .false.                   
      in_sestbl = .false.        
c     read from the sestbl. whether a site-specfic antenna model will be used
      call rdsest( 27,'Use site-specific antenna model',1,epcvmod
     .           , lsess,reqd,ill) 
c     read the ground antenna PCV model for all sites from sestbl.
      call rdsest( 13,'Antenna Model',4,buf4,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      if( buf4(1:1).ne.' ' ) then
        antmod = buf4
        in_sestbl = .true.
      endif
C     override with sittbl. values if present
      reqd = .false.
      call rdsitt( stnnam(1:4),4,'APHS',nout,line,lsite,reqd,ill)
      if( ill.eq.0 ) then
         if( nout.gt.0 ) then
           read( line(1:nout),'(a4)' ) buf4
           if( buf4(1:1).ne.' ' ) then
              if( in_sestbl .and. buf4.ne.antmod ) then 
                write(message,'(a,a4,a,a4,a,a4)')
     .            'Sestbl antenna phase model (',antmod
     .             ,') overridden by sittbl entry (',buf4,') for '
     .             ,stnnam(1:4)   
                call report_stat('WARNING','FIXDRV','mdmake',' '
     .                          , message,0)
               endif
               antmod = buf4
           endif
         endif  
      else
         ierr_sittbl = ill
         return
      endif
c     read sv ant model for all sites from sestbl.
      call rdsest( 16,'SV Antenna Model',4,buf4,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      if( buf4(1:1).ne.' ' ) then
        svantmod = buf4
      endif
      WRITE( 21,'(a1,5x,a4,1x,a4,5x,2a)')  epcvmod,antmod,svantmod
     .  ,'Use site-specific antenna model (Y/N) / antenna model / '
     .  ,' SV antenna model'

c  Write the receiver clock model and y-file line

c     receiver clock poly and epoch by epoch
c        1       0                      poly        [Minimacs]
c        2       poly                   poly        [Atomic clock or good crystal (TI4100)]
c        3       poly + epoch-by-epoch  poly        [All modern receivers]
      write(21,'(i1,1x,a16,2x,a)') klock, yfile,'Clock model / Yaw file'


C  Write the met-model line
                           
      write(21,'(a)') pthline


c  Write the line for zenith delay models and mapping functions
                              
c       no sestbl entry for zenith models, hardwire to Saastamoinen unless UFL for zenith and VMF1 for mapping
      dzen = 'SAAS'
      wzen = 'SAAS'   
      call rdsest(4,'dmap',4,dmap,lsess,reqd,ill)
      if( ill.ne.0 ) ierr_sestbl = ill
      if( dmap(1:1).eq.' '.or.dmap(1:3).eq.'GMF' ) dmap = 'GMFH'
      call rdsest(4,'wmap',4,wmap,lsess,reqd,ill)
      if( ill.ne.0 ) ierr_sestbl = ill
      if( wmap(1:1).eq.' '.or.wmap(1:3).eq.'GMF' ) wmap = 'GMFW'
      if( metmod(1)(1:3).eq.'UFL'.and.dmap.eq.'VMF1' ) then
        dzen = 'VMF1'
        wzen = 'VMF1'
      endif
c       look for sittbl entries    
      sittbl_met = .false.
      call rdsitt( stnnam(1:4),4,'DZEN',nout,buf4,lsite,reqd,ill )
      if( ill.ne.0 ) ierr_sittbl = ill
      if( buf4(1:1).ne.' ') then
        sittbl_met = .true.
        dzen = buf4  
      endif
      call rdsitt( stnnam(1:4),4,'WZEN',nout,buf4,lsite,reqd,ill )  
      if( ill.ne.0 ) ierr_sittbl = ill
      if( buf4(1:1).ne.' ') then
        sittbl_met = .true.
        wzen = buf4  
      endif
      call rdsitt( stnnam(1:4),4,'DMAP',nout,buf4,lsite,reqd,ill ) 
      if( ill.ne.0 ) ierr_sittbl = ill 
      if( buf4(1:1).ne.' ') then
        sittbl_met = .true.
        dmap = buf4  
      endif
      call rdsitt( stnnam(1:4),4,'WMAP',nout,buf4,lsite,reqd,ill )  
      if( ill.ne.0 ) ierr_sittbl = ill 
      if( buf4(1:1).ne.' ') then
        sittbl_met = .true.
        wmap = buf4  
      endif
      if( sittbl_met ) call report_stat('WARNING','FIXDRV','mdmake',' '
     .    , 'One or more sittbl entries for zen/map override sestbl',0)
      write( 21, '(4(a4,1x),a)' )  dzen, wzen, dmap, wmap  
     .         , 'Met models (dryzen wetzen drymap wetmap)'

      close( 21 )
      return     

      end
