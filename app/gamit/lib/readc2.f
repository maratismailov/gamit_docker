      subroutine readc2 (lunit
     .,                sitnam,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,                npart,norb,gnss,nsat,isprn,fl1,fl2
     .,                ndat,dattyp,lambda
     .,                skd,nepoch,inter,ircint,mtime,isessn
     .,                iy,im,id,ihr,min,sec
     .,                offarp,offsl1,offsl2,antdaz, svantdx
     .,                obfiln,tfiln,jfiln
     .,                frame,precmod,nutmod,gravmod,srpmod
     .,                eradmod,antradmod
     .,                ietide,isptide,speopmod
     .,                etidemod,otidemod,atmtide,atmlmod,hydrolmod
     .,                atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap  
     .,                ionsrc,magfield
     .,                antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,                elvcut,nclock,clock
     .,                jdet,tet,jdr,tr
     .,                ut1,ut1dot,xp,xpdot,yp,ypdot
     .,                psi,psidot,eps,epsdot
     .,                avlmet
     .,                nslip,islip,islpst
     .,                niextra,iextra,nrextra,rextra,ncextra,cextra)
                               
c     Read the second part of a C-file which is a header

c     K. Feigl/R. King  July 1989    
c     R. King Sep 1998: Modfied to add 'norb' parameter   
c     R. King Feb 2005  Modified for version 10.2 to add more model names, avg loading 
c                       ocean/atm loading values, and more extra slots    
c     R. King Aug 2010  Modified for release 10.4 to lengthen the names of the
c                       antenna models 
c     R. King Jan 2013  Modified for release 10.41 to add dryzen and wetzen 
c     R. King Mar 2014  Added eradmod,antradmod,ionsrc,magfield
c     R. King Oct 2014  Added gnss, fL1, fL2, replaced iblk by svantbody
c     T. Herring Jan 2020: Converted svantdx into (2,3,NS) to save L1 and L2
c                       values (needed for Galileo). Added antdaz for antenna
c                       azimuth (repr03/itrf2020)

C   HEADER SERIES DATA BLOCK  
C   iflag = 2                                                        I*4
C   sitnam                  : Site name                             C*32 
c   rcvtyp                  : Receiver type                         C*20   
c   rcvnum                  : Receiver serial number                C*20   
c   rcvrsw                  : Receiver software                      C*3    
c   swver                   : Software version number                R*4   
c   anttyp                  : Antenna type                          C*20   
c   antnum                  : Antenna serial number                 C*20  
C   npart                   : Number partials                        I*4 
c   norb                    : Number of orbital partials             I*4      
c   gnss                    : GNSS (G R E C J I)                     C*1
C   nsat                    : Total number of satellites             I*4
C   ischan (nsat)           : PRN in each channel                    I*4
c   fL1 (nsat)              : SV higher frequency each channel       R*8
c   fL2 (nsat)              : SV lower frequency each channel        R*8
C   ndat                    : Number obs types                       I*4
c   dattyp (ndat)           : Data types                             I*4
c     1                           L1: carrier phase (cy)
c     2                           L2: carrier phase (cy)
c     3                           P1: L1 P code pseudo-range (m)
c     4                           P2: L2 P code pseudo-range (m)
c     5                           C1: L1 C/A code pseudo-range (m)
c     6                           D1: L1 doppler frequency (Hz)
c     7                           D2: L2 doppler frequency (Hz)
C   lambda (nsat,ndat)      : wavelength factor                      I*4
c     0                           NO DATA
c     1                           unambiguous, undoubled values
c    -1                             ambiguous, undoubled values
c     2                           unambiguous,   doubled values
c    -2                             ambiguous,   doubled values
c   skd                     : Static (S) Kinematic (K) or Dynamic(D) C*1
C   nepoch                  : Number epochs                          I*4
C   inter                   : Obs interval (sec)                     I*4
c   ircint                  : Original receiver sampling interval(s) I*4
C   mtime                   : Type of time (1=UTC)                   I*4
c   isessn                  : Session number
c   iy                      : Start epoch year                       I*4
c   im                      : Start epoch month                      I*4
c   id                      : Start epoch day                        I*4
c   ihr                     : start epoch hour                       I*4
c   min                     : start epoch minute                     I*4
C   sec                     : start epoch seconds                    R*8
c   offarp                  : ant ref pt offset from mark          3 R*8 
C   offsl1 (3)              : L1 ant offset U,N,E                  3 R*8
C   offsl2 (3)              : L2 ant offset U,N,E                  3 R*8
C MOD TAH 200205: Added antdaz antenna azimuth (from model.h as is offarp)
C   antdaz                   : Antenna azimith (deg)                  R*8
C MOD TAH 200126: Increased to(2,3,nsat) from (3,nsat)   
c   svantdx(2,3,nsat)       : L1 and L2 SV ant offset X, Y, Z        R*8
C   obfiln                  : X-file name                           C*16
C   tfiln                   : T-file name                           C*16
c   jfiln                   : J-file name                           C*16    
c   frame                   : Inertial ref. system (B1950 or J2000)  C*5 
c   precmod                 : Precession model (IAU68 or IAU76)      C*5  
c   nutmod                  : Nutation model (IAU80)                 C*5   
c   gravmod                 : Gravity model for ARC                  C*5  
c   srpmod                  : Radiation-pressure model               C*5  
c   eradmod                 : Earth-radiation model                  C*5
c   antradmod               : Antenna-radiation model                C*5
c   isptide                 : Flag for 12h/24h UT1 pole (binary)     I*4  
c      1 = UT1
c      2 = pole                                                            
c   ietide                  : Flag for tides applied (binary coded)  I*4  
c      1 = solid earth tides
c      2 = K1 frequency dependant earth tide
c      4 = Pole tide
c      8 = Ocean Tide      
c     16 = Atm tide  
c   speopmod                : Short-period EOP model (IERS92 IERS96  C*8 
c   etidemod                : E tide model (IERS92 IERS03)           C*8
c   otidemod                : Ocean tide model (OSO   NAO  )         C*8   
c   atmtide                 : Atmospheric tides (ECMWF)              C*8
c   atmlmod                 : Atmospheric loading model & frame      C*8 
c                              e.g., ECMWF CM
c   hydrolmod               : Hydrological loading model             C*8
c   atmlavg (3)             : Average N E U for atm loading          R*8
c   hydrolavg (3)           : Average N E U for hydrolog. loading    R*8    
c   dryzen                  : Source of P, T for dry zenith delay    C*4
c   wetzen                  : Source of wet zenith delay             C*4 
c   drymap                  : Hydostatic mapping function            C*4
c   wetmap                  : Non-hydrostatis mapping function       C*4  
c   ionsrc                  : Source of ionospheric model            C*4
c   magfield                : Magnetic field for ionosoheric model   C*6 
c   antmod_snx              : SINEX name for antenna PCO/PCV source  C*10
c   antmod                  : Antenna model used (NONE,ELEV,AZEL)    C*4  
c   svantbody(nsat)         : Antenna/body-type                      C*20
c   svantmod_snx(nsat)      : SINEX name for SV antenna PCO/PCV      C*10
c   svantmod (nsat)         : SV antenna model used (NONE,ELEV,AZEL) C*4
c   elvcut                  : Elevation cutoff angle (degrees        R*8    
c   nclock                  : Number of terms in clock polynomial    I*4  
c   clock (nclock)          : Coefficients of clock polynomial       R*8  
c   jdet                    : (PEP) Julian day of t-file IC epoch    I*4
c   tet                     : seconds of day of t-file IC epoch      R*8
c   jdr                     : (PEP) Julian day of earth rotation     I*4
c   tr                      : seconds of day of earth rotationh      R*8
c   UT1                     : UT1 at start (s)                       R*8
c   UT1dot                  : UT1 rate at start (s/day)              R*8
c   xp                      : x pole position at start (arcsec)      R*8
c   xpdot                   : x pole rate at start (arcsec/day)      R*8
c   yp                      : y pole position at start (arcsec)      R*8
c   ypdot                   : y pole rate at start (arcsec/day)      R*8
c   psi                     : nutation in longitude at start (rad)   R*8
c   psidot                  : nutation long  rate at start (rad/s)   R*8
c   eps                     : nutation in obliquity at start (rad)   R*8
c   epsdot                  : nutation obliq rate at start (rad/s)   R*8
c   avlmet                  : availability of met data               I*4
c      Binary coded:              1  pressure
c                                 2  temperature
c                                 4  relative humidity
c                                 8  WVR delay
c   nslip                   : number of extra bias parameter flags   I*4
c   islip  (nslip)          : epoch number of each extra bias flag   I*2
c   islpst (nslip)          : satellite for each extra bias flag     I*2
c   nextra                  : Number extra R*8s (head.rec)           I*4
c   extra (nextra)          : Extra r*8s                             R*8
c   niextra                 : Number extra I*4s (head.rec)           I*4
c   iextra (nextra)         : Extra i*4s                             I*4
c   nrextra                 : Number extra R*8s (head.rec)           I*4
c   rextra (nextra)         : Extra r*8s                             R*8
c   ncextra                 : Number extra c*8s (head.rec)           I*4
c   cextra (nextra)         : Extra r*8s                             C*8



        
      implicit none

      include '../includes/dimpar.h'

c     output variables:
      integer*4          nsat,norb,npart,ndat,dattyp(maxdat)
     .,                  nepoch,inter,mtime,isprn(maxsat)
     .,                  iy,im,id,ihr,min,lambda(maxsat,maxdat)
     .,                  avlmet,nslip,ircint,isessn,ietide,isptide
     .,                  nclock,niextra,iextra(maxext),nrextra,ncextra
      integer*2          islip(maxcsb),islpst(maxcsb)
      integer*4          jdet,jdr
      real*8             sec,offarp(3),offsl1(3),offsl2(3),tet,tr
     .,                  ut1,ut1dot,xp,xpdot,yp,ypdot
     .,                  psi,psidot,eps,epsdot
     .,                  rextra(maxext),elvcut,clock(4)  
C MOD TAH 200126: Changed to svantdx(2,3,maxsat) (added L1/L2 index)
     .,                  atmlavg(3),hydrolavg(3),svantdx(2,3,maxsat)
     .,                  fL1(maxsat),fL2(maxsat)
     .,                  real8
      real*8 antdaz       ! Antenna azimith (deg) (model.h; TAH 200205)
      character*1        skd,gnss
      character*3        rcvrsw
      character*4        antmod,svantmod(maxsat),dryzen,wetzen
     .,                  drymap,wetmap,ionsrc
      character*5        frame,precmod,nutmod,gravmod,srpmod            
     .,                  eradmod,antradmod   
      character*6        magfield
      character*8        etidemod,atmlmod,otidemod,atmtide,hydrolmod
     .,                  speopmod,cextra(maxext)
      character*10       antmod_snx,svantmod_snx(maxsat)
      character*16       tfiln,obfiln,jfiln
      character*20       rctype,rcvnum,anttyp,antnum,svantbody(maxsat)
      character*32       sitnam
      real*4             swver

      integer            iflag,lunit,i,j
      integer*4          ioerr,inqerr,len,rcpar
      character*4        buf4
      character*16       cfname
      character*80       prog_name

      logical            debug
      data debug /.false./

c     get the calling module name for report_stat
      len = rcpar(0,prog_name)

      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag
     .,     sitnam,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,     npart,norb,gnss,nsat,(isprn(i),i=1,nsat)
     .,     (fL1(i),fL2(i),i=1,nsat)
     .,     ndat,(dattyp(i),i=1,ndat)
     .,     ((lambda(i,j),i=1,nsat),j=1,ndat)
     .,     skd,nepoch,inter,ircint,mtime,isessn
     .,     iy,im,id,ihr,min,sec
     .,     (offarp(i),i=1,3),(offsl1(i),i=1,3),(offsl2(i),i=1,3)
C MOD TAH 200205: Add antz
     .,     antdaz
C MOD TAH 200126: Changed size and used F95 bounds
     .,     (svantdx(:,:,j),j=1,nsat),obfiln,tfiln,jfiln
     .,     frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
     .,     ietide,isptide 
     .,     speopmod,etidemod,otidemod,atmtide,atmlmod,hydrolmod
     .,     (atmlavg(i),i=1,3),(hydrolavg(i),i=1,3)
     .,     dryzen,wetzen,drymap,wetmap,ionsrc,magfield
     .,     antmod_snx,antmod
     .,     (svantbody(i),i=1,nsat),(svantmod_snx(i),i=1,nsat)
     .,     (svantmod(i),i=1,nsat) 
c     .,      elvcut,nclock
     .,     elvcut,nclock,(clock(i),i=1,nclock)
     .,     jdet,tet,jdr,tr
     .,     ut1,ut1dot,xp,xpdot,yp,ypdot,psi,psidot,eps,epsdot
     .,     avlmet
     .,     nslip,(islip(i),i=1,nslip),(islpst(i),i=1,nslip)   
     .,     niextra,(iextra(i),i=1,niextra)
     .,     nrextra,(rextra(i),i=1,nrextra)
     .,     ncextra,(cextra(i),i=1,ncextra)

 1000 if (ioerr .ne. 0) then
         inquire ( unit=lunit,name=cfname,iostat=inqerr )
         call report_stat('FATAL',prog_name,'lib/readc2',cfname
     .                   ,'Error reading C-file',ioerr)
      endif

      if (iflag .ne. 2) then
         write(buf4,'(i4)') iflag
         call report_stat('FATAL',prog_name,'lib/readc2',buf4
     .                   ,'Wrong iflag: ',0)
      endif
      
      if ( debug ) then
       print *,'READC2: sitnam= ',sitnam
       print *,'READC2: rctype,rcvnum= ',rctype,rcvnum
       print *,'READC2: rcvrsw,swver= ',rcvrsw,swver
       print *,'READC2: anttyp,antnum= ',anttyp,antnum
       print *,'READC2: npart = ',npart 
       print *,'READC2: norb  = ',norb
       print *,'READC2: gnss = ',gnss 
       print *,'READC2: nsat  = ',nsat
       print *,'READC2: isprn = ',(isprn(i),i=1,nsat)
       print *,'READC2: ndat  = ',ndat
       print *,'READC2: dattyp= ',(dattyp(i),i=1,ndat)
       do  i=1,nsat
          print *,'READC2: lambda ',(lambda(i,j),j=1,ndat)
       enddo
       print *,'READC2: skd   = ',skd
       print *,'READC2: nepoch= ',nepoch
       print *,'READC2: inter = ',inter
       print *,'READC2: ircint= ',ircint
       print *,'READC2: mtime = ',mtime
       print *,'READC2: isessn = ',isessn
       print *,'READC2: iy,im,id,ihr,min,sec ',iy,im,id,ihr,min,sec
       print *,'READC2: offarp= ',(offarp(i),i=1,3)
       print *,'READC2: offsl1= ',(offsl1(i),i=1,3)
       print *,'READC2: offsl2= ',(offsl2(i),i=1,3)
       print *,'READC2: obfiln= ',obfiln
       print *,'READC2: tfiln = ',tfiln
       print *,'READC2: jfiln = ',jfiln
       print *,'READC2: frame,precmod,nutmod= ',frame,precmod,nutmod
       print *,'READC2: gravmod,srpmod,eradmod,antradmod= '
     .                  ,gravmod,srpmod,eradmod,antradmod
       print *,'READC2: ietide,isptide= ',ietide,isptide  
       print *,'READC2: speopmod,etidemod,otidemod,atmtide= '
     .                 ,speopmod,etidemod,otidemod,atmtide
       print *,'READC2: atmlmod,hydrolmod= ',atmlmod,hydrolmod
       print *,'READC2: dryzen,wetzen,drymap,wetmap= '
     .                 ,dryzen,wetzen,drymap,wetmap  
       print *,'READC2: ionsrc,magfield= ',ionsrc,magfield
       print *,'READC2: antmod_snx,antmod,elvcut= '
     .                 ,antmod_snx,antmod,elvcut
       do i=1,nsat
         print *,'READC2: svantbody,svantmod_snx, svantmod= '
     .            ,svantbody(i),svantmod_snx(i),svantmod(i) 
       enddo        
       print *,'READC2: nclock,clock= ',nclock,(clock(i),i=1,nclock)
       print *,'READC2: jdet,tet = ',jdet,tet
       print *,'READC2: jdr,tr = ',jdr,tr
       print *,'READC2: ut1,ut1dot= ',ut1,ut1dot
       print *,'READC2: xp,xpdot= ',xp,xpdot
       print *,'READC2: yp,ypdot= ',yp,ypdot
       print *,'READC2: psi,psidot,eps,epsdot= ',psi,psidot,eps,epsdot
       print *,'READC2: avlmet= ',avlmet
       print *,'READC2: nslip= ',nslip
       print *,'READC2: islip,islpst= ',(islip(i),islpst(i),i=1,nslip)
       print *,'READC2: niextra= ',niextra
       print *,'READC2: iextra = ',(iextra(i),i=1,niextra)     
       print *,'READC2: nrextra= ',nrextra
       print *,'READC2: rextra = ',(rextra(i),i=1,nrextra)
       print *,'READC2: ncextra= ',ncextra
       print *,'READC2: cextra = ',(cextra(i),i=1,ncextra)
      endif
      return
      end
