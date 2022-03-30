c    Common for MODEL.  Last changed by R. King  180320 
              
c Text for C-file header
      integer*4 ntext
      character*80 text(maxtxt)
      common /cftext/ntext,text
      

c Unit numbers and file names: see model/open.f for definitions and assignments
c Global unit number variables now in units.h 
                 
      integer*4  isesfo,istnfo,ihi,ipcv,iuu,iudcb,iuc,iud
     .          ,iuf,iui,iuj,iobs,iuw,iuy,iuz,iueqrn,iepcv
      common /lunits/ isesfo,istnfo,ihi,ipcv,iuu,iudcb,iuc,iud
     .               ,iuf,iui,iuj,iobs,iuw,iuy,iuz,iueqrn,iepcv      
      character*16 pfiln,ifiln,jfiln,lfiln,obfiln,cfiln,tfiln
     .           , ionfiln,metfiln,zfiln,epcvfiln,yfiln
      common /filenames/ pfiln,ifiln,jfiln,lfiln,obfiln,cfiln,tfiln
     .                 , ionfiln,metfiln,zfiln,epcvfiln,yfiln
                  
c Satellite information from the T-file 

c      jdbt, tbt - start PEP JD, sod
c      jdft, tft - finish  PEP JD, sod
c      jdet, tet - JD, sod epoch of initial conditions
c      sdelt     - tabular point interval on t-file (sec0
c      saticst    - initial conditions
c      icsnamt    - names for initial conditions
      integer*4 ntsat,itsat(maxsat),norbprm,nintrs,nepcht
     .        , jdbt,jdft,jdet                                       
      real*8 sdelt,tbt,tft,tet,saticst(maxorb,maxsat)
      character*4 icsnamt(maxorb)
      common /tfcom/tbt,tft,tet,sdelt,saticst,jdbt,jdft,jdet
     .       , ntsat,itsat,norbprm,nintrs,nepcht,icsnamt
                       
                        
c Quantities describing the observations
                                          
c       jd0, t0 - PEP JD, sod start of observations
c       jdend, tend - PEP JD  end of observations
c       
c       inter   - sampling interval on x-file (sec)
c       ircint  - original receiver sampling interval (sec)
      integer*4 jd0,jdend,isessn,mtime,inter,ircint,nepoch,nepoch_use(2)
     .        , ndat,dattyp(maxdat),lambda(maxsat,maxdat)
     .        , nchan,ischan(maxsat),norbpart,npart 
      real*8 t0,tend,fL1(maxsat),fL2(maxsat),elvcut
* MOD TAH 200512:Increased atxfrq from 2 to 3 to allow as secondary
*     low frequeny choice.  Also increased to 4 characters so that
*     4th character can show which one found (* added to name) 
      character*4 atxfrq(3)

      logical fixash,sv_missing(maxsat)
      common /obscom/t0,tend,fL1,fL2,elvcut,mtime,nchan,ischan
     .     , jd0,jdend,isessn,inter,ircint,nepoch,nepoch_use,ndat,dattyp
     .     , lambda,norbpart,npart,fixash,sv_missing,atxfrq
                        
c Qauntities for simulation mode
                  
      real*8 noise(maxdat),simvec0(3,2),dispneu(3)
      logical simulation
      common /simcom/simvec0,noise,dispneu,simulation
      

c L-file coordinates      
c   sitecd  : 4-character site code (used by GAMIT)
c   asite   : 8-character site code (used to document lread return)
c   kepoch0 : apr file epoch ( decimal yrs) 
c   kpos    : Cartesian position from the L-file (m)
c   kvel    : Cartesian velocity from the (apr-style) L-file (m/yr) 
c   kfflg   : flag for type of L-file (spherical = 0, Cartesian = 1)
c   kvflg   : flag for velocities to be used (no=0, yes=1)          
c   jdmidpt : JD of observation midpoint for coordintes
c   tmidpt  : seconds-of-day for observation midpoint for coordinates  
c   

      real*8 kepoch0,kpos,kvel,tmidpt
      integer*4  kfflg, kvflg, jdmidpt
      character*4 sitecd
      character*8 asite

      common /lfcom/ kepoch0, kpos(3), kvel(3), tmidpt
     .             , kfflg, kvflg, jdmidpt
     .             , asite, sitecd
          

c Computed site quantities

c     evec0(3,2)
c     sitepart0(3,3)
c     latr - geodetic latitude in radians
c     lonr - longitude in radians        
c     height - geodetic height in km
c     latr_sph - spherical latitude for q-file adjustments (radians)
c     radius -   spherical radius for q-file adjustments (m)
c     semi - semi-major axis of the ellipsoid
c     finv - inverse flattening of the ellipsoid
c     shft(3) - dx, dy, dz transatlation of the ellipsoid
      real*8 evec0(3,2),sitepart0(3,3),latr,lonr,height
     .       ,semi,finv,shft(3),latr_sph,radius 
      common /sitcom/ evec0,sitepart0,height,latr,lonr,semi,finv,shft
     .       , latr_sph,radius 
         
c Clock model and coefficients
                  
c     klock 
c     clkepc
c     clkrat
c     clkacc
c     clkcub
      integer*4 klock
      real*8 clkepc,clkrat,clkacc,clkcub
      common /clkcom/ clkepc,clkrat,clkacc,clkcub,klock              
                       

c Receiver and antenna information and the current times of validity from 
c station.info, dcb.dat, and antmod.dat 

c     kstarts(5) kstops(5) : yr doy hr min sec of start of current station.info entry
c     kstartr(5) kstopr(5) : yr doy hr min sec of start/stop of current eq/rename   
c     rcvcod               : 6-character GAMIT receiver code    
c     rcvrsn               : 20-character receiver serial number    
c     rcvrswx              : 3-character receiver software code (from makex)
c     swver                : decimal GAMIT code for firmware version
c     rcvers               : 20-character actual firmware version
c     pcncod               : character indicating treatment of diffential code biases ('P' 'C' or 'N') 
c     numdcb               : number of SVs for which we havd DCBs
c     prndcb(maxsat)       : PRNs for DCBs  
c     dcb(maxsat)              : DCBs for each SV 
c     anttyp               : 20-character antenna name including radome 
c     antsn                : 20-character antenna serial number 
c     radome_in            : 5-character radome code for measurements
c     antmod_in            : antenna PCV model requested (AZEL ELEV or NONE)
c     antmod               : antenna model available and used (AZEL ELEV or NONE)
c     antmod_snx           : 10-character SINEX code for antenna model 
c     epcvflg              : Y/N for site-dependent empirical PCV model 
c     offarp(3)            : offset of antenna ARP from monument (U N E) (m)
c     offl1(3) offl2(3)    : mean offset of antenna phase center from monument (U N E) (m) 
c     antdaz               : Antenna deviation from true North (Alignment from True N in log)
c                            Added TAH 200203 for repro3  
c     pcvminelev           : minimum elevation in PCV model from ANTEX file (deg)
c     elvtabl1(maxel),elvtabl2(maxel) : Azimuth-averaged PCVs (90-0 deg) (mm)
c     tablel1(maxel,maxaz),tablel2(maxel,maxaz): PCVs by elevation (90-0 deg) and 
c                                                azimuth (0-360) (mm)
         
      integer*4 kstarts,kstops,kstartr,kstopr,numdcb,prndcb
     .        , nel,naz
      character*1 pcncod,epcvflg                                            
      character*3 rcvrswx 
      character*4 antmod_in,antmod
      character*6 rcvcod,antcod
      character*5 radome_in                  
      character*10 antmod_snx  
      character*20 rctype,rcvers,rcvrsn,anttyp,antsn    
      character*32 sitnam    
c     (lfile has a 12-character site name, station.info and x-file header 
c      16-characters but c-file has 32)
      real*4 swver
      real*8 dcb,offarp,offl1,offl2,pcvminelev,elvtabl1,elvtabl2
     .     , tablel1,tablel2,zen1,zen2,dzen,dazi
      real*8 antdaz  ! Antenna aligment from True N (deg).
      logical newant

      common /rcvantcom/ offarp(3),offl1(3),offl2(3),antdaz
     .                 , pcvminelev, dcb(maxsat)
     .                 , elvtabl1(maxel),elvtabl2(maxel)
     .                 , tablel1(maxel,maxaz),tablel2(maxel,maxaz)
     .                 , zen1,zen2,dzen,dazi,nel,naz
     .                 , kstarts(5),kstops(5),kstartr(5),kstopr(5)  
     .                 , numdcb,prndcb(maxsat),swver,newant,sitnam
     ,                 , rcvcod, pcncod,antcod,anttyp,radome_in
     .                 , antmod_in,antmod,antmod_snx,rctype,rcvers
     .                 , rcvrsn,rcvrswx,antsn,epcvflg 

                 
c Satellite antenna information from svnav.dat and antmod.dat

c     svantbody            : c*20 type of satellite                                      
c     svantmod_in          : antenna model requested (AZEL ELEV or NONE)  
c     svantmod(maxsat)     : antenna model available (AZEL ELEV or NONE)
c     svantmod_snx(maxsat) : 10-character SINEX code for antenna model   
c     svantdx(3,2,maxsat)  : antenna phase center offsets L1 L2 (xyz) (m) 
                                                                      
      character*4 svantmod_in,svantmod(maxsat)
      character*10 svantmod_snx(maxsat)
      character*20 svantbody(maxsat)
      real*8 svantdx(3,2,maxsat),svelvtabl1(maxel,maxsat)
     .     , svelvtabl2(maxel,maxsat),svtabl1(maxel,maxaz,maxsat)
     .     , svtabl2(maxel,maxaz,maxsat),svzen1(maxsat),svzen2(maxsat)
     .     , svdzen(maxsat),svdazi(maxsat)
      integer*4 svnel(maxsat),svnaz(maxsat)
      common/svantcom/ svantdx,svelvtabl1,svelvtabl2,svtabl1,svtabl2
     .               , svzen1,svzen2,svdzen,svdazi,svnel,svnaz
     .               , svantbody,svantmod_in,svantmod,svantmod_snx
                 
c This common stores the model names and controls (see also global.h)

c   ietide    i*4   binary-coded for application of loading tides
* MOD TAH 2020218: Added documentation ietide
*           Bit                                             Value
*             1 = solid earth tides                             1
*             2 = frequency dependant K1 model                  2
*             3 = Pole tide applied to zero mean pole           4
*             4 = ocean tides                                   8
*             5 = Pole tide applied to IERS2010 mean pole      16
*             6 = atmospheric tides                            32
*             7 = Pole tide applied to IERS2020 mean pole      64

c   isptide   i*4   binary-coded for short-period EOP models
*           Bit                                             Value
*             1 = short period pole corrections                 1
*             2 = short period ut1 corrections                  2
*             3 = Ray model (IERS96)                            4
*             4 = IERS 2010 model (IERS10)                      8

c   etidemod  c*8   solid-Earth tides IERS03 etc from sestbl.     
c   otidemod  c*8   ocean-tides from u-file
c   atmlmod   c*8   non-tidal atmospheric loading from u-file
c   atidemod  c*8   atmospheric tidal loading from u-file
c   atmlflg   c*1   apply non-tidal atmospheric loading, Y or N, from ietide from sestbl
c   hydrlflg  c*1   apply hydrologic loading, Y or N, from ietide from sestbl.  
c   hydrolmod c*8   hydrological model, not yet implemented  
c   speopmod  c*8   short-period EOP IERS10 etc, set by isptide read from sestbl. 
c   ionsrc    c*4   ionsosphere source, NONE or GMAP if from IONEX file, from sestbl.
c   magfield  c*6   magnetic field for iono model, IGRF11, from sesbl. 

      integer*4 ietide,isptide
      character*1 atmlflg,hydrlflg
      character*4 ionsrc
      character*6 magfield
      character*8 etidemod,otidemod,atmlmod,atidemod,hydrolmod

      common/modcom/ietide,isptide,etidemod,otidemod,atmlmod,atidemod
     .             ,hydrolmod,ionsrc,magfield,atmlflg,hydrlflg

c This common stores values read from the u-file for ocean tidal loading,
c atmopheric loading (non-tidal and tidal), hydrological loading (not yet 
c implemented), and c meteorological values and mapping function coefficients 
c for the delay c calculation; met values may also come from a RINEX met file.  
c Variables c are described in subroutine gamit/model/readu.f.   rwk 060721
                 
c     max number of ocean tidal constituents from grid
      integer*4 maxotl
      PARAMETER (maxotl=54)
c     max number of values per session for atmospheric loading 
      integer*4 maxatml 
c  MOD TAH 120503: Increased to allow 366 days of 6-hr values (needed in grdtab)
      PARAMETER (maxatml=366*4) 
c     max number of atmospheric tide constituents 
      integer*4 maxatl
      PARAMETER (maxatl=2)
c     max number of values per session for met values
      integer*4 maxmet 
      PARAMETER (maxmet=3000)
c     max number of values per session for mapping function coefficients 
      integer*4 maxmap 
      PARAMETER (maxmap=20)
c     max number of ionospheric maps in a session and max longitude and latitude values 
      integer*4 maxion,maxilon,maxilat
* MOD TAH 141111: Updated maxion to 25 (from 13) to allow 1-hr files introduced 2014/309.
      PARAMETER (maxion=25,maxilon=73,maxilat=73)
      
      character*2 map_name         
      character*3 otlwaves           
      character*8 metmod,mapmod
      integer*4 notl,natml,natl,nmet,nmap  
     .        , ntatml,ntmet,ntmap
      real*4 otides,atml_time,atml_val,atides,met_time,met_val
     .      ,map_time,map_val
      logical lotl,latml,latl,lmet,lmap

      common /ufcom/  
     .   otides(maxotl,6)
     . , atml_time(maxatml), atml_val(maxatml,3)  
     . , atides(maxatl,6)    
     . , met_time(maxmet), met_val(maxmet,3)
     . , map_time(maxmap), map_val(maxmap,9)
     . , notl,natml,natl,nmet,nmap    
     . , ntatml,ntmet,ntmap
     . , lotl,latml,latl,lmet,lmap 
     . , map_name(9),metmod,mapmod,otlwaves(maxotl)

                     
c This common block stores names and values of the meteorlogical models   

c   pres0, temp0, wetvar0 : Constant sea-level values or pressure,
c                           temperature, and water-vapor (RH or WV pressure)
c                           from GPT, STP, or sittbl., used when RINEX 
c                           or u-file values missing         
c   lapse                 : Lapse rate C/km
c   pres, temp, wetvar    : Actual values used at height of station
c   sensor_ht  : height (m) of met sensor (usually station height)
c   sealoc                : Pressure and tempearture local (L) or sealevel (S)
c   ah, aw                : Hydrostatic and wet mapping function coefficients. 
c   dryzen, wetzen, drymap, wetmap:  4-character name of models
c   metsrc                :  Preferred RNX UFL GPT or STP from input file
c   psrc                  :  Pressure source used (RNX GPT STP or input)
c   zsrc                  :  ZHD source used (VMF1 grid)            
c   tsrc                  :  Temperature source used (RNX GPT STP or input)
c   wsrc                  :  Relative humidity (RNX STP or input) or WV pressure (GPT) source used
c   lsrc                  :  Lapse rate source used 
                                                                              
      character*1 sealoc 
      character*3 metsrc,psrc,zsrc,tsrc,wsrc,lsrc,metdef         
      character*4 dryzen,wetzen,drymap,wetmap 
      character*20 metopts 
      real*8 pres0,temp0,wetvar0,zhd0,zendel0,pres,temp,wetvar,lapse
     .     , sensor_ht,ah,aw 
 
      common /metcom/ pres0,temp0,wetvar0,zhd0
     . , zendel0,pres,temp,wetvar,lapse,sensor_ht,ah,aw
     . , metopts,dryzen,wetzen,drymap,wetmap,metsrc,metdef
     . , psrc,zsrc,tsrc,wsrc,lsrc,sealoc 

c This common block stores TEC values from an IONEX file
                                              
c   ilon1                :  First longitude value of ion maps (deg)
c   ilat1                :  First latitude value of ion maps (deg)  
c   dilat                :  Latitude interval (deg)
c   dilon                :  Longitude interval (deg)
c   nilon                :  Number of longitude values per map
c   nilat                :  Number of latitude values per map 
c   ntion                :  Number of ion maps (times)
c   ion_time(ntion)      :  Times of the maps (day-of-year)
c   ion_val(nilon,nilat,ntion) :  Gridded maps 
       
      integer*4 ntion,nilat,nilon
      real*8 ilon1,ilat1,dilat,dilon
      real*4 ion_time,ion_val
      common /ioncom/ ilon1,ilat1,dilon,dilat,nilon,nilat,ntion
     .   ,ion_time(maxion),ion_val(maxilon,maxilat,maxion)    
          
c Other quantities from/for the c-file headers
    
c   psi, psidot          :  Nutation longitude and rate
c   eps, epsdot          :  Nutation obligquity and rate
c   ut1, ut1dot          :  UT1 and rate
c   xp, xpdot            :  Pole x and rate
c   yp, ypdot            :  Pole y and rate 
c   avlmet               :  Binary coded integer for met value availaility
c   atmlavg              :  Average atmospheric loading value over the day (from /solve)
c   hydrolavg            :  Average hydrologic loading value over the day (from /solve)
c   iuttyp               :    
c   nslip.islip,islpst   :  Number cycle slips from input c-file, epochs, SVs 
c   niextra,iextra       :  Spare integer*4 variables 
c   nextra,extra         :  Spare real*8 variables
c   nrextra,rextra       :  Spare character*8 variables
                             
      integer*2 islip(maxcsb),islpst(maxcsb)
      integer*4 avlmet,iuttyp,nparam,nlabel,nslip
     .        , niextra,nrextra,ncextra,iextra(maxext)
      real*8 psi,psidot,eps,epsdot,ut1,ut1dot,xp,xpdot,yp,ypdot
     .     , preval(maxprm),atmlavg(3),hydrolavg(3),rextra(maxext)
      character*8 cextra(maxext)
      character*20 rlabel(maxlab)
      common/chdrecords/psi,psidot,eps,epsdot,ut1,ut1dot,xp,xpdot
     .     , yp,ypdot,atmlavg,hydrolavg,preval
     .     , rextra,avlmet,niextra,nrextra,ncextra,iuttyp,nparam,nlabel
     .     , iextra,nslip,islip,islpst,rlabel,cextra
                    
     .      
c C-file observation records

      integer*4 iepoch,jdobs,okmet,msat,mprn(maxsat)
     .        , ier(maxsat),data_flag(maxsat),isnr(maxdat,maxsat)
     .        , nsave,nspare
      real*8 tobs,rclock0,rclock,svclock(2,maxsat),zendel,atmdel(maxsat)
     .     , elev(maxsat),elvdot(maxsat),azim(maxsat),azmdot(maxsat)
     .     , nadang(maxsat)
     .     , delay(maxdat,maxsat),drate(maxsat),save(maxsav)
     .     , obs(maxdat,maxsat),omc(maxdat,maxsat),tmpart(maxlab,maxsat)
     .     , spare(maxspr)     
      real*4 ampl1(maxsat),ampl2(maxsat),atmlod(3),hydrlod(3)
      common/obsrecords/tobs,rclock0,rclock,svclock,zendel,atmdel
     .     , elev,elvdot, azim,azmdot,nadang,delay,drate
     .     , save,spare,obs,omc,tmpart,ampl1,ampl2
     .     , atmlod,hydrlod
     .     , iepoch,jdobs,msat,mprn,okmet,ier,data_flag,isnr
     .     , nsave,nspare



                                              


