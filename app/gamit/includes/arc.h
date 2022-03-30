c     arc.h  Created 101223/101230; last modified by T. Herrinh 190702

c File unit numbers - arc, arcmrg,eclout, egm96, ephdrd, ertorb,filopn, init, redsat,
c         read_otides,sbfn, shadow, shwprt, wrthed
c     Global unit number variables now in units.h 
      common/arcunits/iarh,iyaw,iyawtmp,iyawtab,isunit,iotide
      integer*4 iarh,iyaw,iyawtmp,iyawtab,isunit,iotide
     
c Binary-coded flag and logical unit for debug and hidden flags - arc, ertorb, shadow, sbfn
      common/debug/idbprt,idebug
      integer*4 idbprt,idebug    

c File names - arc, filopn, read_arc_batch, check_models, wrthed 
      common/filenames/ gfname,ytabfname,afname,tfname,yfname,dbfname
     .                , otidfname
      character*16   gfname,ytabfname,afname,tfname,yfname,dbfname
     .                , otidfname

c Conversion constants - argmnt, dtwopi, hofrad, keplr, prces, redsat, sbfn, shadow, timred
      common/stwopi/twopi,cdr,casr
      real*8 twopi,cdr,casr

c Physical constants - init, keplr,,sbfn, sbfn1, shadow, merit, igs92, wgs84, egm96, egm08 
      common/const/gm(3),ltvel,aunit,sunrad,ertrad,monrad
      real*8 gm,ltvel,aunit,sunrad,ertrad,monrad

      character*16 satnam(maxsat),sname 
c MOD TAH 190922: Added satsvn (unique satellite number).  Need for force
c     modeling of some satellites.  Also added ssvn which is set when get_sat_info
c     is called for satellite isat. 
      integer*4 satsvn(maxsat)   ! Unique number returned as jsat from svnav_read.f
                                 ! or read_svsinex.f.       
      integer*4 iprn,nsats, ssvn

c Satellite names and number - arc, eclout, redsat, shwprt, wrthed
      common/satlnm/satnam,sname,iprn,nsats, ssvn, satsvn
             
c Satellite initial conditions and integration times - arc, adam, init, redsat, start_int, wrthed 
      common/stint/satics(maxyt2),trun0,trunf,tint0,tstop
      real*8 satics,trun0,trunf,tint0,tstop

c Maximum yaw rate and sun sensor bias information - redsat, eclout
      common/yawinfo/yawrate,yaw_entries,bias
      real*8 yawrate
      integer*4 yaw_entries
      character*1 bias
           
c Input time type and 20-character frame name - arc, wrthed
       common/timfrm/time_type,frame_name
       character*4 time_type
       character*20 frame_name

c----------------------------------------------------------------------------------------

c Epochs for integration - arc, eclout, ephdrd, redsat, sbfn, shwprt, wrthed
      common/timxtr/te,tb,tf,delt,tdtoff,jde,jdb,jdf
      real*8 te,tb,tf,delt,tdtoff
      integer*4 jde,jdb,jdf
              
c Integration interval controls - adam,eval,init,start_int,revrec,sbfn,sbout,wrthed
      common/incon/epsi,hc,hmx,hmn,iptr3,l1,diint,stdint,kount,apar     
      real*8 epsi,hc,hmx,hmn,diint,stdint
      integer*4 iptr3,l1,kount        
      character*1 apar

c Flag for forward or backward integration -  adam,init,start_int,sbout
      common/nmstrt/nsign    
      integer*4 nsign

c Integration controls - adam,calcof,init,start_int
      common/adams/dydt(maxyt2,13),tab,dnh,dinthmx,npredt,ncoret,neq,meq    
      real*8 dydt,tab,dnh,dinthmx
      integer*4 npredt,ncoret,neq,meq

c Coefficients of the Adams-Moulton predictor-corrector - adam,calcof
      common/coeff/binom(12,12),cadams(13),dadams(13)     
      real*8 binom,cadams,dadams

c Coefficients for Nordsieck variable-step integerator - eval,start_int
      common/nordcf/anord(maxyt2),bnord(maxyt2),cnord(maxyt2)
     .  ,dnord(maxyt2),enord(maxyt2),l4,l5,n1,n2,fnord(maxyt2,3)
      real*8 anord,bnord,cnord,dnord,enord,fnord
      integer*4 l4,l5,n1,n2

c Integrated coordinates for each of 4 iterations - adam,eval,init,start_int,sbfn,sbout (see also revrec)
      common/output/y(maxyt2,4),k2,nstout,npstep,nrec 
      real*8 y 
      integer*4 k2,nstout,npstep,nrec

c Partial derivative quantities - sbfn,sbfn1
c**       common/paraux/dsbcor(3,15),dadx(3,3)
c**       common/paraux/dsbcor(3,maxorb),dadx(3,3)
      common/paraux/dsbcor(3,maxorb-3),dadx(3,3)
     
      real*8 dsbcor,dadx

c---------------------------------------------------------------------------------------

c Requested degree (and order) for the gravity field, earth tides, and ocean tides
      common/gravdegrees/gravdeg,etidedeg,otidedeg 
      integer*4 gravdeg,etidedeg,otidedeg 

c Stokes harmonic coeffcients for the static gravity field and solid-Earth tides
c     init,igs92,merit,sbfn,sbfn1,sclcof,wgs72,wgs84, egm96, egm08 (increased to degree 12)
      common/harcof/crad,gmhfct,k2mr(3),k2mi(3),k2mp(3),k3m(4)
     .            , czhar(11),cchar(77),cshar(77)
     .            , cztid(11),cctid(77),cstid(77)
     .            , nczone,nctess,nczon1,nctes1,nctes2,zero_tide
      real*8 crad,gmhfct,k2mr,k2mi,k2mp,k3m
     .     , czhar,cchar,cshar,cztid,cctid,cstid
      integer*4 nczone,nctess,nczon1,nctes1,nctes2
      logical zero_tide 
                                          
c Stokes harmonic coefficients for frequency-dependent 2nd-degree zonal 
c and (diurnal) solid-Earth tides
      common/ef2tides/ztidip(6),ztidop(6),dtidip(11),dtidop(11)
     .               , ztid_doodson(6),dtid_doodson(11)
      real*8 ztidip,ztidop,dtidip,dtidop
      character*7 ztid_doodson,dtid_doodson
                                                                                  
c Stokes harmonic coefficients and Doodson numbers for ocean tides 
c     (same degree and order as fixed harmonics but note the addition 
c      of sine zonals 
c     read_otides,sbfn,sbfn1 
      common/otides/ozcp(18,11),ozsp(18,11)
     .            , otcp(18,77),otsp(18,77),otcm(18,77),otsm(18,77)
     .            , czcotid(11),czsotid(11),ctcotid(77),ctsotid(77)
     .            , otid_doodson(18) 
      real*8 ozcp,ozsp,otcp,otsp,otcm,otsm
     .     , czcotid,czsotid,ctcotid,ctsotid
      character*7 otid_doodson

c Lat, long, and Legendre polynomials for gravitational accelerations 
c     sbfn,sbfn1: Increased to 12x12 field 
      common/hrmaux/cntrot(3,3),cntrtd(3,3),cslat,cclat,cslng(12)
     . , cclng(12),cslat1(3),cclat1(3),clng1(3)
     . , cleg(12),cleg1(12),cgleg(77),cgleg1(77),cleg2(12),cgleg2(77)
      real*8 cntrot,cntrtd,cslat,cclat,cslng,cclng,cslat1,cclat1
     . , clng1,cleg,cleg1,cgleg,cgleg1,cleg2,cgleg2

c-----------------------------------------------------------------------------------------
                     
c Planetary ephermeris information
      common/nbody/lbody 
      integer*4 lbody 

c***rwk 180301: Most of commons evcoef, soltab, and luntab  will go away with the new ephemeris code

c Everett interpolation coefficients for luni-solar coordinates - evrtcf,lunred,solred
      common/evcoef/evcf(5,5),fact(9) 
      real*8 evcf,fact

c Earth ephemeris interpolation - ephdrd,solred
      common/soltab/fjdbsn,sdelts,ytrps(5,2,6),yys(10,3)
     .    , nintrss,iendfs,ji0s,jils,iy1s,iy2s,jlasts,jnows 
      real*8 fjdbsn,sdelts,ytrps,yys
      integer*4 nintrss,iendfs,ji0s,jils,iy1s,iy2s,jlasts,jnows
     
c Moon ephereris interpolation - ephdrd,lunred
      common/luntab/fjdbmn,sdeltm,ytrpm(5,2,6),yym(10,3)
     .     , nintrsm,iendfm,ji0m,jilm,iy1m,iy2m,jlastm,jnowm 
      real*8 fjdbmn,sdeltm,ytrpm,yym
      integer*4 nintrsm,iendfm,ji0m,jilm,iy1m,iy2m,jlastm,jnowm

c**** end rwk 180301 

c Sun, Moon, Venus, Jupiter and satellite coordinates - ertorb,sbfn,sbfn1,shadow 
c         (rsbh dimensioned for a 12-degree gravity field)
      common/coraux/sbcor(6),bcor(3),ccor(6),ccor3(3),pccor(6)
     .    , pccor3(3),pbcor(3),vccor(3),vbcor(3),jccor(3),jbcor(3)
     .    , rsb,rsb2,rsb3,rb,rb2,rb3,rc,rc3,rpc,rpc3,rpb2,rpb3
     .    , rvc,rvb,rjc,rjb,rsbh(12)
      real*8 sbcor,bcor,ccor,ccor3,pccor,pccor3,pbcor,vccor,vbcor
     .    , jccor,jbcor,rsb,rsb2,rsb3,rb,rb2,rb3,rc,rc3,rpc,rpc3
     .    , rpb2,rpb3,rvc,rvb,rjc,rjb,rsbh

c Polar motion and UT1 -  sbfn, egm08, wrthed
      common/polmot/xpole,ypole,xpm,ypm,ut1utc 
      real*8 xpole,ypole,xpm,ypm,ut1utc

c---------------------------------------------------------------------------

c Eclipse data - shwprt,eclout 
      common/eclipse/eclipse_start(100),eclipse_end(100)
     .     ,eclipse_beta(100),pjdlast,neclipse,first,eclipse_type(100)
      logical first
      integer neclipse
      real*8 pjdlast,eclipse_start,eclipse_end,eclipse_beta
      character*1 eclipse_type

c Direct and reflected solar radiation quantities - ertorb, ghdred, init, redsat, sbfn, wrthed
c     plus /orbits 
* MOD TAH 190702: Added antpwr to common and added getting this value
*     for the IGS_metadata.snx file (linked to svnav.dat)
      common/radiation/sbmass,antpwr, radprs(3),radcon(13),
     .                snvec(3),yvec(3),
     .                nics,lamprt,modrad,icsnam(maxorb),antbody
      real*8 antpwr  ! Transmission power (W) can be read from IGS_metadata.snx 
      real*8 sbmass,radprs,radcon,snvec,yvec
      integer*4 modrad,nics
      character*4 icsnam
      character*20 antbody
      logical lamprt
c Values for Univeersity College London E-radiation model - filopn,srpfgrd
      common/uclgrd/
     .  busX3xmin,busX3ymin,busX3cellsize,busX3grid(100,120)
     .  ,busY3xmin,busY3ymin,busY3cellsize,busY3grid(100,120)
     .  ,busZ3xmin,busZ3ymin,busZ3cellsize,busZ3grid(100,120)
     .  ,busX4xmin,busX4ymin,busX4cellsize,busX4grid(100,120)
     .  ,busY4xmin,busY4ymin,busY4cellsize,busY4grid(100,120)
     .  ,busZ4xmin,busZ4ymin,busZ4cellsize,busZ4grid(100,120)
     .  ,pnlX3xmin,pnlX3ymin,pnlX3cellsize,pnlX3grid(100,120)
     .  ,pnlY3xmin,pnlY3ymin,pnlY3cellsize,pnlY3grid(100,120)
     .  ,pnlZ3xmin,pnlZ3ymin,pnlZ3cellsize,pnlZ3grid(100,120)
     .  ,pnlX4xmin,pnlX4ymin,pnlX4cellsize,pnlX4grid(100,120)
     .  ,pnlY4xmin,pnlY4ymin,pnlY4cellsize,pnlY4grid(100,120)
     .  ,pnlZ4xmin,pnlZ4ymin,pnlZ4cellsize,pnlZ4grid(100,120)
     .  ,bus3nrows,bus3ncols,pnl3nrows,pnl3ncols
     .  ,bus4nrows,bus4ncols,pnl4nrows,pnl4ncols
C     bus for IIA satellite
      real*8 busX3xmin,busX3ymin,busX3cellsize,busX3grid
     .  ,busY3xmin,busY3ymin,busY3cellsize,busY3grid
     .  ,busZ3xmin,busZ3ymin,busZ3cellsize,busZ3grid
C     bus for IIR-A satellite
      real*8 busX4xmin,busX4ymin,busX4cellsize,busX4grid
     .  ,busY4xmin,busY4ymin,busY4cellsize,busY4grid
     .  ,busZ4xmin,busZ4ymin,busZ4cellsize,busZ4grid
C     solar panel for IIA satellite
      real*8 pnlX3xmin,pnlX3ymin,pnlX3cellsize,pnlX3grid
     .  ,pnlY3xmin,pnlY3ymin,pnlY3cellsize,pnlY3grid
     .  ,pnlZ3xmin,pnlZ3ymin,pnlZ3cellsize,pnlZ3grid
C     solar panel for IIR-A satellite
      real*8 pnlX4xmin,pnlX4ymin,pnlX4cellsize,pnlX4grid
     .  ,pnlY4xmin,pnlY4ymin,pnlY4cellsize,pnlY4grid
     .  ,pnlZ4xmin,pnlZ4ymin,pnlZ4cellsize,pnlZ4grid
C     Number of rows and columns in UCL SRP grid files 
C     assume that the x,y,z for one object will have 
C     the same row column structure.
      real*8 bus3nrows,bus3ncols,pnl3nrows,pnl3ncols
     .  , bus4nrows,bus4ncols,pnl4nrows,pnl4ncols


