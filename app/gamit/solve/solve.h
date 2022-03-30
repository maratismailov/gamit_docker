c     solve.h   - last modified by rwk 190508 

c     common blocks and common parameters for SOLVE package

c     <scratch> -- directory and names for scratch files
      character*128 scratch_dir,ftmp27,ftmp28,ftmp29
      common/scratch/ scratch_dir,ftmp27,ftmp28,ftmp29

c     <nmat>   -- normal matrix, later updated by covariance matrix   
      real*8 a(maxnrm)
      common/nmat/a

c     <umat>   -- right hand term, later updated by ???? 
      real*8 b(maxprm)
      common/umat/b

c     <amat>   -- design matrix for one epoch  
      real*8 c(maxdm)  
      common/amat/c

c     <amatt>  -- transpose of design matrix   
      real*8 ct(maxdm)
      common/amatt/ct

c     <ipnta>  -- index of design matrix in normal matrix 
      integer*4 ipntc(maxdm)
      common/ipnta/ipntc

c     <ipntat> -- index of transpose of design matrix
      integer*4 ipntct(maxdm)
      common/ipntat/ipntct

c     <irowa>  -- row index of design matrix in normal matrix  
      integer*4 irowc(maxnd)
      common/irowa/irowc 

c     <irowat> -- row index of transpose of design matrix  
      integer*4 irowct(maxp2)
      common/irowat/irowct 
                                 
c     <wmat1>  -- one-way phase covariance matrix  
      real*8 cphi(maxwm1)
      common/wmat1/cphi

c     <wmat3>  -- one-way ionosphere covariance matrix  
      real*8 cphik(maxwm1)
      common/wmat3/cphik

c     <wmat2>  -- double-difference phase covariance matrix, working array later 
      real*8 dpdt(maxwm2)
      common/wmat2/dpdt

c     <wmat4>  -- double-difference ionosphere covariance matrix   
      real*8 dpdtk(maxwm2)
      common/wmat4/dpdtk    
                
      integer maxcom,nlcom,ncom21
      parameter (maxcom=maxprm-maxbis,nlcom=(maxcom*(maxcom+1))/2)
      parameter (ncom21=maxcom*maxobs)

c     <lcnom1> -- LC mode normal matrix(N11 part)  
      real*8 alc(nlcom)
      common/lcnom1/alc

c     <lcnom2> -- LC mode right hand term 
      real*8 blc(maxprm)
      common/lcnom2/blc

c     <lcnom3> -- LC mode normal matrix(N21 part) 
      real*8 clc(ncom21)
      common/lcnom3/clc

c     <lcnom4> -- LC mode normal matrix(N22 part)
      real*8 an22(maxwm1),bn22(maxwm1)
      common/lcnom4/an22,bn22
                           
c     <chi2) 
      real*8 r2sum 
      common/chi2/r2sum

c     <point1> -- index array of mapping operator  
      integer*4 ipntd(maxdd)
      common/point1/ipntd

c     <point2> -- index array of mapping operator transpose 
      integer*4 irowd(maxsd)
      common/point3/irowd
                                
c     <point3> -- index array of mapping operator rows
      integer*4 ipntdt(maxdd)
      common/point2/ipntdt

c     <point4> -- index array of mapping operator transpose rows
      integer*4 irowdt
      common/point4/irowdt(maxnd)
                              
c     <wrk1>   -- working array (maxobs) 
      real*8 work(maxobs)
      common/wrk1/work

c     <wrks1>  -- working array (maxobs) 
      real*8 work1(maxobs)
      common/wrks1/work1

c     <wrks2>  -- working array (maxobs) 
      real*8 work2(maxobs)
      common/wrks2/work2
                                    
c     <dmat>   -- working array (maxdd)  
      real*8  d(maxdd)
      common/dmat/d

c     <dmatt>  -- working array (maxdd)
      real*8 dt(maxdd)  
      common/dmatt/dt

c     <constr> -- a priori constraints
      real*8  stat_apr2(maxsit,3),sat_apr2(maxsat,maxorb)
     .      , stat_apr(maxsit,3),sat_apr(maxsat,maxorb)  
     .      , clk_apr(maxsit),clk_apr2(maxsit)
     .      , zen_apr(maxsit),zen_apr2(maxsit)
     .      , zen_mar(maxsit),zen_mar2(maxsit)
     .      , zen_tau(maxsit),zen_tau2(maxsit)
     .      , grad_apr(maxsit,2),grad_apr2(maxsit,2)
     .      , grad_mar(maxsit,2),grad_mar2(maxsit,2)
     .      , grad_tau(maxsit,2),grad_tau2(maxsit,2)
     .      , eop_apr(6),eop_apr2(6),bias_apr              
      integer*4 maxorb_cov  
      parameter (maxorb_cov=maxorb*(maxorb+1)/2)
      real*8  covorbx(maxorb_cov,maxsat),covorbx2(maxorb_cov,maxsat)
      character*3 satwt_type
      common/constr/stat_apr,stat_apr2,sat_apr,sat_apr2
     .         ,clk_apr,clk_apr2
     .         ,zen_apr,zen_apr2,zen_mar,zen_mar2,zen_tau,zen_tau2
     .         ,grad_apr,grad_apr2,grad_mar,grad_mar2,grad_tau,grad_tau2
     .         ,eop_apr,eop_apr2,bias_apr,covorbx,covorbx2,satwt_type

c     <parts>  -- one-way partial derivative array
      real*8 tpart
      integer*4 npartc,npartm(maxsit)
      common/parts/tpart(maxlab,maxsit,maxsat),npartc,npartm

c     <gparts>  -- one-way partials for atmospheric gradients  
      real*8 gpart(maxsit,maxsat,2)
      common/gparts/gpart
            
c     <stime>  -- clock parameters 
      real*8 t00(3),tor(3,maxsit) 
      integer*4 inter,it0(3),itor(3,maxsit)
      common/times/t00,tor,inter,it0,itor

c     <pickwt> -- measurement and ionosphere error model parameters   
      real*8 aphi,bphi,akappa,bkappa
      integer*4 iwght
      common/pickwt/aphi,bphi,akappa,bkappa,iwght

c     <obscnt> -- numbers of observations [DBADST DOPT LSQUAR QHEAD1 QHEAD2]
      integer*4 iuse,nones,nd1obs,nd2obs
      common/obscnt/iuse(maxsit,maxsat),nones,nd1obs,nd2obs
         
c     <block8> -- numbers of non-bias parameters    
      integer*4 lpart
      common/block8/lpart

c     <track>  -- site names in menu  
      character*12 sitnam(maxsit)
      common/track/sitnam 
                                   
c     <global> -- global parameters 
      integer*4 ntpart,nepoch,nobs,nsat,nsite,norb,nlive
      common/global/ntpart,nepoch,nobs,nsat,nsite,norb,nlive

c     <atmprms> -- model and break-points for multiple-zenith-delay parameters
      integer*4 nzen,idtzen(maxatm),ngrad
     .        , idtgrad(maxgrad)
      character*3 zenmod,gradmod
      common/atmprms/nzen,idtzen,ngrad,idtgrad,zenmod,gradmod

c     <files> -- file names  
      character*16 cfiln(maxsit),obfiln(maxsit),tfiln,mfiln
     .           , qfiln,ofiln,hfiln,nfiln
      common/files/cfiln,obfiln,tfiln,mfiln,qfiln,ofiln,hfiln,nfiln
                           
c     <bbii>   -- bias indices and scale factors   
      real*8 bscale(maxbis)
      integer*4 mbias,idxb(maxbis),msig,iband,l1bias,nlres(maxdd)
      logical*4 nlscale
      common/bbii/bscale,mbias,idxb,msig,iband,l1bias,nlres,nlscale

c     <isigs>  -- index of live parameters in total menu + working array 
      integer*4 isigma
      common/isigs/isigma(maxprm)
                                  
c     <bicnt>  -- bias counting information  
      integer*4 ibias,ibcnt,ibcnt1,ibcnt2,ipntb1(maxobs),ipntb2(maxobs)
      common/bicnt/ibias,ibcnt,ibcnt1,ibcnt2,ipntb1,ipntb2  
                 
c     <ajunk>  -- working array (maxdbl*maxobs)  
      real*8 dr(maxdbl*maxobs)
      common/ajunk/dr

c     <lsqvar> -- scale factor 
      real*8 sclerr      
      common/lsqvar/sclerr
   
c     <cutoff> -- epoch coverage information
      real*8 elvcut_solve(maxsit)
      integer*4 minsnr,istart,iend,kepoch,idecim
      common/cutoff/elvcut_solve,minsnr,istart,iend,kepoch,idecim
                                              
c     <stwght> -- effective observation number for each site and satellite
      integer*4 iusest(maxsit),iusesa(maxsat),iatcon
      common/stwght/iusest,iusesa,iatcon
    
c     <bcrit> -- criteria for bias-fixing    
      real*8 wldev,wlsig,wlcut,wlrat,wldmax
     .     , nldev,nlsig,nlcut,nlrat,nldmax 
     .     , prdev,prsig,prcut,bias_rcond
      logical bias_debug
      common/bcrit/wldev,wlsig,wlcut,wlrat,wldmax
     .           , nldev,nlsig,nlcut,nlrat,nldmax
     .           , prdev,prsig,prcut,bias_rcond
     .           , bias_debug
  
c     <constants> -- constants used all over the place   
      real*8 pi,convd
      common/constants/pi,convd

c     <parmflags> -- logical flags for partials, free/fix, and weights
      logical sitpar,satpar,zenpar,gradpar,eoppar
     .       ,sitest,satest,zenest,gradest,eopest,svantest
     .       ,sitwgt,satwgt,zenwgt,gradwgt,eopwgt,svantwgt,clkwgt
      common/parmflags/sitpar,satpar,zenpar,gradpar,eoppar
     .                ,sitest,satest,zenest,gradest,eopest,svantest
     .                ,sitwgt,satwgt,zenwgt,gradwgt,eopwgt,svantwgt
     .                ,clkwgt

c     <flags> -- flags for program flow and output
      integer*4 l2flag,lquick,iseen(maxsat),iqflag,ioflag,ihmode
      logical gloprt,logprt,do_loose 
      real*4 correl_prt,coord_upd_tol
      common/flags/l2flag,lquick,iseen,iqflag,ioflag,gloprt
     .            ,correl_prt,logprt,ihmode,coord_upd_tol,do_loose
                
c     <rxant> -- receiver and antenna characteristics for q-, o- and h-files         
      character*3 rcvr_sw(maxsit)
      character*20 rcvr_type(maxsit),rcvr_sn(maxsit),ant_type(maxsit)
     .           , ant_sn(maxsit)
      real*4 rcvr_swver(maxsit)   
      common/rxant/rcvr_type,rcvr_sn,ant_type,ant_sn,rcvr_swver,rcvr_sw
      
c    <freqs> -- frequency variables       
      real*8 gear
      common/freqs/gear
   
c** --rwk 190508: These added from subroutine-specific commons:    
                               
c     <SOLVE version and login>
      character*3 owner
      character*40 vers
      common/version/vers,owner

      character*80 gline
      common/titles/gline

      character*16 minf,moutf,linf,loutf,ginf,goutf,iinf,ioutf
      common/upfilnam/minf,moutf,linf,loutf,ginf,goutf,iinf,ioutf

      integer*4 ipfil
      common/upfil/ipfil(8)

c     <slvaux> -- original sum of residual square and right hand terms
      real*8 borg(maxprm)
      common/slvaux/borg

      real*8 coords(maxcrd),finv,semi
      integer*4 idatum
      common/sitcrd/coords,finv,semi,idatum
        
c     <l12wl>  -- memory of original adjustments and chi2
      real*8 adorg(maxprm),chi2
      common/l12wl/adorg,chi2
     
      real*8 ddwl(maxbis/2),ddwv(maxbis/2)
      common/wld/ddwl,ddwv
         
c     <dstat>  -- unreliable station index
      integer*4 lbad(maxsit),lbfre(5)
      common/dstat/lbad,lbfre 
             
c     <lhalf> 
      integer*4 lwave(maxsit,maxsat,2),half(maxbis)
      common/lhalf/lwave,half
                          
c     <limit>
      common/limit/limitb
      integer*4 limitb

c     <situse>
      integer*4 jusit(maxsit)
      common/situse/jusit
             
c    <intsam>
      integer*4 intsam(maxsit)   
      common/samint/intsam
      
c      <acbiases> 
      real*8 wlval(maxdd),wlconf(maxdd)
      integer*4 numwl
      character*1 rorx(maxdd) 
      logical miss_bias(maxdd)
      common /acbiases/ wlval,wlconf,numwl,rorx,miss_bias
       
c     <wl> 
      real*8 vwl(maxsit,maxsat),wl0(maxsit,maxsat)
      integer*4 idwl(maxsit,maxsat)
      common/wl/vwl,wl0,idwl

c     <wlmem> 
      real*8 reml1(maxsit,maxsat),reml2(maxsit,maxsat)                           
      common/wlmem/reml1,reml2

c    <wL1> 
      real*8 wl1(maxsit,maxsat)     
      common/wl1/wl1
                    
c     <gaptec> 
      integer*4 igaps(maxsit,maxsat)
      common/gaptec/igaps

c     <dmap> 
      integer*4 nrd,ncd,ipnt2d(maxobs)
      common/dmap/nrd,ncd,ipnt2d 

c     <sats> -- satellite PRNs from c-files
      integer*4 isprn(maxsat) 
      common/sats/isprn 
                       
c     <bbii2> 
      real*8 bdev(maxbis)
      integer*4 nfix(maxbis)
      common/bbii2/bdev,nfix
            
c     <ambopt> -- options for ambiguity resolution
      common/ambopt/ noptin
      integer*4 noptin
            
c     <erp1> <erp2> <ante> -- rotaton and antenna parameters for H-file
      integer*4 jde,jdr
      real*8 te,tr,ut1,ut1dot,xp,xpdot,yp,ypdot
      common/erp1/jde,jdr,te,tr,ut1,ut1dot,xp,xpdot,yp,ypdot
      real*8 psi,psidot,eps,epsdot     
      common/erp2/psi,psidot,eps,epsdot        
      real*8 anoff,svantdx
C MOD TAH 200126: Changed to svantdx(3,2,maxsat) (added L1/L2 index)
C MOD TAH 200205: Added s_antdaz(maxsit) to save the antenna offset
C     from true North (deg).
      real*8 s_antdaz(maxsit)
      common/ante/anoff(maxsit,3,3),svantdx(3,2,maxsat), s_antdaz
                  
c     <data_noise> -- a priori data noise and models
      real*8 sit_err(maxsit),sit_elv(maxsit),sat_err(maxsat)
      character*10 err_mod(maxsit)
      common /data_noise/sit_err,sit_elv,sat_err,err_mod
 
c     <loading> --values for estimating and recording average loading
      real*8 atmload(3,maxsit),hydload(3,maxsit),amatl(5,5,maxsit)
     .      ,batm(5,maxsit),bhyd(5,maxsit),sumweight(maxsit)
     .      ,count(maxsit),sumatmload(3,maxsit),sum1atmload(3,maxsit)
     .      ,sumhydload(3,maxsit),sum1hydload(3,maxsit)
      logical arithavg(maxsit)
      common /loading/atmload,hydload,amatl,batm,bhyd,sumweight,count
     .      , sumatmload,sum1atmload,sumhydload,sum1hydload,arithavg
     
c     <biasod> 
      common/biasod/ibod
      integer ibod

c     <keycom>
      character*5  keyword(15)
      common/keycom/keyword   
      
c**  rwk 190516 -- add this to faciliate bias removal print and debugging

c     <current epoch>    
      integer*4 iepoch 
      common/epcom/iepoch
