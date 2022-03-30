c     orbfit.h --  common variables for orbfit  -- last modied rwk 140930

c       Dimensions
       
c       dimension limits: maxsat and maxorb defined in ../includes/dimpar.h
c                         maxtfil and maxglb defined in ../includes/orbits.h
      integer*4   maxepc,mxoprm
      parameter   (maxepc=(maxtfil-1)*96+28)
      parameter   (mxoprm=maxglb+maxorb*maxsat)
       
c       Unit numbers
      integer*4   lucmd,lut(maxtfil),luplt(maxsat),lurms,lufit,lug,lusvs
     .           ,iscrn,iprnt,inut,iut1,ipole
      common/units/lucmd,lut,luplt,lurms,lufit,lug,lusvs,iscrn,iprnt
     .            ,inut,iut1,ipole

c       T-file header information 
      integer*4   nics(maxtfil),nintrs(maxtfil),nepcht(maxtfil)
     .        ,   nprn(maxtfil),itsat(maxsat,maxtfil)
     .        ,   jdb(maxtfil),jdstp(maxtfil),jde
      real*8      tb(maxtfil),tstp(maxtfil),te,delt(maxtfil)
     .        ,   satics(maxorb,maxsat) 
      character*1 gnss(maxsat,maxtfil)  
      character*4 icsnam(maxorb)
      character*5 precmod,frame,srpmod,nutmod,gravmod,eradmod,antradmod
      character*16 satnamt(maxsat,maxtfil)
      common /tfhdrs/tb,tstp,te,delt,satics,nics,nintrs,nepcht,nprn
     .             ,itsat,jdb,jdstp,jde,satnamt,icsnam
     .             ,precmod,frame,srpmod,nutmod,gravmod
     .             ,eradmod,antradmod,gnss

c       SVs and epochs actually used
      integer*4 nsat,isat(maxsat),nepoch,iepstart,iepstop 
      character*16 satnam(maxsat)
      common /svepochs/ nsat,isat,nepoch,iepstart,iepstop,satnam

c       Reference vector, O-Cs and partials for every epoch (save to calculate residuals)
      real*8      ref_pos_vel(6,maxsat,maxepc),omc(3,maxsat,maxepc)
     .           ,part(3,mxoprm,maxsat,maxepc),drac(3,maxsat,maxepc)
     .           ,prefit_sum2
      integer*4 nobs(maxsat)
      logical   lobs(maxsat,maxepc)
      common/omcpart/ ref_pos_vel,omc,part,drac,prefit_sum2,nobs,lobs
                  
c       Parameter values, adjustments, and normal equations
      integer*4   nparam,islot(mxoprm) 
      real*8      apr_prm(mxoprm),adjust(mxoprm)
      real*8      amat(mxoprm,mxoprm),bvec(mxoprm)
      character*30 prmnam(mxoprm)
      common/nrmcom/nparam,islot,amat,bvec,apr_prm,adjust,prmnam

c       Max orbit misfit tolerance
      integer     nbad_sat, ibad_sat(maxsat)
      real*8      max_fit_tol, bad_rmstot
      common/fittol/bad_rmstot,max_fit_tol,nbad_sat,ibad_sat
       
