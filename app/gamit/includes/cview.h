c     CVIEW.FTI
c     for larger versions 

c     maximum number of sites in cview
      integer ncvsit
      parameter(ncvsit=80)

c     all the variables which occur in common blocks 
c     should be declared here

c     Time tags, expressed as seconds of day
      real*8 tx(maxepc)
      real*8 tag(maxepc,ncvsit)
      common/timetg/tx,tag

c     phases
      real*8 yl1(maxepc,maxsat,ncvsit)
     .,      yl2(maxepc,maxsat,ncvsit)
      common/rdata/yl1,yl2

c     pseudoranges in cycles
      real*8 pr1(maxepc,maxsat,ncvsit)
     .,      pr2(maxepc,maxsat,ncvsit)
      common/prange/pr1,pr2

c     multipath delay 
c     comment this out to save time and space
ckf   real*4   rd1(maxepc,maxsat,ncvsit)
ckf   .,       rd2(maxepc,maxsat,ncvsit)        
ckf   .,       rw1(maxepc),rw2(maxepc)
ckf    common /multip/rd1,rd2,rw1,rw2

c     GAMIT error flags                  
c     readc,writec,makobs,readx,readc
      integer*4 ierr(maxepc,maxsat,ncvsit)
      common/errflg/ierr

c     satellite elevations in radians 
c     careful of real*4!
      real*4 ela(maxepc,maxsat,ncvsit)
      real*8 elw(maxepc)
      common /elevat/elw,ela

c     satellite azimuths in radians clockwise from North
c     careful of real*4!
      real*4 aza(maxepc,maxsat,ncvsit)
      real*8 azw(maxepc)
      common /azimut/azw,aza

c     clock offset terms in seconds
c     careful of real*4!
      real*4  clk(maxepc,maxsat,ncvsit)
      real*8  clw(maxepc)
      common /clocks/clk,clw

c     readm,addadj,readc,pointr,writev
      real*8 adjust(maxprm)
      integer*4 idms(maxprm),mtpart,islot(maxprm)
      common/postft/adjust,idms,mtpart,islot 

c     multiple zenith-delay and gradient parameters
c      readc,writec,addadj,finish,getfil,readm,write_vciew
      integer*4 numzen,numgrad,idtzen(maxatm),idtgrad(maxgrad/2)
      character*3 zenmod,gradmod
      common/zendel/numzen,numgrad,idtzen,idtgrad,zenmod,gradmod

c     data interval
      integer*4 inter,ircint,inext,jnext(ncvsit)
           
c     start time in year,month,day,hours,minutes,seconds
      integer*4 iit0(3)
      real*8 tt00(3)

c     ireadm,readc,ptitle
      common/times/tt00,inter,ircint,inext,jnext,iit0

c     values of L1 and L2 from and for C-files and their quality flags
      real*8  cl1(maxepc),cl2(maxepc)
      integer*4 kcc(maxepc)
      common /holdc/cl1,cl2,kcc

c     working values of L1 and L2 and their quality flags
      real*8  wl1(maxepc),wl2(maxepc)
      integer*4 kww(maxepc)
      common /workin/wl1,wl2,kww

c     buffered values of L1 and L2 and their quality flags
      real*8  bl1(maxepc),bl2(maxepc)
      integer*4 kbb(maxepc)
      common /buffin/bl1,bl2,kbb

c     values of P1 and P2 from C-files
c     only display these, don't edit them
c     assume that the error flags in kww and kww apply to them, too
      real*8  pc1(maxepc),pc2(maxepc)
      common /rangin/pc1,pc2

c     set up lambda arrays containing wavelength factors
      integer*2 lambds(ncvsit,maxsat,maxdat)         
c        lambda = 1 for unambiguous, undoubled values
c        lambda =-1 for   ambiguous, undoubled values
c        lambda = 2 for unambiguous,   doubled values
c        lambda =-2 for   ambiguous,   doubled values
      common /lambd/lambds   

c     first and last epochs of marginal data
      integer*2 kk0(maxsat,ncvsit),kk1(maxsat,maxsit)
      common /kk01/kk0,kk1

c     GNSS tyep and array of PRN (satellite ID) numbers
      integer*4 isprn(maxsat)
      character*1 gnss
      common /prns/isprn,gnss

c     primary L-band frequencies and their combinations assigned in 
c     getfil from values in libraries/includes/freq_def.h
      real*8 fL1(maxsat),fL2(maxsat),gear(maxsat),faclr(maxsat)
     .     , facl2(maxsat),facwl(maxsat)
      common/freqs/fL1,fL2,gear,faclr,facl2,facwl
ccc      parameter (faclr = 1.983683984543d0)
ccc      parameter (facl2 = 0.220779220779d0)
ccc     60/77
ccc      parameter (gear  = 60.0d0/77.0d0)
ccc     17/137
ccc      parameter (facwl = 17.0d0/137.0d0)
ccc      parameter (frq1  = 154.0d0 * 10.23d6)
ccc      parameter (frq2  = 120.0d0 * 10.23d6)   


c     speed of light in meters/second              
      real*8 clight
      parameter (clight = 299792458.d0)

c     number of possible menus (should be an even number, because there 2 rows)
      integer nmenus
      parameter (nmenus=20)

c     number of possible data types
      integer ntypes
      parameter (ntypes=18)              

c     maximum order of polynomial fit
      integer maxply
      parameter (maxply=19)

c     maximum number of plotting windows
      integer   maxwin                  
      parameter (maxwin = 5)
           
c     maximum number of graph types
c     time series, skymap, and spectrum aand allan
      integer   maxgph 
      parameter (maxgph = 4)
           
c     length of a big GAP, in 1-way and dd
      integer   biggap1,biggap
      parameter (biggap1 = 60)
      parameter (biggap  = 15)
           
c     maximum size of the list KLIST
      integer   maxlist,klist
      parameter (maxlist = 500)
      common /cviewlist/klist(maxlist,5)
           
c     elevation angle cutoff
      real*8 elvcut(maxsit)
      common /cutoff/elvcut

c***  KLUGE warning: With the current structure of forming 
c     site/sat differences first and then frequency combinations 
c     prevents forming between-satellite (or double) differences
c     for the Glonass FDMA satellites.  By putting into common
c     the pointer to the first satellite (channel), it will be
c     possible to reference correctly the frequency quantities for
c     one-way and between-site (single) differences.  This variable
c     is set from ichan1 in subroutine editor. rwk 151230
      integer*4 jsat1
      common /freq_index/jsat1 








