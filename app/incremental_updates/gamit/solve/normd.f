c
      Subroutine NORMD
c
c     open simultaneously observed c-files and fill normal equations
c
c flow chart :         `/
c
c     NORMD
c        |
c        +---> (initialization)
c        |
c        +---> (epoch loop)
c           |
c           +---> EPTIME                        } get epoch time
c           |
c           +---> (c-file loop)
c              |
c              +---> GETCDA                     } read data records
c              |
c              +---> FILOMC                     } fill o-c data arrays
c           |  (if l2flag>0)
c           +---> PSEUWL                        } derive one-way WL data
c           |       (single-difference)
c           |    +---> BSTAT -------+           } remove useless one-way data
c           +--->|  (double-difference)
c                +---> DOUBLE       |           } remove useless one-way data
c           |           (yes)       |
c           +---> DCHECK-->skip DOUBLE,OPERA    } compare effective useage with last epoch
c           |                       |
c           +---> STACK             |           } determine if need to stack data
c              |  (yes)             |
c              +---> FORMN2         |           } form stackable part of normal matrix
c           |                       |
c           +---> LRGINT            |           } take off large integers
c           |                       |
c           +---> DOUBLE            |           } form double-difference operator
c           |                       |
c           +---> FILOBS    <-------+           } fill observation matrix
c           |
c           +---> OPERA ---> OPERA2             } fill weighting matrix
c           |
c           +---> FXBIAS                        } fix biases for gap and flag
c           |
c           +---> FILPAR                        } fill partials into a submatrix
c           |
c           +---> FORMN1                        } form normal matrix (non-stackable part)
c           |
c           +---> FORMN2                        } form normal matrix (stackable part)
c           |
c           +---> (chi2)                        } sum up residual squares
c        | (if l2flag>0)
c        +---> PSEUWL                           } take off large integers from WL data
c        |
c        +---> QHEAD2                           } output number of observables
c        |
c        +---> FXBIAS                           } fix biases for gap and flag
c        |
c        +---> APPLY                            } atmos. and ionos. constraints
c
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'        
      include 'parameters.h' 
      include 'models.h'

      character*1 frempc(12)
      character*256 message    

      logical istack 

      character*3 upperc
             
      integer*4 ierrfl(maxsat)
     .        , iusebk(maxsit,maxsat)
     .        , iphi(maxsit,maxsat),iphibk(maxsit,maxsat)
     .        , iquick(maxsit,maxsat),iertag(maxsit,maxsat)   
     .        , iphilv(maxobs),iwork(maxprm),isnr(maxdat,maxsat)
     .        , ij4,ibias1,ibias2,ibias3
     .        , isame,ifile,istat,igood,ipflg,iobs,ieflg,iones
     .        , ifac,iobsc,iobs2,iones2,irow,ipcnt,jsat,ietag,is,i1,key
     .        , j
cd    integer*4 ii,jj

      real*8  first(maxsit,maxsat),omc(maxdat,maxsat)
     .      , phi(maxsit,maxsat),phi2(maxsit,maxsat)    
     .      , first2(maxsit,maxsat),al1(maxobs),al2(maxobs)   
     .      , bl1(maxobs),bl2(maxobs),scpu(maxdm)
     .      , el(maxsat),elev(maxsit,maxsat),small
     .      , g1,g2,gearf,gsq,r2lc,dinter,r2sum1,r2sum2,dum

* MOD TAH 190612: Added data statement for iones variable which is used 
*     before it is initialized.  iones is first seen in the code in formn2
*     but is initialized later (in the code) in filobs.  There are goto 
*     statements so it is possible filobs is called before formn2 but this
*     does not seem to be case.  Using -finit-integer=100000 gfortran allows
*     this to be tested.  Using the data statement best emmulates the 
*     assumption on zero memory and the local routines remember there values.
*     Actual bug found so data statement not needed.
*     data iones / 0 /


c ***DEBUG 
c      logical debug1/.true./,debug2/.true./
      logical debug1/.false./,debug2/.false./
      integer*4 i

      equivalence  (scpu,dr)
      equivalence  (isigma,iwork)

c  gear  now set in gethed.f from fL1, fL2 and stored in solve.h
c*      data gear/0.779220779220779d+00/
c*       1/(77/60-60/77)    
c  gearf is either gear or, if l2-only, =0, and used to call 
c       formn1, formn2, opera, opera2.
c        
c     data facl2/1.983683984542722d0/
      data small/1.0d-12/

cd     if(debug1) then 
cd        print *,'NORMD start islot1(1) islot2(1,1) r2sum coords(1-3) '
cd     .    ,islot1(1),islot2(1,1),r2sum,coords(1),coords(2),coords(3)  
cd        print *,'  maxprm ntpart rlabel(1) rlabel(456)  '
cd     .         ,   maxprm,ntpart,rlabel(1),rlabel(456)
cd      endif 

c     initialization to avoid compiler warning
      g1 = 0.d0
      g2 = 0.d0
c
c     for L1-only or L2-only solution
      if(l2flag.le.0) gearf = 0.d0
      if(l2flag.gt.0) gearf = gear
      gsq = gearf*gearf
c
c     key = 1: record LC mode normal matrix
      key = 0
      if (ibcnt.gt.ibcnt1) key = 1
c     write(*,*) 'ibcnt,ibcnt1,key',ibcnt,ibcnt1,key
c
      r2lc = 0.0d0
c     number of biases this session
      ibias1 = nsat*nsite
      ibias2 = (ibias1*(ibias1-1))/2
      ibias3 = ibias1*(ibias1+1)/2
      call zero1d(1,ibias1,wl1)
      call zero1d(1,ibias2,cphi)
      call zero1d(1,ibias2,cphik)
      call zero1d(1,ibias3,an22)
      call zero1d(1,ibias3,bn22) 
c      write(*,'(a,i7,1x,/,(10e9.2))') 'NORMD ibias3 BN22 '
c     .        ,ibias3,(bn22(i),i=1,ibias3+1)
      call zero1d(1,nsat,el)
      ij4 = lpart*ibias1
      call zero1d(1,ij4,clc)
c
      call zero2d(maxsit,maxsat,nsite,nsat,first)
      call zero2d(maxsit,maxsat,nsite,nsat,first2)
      call zero2d(maxsit,maxsat,nsite,nsat,wl0)
      call zero2d(maxdat,maxsat,maxdat,nsat,omc)
      call zero2d(maxsit,maxsat,nsite,nsat,elev)
      call zero2i(maxsit,maxsat,nsite,nsat,iquick)
      call zero2i(maxsit,maxsat,nsite,nsat,iphibk)
      call zero2i(maxsit,maxsat,nsite,nsat,iusebk)
      call zero2i(maxsit,maxsat,nsite,nsat,idwl)
      call zero2i(maxsit,maxsat,nsite,nsat,iertag)
c
      dinter = dble(inter)
* MOD TAH 190612: Added back code from active normd.f to initialize
*     the isame variable.  A non-zero value can cause istack to be
*     set to false is stack which in turn will cause formn2 to be 
*     called before iones is initialized.  No comments in code as to
*     why orginal lines were removed.
      isame = 0
      if (lquick.ge.1) isame = 1
* END of restored lines.

      istack = .true.
      call report_stat('STATUS','SOLVE','normd',' '
     .          , 'Reading data and forming normal equations',0)

c** rwk 190813   Zero out the operator for remaining biases
        do i=1,maxobs
           iphilv(i) = 0 
         enddo


c Begin epoch loop
      do 100 iepoch = 1,iend

       if(debug1) print *,'NORMD epoch ',iepoch  
cd       if(debug.and.iepoch.gt.573 ) stop 1   
     

c*     this for debug, put epoch number in common
         kepoch = iepoch 
         if( mod(iepoch,200).eq.0 ) then
           call eptime( frempc )
           if( iepoch.lt.istart ) then
              write( message,'(a,i4,a,4x,12a1,a,i4)') ' Epoch <'
     .              ,iepoch,' >', frempc,' < start epoch ',istart
               call report_stat('STATUS','SOLVE','normd',' ',message,0)
           else
             write( message,'(a,i4,a,4x,12a1)') ' Epoch <',iepoch,' >'
     .                                         , frempc
             call report_stat('STATUS','SOLVE','normd',' ',message,0)
           endif
         endif
                      
c** rwk 190813: Need to zero-out the mapping operator completely since
c**             it is checked to set iphilv.
         do j=1,maxsat
           do i=1,maxsit 
             iphi(i,j) = 0 
           enddo
         enddo         
 
      ifile = 30      
             
      if(debug2) print *,'NORMD calling getcda nsat  ',nsat
      do istat = 1,nsite
         ifile = ifile+1
c         read data records
cd         if(iepoch.eq.952.or.iepoch.eq.953)  
cd      .       print *,'NORMD reading iepoch ifile istat '
cd     .                              ,iepoch,ifile,istat
         call getcda(ifile,istat,ierrfl,omc,isnr,el)
        if(debug2) print *,'NORMD iepoch istat aft getcda A '
     .      ,iepoch,istat,(a(i),i=1,3) 
c         fill phase and wl arrays with omc data
         call filomc(istat,phi,phi2,iphi,omc,
     .        isnr,ierrfl,iertag,el,elev)           
         if(debug2) print *,'NORMD iepoch istat aft filomc A '
     .       ,  iepoch,istat,(a(i),i=1,3) 
      enddo                                         
      if(debug1) then 
        print *,'NORMD aft filomc iphi(4,,i) ',iepoch,(iphi(4,i),i=1,32)
        print *,'NORMD aft filomc iphilv(97) ',iphilv(97)
      endif 

      if (iepoch.lt.istart) then
c*        comment this out.  why needed?   rwk 950812
         goto 100
      endif

      if (l2flag.gt.0) call pseuwl(iphi,1)

c       remove useless double difference observations
       call uselss(phi,iphi,igood)
       if(l2flag.gt.0) call uselss(phi2,iphi,igood)
       if(debug2) print *,'NORMD aft uselss igood nsite nsat nlive '
     .                                      ,igood,nsite,nsat,nlive 
      if (igood.eq.0.and.iepoch.lt.iend) go to 100
cd      if(debug2) then 
cd        do ii= 1,nsite
cd          print*,'NORMD aft uselss site,iepoch,iuse: ',cfiln(ii),iepoch
cd     .        ,(iuse(ii,jj),jj=1,nsat)
cd        enddo
cd      endif     
      if(debug1) then 
        print *,'NORMD aft uselss iphi(4,i) ',(iphi(4,i),i=1,32)  
      endif 

c      stop 'In epoch loop'
c      check if operator changed from last epoch
      call dcheck(iphi,iphibk,iertag,ipflg,ieflg)
      if(debug1) print *,'NORMD aft dcheck nsat igood ',nsat,igood 
c      determine if the data should be stacked

      if (iepoch.lt.iend) call stack(istack,isame,ipflg,ieflg)    
      if(debug1) print *,'NORMD aft stack nsat igood ',nsat,igood 
      if (iepoch.eq.iend) then
         istack = .false.
         ipflg = 1
      endif

      if (.not.istack) then
         call formn2(gearf,iphilv,iones,isame,key,bl1,bl2)
        if(debug1) print *,'NORMD aft formn2 nsat,igood ',nsat,igood
        if(debug1)  print *,'NORMD aft formn2 iphilv(97) ',iphilv(97)

         if (iepoch.eq.iend.and.igood.eq.0) go to 100
         isame = 1
         istack = .true.
      endif
c       take off large integer from carrier-beat-phases
      call lrgint( phi,iphi,first,iertag )
      if(l2flag.gt.0) call lrgint(phi2,iphi,first2,iertag )
      if(debug1) then 
        print *,'NORMD aft lrgint nsat igood ipflg '
     .     ,nsat,igood,ipflg
        print *,'  iphi(4,i) ',(iphi(4,i),i=1,32)  
      endif 
      if(ipflg.eq.1)
     .   call double(iphi,iobs,iusebk)  

      call add2i(maxsit,maxsat,nsite,nsat,iusebk,iuse)

      if(iobs.eq.0.and.lquick.ge.1) go to 60
      if(iobs.eq.0) go to 100
      nd2obs = nd2obs+iobs
c       fill observation matrix (by station)
      call filobs(gearf,iphi,phi,phi2,iones,al1,al2)
      nones = nones+iones
      if(iones.eq.0.or.ipflg.eq.0) go to 60
c       form w-covariance matrix                         
      if(debug2) print *,'NORMD calling opera nsat nsat ',nsat 
      call opera(phi,iobs,iones,gearf,elev) 

c       form we-covariance matrix for ionospheric constraint
c       skip if single-differences
      if(l2flag.ge.2)
     .   call opera2(phi,iobs,iones,gearf,elev)   
c       4 oneways per double difference, 2 oneways per single difference

      ifac = 4
      iobsc = ifac*iobs
      iobs2 = iobs+2
      iones2 = iones+2
      call zero1d(1,iobsc,d)
      call zero1d(1,iobsc,dt)
      call zero1i(1,iobsc,ipntd)
      call zero1i(1,iobsc,ipntdt)
      call zero1i(1,iobs2,irowd)
      call zero1i(1,iones2,irowdt)
 60   continue
      irow = 0
      ipcnt = 0      

c   Loop over stations
    
cd      if(debug2) print *,'NORMD nsite nsat ',nsite,nsat 
      do  istat = 1,nsite
              
c        zenith model now set in CFMRG and read from M-file.
c        set the number of zenith delays for this station and the break epochs
         if( upperc(zenmod).ne.'CON' .and. upperc(zenmod).ne.'PWL' )
     .      call report_stat('FATAL','SOLVE','normd',' '
     .                      ,'Zenmod not correct from M-file',0)
                  
c    Loop over satellites
         
         do j=1,3            
           atmlavg(istat,j) = 0.d0 
           hydrolavg(istat,j) = 0.d0  
         enddo
         do  jsat = 1,nsat
           ietag = iertag(istat,jsat)
           is = iphi(istat,jsat)
           if(debug1) print *,'NORMD istat,jsat ietag is '
     .       ,istat,jsat,ietag,is 
c          skip if no observations
           if(is.ne.0.0) then
c            fix biases in the following two cases:
c            1. quick solution mode after a gap (put virtual bias in filomc)
c            2. regular solution mode with a bias flag (ietag = 10)
c
             if(ietag.eq.10) then
               call fxbias(-jsat,istat,ietag,r2lc)

               iertag(istat,jsat) = 0
               if(debug1) print *,'NORMD ietag jsat istat r2sum r2lc '
     .                                , ietag,jsat,istat,r2sum,r2lc 
             endif

             iseen(jsat) = iseen(jsat)+1
             if(debug1) print *,'NORMD iepoch jsat iseen '
     .         ,iepoch,jsat,iseen 
* MOD TAH 190601: Stop below seems to be debug.  Removed to test solve.
C            if(iepoch.gt.949) stop 1 
             irow = irow+1
c              fill partials into a submatrix {c}
cd              if(debug1) then
cd                print *,'NORMD iepoch istat jsat islot1 '
cd     .                       , iepoch,istat,jsat,islot1 
cd                print *,'NORMD ipcnt npart ',ipcnt,ntpart
cd                do ii=1,ntpart
cd                  print *,rlabel(ii)
cd                enddo
cd              endif                                 
              call filpar( istat, jsat, ipcnt )
              irowc(irow+2) = ipcnt 
c             if loading used, increment the normal equation for estimating the average
              if( atmlmod.ne.'       '.or.hydrolmod.ne.'        ' ) then
                call avgload( 1,istat,jsat,elev,iphi )
              endif            
           endif 
c        end loop on satellites
         enddo 
c     end loop over stations
      enddo   

      if(irow.eq.0) go to 100
      irowc(1) = ipcnt
      irowc(2) = 0
c       transpose design (partials) matrix
cd      print *,'NORMD xtra iepoch ',iepoch 
cd      if( iepoch.ge.940.and.iepoch.le.980)  then                 
cd         print *,'maxprm maxobs maxnd maxdm maxwm1 '
cd     .          , maxprm,maxobs,maxnd,maxdm,maxwm1 
cd          jj= 0 
cd          print *,'NORMD iepoch iones lpart ',iepoch,iones,lpart
cd          do ii=1,100               
cd           write(*,'(30i5)') (ipntc(jj+j),j=1,30)
cd           jj= jj+30
cd          enddo
cd      endif           
      call transx(c,ct,ipntc,ipntct,irowc,irowct,iones,lpart,iwork)     


c     copy ct to scpu.
      if (istack.and.isame.gt.0) then
         if (isame.eq.1.or.lquick.ge.1) then
            call copy1d(1,ipcnt,0,ct,scpu)
            call copy1d(1,iones,0,al1,bl1)
            call copy1d(1,iones,0,al2,bl2) 

         else
            call add1(1,ipcnt,ct,scpu)
            call add1(1,iones,al1,bl1)
            call add1(1,iones,al2,bl2)

         endif
      endif
c       form first part of normal matrix (not stackable)    
      if(debug1)  print *,'NORMD bef formn1 iphilv(97) ',iphilv(97)
      call formn1(gearf,iphi,iphilv,iones,0,key,al1,al2)             
      if(debug1)  print *,'NORMD aft formn1 iphilv(97) ',iphilv(97)

c     add contribution to chi-squared
      r2sum1 = 0.d0
      r2sum2 = 0.d0 
                                           
      call gtpgx(al1,cphi,iones,1,r2sum1,work,dum)

      if(debug2) print *,'NORMD iones r2sum1 ',iones,r2sum1
      if(l2flag.gt.1) call gtpgx(al2,cphik,iones,1,r2sum2,work,dum) 
      if(debug2) print *,'NORMD iones r2sum2 ',iones,r2sum2
      if(l2flag.ge.-1) then
      g1 = 1.d0-gsq
      g1 = g1*g1
      g2 = 4.d0*gsq
      elseif(l2flag.eq.-2) then
      g1=1.d0
      g2=0.d0
      endif   
      if(debug1) print *,'NORMD epoch bef incre r2sum ',iepoch,r2sum 
      r2sum = r2sum+(g1*r2sum1+g2*r2sum2)
      r2lc = r2lc+g1*r2sum1
      if(debug1)  print *,'NORMD g1 r2sum1 g2 r2sum2 r2sum,r2lc '
     .   ,g1,r2sum1,g2,r2sum2,r2sum,r2lc    
cd     if( kepoch.eq.281 ) stop
      call zero1d(1,iones,al1)
      call zero1d(1,iones,al2)
      call zero1d(1,ipcnt,c)
      call zero1i(1,ipcnt,ipntc)
      i1 = iones+2
      call zero1i(1,i1,irowc)  
            
      do jsat = 1,nsat
        do istat = 1,nsite
           iphibk(istat,jsat) = iphi(istat,jsat)
        enddo
      enddo                       
  
cd      if(debug1) print *,'NORMD end epoch loop nsat  ',nsat 
c     end epoch loop  

  100 continue

      do jsat = 1,nsat
        do istat = 1,nsite
          wl1(istat,jsat) = first(istat,jsat)-first2(istat,jsat)
          if (dabs(first(istat,jsat)).lt.small.and.
     .       dabs(first2(istat,jsat)).lt.small)
     .    wl1(istat,jsat) = reml1(istat,jsat)-reml2(istat,jsat)
        enddo
      enddo                       
c
      if (l2flag.gt.0) call pseuwl(iphi,2)
c
c     print number of double difference observations by satellite
c
      call qhead2(phi,phi2,2)              
c
c     take out implicit biases, for quick algorithm, single
c     differences, and for implicitly declared biases
      if(debug1) print *,'NORMD before fxbias ',nsat 
      do istat = 1,nsite
         call fxbias (nsat,istat,0,r2lc)
         if(debug1) print *,'NORMD nsat istat r2lc '
     .                        ,    nsat,istat,r2lc 
      enddo
      bn22(1) = r2lc

c Calculate atmospheric-delay constraint matrix for  double differences

      if( iatcon.eq.1 )  call apply

c Calculate average loading corrections for the h-file

      if( atmlmod.ne.'       '.or.hydrolmod.ne.'        ' ) then    
c          all/ arguments dummy here except the first
           call avgload( 2,istat,jsat,elev,iphi )
      endif 

      return
      end
