
      Subroutine DBIAS( mbias0,last_nonbias,last_bias )

C---- Calculate new N12,N21,N22,U2 which correspond independent bias
C        parameters.
C     Formula:
C        D = D-optimal,  D~ = transpose(D)
C        D^ = inverse(DD~)*D
C        N12 = n12*D^~ = transpose(N21)
C        N21 = D^*n21
C        N22 = D^*n22*D^~
C        U2 = D^*u2
C           -DND- 871230
C---- Divide D into several small matrixes to save space and time.
C            | D(L1)     0   |         | D^(L1)    0    |
C        D = |               | ,  D^ = |                |
C            |   0     D(L1) |         |   0    D^(L1)  |
C        D^(L1)=inverse(D(L1)D(L1)~)*D(L1)
C                | n11 n12 n13 |              | N11 N12 N13 0  |
C       old  n = | n21 n22 n23 | ,   new  N = | N21 N22 N23 0  |
C                | n31 n32 n33 |              | N31 N32 N33 0  |
C                                             |  0   0   0  0  |
C        N11 = n11
C        N21 = D^(L1)*n21     N31 = D^(L1)*n31
C        N22 = D^(L1)*n22*D^(L1)~
C        N32 = D^(L1)*n32*D^(L1)~
C        N33 = D^(L1)*n33*D^(L1)~
C        U2  = D^(L1)*u2      U3 = D^(L1)*u3
C     For simplicity, using D^ represents D^(L1).
C        N21 = | N21 |        U2 = | U2 |
C              | N31 |             | U3 |
C              | N22  N23 |
C        N22 = | N32  N33 |
c           -DND- 880127
C
C     Revised to fit multi-session mode.
C                                 -DND-   880227

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer maxbas
      parameter (maxbas=maxsit*(maxsit-1)/2)

      integer lm2,lm3,iel21,ncomon,lone,ier,nxa,jc,ic,ir,i1,lm1
     .      , mbias0,last_nonbias,last_bias,jc1,ibnd,irow,le22
     .      , i2,l2,ir1,ifirst,i,j,k

      real*8 dl1,dl2,dlc,dll,small,ddr,rcond

* MOD TAH 980121: Added new function jelf(i) which returns
*      jelf = i*(i-1)/2.  Routine also checks that the calculation
*     done correctly.  (Check of HPUX 10.2 +O3 Bug)
      integer jelf
c
      real*8 temp12(maxwm2*2)
      real*8 worklc(maxobs) 
        
c     initialization to avoid compiler warnings, but traps may be needed
      dl2 = 0.d0
      dl1 = 0.d0
      dll = 0.d0
      jc1 = 0


c      print *,'DBIAS irowd ',irowd(1)
      small = 1.0d-11
      ifirst = idxb(mbias0+1)
      if (ifirst.lt.0) ifirst = -ifirst
      ncomon = idxb(1)
      if (ncomon.lt.0) ncomon = -ncomon
      ncomon = ncomon-1
      lone = nsat*nsite
      ibnd = iband
      if (l2flag.eq.1) ibnd = 1
c
c---- find number of rows and number of columns for d
      nrd = irowd(1)
      ncd = irowdt(1)
c
c---- calculate dd~  (nrd x nrd)
      call atpdal(d,ipntd,irowd,nrd,ncd,1.0d0,dpdt,work2)
c
c---- calculate inverse(dd~)
c      call nrmsc2(dpdt,work2,work3,nrd,0,0)   
      call inver2 (dpdt,work2,1,nrd,rcond,ier)
      if (ier.eq.130)
     .  call report_stat('WARNING','SOLVE','dbias',' '
     .          , 'Bad inverse DPDT',0)
c      call nrmsc2(dpdt,work2,work3,nrd,1,0)
c---- calculate d^
      call dbar( dr ) 

c------------------------- calculate n21 (nrd x nxa)
       nxa = ncomon
      do 120 j = 1,nxa
c        mode=1 : enlarge d^ to fit n21 matrix.
         call enlarg(ifirst,j,ncomon,lone,1)
c        update n21 by column
         do 110 i = 1,nrd
            ir = i+last_nonbias
* MOD TAH use the jelf function
            lm3 = j+jelf(ir)
            a(lm3) = dpdt(i)
            iel21 = (i-1)*ncomon+j
            clc(iel21) = bn22(i)
            if (ibnd.eq.1) goto 110
            ir = ir+nrd
C           lm3 = j+ir*(ir-1)/2
            lm3 = j+jelf(ir)
            a(lm3) = dpdtk(i)
 110     continue
 120  continue    
c*** rwk: By my reckoning, BN22 has never been filled before DIBAS is called.
c         Print out to see
c          write(*,'(a,i5,/,(10e9.2))') 'DBIAS nrd BN22 '
c     .        ,nrd,(bn22(i),i=1,nrd+1)

c---- for multi-session mode, biases are uncorrelated.
      if (mbias0.gt.0) then
         i1 = last_nonbias+1
         i2 = last_nonbias+nrd*iband
         do 144 i = i1,i2
C           lm3 = i*(i-1)/2
            lm3 = jelf(i) 
            call zero1d(lm3+nxa+1,lm3+last_nonbias,a)
 144     continue
      endif
c
c--------------------------------- calculate new u2
c     enlarge d^ to fit old u2
      call enlarg(ifirst,j,ncomon,lone,2)
c     update u2
      call copy1d(1,nrd,last_nonbias,dpdt,b)
      call copy1d(1,nrd,last_nonbias,bn22,blc)  
      ir = last_nonbias+nrd
      if (l2flag.gt.1) call copy1d(1,nrd,ir,dpdtk,b) 
      i1 = last_nonbias+nrd*ibnd+1
      i2 = mbias+nxa
      call zero1d(i1,i2,b)
      call zero1d(i1,i2,blc)
c
c---- calculate new n22
c       compute d^*n22
      do 320 i = 1,nrd
         do 330 j = 1,lone
            dl2 = 0.0d0
            dlc = 0.0d0
            jc = j+ifirst-1
            if (l2flag.gt.1) then
               dl1 = 0.0d0
               dll = 0.0d0
               jc1 = jc+lone
            endif
            do 300 ic = 1,ncd
c              enlarge d^ to fit n22 matrix.
               lm1 = (ic-1)*nrd+i
               ddr = dr(lm1)
               if (dabs(ddr).lt.small) goto 300
               i1 = ipnt2d(ic)
               irow = i1+ifirst-1
               lm2 = jc+jelf(irow)
               le22 = j+jelf(i1)
               if(irow.lt.jc) lm2 = irow+jelf(jc)
               if(irow.lt.jc) le22 = i1+jelf(j)
               dl2 = dl2+ddr*a(lm2)
               dlc = dlc+ddr*an22(le22)
               if (l2flag.le.1) goto 300
               irow = irow+lone
               lm2 = jc+jelf(irow)
               if(irow.lt.jc) lm2 = irow+jelf(jc)
               dl1 = dl1+ddr*a(lm2)
               lm3 = jc1+jelf(irow)
               if(irow.lt.jc1) lm3 = irow+jelf(jc1)
               dll = dll+ddr*a(lm3)
 300        continue
            work2(j) = dl2
            worklc(j) = dlc
            if (l2flag.le.1) goto 330
            work1(j) = dl1
            work(j) = dll
 330     continue
c       compute d^*n22*d^~ (column i)
         do 370 k = 1,i
            lm1 = jelf(i)+k
            dl1 = 0.0d0
            dlc = 0.0d0
            if (l2flag.gt.1) dl2 = 0.0d0
            do 380 ic = 1,ncd
               i1 = ipnt2d(ic)
               l2 = (ic-1)*nrd+k
               ddr = dr(l2)
               if (dabs(ddr).lt.small) goto 380
               dl1 = dl1+work2(i1)*ddr
               dlc = dlc+worklc(i1)*ddr
               if (l2flag.le.1) goto 380
               dl2 = dl2+work(i1)*ddr
 380        continue
            dpdt(lm1) = dl1
            bn22(lm1) = dlc 
            if (l2flag.ge.2) dpdtk(lm1) = dl2
 370     continue
         if (l2flag.le.1) goto 320
         do 372 k = 1,nrd
            lm1 = (i-1)*nrd+k
            dl1 = 0.0d0
            do 374 ic = 1,ncd
               i1 = ipnt2d(ic)
               l2 = (ic-1)*nrd+k
               dl1 = dl1+work1(i1)*dr(l2)
 374        continue
         temp12(lm1) = dl1
 372     continue
 320  continue
c
c     update n22
      do 310 i = 1,nrd
         ir = i+last_nonbias
         ir1 = ir+nrd
         lm1 = jelf(i)
         lm3 = last_nonbias+jelf(ir)
         jc = lm3-lm1  
         call copy1d(lm1+1,lm1+i,jc,dpdt,a) 
         call copy1d(lm1+1,lm1+i,0,bn22,an22)    
         if (l2flag.gt.1) then
            lm3 = last_nonbias+nrd+jelf(ir1)
            jc = lm3-lm1
            call copy1d(lm1+1,lm1+i,jc,dpdtk,a)
            lm3 = last_nonbias+jelf(ir1)
            lm1 = (i-1)*nrd
            jc = lm3-lm1
            call copy1d(lm1+1,lm1+nrd,jc,temp12,a)
         endif
 310  continue      
      i1 = last_nonbias+nrd*ibnd+1
      i2 = mbias+nxa
      ir = jelf(i1)+1
      jc = jelf(i2)
      call zero1d(ir,jc,a)
c
c---- update index of last live bias parameter
      last_bias = last_nonbias+nrd*iband       
c      print *,'DBIAS iband nrd last_nonbias last_bia '
c     .      , iband,nrd,last_nonbias,last_bias

c
      return
      end
c---------------------------------------------------------
      Subroutine ENLARG(ifirst,j,ncomon,lone,mode)
c
c     mode=1 : enlarge D^ to fit n21 matrix.
C     mode=2 : enlarge D^ to fit old U2
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer jelf,ifirst,ncomon,i1,lone,mode,ic,lm1,irow,lm2,iel21,i,j
             

      real*8 small,dl1,dl2,dlc,ddr
c
      small = 1.0d-11
      do 130 i = 1,nrd
         dl1 = 0.0d0
         dl2 = 0.0d0
         dlc = 0.0d0
         do 100 ic = 1,ncd
            lm1 = (ic-1)*nrd+i
            ddr = dr(lm1)
            if (dabs(ddr).lt.small) goto 100
            i1 = ipnt2d(ic)
            irow = i1+ifirst-1
            if (mode.eq.1) then
               lm2 = j+jelf(irow)
               dl1 = dl1+ddr*a(lm2)
               iel21 = (i1-1)*ncomon+j
               dlc = dlc+ddr*clc(iel21)
               if (l2flag.le.1) goto 100
               irow = irow+lone
               lm2 = j+jelf(irow)
               dl2 = dl2+ddr*a(lm2)
            endif
            if (mode.eq.2) then
               dl1 = dl1+ddr*b(irow)
               dlc = dlc+ddr*blc(irow)
               if (l2flag.le.1) goto 100
               irow = irow+lone
               dl2 = dl2+ddr*b(irow)
            endif
 100     continue
         dpdt(i) = dl1
         dpdtk(i) = dl2
         bn22(i) = dlc
 130  continue
c
      return
      end

