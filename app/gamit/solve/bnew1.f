Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      Subroutine BNEW1(ilive0,ilive,ifast)
c
C---- Calculate updated CHI-square only.
c----                                 -DND-    880408
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer nb12,nb22,ilive,ilive0,ifast,ifix,ib,j1,ier,ij,i,j
      real*8 bdev0(20),temp(maxprm),tt,rcond
c
C---- Algorithm introduction :
C     The nomal equation ---
C        | N11 N12 | | X(old) |   | U1 |
C        |         | |        | = |    |                 (1)
C        | N21 N22 | |   B*   |   | U2 |
C     Where B* is a subset of bias parameters which will be fixed.
C           X(old) includes all live parameters except B*.
C              Notice : X(old) includes bias parameters which will
C                       keep free.
C     Before bias-fixing,
C        | X(old) |   | Q11 Q12 | | U1 |
C        |        | = |         | |    |                 (2)
C        |   B*   |   | Q21 Q22 | | U2 |
C     Where
C        | Q11 Q12 |   | N11 N12 |^
C        |         | = |         |                       (3)
C        | Q21 Q22 |   | N21 N22 |
C       symbol ^ means inverse.
C     According to the matrix partition theory,
C        Q11 = N11^ + N11^ N12 Q22 N21 N11^
C        Q12 = -N11^ N12 Q22
C        Q21 = -Q22 N21 N11^                             (4)
C        Q22 = ( N22 - N21 N11^ N12 )^
C     After bias-fixing, all B* were fixed as B(fix). Meanwhile, X(old)
C     were updated to X(new).
C     Now, the normal equation becomes :
C        | N11 N12 | | X(new) |   | U1 |
C        |         | |        | = |    |                 (5)
C        | N21 N22 | | B(fix) |   | U2 |
C     and
C        X(new) = N11^ U1 - N11^ N12 B(fix)              (6)
C     It is easy to get the expressions from (4)
C        N11^ = Q11 - Q12 Q22^ Q21
C        -N11^ N12 = Q12 Q22^                            (7)
C     Put (7) into (6)
C        X(new) = Q11 U1 - Q12 Q22^ ( Q21 U1 - B(fix))   (8)
C     From (2), we get another relation :
C        X(old) = Q11 U1 - Q12 Q22^ ( Q21 U1 - B* )      (9)
C     Comparing (8) and (9)
C        X(new) = X(old) + Q12 Q22^ ( B(fix) - B* )      (10)
C     The updated matrix is :
C        Q(new) = N11^ = Q11 - Q12 Q22^ Q21              (11)
C     The Q11, Q12, Q22 were already calculated from pre-bias-
C     fixing solution and stored in array A.  In old approach, we
C     have to inverse N11 . Now, we need inverse only Q22 to get
C     X(new) and Q(new).
C     Remember the dimension of N11 >> dimention of Q22.  We
C     never worry about the computational nightmare of
C     bias-fixing !
C     The old CHI2
C                            | U1 |~ | X(old) |
C        CHI2(old) = R2sum - |    |  |        |          (12)
C                            | U2 |  |  B*    |
C     Where symbol ~ means transpose.
C     After some manipulations, we get the formula :
C
C        CHI2(new) = CHI2(old) + dB~ Q22^ dB             (13)
C     Where
C        dB = B(fix) - B*
c--------------------------------------------------------------
      ifix = ilive0-ilive
      if (ifix.eq.0) then
         write (6,'(1x,''No new bias has been fixed.'')')
         goto 200
      endif

c      write (6,'(5x,''BNEW1 is working ......'')')
c---- ifast = 1 : using stored inverse(Q22), calculate chi2 directly.
      if (ifast.eq.1) then
c---- calculate inverse(Q22)*dB
         do 40 i = 1,ifix
            temp(i) = 0.0d0
            do 50 j = 1,ifix  
               if( i.lt.j ) then
                 ij = j*(j-1)/2+i  
               else
                 ij = i*(i-1)/2+j 
               endif
               temp(i) = temp(i)+dpdt(ij)*bdev(j)
 50         continue
            bdev0(i) = temp(i)
 40      continue
         goto 110
      endif
      do 100 i = 1,ifix
         bdev0(i) = bdev(i)
         ib = nfix(i)
         nb12 = ib*(ib-1)/2
c
c----    form Q22 submatrix.
         do 20 j = 1,i
            ib = nfix(j)
            nb22 = nb12+ib
            j1 = i*(i-1)/2+j
            dpdt(j1) = a(nb22)
 20      continue
 100  continue
c---- calculate invers(Q22) and invers(Q22)*dB
      if (ifix.eq.1) then
         dpdt(1) = 1.0d0/dpdt(1)
         bdev0(1) = dpdt(1)*bdev0(1)
         goto 110
      endif
c      call invers(dpdt,bdev0,3,ifix,ier)
      call inver2(dpdt,bdev0,3,ifix,rcond,ier)
c---- calculate dB~*invers(Q22)*dB
 110  continue
      tt = 0.0d0
      do 170 i = 1,ifix
         tt = tt+bdev(i)*bdev0(i)
 170  continue
c---- calculate new CHI-square.
      chi2 = sclerr**2*dble(nobs-ilive0)+tt

 200  continue
c
      return
      end
