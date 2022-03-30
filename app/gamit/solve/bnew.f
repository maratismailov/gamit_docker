c
      Subroutine BNEW( nlive0,ifix )
c
c---- sequentially update solution and matrix.
c---- it is more efficient.
c----                                 -DND-    880407
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      character*256 message

      integer nb12,nb22,ij22,ij11
     .      , nlive0,ifix,ib,ii,i1,ier,j1,k1,l1,i,j,k,l

      real*8 temp(maxprm),bdev0(maxbis),tt,tw,co,rcond

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
C--------------------------------------------------------------
C---- Reorder array A put all fixed bias rows at the bottom.
c
      if (ifix.eq.0) then
         write (6,'(1x,''No new bias has been fixed.'')')
         goto 200
      endif
c      write(6,'(5x,''BNEW is working ...... (n.rms)='',f12.5)') sclerr
         
      nb12 = 0
      do 100 i = 1,ifix
         ib = nfix(i)-i+1
         bdev0(i) = bdev(i)
c
c----    reorder inversed normal matrix.
         if (ib.ge.nlive0) goto 100
         nb12 = ib*(ib-1)/2
         nb22 = nb12+ib
         call copy1d(nb12+1,nb22,-nb12,a,temp)
         temp(nlive0) = temp(ib)
         do 30 j = ib+1,nlive0
            nb22 = j*(j-1)/2
            ij11 = nb22+ib
            temp(j-1) = a(ij11)
            do 40 k = 1,j
               if (k.eq.ib) goto 40
               k1 = k
               if (k.gt.ib) k1 = k-1
               ij11 = nb12+k1
               ij22 = nb22+k
               a(ij11) = a(ij22)
 40         continue
            nb12 = nb22
 30      continue
         call copy1d(1,nlive0,nb22,temp,a)
 100  continue
c
c---- calculate invers(q22) and invers(q22)*db
      if (ifix.eq.1) then
         nb12 = (nlive+1)*nlive/2
         ij11 = nb12+nlive+1
         dpdt(1) = 1.0d0/a(ij11)
         bdev(1) = dpdt(1)*bdev(1)
         goto 110
      endif
      ii = 0
      do 60 i = 1,ifix
         i1 = nlive+i
         nb12 = i1*(i1-1)/2
         nb22 = nb12+nlive
         do 70 j = 1,i
            ii = ii+1
            ij11 = nb22+j
            dpdt(ii) = a(ij11)
 70      continue
 60   continue
c      call invers(dpdt,bdev,3,ifix,ier)
      call inver2(dpdt,bdev,3,ifix,rcond,ier ) 
      if( ier.ne.0 ) then
        write(message,'(a,d12.5,a)') 
     .    'Non-zero error from inver2, rcond=',rcond,' ier='   
        call report_stat('WARNING','SOLVE','bnew',' ',message,ier)
      endif

c---- calculate q12*invers(q22)*db
 110  ii = 0
      do 80 i = 1,ntpart
         if (free(i).eq.0) goto 80
         ii = ii+1
         if (ifix.eq.1) then
            ij11 = nb12+ii
            adjust(i) = adjust(i)+a(ij11)*bdev(1)
            goto 80
         endif
         tw = 0.0d0
         do 90 j = 1,ifix
            j1 = nlive+j
            ij11 = j1*(j1-1)/2+ii
            tw = tw+a(ij11)*bdev(j)
 90      continue
         adjust(i) = adjust(i)+tw
 80   continue
c---- calculate invers(n) and q12*invers(q22)*q21
      do 140 j = 1,nlive
         do 150 i = j,nlive
            nb22 = i*(i-1)/2+j
            if (ifix.eq.1) then
               ij11 = nb12+j
               temp(1) = dpdt(1)*a(ij11)
               ij22 = nb12+i
               a(nb22) = a(nb22)-a(ij22)*temp(1)
               goto 150
            endif
            tt = 0.0d0
            do 130 k = 1,ifix
               tw = 0.0d0
               do 120 l = 1,ifix
                  l1 = nlive+l
                  ij11 = l1*(l1-1)/2+j
                  if (k-l) 124,122,122
 122              i1 = k*(k-1)/2+l
                  goto 126
 124              i1 = l*(l-1)/2+k
 126              continue
                  if( i1.gt.maxwm2.or.ij11.gt.maxnrm ) then
                     write(message,'(a,4i4)') 'i1 maxwm2 ij11 maxnrm '
     .                                      ,  i1,maxwm2,ij11,maxnrm
                     call report_stat('WARNING','SOLVE','bnew',' '
     .                               , message,0)
                  endif
                  tw = tw+dpdt(i1)*a(ij11)
 120           continue
               temp(k) = tw
 130        continue
            do 160 k = 1,ifix
               i1 = nlive+k
               ij22 = i1*(i1-1)/2+i
               tt = tt+a(ij22)*temp(k)
 160        continue
            a(nb22) = a(nb22)-tt
 150     continue
 140  continue
c---- calculate new chi-square.
      tt = 0.0d0
      do 170 i = 1,ifix
         tt = tt+bdev0(i)*bdev(i)
 170  continue
      co = sclerr**2*dble(nobs-nlive0)+tt
      tt = dabs(co)/dble(nobs-nlive)
      sclerr = dsqrt(tt)

 200  continue
c
      return
      end
