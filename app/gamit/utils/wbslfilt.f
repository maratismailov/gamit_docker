      program wblsfilt

c     filter baseline (.bsl) files to compute RMS and other goodies
c     modified to handle baseline correlations: mhm 871201

      implicit none

      integer*4 kk,in,mode,ndof,idoy,nset,nsiz,nrow
     .        , ntotal,maxprm,ngood,iyear
     .        , ios,i1,i2,i,j,k

      real*8 damp,amat,bmat,cond,sqrtcov,bcov,comp,corr
     .     , temp,chisq,redchisq,dev,cov,rms,sum,wdata,chi2
     .     , work1,work2,work3,work4,work5,work6

c     baseline name is 9 letters long eg. VNDN_BLHL
      integer*4 maxbsl
      parameter (maxbsl = 5000,maxprm = 200)
      character*9    bslnam(maxbsl),buf9
      character*80   file1

c     North, East, Up, Length components, with sigmas
      dimension comp(maxbsl,4,2)
      dimension rms(maxbsl,4,2)
      dimension corr(maxbsl,3),sqrtcov(maxbsl,3,3)
      dimension bcov(3,3),work1(3,3),work2(3,3),work3(3)
      dimension amat(maxprm,3),bmat(3,maxprm),work4(maxprm,3)
      dimension work5(3,3),work6(3),wdata(maxprm)
      dimension cov(3,3),temp(4)
      dimension chi2(4)
      real*8    mmat(3),mean(maxbsl,7)
c     mean N,E,U
      dimension sum(maxbsl,4)
c     deviation from mean N,E,U,L
      dimension dev(maxbsl,4)
c
c     Year and day number
      dimension iyear(maxbsl)
      dimension idoy (maxbsl)

      character*4    buf4a1,buf4b1,buf4a2,buf4b2,buf4a,buf4b
      logical        ldrop(maxbsl)
      data damp/5.0d-08/

c     initialize everything
      ntotal = 0
      do i=1,maxbsl
         do j=1,4
            do k=1,2
               comp(i,j,k) = 0.
               rms(i,j,k)  = 0.
            enddo
            sum(i,j)  = 0.
            dev(i,j)  = 0.
         enddo
         do j=1,3
            corr(i,j) = 0.
         enddo
      enddo
      do i=1,4
         chi2(i) = 0.
      enddo

      in = 0
      file1 = 'standard input'
      do i = 1,maxbsl

         read (unit = 5,
     .      fmt     = 2897,
     .      iostat  = ios,
     .      err     = 1010,
     .      end     = 1020) bslnam(i),iyear(i),idoy(i),
     .                     ((comp(i,j,k),k=1,2),j=1,4),(corr(i,j),j=1,3)

         in = in + 1

c        the critical format statement has the same number in LSQX/LSQDO1.FTN
 2897    format (a9,1x,i4,1x,i3,2x,4(2x,f13.4,4x,f8.4,1x),30x,3f10.5)

c        print 2897, bslnam(i),iyear(i),idoy(i),
c     .                 ((comp(i,j,k),k=1,2),j=1,4),(corr(i,j),j=1,3)

      enddo

 1010 continue
      if (ios .ne. 0) then
         write (6,*) 'WBSLFILT: file error'
         call ferror (ios,6)
      endif

 1020 continue

      i1 = 1
      i2 = 0
      nset = 0
      do i = 1,in

c        if a baseline should be excluded from the stats, it does not
c        have a "_' between the two stations
         buf9 = bslnam(i)
         if (buf9(5:5) .ne. '_') then
            ldrop(i) = .true.
            buf9(5:5) = '*'
            bslnam(i) = buf9
         else
            ldrop(i) = .false.
         endif

c        station names of first baseline in set
         buf9 = bslnam(i1)
         buf4a1 = buf9(1:4)
         buf4b1 = buf9(6:9)

c        current station names
         buf9 = bslnam(i+1)
         buf4a2 = buf9(1:4)
         buf4b2 = buf9(6:9)

c        assume sorted input
c        if (bslnam(i+1) .ne. bslnam(i1)) then
         if (.not. (buf4a1 .eq. buf4a2 .and. buf4b1. eq. buf4b2)) then
            i2 = i
            nsiz = i2 - i1 + 1
            nset = nset + 1
            ngood = 0
            nrow = 0

c           count the number of good days
            do j = i1,i2
               if (.not. ldrop(j)) ngood = ngood + 1
            enddo

            if (ngood .gt. 1) then

               do j = i1,i2

                if (.not. ldrop(j)) then
                  nrow = nrow+1

c                 the covariance matrix for the particular baseline
                  do k = 1,3
                  do kk = 1,3
                     if (k.eq.kk) then
c                       No longer need the factor of 3 kf 910109
c                       cov(k,kk) = (comp(j,k,2)/3.d0)**2
                        cov(k,kk) = (comp(j,k,2))**2
                     else
c                       No longer need the factor of 3 kf 910109
                        cov(k,kk) =  corr(j,k+kk-2)*
c    .                   (comp(j,k,2)*comp(j,kk,2)/9.d0)
     .                   (comp(j,k,2)*comp(j,kk,2))
                     endif
                  enddo
                  enddo

c                 find inverse 'square root' of covariance
                  mode = 2
                  call geninv(cov,work5,work1,work2,work3,
     .                 3,3,3,3,3,3,damp,cond,mode)
                  if(cond.gt.500)
     .              print*,'high condition number on covariance:',cond

c                 fill weighted partial matrix
                  do k = 1,3
                  do kk = 1,3
                     sqrtcov(j,k,kk) = work5(k,kk)
                     amat((nrow-1)*3+k,kk) = sqrtcov(j,k,kk)
                  enddo
                  enddo

c                 weight the data
                  do k = 1,3
                     wdata((nrow-1)*3+k) = sqrtcov(j,k,1)*comp(j,1,1) +
     .            sqrtcov(j,k,2)*comp(j,2,1)+sqrtcov(j,k,3)*comp(j,3,1)
                  enddo

                endif
               enddo

c              generalize invert the weighted partial matrix
               mode = 1
               call geninv(amat,bmat,work4,work5,work6,
     .           maxprm,3,3,maxprm,ngood*3,3,damp,cond,mode)

c              multiply result by weighted data to find estimate of baseline vector
               call dotd(bmat,wdata,mmat,3,maxprm,maxprm,1,3,1,
     .                   3,ngood*3,1,1)

c              find unscaled covariance of baseline vector
               call dotd(bmat,bmat,bcov,3,maxprm,3,maxprm,3,3,
     .                   3,ngood*3,3,3)

c              find length and unscaled sigma
               do k = 1,3
                  sum(nset,k) = mmat(k)
               enddo
               sum(nset,4) = mmat(1)**2 + mmat(2)**2 + mmat(3)**2
               sum(nset,4) = dsqrt(sum(nset,4))

c              now calculate the squares of the deviations from the mean
               do k = 1,4
                  temp(k) = 0.d0
               enddo
               do j = i1,i2
                  if (.not. ldrop(j)) then
                     ntotal = ntotal + 1
                     do k=1,4
                        dev(j,k)      = comp(j,k,1) - sum(nset,k)
                        rms(nset,k,1) = rms(nset,k,1)
     .                                   + (dev(j,k)/comp(j,k,2))**2
                        temp(k) = temp(k) + comp(j,k,2)**(-2)

c                       accumulate total normalized variance
c                       for eventual F-test statistics.
                        chi2(k) = chi2(k) + (dev(j,k)/comp(j,k,2))**2
c                       use unnormalized variance kurt 910216
c                       chi2(k) = chi2(k) + (dev(j,k))**2
                     enddo
                  else
c                    calculate the deviation from the mean even if the baseline
c                    is not to be included in the RMS calculations.
                     do k=1,4
                        dev(j,k)      = comp(j,k,1) - sum(nset,k)
                     enddo
                  endif
               enddo

c              now take the square root
               do k=1,4
                  rms(nset,k,1) =
     .                dsqrt(ngood/(ngood-1)*rms(nset,k,1)/temp(k))
               enddo

c              find chi-square
               chisq = 0.d0
               do j = i1,i2
                  if (.not. ldrop(j)) then
                     do k = 1,3
                        work3(k) = dev(j,1)*sqrtcov(j,k,1) +
     .                  dev(j,2)*sqrtcov(j,k,2)+dev(j,3)*sqrtcov(j,k,3)
                     enddo
                     chisq = chisq + work3(1)**2 +
     .                  work3(2)**2 + work3(3)**2
                  endif
               enddo
               redchisq = dsqrt(chisq/((ngood-1)*3))

c              scale diagonals of covariance and find the scaled sigma
c                  of the length:  provides sigmas of the weighted mean
               mean(nset,4) = 0.d0
               do k = 1,3
                  mean(nset,k) = redchisq*dsqrt(bcov(k,k))
                  mean(nset,4) = mean(nset,4) +bcov(k,k)*mmat(k)**2
               enddo
               mean(nset,4) = mean(nset,4) +
     .           2.d0*( bcov(1,2)*mmat(1)*mmat(2) +
     .                  bcov(1,3)*mmat(1)*mmat(3) +
     .                  bcov(2,3)*mmat(2)*mmat(3) )
               mean(nset,4) = redchisq*dsqrt(mean(nset,4))/sum(nset,4)

c              correlations
               mean(nset,5) = bcov(1,2)/dsqrt(bcov(1,1)*bcov(2,2))
               mean(nset,6) = bcov(1,3)/dsqrt(bcov(1,1)*bcov(3,3))
               mean(nset,7) = bcov(2,3)/dsqrt(bcov(2,2)*bcov(3,3))

            else
c              dummy it if we do not have enough observations
               do j = i1,i2
                  do k=1,4
                     sum(nset,k)   = comp(j,k,1)
                     rms(nset,k,1) = comp(j,k,2)
                     mean(nset,k)  = 0.d0
c                    pass on the value instead of its deviation from the mean
                     dev(j,k)   = comp(j,k,1)
                  enddo
                  do k = 5,7
                     mean(nset,k) = 0.d0
                  enddo
               enddo
            endif

c           second element of rms contains ppm calculation as parts per billion
            do k=1,4
               rms(nset,k,2) = dabs(1.0e9*rms(nset,k,1)/sum(nset,4))
            enddo

c           now print everything out
            if (ngood .gt. 1) then
               write (buf4a,22) ngood
 22            format (1x,i3)
               write (buf9,23) buf4a1,buf4b1
 23            format (a4,'_',a4)
               print 2899, buf9,'RMS ',buf4a,
     .                  (sum(nset,j),rms(nset,j,1),j=1,4)
               print 2899, buf9,'PPB ',buf4a,
     .                  (sum(nset,j),rms(nset,j,2),j=1,4)
               print 2899, buf9,'SIG ',buf4a,
     .                  (sum(nset,j),mean(nset,j),j=1,4)
               print 2898, buf9,'CORR',buf4a,
     .                  (mean(nset,j),j=5,7),redchisq
            endif

            do j = i1,i2
               write (buf4a,20) iyear(j)
 20            format (i4.4)
               write (buf4b,21) idoy(j)
 21            format ('.',i3.3)

               if (ngood .eq. 1 .and. .not. ldrop(j)) then
c                 replace underscore with hyphen if only 1 good obs
                  write (buf9,24) buf4a1,buf4b1
 24               format (a4,'-',a4)
               else
                  buf9 = bslnam(j)
               endif

               print 2899, buf9,buf4a,buf4b,
     .                 (dev(j,k),redchisq*comp(j,k,2)/3.d0
     .                  ,k=1,4)
            enddo
            print 30
 30         format ('---------')


            i1 = i2 + 1
         endif

 2899    format (1x,a9,1x,a4,a4,1x,
     .         'N',1x,f13.4,1x,'+-',1x,f8.4,1x,
     .         'E',1x,f13.4,1x,'+-',1x,f8.4,1x,
     .         'U',1x,f13.4,1x,'+-',1x,f8.4,1x,
     .         'L',1x,f13.4,1x,'+-',1x,f8.4,1x)
 2898    format (1x,a9,1x,a4,a4,1x,
     .         'N-W',6x,f6.4,13x,'N-U',6x,f6.4,13x,'W-U',6x,f6.4,13x,
     .         'ST DEV OF UNIT WGT',1x,f8.4)

      enddo

c     print chi2 per degree of freedom
c     degrees of freedom
      ndof = ntotal - nset
      write (buf4a,22) ndof
      print 2899, buf9,'CHI ',buf4a,(0.,chi2(j)/ndof,j=1,4)
c     print 2899, buf9,'CHI ',buf4a,(0.,chi2(j),j=1,4)

      end



      subroutine geninv(a,b,u1,u2,l1,ia,ja,ib,jb,ma,na,
     +  damp,cond,mode)
c
c  Moore-Penrose generalized inverse
c     mode = 1 :  generalized inverse a^
c     mode = 2 :  positive definite 'square root' of the positive
c                 definite symetric matrix a
c     mode = 3 :  The orthogonal projection operator which when pre-
c                 multiplying a matrix m, projects the m-space into the
c                 null space of a (i.e. qa = I - a.a^). The so-called
c                 denuisance operator
c
c  Given:  a(maxna) which dimensioned (iaxja),with singular
c          value cutoff damp
c
c  Returns: b(naxma) dimensioned (ibxjb), condition number=
c           max/min non-zero singular value
c
c  Method:  svd decomposition   a=u(maxna).sv(naxna).v(naxna)
c
c   8/12/86 MHM  Moore-Penrose implemented
c  11/30/88 MHM  Square Root and Denuisance Operators implemented
c
ccccccccccccccccccccccccccccccccccccccccccccc 

      implicit none

      integer i,j,k,ia,ja,ib,jb,na,ma,mode
      real*8 a(ia,ja),b(ib,jb),u1(ia,ja),u2(ja,ja),l1(ja)
      real*8 cond,lmax,lmin,damp,sv1
c
      goto (10,20,30), mode
         print*,' This geninv mode does not exist: ',mode
         stop
c
c  Moore-Penrose
c
10     call svd1(ia,ja,ma,na,a,u1,u2,l1,1)
       call condnum(ia,ma,l1,lmax,lmin,cond)
c
        do 102 j=1,ma
        do 102 i=1,na
102       b(i,j)=0.d0
        do 101 k=1,na
c         print*,l1(k)
             if(l1(k).ge.damp) then
          do 100 j=1,ma
          do 100 i=1,na
            b(i,j)=b(i,j) + u2(i,k)*u1(j,k)/l1(k)
100       continue
            endif
101     continue
      return
c
c  Square Root
c
20       continue
       call svd1(ia,ja,ma,na,a,u1,u2,l1,2)
       call condnum(ia,ma,l1,lmax,lmin,cond)
c
        do 202 j=1,ma
        do 202 i=1,na
202       b(i,j)=0.d0
        do 201 k=1,na
c       print*,l1(k)
             if(l1(k).ge.damp) then
             sv1=dsqrt(l1(k))
          do 200 j=1,ma
          do 200 i=1,na
             b(i,j)=b(i,j) + u1(i,k)*u1(j,k)/sv1
200       continue
             endif
201     continue
      return
c
c  Denuisance Operator
c
30     call svd1(ia,ja,ma,na,a,u1,u2,l1,2)
       call condnum(ia,ma,l1,lmax,lmin,cond)
c
        do 302 j=1,ma
        do 302 i=1,ma
302       b(i,j)=0.d0
        do 301 k=1,na
          do 300 j=1,ma
          do 300 i=1,ma
             b(i,j)=b(i,j) - u1(i,k)*u1(j,k)
300       continue
301     continue
        do 303 i = 1,ma
303         b(i,i) = b(i,i) + 1.d0
      return
c
c
      end
      SUBROUTINE DOTD(A,B,C,MA,NA,MB,NB,MC,NC,L,M,N,IFLAG)
C
C     ALL MATRACIES ARE DOUBLE PRECESSION
C     DIMENSIONS OF ARRAYS IN CALLING PROGRAM ARE :
C     A(MA,NA) , B(MB,NB) , C(MC,NC)
C     IE. IF B IS A VECTOR THEN IT IS DIMENSIONED B(MB,1)
C
C     THE SIZES OF THE ARRAYS TO BE MULTIPLIED ARE:
C     A(L,M)  *  B(M,N)  =  C(L,N)
C
C     IFLAG  =  1  A * B = C
C            =  2  A**TR * B = C
C            =  3  A * B**TR = C
C            =  4  C  -  A * B = C
C            =  5  C  -  A**TR * B = C
C            =  6  C  -  A * B**TR = C
C  
      implicit none
  
      integer*4 ma, na, mb, nb, mc, nc, l, m, n, iflag, i, j, k
      real*8  A(MA,NA) , B(MB,NB) , C(MC,NC)

      if (iflag .le. 3) then
         DO 1 K=1,N
         DO 1 I=1,L
    1    C(I,K)=0.0D0
      end if

      GO TO (10,20,30,40,50,60) , IFLAG
C
C     A * B = C
C
   10 DO 11 K=1,N
      DO 11 J=1,M
      DO 11 I=1,L
   11 C(I,K)=C(I,K) + A(I,J)*B(J,K)
      RETURN
C
C     A**TR * B = C
C
   20 DO 21 K=1,N
      DO 21 I=1,L
      DO 21 J=1,M
   21 C(I,K)=C(I,K) + A(J,I)*B(J,K)
      RETURN
C
C     A * B**TR = C
C
   30 DO 31 J=1,M
      DO 31 K=1,N
      DO 31 I=1,L
   31 C(I,K)=C(I,K) + A(I,J)*B(K,J)
      RETURN

C
C     C = C - A * B
C
   40 DO 41 K=1,N
      DO 41 J=1,M
      DO 41 I=1,L
   41 C(I,K)=C(I,K) - A(I,J)*B(J,K)
      RETURN
C
C     C = C - A**TR * B
C
   50 DO 51 K=1,N
      DO 51 I=1,L
      DO 51 J=1,M
   51 C(I,K)=C(I,K) - A(J,I)*B(J,K)
      RETURN
C
C     C = C - A * B**TR
C
   60 DO 61 J=1,M
      DO 61 K=1,N
      DO 61 I=1,L
   61 C(I,K)=C(I,K) - A(I,J)*B(K,J)
      RETURN

      END
      subroutine condnum(ipar,n,sv,lmax,lmin,cond)
c   Calculates the condition number equaled to the ratio
c of the largest to the smallest non-zero singular value.
c 
      implicit none

      integer*4 ipar,n,i
      double precision lmax,lmin,sv(ipar),cond
c
      lmax=0.d0
      lmin=1.d50
      do 10 i=1,n
       if(sv(i).ge.5.0d-8) then
       lmax=dmax1(lmax,sv(i))
       lmin=dmin1(lmin,sv(i))
      endif
10    continue
      cond=lmax/lmin
      return
      end
      SUBROUTINE SVD1(M1,N1,M,N,A,U,V,Q,INDEX)
C$$$$$ CALLS NO OTHER ROUTINES
C  SINGULAR VALUE DECOMPOSITION)  FOR ALGO PROGRAM SEE WILKINSON+REINSCH
C  HANDBOOK FOR AUTOMATIC COMPUTATION VOL 2 - LINEAR ALGEBRA, PP140-144
C  TRANSLATED FROM ALGOL BY R.L.PARKER
C  THE MATRIX A(M,N) IS DECOMPOSED.  SINGULAR VALUES IN Q, PRE-MATRIX IN U,
C  POST-MATRIX IN V.   INDEX MAY BE 1,2,3 OR 4.  IF 1, FIND U,V. IF 2, FIND
C  ONLY U. IF 3, FIND ONLY V. IF 4, FIND NEITHER. IN ALL CASES, THE ARRAY  U
C  MUST BE SUPPLIED AS IT IS USED AS WORKING SPACE FOR THE ROUTINE.
C  PROGRAM ALTERED BY PAUL SILVER 4/15 TO HANDLE UNPACKED ARRAYS.
C  M1,N1 ARE DIMENSIONS IN MAIN ROUTINE.M,N ARE ACTUAL DIMENSIONS TO
C  BE USED IN THE SUBROUTINE. 

      implicit none

      integer*4 lplus,iback,kback,lback,l1,m1,n1,index
     .        , i,j,k,l,m,n

      real*8 a,c,e,f,g,h,q,s,u,v,x,y,z,eps,tol

      DIMENSION A(M1,N1),U(M1,N1),V(N1,N1),Q(N1)
      DIMENSION E(1000) 
C
      l = 0
C
      EPS=1.0E-10
      TOL=1.0E-35
      DO 1100 I=1,M
      DO 1100 J=1,N
 1100 U(I,J)=A(I,J)
C  HOUSEHOLDER REDUCTION TO BI-DIAGONAL FORM
      G=0.0
      X=0.0
      DO 2900 I=1,N
      E(I)=G
      S=0.0
      L=I+1
      DO 2100 J=I,M
 2100 S=U(J,I)**2 + S
      IF (S .LT. TOL) GO TO 2500
      F=U(I,I)
      G=-DSIGN(SQRT(S),F)
      H=F*G - S
      U(I,I)=F - G
      IF (L.GT.N) GO TO 2501
      DO 2400 J=L,N
      S=0.0
      DO 2200 K=I,M
 2200 S=U(K,I)*U(K,J) + S
      F=S/H
      DO 2300 K=I,M
 2300 U(K,J)=U(K,J) + F*U(K,I)
 2400 CONTINUE
      GO TO 2501
 2500 G=0.0
C
 2501 CONTINUE
      Q(I)=G
      S=0.0
      IF (L.GT.N) GO TO 2601
      DO 2600 J=L,N
 2600 S=U(I,J)**2 + S
 2601 IF (S.LT.TOL) GO TO 2800
      F=U(I,I+1)
      G=-DSIGN(SQRT(S),F)
      H=F*G - S
      U(I,I+1)=F - G
      IF (L.GT.N) GO TO 2651
      DO 2650 J=L,N
 2650 E(J)=U(I,J)/H
 2651 CONTINUE
      IF (L.GT.M) GO TO 2850
      DO 2700 J=L,M
      S=0.0
      IF (L.GT.N) GO TO 2700
      DO 2670 K=L,N
 2670 S=U(J,K)*U(I,K) + S
      DO 2690 K=L,N
 2690 U(J,K)=U(J,K) + S*E(K)
 2700 CONTINUE
      GO TO 2850
 2800 G=0.0
 2850 Y=DABS(Q(I)) + DABS(E(I))
      IF (Y .GT. X) X=Y
 2900 CONTINUE
C
C  ACCUMULATION OF RIGHT-HAND TRANSFORMS (V)
C
      GO TO (3000,3701,3000,3701       ),INDEX
 3000 CONTINUE
      DO 3700 IBACK=1,N
      I=N+1-IBACK
      IF (G .EQ. 0.0) GO TO 3500
      H=U(I,I+1)*G
      IF (L.GT.N) GO TO 3500
      DO 3100 J=L,N
 3100 V(J,I)=U(I,J)/H
      DO 3400 J=L,N
      S=0.0D0
      DO 3200 K=L,N
 3200 S=U(I,K)*V(K,J) + S
      DO 3300 K=L,N
 3300 V(K,J)=V(K,J) + S*V(K,I)
 3400 CONTINUE
 3500 CONTINUE
      IF (L.GT.N) GO TO 3601
      DO 3600 J=L,N
      V(J,I)=0.0D0
 3600 V(I,J)=0.0D0
 3601 V(I,I)=1.0D0
      G=E(I)
      L=I
 3700 CONTINUE
 3701 CONTINUE
C
C  ACCUMULATION OF LEFT-HAND TRANSFORMS
      GO TO (4000,4000,4701,4701       ),INDEX
 4000 CONTINUE
      DO 4700 IBACK=1,N
      I=N+1-IBACK
      L=I+1
      G=Q(I)
      IF (L.GT.N) GO TO 4101
      DO 4100 J=L,N
 4100 U(I,J)=0.0D0
 4101 IF (G.EQ. 0.0D0) GO TO  4500
      H=U(I,I)*G
      IF (L.GT.N) GO TO 4401
      DO 4400 J=L,N
      S=0.0D0
      DO 4200 K=L,M
 4200 S=U(K,I)*U(K,J) + S
      F=S/H
      DO 4300 K=I,M
 4300 U(K,J)=U(K,J) + F*U(K,I)
 4400 CONTINUE
 4401 CONTINUE
      DO 4550 J=I,M
 4550 U(J,I)=U(J,I)/G
      GO TO 4700
 4500 CONTINUE
      DO 4600 J=I,M
 4600 U(J,I)=0.0D0
 4700 U(I,I)=U(I,I) + 1.0D0
C
C  DIAGONALIZATION OF BI-DIAGONAL FORM
 4701 EPS=EPS*X
      DO 9000 KBACK=1,N
      K=N+1-KBACK
C  TEST F-SPLITTING
 5000 CONTINUE
      DO 5100 LBACK=1,K
      L=K+1-LBACK
      IF (DABS(E(L)).LE. EPS) GO TO 6500
      IF (DABS(Q(L-1)) .LE. EPS) GO TO 6000
 5100 CONTINUE
C  CANCELLATION OF E(L), IF L.GT. 1
 6000 C=0.0D0
      S=1.0D0
      L1=L - 1
      DO 6200 I=L,K
      F=S*E(I)
                   E(I)=C*E(I)
      IF (DABS(F) .LE. EPS) GO TO 6500
      G=Q(I)
      Q(I)=DSQRT(F*F + G*G)
      H=Q(I)
      C=G/H
      S=-F/H
      GO TO (6050,6050,6200,6200       ),INDEX
 6050 CONTINUE
      DO 6100 J=1,M
      Y=U(J,L1)
      Z=U(J,I)
      U(J,L1)=Y*C + Z*S
      U(J,I)=-Y*S + Z*C
 6100 CONTINUE
 6200 CONTINUE
C  TEST F-CONVERGENCE
 6500 Z=Q(K)
      IF (L .EQ. K) GO TO  8000
C  SHIFT FROM BOTTOM 2 X 2 MINOR
      X=Q(L)
      Y=Q(K-1)
      G=E(K-1)
      H=E(K)
      F=((Y-Z)*(Y+Z) + (G-H)*(G+H))/(2.0D0*H*Y)
      G=DSQRT(F*F + 1.0D0)
      F=((X-Z)*(X+Z) + H*(Y/(F + DSIGN(G,F))-H))/X
C  NEXT Q-R TRANSFORMATION
      C=1.0D0
      S=1.0D0
      LPLUS=L + 1
      DO 7500 I=LPLUS,K
      G=E(I)
      Y=Q(I)
      H=S*G
      G=C*G
      Z=DSQRT(F*F + H*H)
      E(I-1)=Z
      C=F/Z
      S=H/Z
      F=X*C + G*S
      G=-X*S + G*C
      H=Y*S
      Y=Y*C
      GO TO (7100,7201,7100,7201       ),INDEX
 7100 DO 7200 J=1,N
      X=V(J,I-1)
      Z=V(J,I)
      V(J,I-1)=X*C + Z*S
      V(J,I)=-X*S + Z*C
 7200 CONTINUE
 7201 Z=DSQRT(F*F + H*H)
      Q(I-1)=Z
      C=F/Z
      S=H/Z
      F=C*G + S*Y
      X=-S*G + C*Y
      GO TO (7300,7300,7500,7500       ),INDEX
 7300 DO 7400 J=1,M
      Y=U(J,I-1)
      Z=U(J,I)
      U(J,I-1)=Y*C + Z*S
      U(J,I)=-Y*S + Z*C
 7400 CONTINUE
 7500 CONTINUE
      E(L)=0.0
      E(K)=F
      Q(K)=X
      GO TO  5000
C  CONVERGENCE
 8000 IF (Z .GE. 0.0D0) GO TO 9000
C  Q IS MADE NON-NEGATIVE
      Q(K)=-Z
      GO TO (8100,9000,8100,9000       ),INDEX
 8100 DO 8200 J=1,N
 8200 V(J,K)=-V(J,K)
 9000 CONTINUE
      RETURN
      END                                                               
