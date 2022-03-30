C     SUBROUTINE cholsk (A,B,IJOB,N,IER)
C
C
C   FUNCTION            - IN-PLACE MATRIX INVERSION AND LINEAR
C                           EQUATION SOLUTION - POSITIVE DEFINITE
C                           MATRIX - SYMMETRIC STORAGE MODE.
C   PARAMETERS   A      - INPUT/OUTPUT VECTOR OF LENGTH N(N+1)/2. ON
C                           INPUT, A CONTAINS THE N BY N POSITIVE
C                           DEFINITE SYMMETRIC MATRIX STORED IN
C                           SYMMETRIC STORAGE MODE. FOR OUTPUT VECTOR A,
C                           SEE PARAMETER IJOB.
C                B      - INPUT/OUTPUT VECTOR OF LENGTH N WHEN IJOB =
C                           2 OR 3. ON INPUT, B CONTAINS THE RIGHT-HAND
C                           SIDE OF THE EQUATION AX = B. ON OUTPUT THE
C                           SOLUTION X REPLACES B WHEN IJOB = 2 OR 3.
C                           OTHERWISE, B IS NOT USED.
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES WHEN
C                           I = 1, INVERT MATRIX A. A IS REPLACED BY ITS
C                             INVERSE STORED IN SYMMETRIC STORAGE MODE.
C                           I = 2, SOLVE THE EQUATION AX = B. A IS
C                             REPLACED BY THE DECOMPOSED MATRIX L SUCH
C                             THAT A = L*L-TRANSPOSE. L IS STORED IN
C                             SYMMETRIC STORAGE MODE. THE DIAGONAL OF L
C                             CONTAINS THE RECIPROCALS OF THE ACTUAL
C                             DIAGONAL ELEMENTS OF L.
C                           I = 3, SOLVE AX = B AND INVERT MATRIX A.
C                             A IS REPLACED BY ITS INVERSE AND THE SOLU-
C                             TION X REPLACES B.
C                           I = 4, CHECK THE CHARACTER OF MATRIX A ONLY
C                N      - ORDER OF A. (INPUT)
C                IER    - ERROR PARAMETER.
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT IJOB WAS LESS THAN
C                             1 OR GREATER THAN 3.
C                           IER = 130 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE
C
C    NOTE : CANNOT BE USED FOR INNER CONSTRAINT SOLUTIONS
C           WHEN CONSTRAINTS ARE ADDED AS A BLOCK - MUST BE
C           COLLAPSED INTO NORMAL EQUATIONS BEFORE INVERSION
C
C-----------------------------------------------------------------------
C
      SUBROUTINE cholsk (A,B,IJOB,N,IER)
C
      DOUBLE PRECISION   A(n*(n+1)/2),B(n),ONE,Q,RN,S,SIXTN,DSQRT,T,X
      real*8 zero,eps
      integer ijob,n,ier,ip,i,j,iq,ir,isw,nm1,ii,jj,ip1,jm1
      integer jl,ji,li,kk,k,lj,is,im1,l,iw,n1,numberi,numberj
      DATA               ZERO/0.0D0/,ONE/1.0D0/,SIXTN/16.0D0/   
c.....Machine precision = 2.5D-14, cumulated error -- D-10
      DATA               EPS/2.5D-14/
c     print*, 'Machine precision =0.1D-14'
c     DATA               EPS/0.1D-14/
      IF (IJOB.LT.1.OR.IJOB.GT.4) GO TO 115
      ISW = 0
C                                  CHOLESKY DECOMPOSITION OF A
C                                     A = L*L-TRANSPOSE
c     print *,'n',n
c     print *,'A=', A
c     print *,'B=', B
c     print *,'ijob= ',ijob
c     print *,' in CHOLSK   -----------  ********'
      RN = ONE/(N*SIXTN)
      IP = 1
      IER = 0
      DO 30 I=1,N
         IQ = IP
         IR = 1
         DO 25 J=1,I
            X = A(IP)
            IF (J.EQ.1) GO TO 10
            DO 5 K=IQ,IP1
               X = X-A(K)*A(IR)
               IR = IR+1
    5       CONTINUE
   10       IF (I.NE.J) GO TO 15
            Q = A(IP)+X*RN
c            IF (Q-A(IP).LE.EPS) GO TO 120
            IF (x.LE.EPS) then
              print*, 'x= ',x
              print*, 'ip= ',ip
              GO TO 120
            endif
      A(IP)=ONE/DSQRT(X)
      GO TO 20
   15       A(IP) = X*A(IR)
   20       IP1 = IP
            IP = IP+1
            IR = IR+1
   25    CONTINUE
   30 CONTINUE    
      IF (IJOB.EQ.4) GOTO  9005
      IF (ISW.EQ.0.AND.IJOB.NE.1) GO TO 75
C                                  FORM A-INVERSE
   35 NM1 = N-1
      IF (N.EQ.1) GO TO 55
C                                  COMPUTE L-INVERSE FIRST
      II = 1
      DO 50 I=1,NM1
         IP1 = I+1
         JM1 = I
         JJ = II
         DO 45 J=IP1,N
            S = ZERO
            LI = II
            JI = JJ+I
            JL = JI
            DO 40 L=I,JM1
               S = S+A(LI)*A(JL)
               JL = JL+1
               LI = LI+L
   40       CONTINUE
            JJ = JJ+J
            A(JI) = -A(JJ)*S
            JM1 = J
   45    CONTINUE
         II = II+IP1
   50 CONTINUE
   55 II = 0
C                                  NOW FORM A-INVERSE
      DO 70 I=1,N
         JJ = II
         DO 65 J=I,N
            S = ZERO
            JI = JJ+I
            LI = JI
            JJ = JJ+J
            LJ = JJ
            DO 60 L=J,N
               S = S+A(LI)*A(LJ)
               LI = LI+L
               LJ = LJ+L
   60       CONTINUE
            A(JI) = S
   65    CONTINUE
         II = II+I
   70 CONTINUE
      IF (IJOB.EQ.1.OR.ISW.EQ.1) GO TO 9005
C                                  SOLVE AX = B
   75 ISW = 1
      IP = 1
      IW = 0
C                                  SOLUTION OF LY = B
      IM1 = 0
      DO 95 I=1,N
         T = B(I)
         IF (IW.EQ.0) GO TO 85
         IP = IP+IW-1
         DO 80 K=IW,IM1
            T = T-A(IP)*B(K)
            IP = IP+1
   80    CONTINUE
         GO TO 90
   85    IF (T.NE.ZERO) IW = I
         IP = IP+IM1
   90    B(I) = T*A(IP)
         IP = IP+1
         IM1 = I
   95 CONTINUE
C                                  SOLUTION OF UX = Y
C                                  WHERE U = L-TRANSPOSE
      N1 = N+1
      DO 110 I=1,N
         II = N1-I
         IP = IP-1
         IS = IP
         IQ = II+1
         T = B(II)
         IF (N.LT.IQ) GO TO 105
         KK = N
         DO 100 K=IQ,N
            T = T-A(IS)*B(KK)
            KK = KK-1
            IS = IS-KK
  100    CONTINUE
  105    B(II) = T*A(IS)
  110 CONTINUE
      IF (IJOB.EQ.3) GO TO 35
      GO TO 9005
C                                  IJOB OUT OF RANGE
  115 IER = 129
      WRITE(7,200)
  200 FORMAT(//,5X,'IJOB INPUT PARAMETER OUT OF RANGE')
      GO TO 9005  
  120 IER=130
      IF (IJOB.EQ.4) GOTO 9005
      WRITE(*,300) I,J,x
  300 FORMAT(//,5X,'Matrix is algorithmically not positive definite',
     *   /, 8X,'at row',i4,' and column',i4,' X=',e12.6)
      numberi=aint(float(I)/6.) 
      numberj=aint(float(J)/6.) 
      print*, 'isit= ',numberi+1,'xyz dot =',
     .         aint( ( float(I)/6.-numberi) *6 )
      print*, 'jsit= ',numberj+1,'xyz dot =',
     .         aint( ( float(J)/6.-numberj) *6 )
 9005 RETURN
      END
