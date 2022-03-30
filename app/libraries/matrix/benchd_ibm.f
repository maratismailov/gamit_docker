      PROGRAM BENCHD
                
      Implicit none

      integer*4 mxp
c     integer*4 lda,ldaa,n,info,i,ntimes,lda,i
      integer*4 lda,n,info,i,ntimes

      parameter ( mxp =  4800 )
      
      DOUBLE PRECISION A(mxp,mxp),B(mxp), X(mxp)
      DOUBLE PRECISION TIME(8,6),CRAY,OPS,TOTAL,NORMA,NORMX
      DOUBLE PRECISION RESID,RESIDN,EPS,EPSLON

      double precision rn

*     Put A in common so that we donot get a stack overflow.
      common A
 
      real*4 TSECOND,t1,tm,tm2
      INTEGER IPVT(mxp), offset,ierr, Npassed
      
      character*64 string
      
*     Get the size from the runstring
      print *,' '
      print *,' BENCHD: Benchmark program using matrix factorization'
      print *,' ----------------------------------------------------'
      print *,' '
      
*     See if user has passed dimensions (need to check if getarg
*     count starts at 0 or 1.      
      call getarg(0, string)
      
      if( string(1:1).eq.' ' ) then
         offset = 1
      else
         offset = 0
      end if

*     Set the default size and see if user has passed value      
      N = 500
      
      call getarg(1+offset,string)

      read(string,*,iostat=ierr) Npassed
      if( ierr.eq.0 ) then
         if( Npassed.le.1 .or.Npassed.gt.mxp ) then
	     print *,' BENCHD: Size ',Npassed,
     .               ' not correct'
             print *,' Usage: % benchd <dim>'
  
         else
	     N = Npassed
	 end if
      else
       	 print *,' Usage: % benchd <dim>.  Using defaults'
      end if

*     Set lead dimension to max size      	   
      LDA = mxp
C
      CRAY = .056
c      WRITE(6,1)
c    1 FORMAT(' PLEASE SEND THE RESULTS OF THIS RUN TO:'//
c     $       ' JACK J. DONGARRA'/
c     $       ' MATHEMATICS AND COMPUTER SCIENCE DIVISION'/
c     $       ' ARGONNE NATIONAL LABORATORY'/
c     $       ' ARGONNE, ILLINOIS 60439'//
c     $       ' TELEPHONE: 312-972-7246'//
c     $       ' ARPANET: DONGARRA@ANL-MCS'/)

      rn = n
      OPS = (2.0D0*rN**3)/3.0D0 + 2.0D0*rN**2

C
         CALL MATGEN(A,LDA,N,B,NORMA)
         T1 = TSECOND()
         CALL DGEFA(A,LDA,N,IPVT,INFO)
         TIME(1,1) = TSECOND() - T1
         T1 = TSECOND()
         CALL DGESL(A,LDA,N,IPVT,B,0)
         TIME(1,2) = TSECOND() - T1
         TOTAL = TIME(1,1) + TIME(1,2)
C
C     COMPUTE A RESIDUAL TO VERIFY RESULTS.
C
         DO 10 I = 1,N
            X(I) = B(I)
   10    CONTINUE
         CALL MATGEN(A,LDA,N,B,NORMA)
         DO 20 I = 1,N
            B(I) = -B(I)
   20    CONTINUE
         CALL DMXPY(N,B,N,LDA,X,A)
         RESID = 0.0
         NORMX = 0.0
         DO 30 I = 1,N
            RESID = DMAX1( RESID, DABS(B(I)) )
            NORMX = DMAX1( NORMX, DABS(X(I)) )
   30    CONTINUE
         EPS = EPSLON(1.0D0)
         RESIDN = RESID/( N*NORMA*NORMX*EPS )
         WRITE(6,40)
   40    FORMAT('     NORM. RESID      RESID           MACHEP',
     $          '         X(1)          X(N)')
         WRITE(6,50) RESIDN,RESID,EPS,X(1),X(N)
   50    FORMAT(1P5E16.8)
C
         WRITE(6,60) N
   60    FORMAT(//'    TIMES ARE REPORTED FOR MATRICES OF ORDER ',I6)
         WRITE(6,70)
   70    FORMAT(6X,'DGEFA',6X,'DGESL',6X,'TOTAL',5X,'MFLOPS',7X,'UNIT',
     $         6X,'RATIO')
C
         TIME(1,3) = TOTAL
         TIME(1,4) = OPS/(1.0D6*TOTAL)
         TIME(1,5) = 2.0D0/TIME(1,4)
         TIME(1,6) = TOTAL/CRAY
         WRITE(6,80) LDA
   80    FORMAT(' TIMES FOR ARRAY WITH LEADING DIMENSION OF',I6)
         WRITE(6,110) (TIME(1,I),I=1,6)
C
         CALL MATGEN(A,LDA,N,B,NORMA)
         T1 = TSECOND()
         CALL DGEFA(A,LDA,N,IPVT,INFO)
         TIME(2,1) = TSECOND() - T1
         T1 = TSECOND()
         CALL DGESL(A,LDA,N,IPVT,B,0)
         TIME(2,2) = TSECOND() - T1
         TOTAL = TIME(2,1) + TIME(2,2)
         TIME(2,3) = TOTAL
         TIME(2,4) = OPS/(1.0D6*TOTAL)
         TIME(2,5) = 2.0D0/TIME(2,4)
         TIME(2,6) = TOTAL/CRAY
C
         CALL MATGEN(A,LDA,N,B,NORMA)
         T1 = TSECOND()
         CALL DGEFA(A,LDA,N,IPVT,INFO)
         TIME(3,1) = TSECOND() - T1
         T1 = TSECOND()
         CALL DGESL(A,LDA,N,IPVT,B,0)
         TIME(3,2) = TSECOND() - T1
         TOTAL = TIME(3,1) + TIME(3,2)
         TIME(3,3) = TOTAL
         TIME(3,4) = OPS/(1.0D6*TOTAL)
         TIME(3,5) = 2.0D0/TIME(3,4)
         TIME(3,6) = TOTAL/CRAY
C
         NTIMES = 10
         TM2 = 0
         T1 = TSECOND()
         DO 90 I = 1,NTIMES
            TM = TSECOND()
            CALL MATGEN(A,LDA,N,B,NORMA)
            TM2 = TM2 + TSECOND() - TM
            CALL DGEFA(A,LDA,N,IPVT,INFO)
   90    CONTINUE
         TIME(4,1) = (TSECOND() - T1 - TM2)/NTIMES
         T1 = TSECOND()
         DO 100 I = 1,NTIMES
            CALL DGESL(A,LDA,N,IPVT,B,0)
  100    CONTINUE
         TIME(4,2) = (TSECOND() - T1)/NTIMES
         TOTAL = TIME(4,1) + TIME(4,2)
         TIME(4,3) = TOTAL
         TIME(4,4) = OPS/(1.0D6*TOTAL)
         TIME(4,5) = 2.0D0/TIME(4,4)
         TIME(4,6) = TOTAL/CRAY
C
         WRITE(6,110) (TIME(2,I),I=1,6)
         WRITE(6,110) (TIME(3,I),I=1,6)
         WRITE(6,110) (TIME(4,I),I=1,6)
C 110    FORMAT(6(1PE11.3))
  110    format(6(f10.4,1x))
C
      STOP
      END
      SUBROUTINE MATGEN(A,LDA,N,B,NORMA)

      implicit none
     
      integer i,j,n,lda,init       
      DOUBLE PRECISION A(LDA,1),B(1),NORMA
 
c
C	the above two statements were modified by HBPapo 1999Dec9
c
      INIT = 1325
      NORMA = 0.0
      DO 30 J = 1,N
         DO 20 I = 1,N
            INIT = MOD(3125*INIT,65536)
            A(I,J) = (INIT - 32768.0)/16384.0
            NORMA = DMAX1(A(I,J), NORMA)
   20    CONTINUE
   30 CONTINUE
      DO 35 I = 1,N
          B(I) = 0.0
   35 CONTINUE
      DO 50 J = 1,N
         DO 40 I = 1,N
            B(I) = B(I) + A(I,J)
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO,INCX,INCY
      DOUBLE PRECISION A(LDA,1)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
	INCX=1
	INCY=1
               CALL DAXPY(N-K,T,A(K+1,K),INCX,A(K+1,J),INCY)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(1),JOB,INCX,INCY
      DOUBLE PRECISION A(LDA,1),B(1)
C
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
C        OR DGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
	INCX=1
	INCY=1
            CALL DAXPY(N-K,T,A(K+1,K),INCX,B(K+1),INCY)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
	INCX=1
	INCY=1
            CALL DAXPY(K-1,T,A(1,K),INCX,B(1),INCY)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
  
      implicit none

c     integer n
C	the above statement was blocked HBPapo 1999Dec9
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 CONTINUE
      DO 30 I = 1,N
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)

      implicit none

C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 CONTINUE
      DO 30 I = 1,N
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      DDOT = DTEMP
      RETURN
      END
c
      SUBROUTINE  DSCAL(N,DA,DX,INCX)

      implicit none

c     integer nincx
C
C     SCALES A VECTOR BY A CONSTANT.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 CONTINUE
      DO 30 I = 1,N
        DX(I) = DA*DX(I)
   30 CONTINUE
      RETURN
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. DABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
C
      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION EPSLON (X)
      DOUBLE PRECISION X
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
      DOUBLE PRECISION A,B,C,EPS
C
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C        1.  THE BASE USED IN REPRESENTING DFLOATING POINT
C            NUMBERS IS NOT A POWER OF THREE.
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C            THE ACCURACY USED IN DFLOATING POINT VARIABLES
C            THAT ARE STORED IN MEMORY.
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C     ASSUMPTION 2.
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C            C  IS NOT EXACTLY EQUAL TO ONE,
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM
C                 THE NEXT LARGER DFLOATING POINT NUMBER.
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
C
C     *****************************************************************
C     THIS ROUTINE IS ONE OF THE AUXILIARY ROUTINES USED BY EISPACK III
C     TO AVOID MACHINE DEPENDENCIES.
C     *****************************************************************
C
C     THIS VERSION DATED 4/6/83.
C
      A = 4.0D0/3.0D0
   10 B = A - 1.0D0
      C = B + B + B
      EPS = DABS(C-1.0D0)
      IF (EPS .EQ. 0.0D0) GO TO 10
      EPSLON = EPS*DABS(X)
      RETURN
      END 

      SUBROUTINE MM (A, LDA, N1, N3, B, LDB, N2, C, LDC)

      implicit none      

      integer lda,n1,n3,ldb,n2,ldc,i,j

      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*)
C
C   PURPOSE:
C     MULTIPLY MATRIX B TIMES MATRIX C AND STORE THE RESULT IN MATRIX A.
C
C   PARAMETERS:
C
C     A DOUBLE PRECISION(LDA,N3), MATRIX OF N1 ROWS AND N3 COLUMNS
C
C     LDA INTEGER, LEADING DIMENSION OF ARRAY A
C
C     N1 INTEGER, NUMBER OF ROWS IN MATRICES A AND B
C
C     N3 INTEGER, NUMBER OF COLUMNS IN MATRICES A AND C
C
C     B DOUBLE PRECISION(LDB,N2), MATRIX OF N1 ROWS AND N2 COLUMNS
C
C     LDB INTEGER, LEADING DIMENSION OF ARRAY B
C
C     N2 INTEGER, NUMBER OF COLUMNS IN MATRIX B, AND NUMBER OF ROWS IN
C         MATRIX C
C
C     C DOUBLE PRECISION(LDC,N3), MATRIX OF N2 ROWS AND N3 COLUMNS
C
C     LDC INTEGER, LEADING DIMENSION OF ARRAY C
C
C ----------------------------------------------------------------------
C
      DO 20 J = 1, N3
         DO 10 I = 1, N1
            A(I,J) = 0.0
   10    CONTINUE
         CALL DMXPY (N2,A(1,J),N1,LDB,C(1,J),B)
   20 CONTINUE
C
      RETURN
      END                             
c
      SUBROUTINE DMXPY (N1, Y, N2, LDM, X, M)

      implicit none
      
      integer n1,n2,ldm,i,j,jmin

      DOUBLE PRECISION Y(*), X(*), M(LDM,*)
C
C   PURPOSE:
C     MULTIPLY MATRIX M TIMES VECTOR X AND ADD THE RESULT TO VECTOR Y.
C
C   PARAMETERS:
C
C     N1 INTEGER, NUMBER OF ELEMENTS IN VECTOR Y, AND NUMBER OF ROWS IN
C         MATRIX M
C
C     Y DOUBLE PRECISION(N1), VECTOR OF LENGTH N1 TO WHICH IS ADDED
C         THE PRODUCT M*X
C
C     N2 INTEGER, NUMBER OF ELEMENTS IN VECTOR X, AND NUMBER OF COLUMNS
C         IN MATRIX M
C
C     LDM INTEGER, LEADING DIMENSION OF ARRAY M
C
C     X DOUBLE PRECISION(N2), VECTOR OF LENGTH N2
C
C     M DOUBLE PRECISION(LDM,N2), MATRIX OF N1 ROWS AND N2 COLUMNS
C
C ----------------------------------------------------------------------
C
C   CLEANUP ODD VECTOR
C
      J = MOD(N2,2)
      IF (J .GE. 1) THEN
         DO 10 I = 1, N1
            Y(I) = (Y(I)) + X(J)*M(I,J)
   10    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF TWO VECTORS
C
      J = MOD(N2,4)
      IF (J .GE. 2) THEN
         DO 20 I = 1, N1
            Y(I) = ( (Y(I))
     $             + X(J-1)*M(I,J-1)) + X(J)*M(I,J)
   20    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF FOUR VECTORS
C
      J = MOD(N2,8)
      IF (J .GE. 4) THEN
         DO 30 I = 1, N1
            Y(I) = ((( (Y(I))
     $             + X(J-3)*M(I,J-3)) + X(J-2)*M(I,J-2))
     $             + X(J-1)*M(I,J-1)) + X(J)  *M(I,J)
   30    CONTINUE
      ENDIF
C
C   CLEANUP ODD GROUP OF EIGHT VECTORS
C
      J = MOD(N2,16)
      IF (J .GE. 8) THEN
         DO 40 I = 1, N1
            Y(I) = ((((((( (Y(I))
     $             + X(J-7)*M(I,J-7)) + X(J-6)*M(I,J-6))
     $             + X(J-5)*M(I,J-5)) + X(J-4)*M(I,J-4))
     $             + X(J-3)*M(I,J-3)) + X(J-2)*M(I,J-2))
     $             + X(J-1)*M(I,J-1)) + X(J)  *M(I,J)
   40    CONTINUE
      ENDIF
C
C   MAIN LOOP - GROUPS OF SIXTEEN VECTORS
C
      JMIN = J+16
      DO 60 J = JMIN, N2, 16
         DO 50 I = 1, N1
            Y(I) = ((((((((((((((( (Y(I))
     $             + X(J-15)*M(I,J-15)) + X(J-14)*M(I,J-14))
     $             + X(J-13)*M(I,J-13)) + X(J-12)*M(I,J-12))
     $             + X(J-11)*M(I,J-11)) + X(J-10)*M(I,J-10))
     $             + X(J- 9)*M(I,J- 9)) + X(J- 8)*M(I,J- 8))
     $             + X(J- 7)*M(I,J- 7)) + X(J- 6)*M(I,J- 6))
     $             + X(J- 5)*M(I,J- 5)) + X(J- 4)*M(I,J- 4))
     $             + X(J- 3)*M(I,J- 3)) + X(J- 2)*M(I,J- 2))
     $             + X(J- 1)*M(I,J- 1)) + X(J)   *M(I,J)
   50    CONTINUE
   60 CONTINUE
      RETURN
      END     

      REAL FUNCTION TSECOND()

c      INTEGER IT(4) 
      REAL eTIME, TARG(2)
C
c     CALL gettimeofday(IT,IT(3))
c     IT(1)=MOD(IT(1),86400)
c     TSECOND=IT(1)*1.0D0+IT(2)*1.0D-6 
c     TSECOND = secnds(v1)
C     dummy = dtime(targ)
      TSECOND = etime_(targ)
      RETURN
      END
