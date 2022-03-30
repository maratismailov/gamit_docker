cC*****************************************************************************
C*****                                                                   *****
C*****                             RAND5.                                *****
C*****                                                                   *****
C*****               Pseudo Random Number Generator Package              *****
C*****                            Version 5                              *****
C*****               David M. Krowitz May 24, 1989.                      *****
C*****                                                                   *****
C*****      Copyright (c) 1989                                           *****
C*****      David M. Krowitz                                             *****
C*****      Massachusetts Institute of Technology                        *****
C*****      Department of Earth, Atmosheric, and Planetary Sciences      *****
C*****************************************************************************



C****************************************************************************
C
C           Routine to seed the random number generator using
C           the system's time of day. 
C                         
C****************************************************************************
                           
      Block Data RANDIT
      COMMON /RANDCOM/ INITFLAG,RANTAB,INDEX1,INDEX2,INDEX3
      LOGICAL INITFLAG
      INTEGER*4 RANTAB(0:54)
      INTEGER*4 INDEX1,INDEX2,INDEX3
      DATA INITFLAG /.FALSE./
      End

      SUBROUTINE RAND5_INIT 

      implicit none

      INTEGER*4 I,J,K
c     temporarily commented
c      INTEGER*2 CLOCK(3)

      COMMON /RANDCOM/ INITFLAG,RANTAB,INDEX1,INDEX2,INDEX3
      LOGICAL INITFLAG
      INTEGER*4 RANTAB(0:54)
      INTEGER*4 INDEX1,INDEX2,INDEX3

C
C   Get the 48-bit long time of day from the system and turn it
C   into a 32-bit long integer.
C          
c****TEMP:  Comment out the following for Sun testing; then replace later
c      CALL TIME_$CLOCK (CLOCK)
c      I = INT(CLOCK(1))
c      J = INT(CLOCK(2))*2**4
c      K = INT(CLOCK(3))*2**8  
c**added to provide a (constant) seed if this routine is called
      i=1
      j=0
      k=0
      RANTAB(0) = I+J+K
C
C   Use a linear congruential pseudo random number generator to
C   initialize the table for the RAND5_UNIF function. Note that
C   Knuth's rules for 'good' generators requires a multiplier of
C   the form ending in ...x21 where 'x' is an even digit. The
C   choice of 31415821 is taken from one of his examples.
C
      DO 100 I = 1,54
        RANTAB(I) = RANTAB(I-1)*31415821+1
100   CONTINUE

      INITFLAG = .TRUE.
      INDEX1 = 0
      INDEX2 = 23
      INDEX3 = 54

      RETURN
      END




C****************************************************************************
C
C           Uniform Random Number Generator.
C           Generates a pseudo-random number in the range 0.0 to 1.0 with
C           a uniform distribution. The first time RAND5_UNIF is called the
C           argument is used to seed the random number generator if the
C           random number generator has not already been seeded by a call
C           to the RAND5_INIT routine. If the seed passed on the first call
C           is greater than 0, then the seed is used to initialize the
C           random number generator, otherwise a built in seed number is
C           used to set up the initial table. Since IEEE standard floating
C           point numbers have a 24 bit mantissa (one bit being the hidden bit)
C           we generate random integers in the range 0 to 2**24-1 and then
C           divide them by 2**24-1 to get a real number in the range 0 to
C           1.0, inclusive.
C
C****************************************************************************
          
      REAL*4 FUNCTION RAND5_UNIF (ISEED)
       
      implicit none

      COMMON /RANDCOM/ INITFLAG,RANTAB,INDEX1,INDEX2,INDEX3
      LOGICAL INITFLAG
      INTEGER*4 RANTAB(0:54)
      INTEGER*4 INDEX1,INDEX2,INDEX3,cand
      integer*4 i

      INTEGER*4 ISEED

C
C   If random number generator has not already been seeded
C   by calling the RAND5_INIT routine, then set up the initial
C   table for the additive congruential method used here. If
C   the seed number supplied by the caller is 0, then use our
C   own seed number.
C
      IF (INITFLAG) GOTO 100

      IF (ISEED.GT.0) THEN
        RANTAB(0) = ISEED
      ELSE
        RANTAB(0) = 31415926
      ENDIF

C
C   Use a linear congruential pseudo random number generator to
C   initialize the table for the RAND5_UNIF function. Note that
C   Knuth's rules for 'good' generators requires a multiplier of
C   the form ending in ...x21 where 'x' is an even digit. The
C   choice of 31415821 is taken from one of his examples.
C

      DO 10 I = 1,54
        RANTAB(I) = RANTAB(I-1)*31415821+1
10    CONTINUE

      INITFLAG = .TRUE.
      INDEX1 = 0
      INDEX2 = 23
      INDEX3 = 54

C
C   Generate the next random integer in the range 0 to 2**24-1
C   and then convert it into an IEEE format real number in the
C   range 0.0 to 1.0
C
100   INDEX1 = INDEX1+1
      IF (INDEX1.GE.55) INDEX1 = 0
      INDEX2 = INDEX2+1
      IF (INDEX2.GE.55) INDEX2 = 0
      INDEX3 = INDEX3+1
      IF (INDEX3.GE.55) INDEX3 = 0

      RANTAB(INDEX1) = RANTAB(INDEX2)+RANTAB(INDEX3)

      RAND5_UNIF = FLOAT(cand(RANTAB(INDEX1),2**24-1))/FLOAT(2**24-1)

      RETURN
      END




C****************************************************************************
C
C           Vector Uniform Random Number Generator.
C           Same as RAND5_UNIF function except it fills a complete table
C           full of uniformly distributed pseudo-random numbers in the
C           range 0.0 to 1.0. This routine avoids the overhead of making
C           repeated functions calls when you need a lot of random numbers
C           at once.
C
C****************************************************************************

      SUBROUTINE RAND5_VEC_UNIF (ISEED,TABLE,NUM)

      implicit none

      COMMON /RANDCOM/ INITFLAG,RANTAB,INDEX1,INDEX2,INDEX3
      LOGICAL INITFLAG
      INTEGER*4 RANTAB(0:54)
      INTEGER*4 INDEX1,INDEX2,INDEX3
      integer*4 i

      INTEGER*4 ISEED,NUM,cand
      REAL*4 TABLE(NUM)

C
C   If random number generator has not already been seeded
C   by calling the RAND5_INIT routine, then set up the initial
C   table for the additive congruential method used here. If
C   the seed number supplied by the caller is 0, then use our
C   own seed number.
C
      IF (INITFLAG) GOTO 100

      IF (ISEED.GT.0) THEN
        RANTAB(0) = ISEED
      ELSE
        RANTAB(0) = 31415926
      ENDIF

C
C   Use a linear congruential pseudo random number generator to
C   initialize the table for the RAND5_UNIF function. Note that
C   Knuth's rules for 'good' generators requires a multiplier of
C   the form ending in ...x21 where 'x' is an even digit. The
C   choice of 31415821 is taken from one of his examples.
C

      DO 10 I = 1,54
        RANTAB(I) = RANTAB(I-1)*31415821+1
10    CONTINUE

      INITFLAG = .TRUE.
      INDEX1 = 0
      INDEX2 = 23
      INDEX3 = 54

C
C   Generate the next random integer in the range 0 to 2**24-1
C   and then convert it into an IEEE format real number in the
C   range 0.0 to 1.0
C
100   DO 200 I = 1,NUM
        INDEX1 = INDEX1+1
        IF (INDEX1.GE.55) INDEX1 = 0
        INDEX2 = INDEX2+1
        IF (INDEX2.GE.55) INDEX2 = 0
        INDEX3 = INDEX3+1
        IF (INDEX3.GE.55) INDEX3 = 0

        RANTAB(INDEX1) = RANTAB(INDEX2)+RANTAB(INDEX3)

        TABLE(I) = FLOAT(cand(RANTAB(INDEX1),2**24-1))/FLOAT(2**24-1)

200   CONTINUE

      RETURN
      END




C****************************************************************************
C
C           Decaying Exponential Random Number Generator.
C           Generate random numbers in the range 0.0 to 1.0E+38 (largest
C           positive single precision real number) with a distribution of
C           the form:
C                               1    -X/TAU
C                              ---  E               , TAU > 0
C                              TAU
C
C****************************************************************************

      REAL FUNCTION RAND5_EXP(TAU)
                   
      implicit none
      REAL TAU
      REAL RAND5_UNIF

      RAND5_EXP = -TAU*ALOG(RAND5_UNIF(0))

      RETURN
      END





C****************************************************************************
C
C           Vector Decaying Exponential Random Number Generator.
C           Same as the RAND5_EXP function except that it fills a
C           complete table of pseudo-random numbers with a distribution
C           of a decaying exponential. This routine avoids the overhead of
C           making repeated function calls when you need a lot of values
C           all at once.
C
C****************************************************************************

      SUBROUTINE RAND5_VEC_EXP(TAU,TABLE,NUM)
                   
      implicit none
      REAL TAU
      REAL RAND5_UNIF
      INTEGER*4 NUM,i
      REAL TABLE(NUM)

      DO 100 I = 1,NUM
        TABLE(I) = -TAU*ALOG(RAND5_UNIF(0))
100   CONTINUE

      RETURN
      END





C****************************************************************************
C
C           Lorentzian Random Number Generator.
C           Generate random numbers in the range -1.0E+38 to 1.0E+38 (largest
C           negative single precision real number to largest positive single
C           precision real number) with a distribution of the form:
C
C             1          GAMMA/2
C           ----   ---------------------     , GAMMA > 0 == WIDTH AT HALF MAX
C                        2             2
C            PI    (X-MU)  +  (GAMMA/2)      , MU == MEAN VALUE
C
C****************************************************************************

      REAL FUNCTION RAND5_LOREN (GAMMA,MU)
                   
      implicit none

      REAL PI
      PARAMETER (PI = 3.1415926)

      REAL GAMMA,MU
      REAL RAND5_UNIF
      REAL OLD_GAMMA,OLD_MU
      REAL CONST1,CONST2
      SAVE OLD_GAMMA,OLD_MU
      SAVE CONST1,CONST2

      DATA OLD_GAMMA,OLD_MU/-1.0E38,-1.0E38/

C
C           Avoid overhead of recalculating the ATAN function
C           if GAMMA and MU haven't change since the last call.
C
      IF ((GAMMA.EQ.OLD_GAMMA).AND.(MU.EQ.OLD_MU)) GOTO 100
      CONST1 = GAMMA/2.0
      CONST2 = ATAN(-2.0*MU/GAMMA)
      OLD_MU = MU
      OLD_GAMMA = GAMMA

100   RAND5_LOREN = CONST1*TAN(PI*RAND5_UNIF(0)+CONST2)+MU

      RETURN
      END





C****************************************************************************
C
C           Vector Lorentzian Random Number Generator.
C           Same as the RAND5_LOREN function except that it fills a
C           complete table of pseudo-random numbers with a distribution
C           of a Lorentzian function. This routine avoids the overhead of
C           making repeated function calls when you need a lot of values
C           all at once.
C
C****************************************************************************

      SUBROUTINE RAND5_VEC_LOREN (GAMMA,MU,TABLE,NUM)
  
      implicit none

      REAL PI
      PARAMETER (PI = 3.1415926)

      REAL GAMMA,MU
      INTEGER*4 NUM,i
      REAL TABLE(NUM)
      REAL RAND5_UNIF
      REAL OLD_GAMMA,OLD_MU
      REAL CONST1,CONST2
      SAVE OLD_GAMMA,OLD_MU
      SAVE CONST1,CONST2

      DATA OLD_GAMMA,OLD_MU/-1.0E38,-1.0E38/

C
C           Avoid overhead of recalculating the ATAN function
C           if GAMMA and MU haven't change since the last call.
C
      IF ((GAMMA.EQ.OLD_GAMMA).AND.(MU.EQ.OLD_MU)) GOTO 100
      CONST1 = GAMMA/2.0
      CONST2 = ATAN(-2.0*MU/GAMMA)
      OLD_MU = MU
      OLD_GAMMA = GAMMA

100   DO 200 I = 1,NUM
        TABLE(I) = CONST1*TAN(PI*RAND5_UNIF(0)+CONST2)+MU
200   CONTINUE

      RETURN
      END





C****************************************************************************
C
C           Gaussian Random Number Generator.
C           Generate random numbers in the range -1.0E+38 to 1.0E+38 (largest
C           negative single precision real number to largest positive single
C           precision real number) with a distribution of the form:
C
C                                              2
C                  1             -1  ( X - MU )
C           _______________      --- (--------)     , SIGMA > 0 == STD. DEVIATION
C              _____              2  ( SIGMA  )
C             /                E                    , MU == MEAN VALUE
C           \/ 2 PI   SIGMA
C                         
C****************************************************************************
                          
      REAL FUNCTION RAND5_GAUSS (SIGMA,MU)
                   
      implicit none
        
      LOGICAL INIT        

      integer intv,samples,interval,i
                                       
      real rand5_unif  

      double precision pi,delta,sigma,mu,r,xx
      double precision x,erf,gauss,ierf,dx,y

      PARAMETER (PI = 3.1415926d0)
      PARAMETER (DELTA = 6.0d0)
      PARAMETER (SAMPLES = 200)

      dimension x(samples),erf(samples),intv(samples)


      SAVE INIT           
      SAVE INTV
      SAVE X,ERF

      DATA INIT/.FALSE./
C
C           Define the unit Gaussian function.
C
      GAUSS(y) = 1.0d0/(dSQRT(2.0E0*PI))*DEXP(-1.0d0/2.0d0*y**2)

C
C           If not already initialized, tabulate the integral of the
C           unit gaussian funtion using Simpson's method.
C
      IF (INIT) GOTO 1000
      DX = 2.0E0*DELTA/FLOAT(SAMPLES-1)
      X(1) = -DELTA
      ERF(1) = 0.0E0
      DO 10 I = 2,SAMPLES
        X(I) = X(I-1)+DX
        ERF(I) = ERF(I-1)+(GAUSS(X(I-1))+4.0E0*GAUSS(X(I)-DX/2.0E0)
     &                     +GAUSS(X(I)))*DX/6.0E0
10    CONTINUE

C
C           Set up a table for looking up what interval of
C           the inverse ERF function that a particular value
C           would fall into (since the ERF function was
C           evaluated at equally spaced intervals, its
C           inverse does not have equally space points).
C

      INTERVAL =1
      DO 30 I = 1,SAMPLES
        IERF = (I-1)*(1.0E0/FLOAT(SAMPLES-1))
15      IF (IERF.LT.ERF(INTERVAL)) GOTO 20
        INTERVAL = INTERVAL+1
        GOTO 15
20      INTV(I) = INTERVAL
30    CONTINUE

C
C           Initialization done.
C
      INIT = .TRUE.
    
C
C           Look up the inverse ERF funtion for a given random number
C           between 0.0 and 1.0 and interpolate result.
C

1000  R = RAND5_UNIF(0)
      INTERVAL = INTV(INT(R*(SAMPLES-1)+1))
1010  IF (R.LT.ERF(INTERVAL)) GOTO 1100
      INTERVAL = INTERVAL+1
      GOTO 1010

1100  XX = X(INTERVAL-1)+(X(INTERVAL)-X(INTERVAL-1))/
     &     (ERF(INTERVAL)-ERF(INTERVAL-1))*(R-ERF(INTERVAL-1))

      rand5_gauss = SIGMA*XX+MU

      RETURN
      END





C****************************************************************************
C
C           Vector Gaussian Random Number Generator.
C           Same as the RAND5_GAUSS function except that it fills a
C           complete table of pseudo-random numbers with a distribution
C           of a Gaussian function. This routine avoids the overhead of
C           making repeated function calls when you need a lot of values
C           all at once.
C                         
C****************************************************************************
                          
      SUBROUTINE RAND5_VEC_GAUSS (SIGMA,MU,TABLE,NUM)

      implicit none

      INTEGER SAMPLES,interval,i
      double precision pi,delta,r,xx,gauss,ierf,dx,y
      PARAMETER (PI = 3.1415926d0)
      PARAMETER (DELTA = 6.0d0)
      PARAMETER (SAMPLES = 200)

      REAL SIGMA,MU
      INTEGER*4 NUM
      REAL TABLE(NUM)

      LOGICAL INIT        
      INTEGER*2 INTV(SAMPLES)
      REAL RAND5_UNIF     
      DOUBLE PRECISION X(SAMPLES),ERF(SAMPLES)
      SAVE INIT           
      SAVE INTV
      SAVE X,ERF

      DATA INIT/.FALSE./
C
C           Define the unit Gaussian function.
C
      GAUSS(y) = 1.0E0/(dSQRT(2.0E0*PI))*DEXP(-1.0E0/2.0E0*y**2)

C
C           If not already initialized, tabulate the integral of the
C           unit gaussian funtion using Simpson's method.
C
      IF (INIT) GOTO 1000
      DX = 2.0E0*DELTA/FLOAT(SAMPLES-1)
      X(1) = -DELTA
      ERF(1) = 0.0E0
      DO 10 I = 2,SAMPLES
        X(I) = X(I-1)+DX
        ERF(I) = ERF(I-1)+(GAUSS(X(I-1))+4.0E0*GAUSS(X(I)-DX/2.0E0)
     &                     +GAUSS(X(I)))*DX/6.0E0
10    CONTINUE

C
C           Set up a table for looking up what interval of
C           the inverse ERF function that a particular value
C           would fall into (since the ERF function was
C           evaluated at equally spaced intervals, its
C           inverse does not have equally space points).
C

      INTERVAL =1
      DO 30 I = 1,SAMPLES
        IERF = (I-1)*(1.0E0/FLOAT(SAMPLES-1))
15      IF (IERF.LT.ERF(INTERVAL)) GOTO 20
        INTERVAL = INTERVAL+1
        GOTO 15
20      INTV(I) = INTERVAL
30    CONTINUE

C
C           Initialization done.
C
      INIT = .TRUE.
    
C
C           Look up the inverse ERF funtion for a given random number
C           between 0.0 and 1.0 and interpolate result.
C

1000  DO 2000 I = 1,NUM
        R = RAND5_UNIF(0)
        INTERVAL = INTV(INT(R*(SAMPLES-1)+1))
1010    IF (R.LT.ERF(INTERVAL)) GOTO 1100
        INTERVAL = INTERVAL+1
        GOTO 1010

1100    XX = X(INTERVAL-1)+(X(INTERVAL)-X(INTERVAL-1))/
     &       (ERF(INTERVAL)-ERF(INTERVAL-1))*(R-ERF(INTERVAL-1))

        TABLE(I) = SIGMA*XX+MU
2000  CONTINUE

      RETURN
      END
