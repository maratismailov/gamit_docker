Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      FUNCTION FUNCOF(T,A,INTB,B,C,D)
C
C       M.E.ASH  SEPT 1967   FUNCTION SUBROUTINE FUNCOF
C       THIRD ORDER POLYNOMIALS FOR FUNDAMENTAL ARGUMENTS OF BROWN LUNAR
C       THEORY ARE EVALUATED

      implicit none

      real*8 funcof,ab,bb,tt,poly,a,b,c,d,t
      integer*4 intb
C
C           REDUCE NUMBER OF REVOLUTIONS
C*      INT   =  T
C*      TINT  =  INT
C*      TT    =  T-TINT
C*      INT   =  INT*INTB
C*      AB    =  INT-((INT/10000)*10000)
C*      BB    =  INTB
C       REPLACE THE ABOVE 5 C* STATMENTS WITH THE FOLLOWING 4 BECAUSE OF LARGE
C       INTEGER PROBLEM
      TT    =  DMOD(T,1.D0)
      BB    =  INTB
      AB    =  DINT(T)*BB
      AB    =  AB-(DINT(AB/10000.D0)*10000.D0)
C
C       EVALUATE POLYNOMIAL
      POLY  =  A+(AB+BB*TT)*1.0D-4+T*(B+T*(C+T*D))
C
C       GET BETWEEN -1 AND 1 REVOLUTIONS
C*      INT   =  POLY
C*      TINT  =  INT
C*      FUNCOF=  POLY-TINT
C       THE FOLLOWING STATEMENT REPLACES THE ABOVE 3 C* STATEMENTS
      FUNCOF=  DMOD(POLY,1.D0)
C
      RETURN
      END
