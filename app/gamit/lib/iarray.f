      FUNCTION IARRAY( JTEST,JARRAY,N )
C
C     Determine the position in the array JARRAY of the number JTEST
C
C        R. King 8 June 1987
C
C
      INTEGER JARRAY(*),jtest,iarray,n,i
C
        IARRAY = 0
        DO 1 I = 1,N
           IF(JTEST.EQ.JARRAY(I)) THEN
              IARRAY = I
           ENDIF
1       CONTINUE

        RETURN
        END
