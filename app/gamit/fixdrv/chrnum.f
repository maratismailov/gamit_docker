      CHARACTER * 3 FUNCTION  CHRNUM( N )
C     function to convert number N (<1000) to 3-digit character
      integer*4 n
      CHARACTER * 3  LIST
      WRITE( LIST, '(I3)' )  N
      CHRNUM = LIST
      IF( CHRNUM(1:1) .EQ. ' ' )  CHRNUM(1:1) = '0'
      IF( CHRNUM(2:2) .EQ. ' ' )  CHRNUM(2:2) = '0'
      RETURN
      END
