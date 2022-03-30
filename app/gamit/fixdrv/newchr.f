      SUBROUTINE  NEWCHR( CHR )
C
C     Update year/iteration code
C
      CHARACTER*1  CHR

      IF( ICHAR(CHR).GE.ICHAR('0') .AND.
     +    ICHAR(CHR).LE.ICHAR('9') ) THEN
         CHR = 'a'
      ELSEIF( ICHAR(CHR).GE.ICHAR('A') .AND.
     +        ICHAR(CHR).LT.ICHAR('Z') ) THEN
         CHR = CHAR( ICHAR(CHR)+1+ICHAR('a')-ICHAR('A') )
      ELSEIF(ICHAR(CHR).GE.ICHAR('a') .AND.
     +         ICHAR(CHR).LT.ICHAR('z') ) THEN
         CHR = CHAR( ICHAR(CHR)+1 )
      ELSE
         CHR = 'a'
      ENDIF
C
      RETURN
      END

