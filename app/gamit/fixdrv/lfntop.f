      INTEGER* 4  FUNCTION  LFNTOP( FLNAME )
C     Function to return the location of first character of file name
C     in the file path name
C        S.Shimada        05/09/90  at IGPP
C
      INTEGER* 4  I
      CHARACTER   FLNAME*(*)
C
      IF( LEN(FLNAME) .EQ. 0 )
     .   call report_stat('FATAL','FIXDRV','lfntop','flname'
     .                   ,'Improper file name:',0)
      DO  10  I = LEN(FLNAME), 1, -1
         IF( FLNAME(I:I).EQ.'/' )   THEN
            LFNTOP = I + 1
            RETURN
         ENDIF
   10 CONTINUE
      LFNTOP = 1
      RETURN
      END
