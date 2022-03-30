      SUBROUTINE GETCOR(ITYP,IDMS,ILOC,NUMSIT,X,
     1                  SITNAM,IFILE)
C
C  Gets coordinates from terminal or file and converts
C  to Cartesian

      implicit none

      include '../includes/tform.h'

      integer*4 numsit,ifile,ityp,idms,iloc

      real*8 x(3,nsdim)

      character*16 sitnam(nsdim)
    

      IF( ITYP.EQ.1 ) THEN
         CALL GETCAR(NUMSIT,X,SITNAM,IFILE)
C
      ELSE IF( ITYP.EQ.2 ) THEN
         CALL GETSPH(IDMS,NUMSIT,X,SITNAM,IFILE)
C
      ELSE IF( ITYP.EQ.3 ) THEN
         CALL GETGEO(IDMS,NUMSIT,X,SITNAM,IFILE) 
C
      ELSE IF( ITYP.EQ.4 ) THEN
         CALL GETLOC(ILOC,NUMSIT,X,SITNAM,IFILE)
C
      ELSE IF( ITYP.EQ.5 ) THEN
         CALL GETCYL(NUMSIT,X,SITNAM,IFILE)
C
      ENDIF
C             
      RETURN
      END
