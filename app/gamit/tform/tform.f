      Program TFORM
C
C     R. King 3 February 1987
C     MHMurray 29 Oct 1987   Modified for Apollo
c     R. King 15 October 2002  - modified for limited formats and addition of 
c                                command-line program convgeo
C
C     Compute transformations between sets of coordinates.  TFORM
C     uses Cartesian coordinates internally, but can input or output
C     coordinates in one of several representations:
C
C
C       1  Cartesian -- X,Y,Z or DX,DY,DZ in meters
C       2  Spherical -- LAT,LONG(E.),RAD in dec. deg or deg,min,sec, m
C       3  Geodetic  -- LAT,LONG(E.),HT in decimal deg or deg,min,sec
C                       and meters, plus A,FINV,DX,DY,DZ in meters
C       4  Cylindrical- RHO,LONG(E.),Z in m, decimal deg. or deg,min,sec
C       5  Local     -- U,V,W in meters, plus geodetic coodinates
C
C     TFORM can estimate the transformation parameters between sets
C     of coordinates in two reference system, OR it can transform
C     coordinates from one system to another from a set of input
C     parameters.
C                      
c     Coordinates are input and output either via the terminal/screen or
c     files, with the following allowed formats:
c
c       Terminal/screen:  a) free-format, with no station name (meters or decimal degrees)
c                         b) GAMIT L-file with site id in columns 1-4, name in columns 5-16 (deg/min/sec and meters)
c
c       File:  a) GLOBK apr format (8-character station name in column 2, 
c                    free-format coordinates  (meters or decimal degrees)
c              b) GAMIT L-file with site id in column 1, name in columns 5-16 (deg/min/sec and meters)
c
c       NOTE:  TFORM no longer expects/supports files written with headers
c              and format statements


      implicit none 

c       tform.h has dimensions of station and coordinate arrays and 
c       a common for global unit numbers (terminal, screen, gdetic.dat)

      include '../includes/tform.h'

      integer*4 numsit,ityp1,idms1,iloc1,iop,ityp2
     .        , idms2,iloc2,ityp3,idms3,iloc3
     .        , ifile1, ifile2, ifile3, i,j     
                                
      real*8 x1(3,nsdim), x2(3,nsdim), x2n(3,nsdim)
     .      ,omc(ndim), error(ndim), deriv(ndim,7)
     .      ,t(3), r(3), xlen, scale, epoch
            
      character*5 datum
      character*16 sitnam(nsdim)

c  Output program header

      CALL TVERSN

c  Assign unit numbers  (all in common in tform.h)
              
c     these four global, stored in common in tform.h
      iterm = 5
      iscrn = 6
      iprnt = 6   
c     gdetic.dat
      idatum = 25
c     Input data files (unit number is used to detect whether input (1 or 2) or output (3)
      ifile1 = 0
      ifile2 = 0
      ifile3 = 0   
c
      iloc1 = 0
      iloc2 = 0

c  Input initial coordinates

    1 CONTINUE
      NUMSIT = 0
      WRITE(ISCRN,5)
    5 FORMAT(/,1X,'Specify initial coordinate characteristics:'/)
C
      IFILE1 = 1
      CALL FORSYS(ITYP1,IDMS1,ILOC1,IFILE1)
C
C
C   Get Input Coordinates
C
      CALL GETCOR(ITYP1,IDMS1,ILOC1,NUMSIT,X1,
     1                  SITNAM,IFILE1) 
         
C   Set up Operation

    9 WRITE(ISCRN,10)
C   10 FORMAT(/,1X,'Specify desired operation:',/,3X,
   10 FORMAT(/,1X,
     1 '1=Add',/,1X,
     1 '2=Subtract',/,1X,
     1 '3=Cartesian length',/,1X,
     2 '4=Convert coordinate type (e.g. cartesian=>spherical)',/,1X,
     2 '5=Transform coordinates to another system (e.g. WGS=>SV6)',/,1X,
     2 '6=Estimate transformation',/,1X,
     3 'Specify operation: ',$)
      READ(ITERM,*) IOP
C
      IF( IOP.LT.1 .OR. IOP.GT.6) THEN
         WRITE(ISCRN,1000)
         GOTO 9
      ENDIF
C
C   Add
C
      IF( IOP.EQ.1 ) THEN
         WRITE(ISCRN,50)
   50    FORMAT(/,1X,
     1    'Specify characteristics of coordinates to be added:'/)
C
         IFILE2 = 2
         CALL FORSYS(ITYP2,IDMS2,ILOC2,IFILE2)
C
         CALL GETCOR(ITYP2,IDMS2,ILOC2,NUMSIT,X2,
     1                  SITNAM,IFILE2)
         DO J = 1, NUMSIT
         DO I = 1, 3
            X2(I,J) = X1(I,J) + X2(I,J)
         ENDDO
         ENDDO
C
C   Subtract
C
      ELSEIF( IOP.EQ.2 ) THEN
         WRITE(ISCRN,60)
   60    FORMAT(/,1X,
     1    'Specify characteristics of coordinates to be subtracted:'/)
C
         IFILE2 = 2
         CALL FORSYS(ITYP2,IDMS2,ILOC2,IFILE2)
C
         CALL GETCOR(ITYP2,IDMS2,ILOC2,NUMSIT,X2,
     1                  SITNAM,IFILE2)
         DO J = 1, NUMSIT
         DO I = 1, 3
            X2(I,J) = X1(I,J) - X2(I,J)
         ENDDO
         ENDDO
C
C   Length
C
      ELSEIF( IOP.EQ.3 ) THEN
         WRITE(ISCRN,70)
   70    FORMAT(/,1X,'Length of current Cartesian vectors'/)
C
         DO J = 1, NUMSIT
            XLEN = 0.D0
            DO I = 1, 3
               XLEN = XLEN + X1(I,J)*X1(I,J)
               X2(I,J) = X1(I,J)
            ENDDO
            XLEN = DSQRT(XLEN)
            WRITE(ISCRN,75) SITNAM(J),XLEN
   75       FORMAT(1X,A16,3F15.4)
         ENDDO
         GOTO 199
C
C   Convert
C
      ELSEIF( IOP.EQ.4 ) THEN
         DO J = 1, NUMSIT
         DO I = 1, 3
            X2(I,J) = X1(I,J)
         ENDDO
         ENDDO
C
C   Transform
C
      ELSEIF( IOP.EQ.5 ) THEN
         CALL GETTFR( SCALE,T,R )
         DO J=1,NUMSIT
            CALL TF2REF( X1(1,J),X2(1,J),SCALE,T,R )
         ENDDO
C
C    Estimate Transformation
C
      ELSEIF( IOP.EQ.6 ) THEN
         WRITE(ISCRN,80)
   80    FORMAT(/,1X,'Specify System 2 characteristics:',/,3x,
     1      'Note:  System 2 is rotated into System 1'/)
C
         IFILE2 = 2
         CALL FORSYS(ITYP2,IDMS2,ILOC2,IFILE2)
C
         CALL GETCOR(ITYP2,IDMS2,ILOC2,NUMSIT,X2,
     1                  SITNAM,IFILE2)
C
         CALL ESTREF( NUMSIT,SITNAM,X1,X2,X2N,OMC,ERROR,DERIV )
C
         DO J = 1, NUMSIT
         DO I = 1, 3
            X2(I,J) = X2N(I,J)
         ENDDO
         ENDDO
      ENDIF
C
C    Output results
C
         WRITE(ISCRN,100)
  100    FORMAT(/,1X,'Specify output coordinate characteristics:'/)
C
         IFILE3 = 3
         CALL FORSYS(ITYP3,IDMS3,ILOC3,IFILE3)
C
         IF(ITYP3.EQ.1) THEN
            CALL PUTCAR( NUMSIT,X2,SITNAM,IFILE3 )
C
         ELSE IF(ITYP3.EQ.2) THEN
            CALL PUTSPH( IDMS3,NUMSIT,X2,SITNAM,IFILE3)
C
         ELSE IF(ITYP3.EQ.3) THEN  
c           these entries read interactively for TFORM but passed by CONVGEO
            datum = '     '     
            epoch = 0.d0
            CALL PUTGEO( IDMS3,NUMSIT,X2,SITNAM,datum,epoch,IFILE3)
C
         ELSE IF(ITYP3.EQ.4) THEN
            CALL PUTLOC( ILOC3,NUMSIT,X2,SITNAM,IFILE3)
C
         ELSE IF(ITYP3.EQ.5) THEN
            CALL PUTCYL( NUMSIT,X2,SITNAM,IFILE3 )
C
         ENDIF
C
C   Reset, redo options
C
C
  199 CONTINUE
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
C
      WRITE(ISCRN,200)
C  200 FORMAT(/,1X,'Specify action:',/,3x,
  200 FORMAT(/,1X,
     1   '1=Continue    2=Start over   3=Quit: ',$)
      READ(ITERM,*) IOP
C
      IF( IOP.EQ.1 ) THEN
         DO J = 1, NUMSIT
         DO I = 1, 3
            X1(I,J) = X2(I,J)
         ENDDO
         ENDDO
         GOTO 9
      ELSEIF( IOP.EQ.2 ) THEN
         GOTO 1
      ELSEIF( IOP.EQ.3 ) THEN
         GOTO 999
      ENDIF
C
C    Finished
C
 1000 FORMAT(/,1X,'Not a valid option')
  999 CONTINUE
      STOP
      END
