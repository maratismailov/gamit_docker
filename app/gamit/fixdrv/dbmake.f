      SUBROUTINE  DBMAKE( MFILE, MSESSN, DOY, BFIL2, deldbl )
C     subroutine to make DBLCLN batch file
C
C     Y. Bock     03/17/91
c                 Modified 91/07/17 YB
C
C     input arguments
C        MFILE  : M-file name                        input
C        MSESSN : Number of session                  input
C        DOY    : Day of year                        input
C        BFIL2  : secondary batch file name          input
C
      include '../includes/dimpar.h'
C
      INTEGER msessn

      CHARACTER* 1  YES, ans, deldbl
      CHARACTER* 3  DOY
      CHARACTER*16  BFIL2, MFILE
C
C     write "call to DBLCLN" in main batch file
C
      WRITE( 17, '(A,A16)' )  'dblcln < ', BFIL2
C
C     open secondary batch file and select mode
C
      OPEN( 21, FILE=BFIL2, STATUS='UNKNOWN' )
C
C     example of DBLCLN batch file:
C     msuma1.212     : M-file
c     y              : Use post-fit residuals
c     y              : Use all stations
c     y              : Use this series of C-files
c     1              : Increment C-file name and do not delete input C-file
c                              2 = increment and delete
c                              3 = overwrite (delete input, no increment)
C                              4 = do not write any C-files

C        M-file name
C
         WRITE( 21, '(A16)' )  MFILE
C
         IF(MSESSN.GT.1) WRITE( 21, '(A3)' ) DOY
         yes='y'
         WRITE( 21, '(A1,10x,A)' ) yes, 'Use post-fit residuals?'
         WRITE( 21, '(A1,10x,A)' ) yes, 'Use all these stations?'
         WRITE( 21, '(A1,10x,A)' ) yes, 'Use this series of C-files?'
c        set default to delete C-file
         ans = '0'
         if( deldbl.eq.'Y' ) then
           ans = '2'
           WRITE( 21, '(A1,10x,A)' )
     .          ans, 'Increment name and delete old C-file (1=keep)'
         else
           ans = '1'
           WRITE( 21, '(A1,10x,A)' )
     .          ans, 'Increment name and keep old C-file (0=delete)'

         endif

      CLOSE( 21 )
      RETURN
      END
