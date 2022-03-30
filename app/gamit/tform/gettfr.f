      SUBROUTINE GETTFR( SCALE,T,R )
C
C      Read transformation parameters from the terminal or a file.

      implicit none

      include '../includes/tform.h'

      integer*4 tfile,index,nskip,index1,isign,i,j
      real*8 scale,t,r

      CHARACTER*80 RECORD
      CHARACTER*8 NAME1,NAME2,BLANK8
      DIMENSION T(3),R(3)

* MOD TAH 921002: Added home_dir string and used SUN system getenv
*     subroutine to return the user home directory.  The opening
*     of tform.dat below assumes that the user has gu link in
*     directory and that this link points to correct location.  This
*     is the same assumption as in the original code.  The mod makes
*     the feature work.

* Home_dir  - String to hold the name of the users home directory
* tform_file - Full name and path to tform.dat file
* nblen     - GAMIT lib function to return of string

c*    Temporarily comment out to make this work without a gu.
c*    Note also the name change to tform.dat (vice tform.data) to
c*    be consistent with other GAMIT tables.   rwk 921005
c*      character*128 home_dir
        character*128 tform_file
c*      integer*4 nblen

      DATA BLANK8/'       '/
C
C      Read and display the available transformations
C
       TFILE=4

* MOD TAH 921002: Generate name of the tform.data file.
*     The following code replaces the line below.
c      OPEN (UNIT=TFILE, FILE='~/gu/tables/tform.data',status='unknown')

*     Get the users home directory
c*      call getenv('HOME', home_dir )     comment for now - rwk

*     Generate tform.data file name
c*      tform_file = home_dir(1:max(1,nblen(home_dir))) //
c*     .              '/gu/tables/tform.data'
      tform_file = 'tform.dat'
      OPEN (UNIT=TFILE, FILE=tform_file,status='unknown')

* END OF MOD TAH 921002 and RWK 921005:

         REWIND (TFILE)
C
         DO 75 I=1,100
            READ(TFILE,70,END=80) RECORD
            WRITE(ISCRN,70) RECORD
   75    CONTINUE
   70    FORMAT(a80)
C
   80    CONTINUE
         WRITE(ISCRN,85)
   85    FORMAT(/,1X,'Enter index of transformation desired',
     .               ' (negative if inverse):')
         READ(ITERM,*) INDEX
         REWIND TFILE
C
         IF(INDEX.NE.0) THEN
C
         NSKIP=ABS(INDEX)+4
         DO 90 I=1,NSKIP
   90       READ(TFILE,91)
   91    FORMAT()
C
         READ(TFILE,92) INDEX1,NAME1,NAME2,SCALE,
     1        (T(I),I=1,3),(R(J),J=1,3)
   92    FORMAT(I2,2X,A8,4X,A8,F8.5,2X,3F8.4,1X,3F8.5)
         WRITE(IPRNT,96) NAME1,NAME2,SCALE,
     1        (T(I),I=1,3),(R(J),J=1,3)
   96    FORMAT(/,1X,'Coordinates transformed from System ',A8,
     1    ' to System ',A8,/,3X,'Scale=',F8.5,' ppm    DX=',F9.4,
     2    '  DY=',F9.4,'  DZ=',F9.4,' M',/,3X,'R1 =',F8.5,'   R2 =',
     3    F8.5,'    R3=',F8.5)
C
      ELSE
C
C
C           Obtain the transformations parameters interactively
C
         WRITE(ISCRN,30)
   30    FORMAT(/,1X,'Enter scale difference (ppm) (New-Old)')
         READ(ITERM,*) SCALE
         WRITE(ISCRN,40)
   40    FORMAT(/,1X,'Enter translation: DX DY DZ  (New-Old)')
         READ(ITERM,*) T
         WRITE(ISCRN,50)
   50    FORMAT(/,1X,'Enter orientation angles: R1  R2  R3  (arcsec)')
         READ(ITERM,*) R
         NAME1= BLANK8
         NAME2= BLANK8
         INDEX= 1
      ENDIF
C
C
C       Change units of scale
C
      SCALE= DBLE(ISIGN(1,INDEX))*SCALE*1.D-6
      DO I = 1, 3
        T(I) = DBLE(ISIGN(1,INDEX))*T(I)
        R(I) = DBLE(ISIGN(1,INDEX))*R(I)
      ENDDO
C
      RETURN
      END
