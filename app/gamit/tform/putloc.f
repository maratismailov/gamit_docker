      SUBROUTINE PUTLOC( ILOC,NUMSIT,X,SITNAM,OFILE)
C
C       Write out a set or file of local coordinates.

      implicit none

      include '../includes/tform.h'

      integer*4 ofile,iloc,numsit,ifile4,ityp4,idms4,iloc4
     .        , iflag,ioerr,i

      real*8 radcon,x(3,nsdim),x4(3,nsdim)
     .      ,shft(3),semi,finv,tx,ty,tz,slat,slon,clat,clon
     .     , north,east,up,blat,blong,bh

      character*10 tmpnam
      character*16 sitnam(nsdim)

      logical fcheck

      RADCON=DATAN(1.D0)/45.D0

C
      WRITE(ISCRN,200)
  200 FORMAT(/,1X,
     1  'Specify coordinate input for the local vector location:',/,3x,
     2  'Note:  All system types, except Geodetic (3), assume ',/,3x,
     2  '       the Up vector is parallel to the geocentric ',/,3x,
     3  '       location vector.  Geodetic assumes the up vector',/,3x,
     4  '       is parallel to the local ellipsoidal normal vector.'/)
C
      IFILE4 = 4
  205 CALL FORSYS(ITYP4,IDMS4,ILOC4,IFILE4 )
      CALL GETCOR(ITYP4,IDMS4,ILOC4,NUMSIT,X4,
     4                  SITNAM,IFILE4)
      IF( ITYP4.EQ.4 ) THEN
         WRITE(ISCRN,210)
  210    FORMAT(/,1X,'Do not be redundant, try another type.'/)
         GOTO 205
      ENDIF


C  Convert to decimal degrees

      IF( IPRNT.GT.0 ) WRITE(IPRNT,'(/,A)') ' Result:'

c-----Loop over sites
      DO I = 1, NUMSIT

c     If type is geodetic, get the transformation parameters
      IF( ITYP4.EQ.3) THEN

c        do get the transformations only once
         IF( I.EQ.1 ) then
c**         old subroutine: revive eventually but use standard tables
c**         CALL GEOTAB(ISYS4,SEMI,FINV,TX,TY,TZ)
            tmpnam = 'gdetic.dat'
            if( fcheck(tmpnam)) then
             open(unit=idatum,file=tmpnam,status='old',iostat=ioerr
     .                   ,err=5)
    5        if (ioerr .ne. 0) then
              write (iscrn,*) 'OPEN: error opening datum file: ',tmpnam
             endif
c**             call gdatum(idatum,iterm,iscrn,iprnt,semi,finv,shft)
c            gdatum returns kilometers
             close(idatum)
            else
             write(iscrn,*) 'Missing geodetic datum file (gdetic.dat)'
             write(iscrn,*) 'Using default WGS84'
               SEMI= 6378.137D0
               FINV= 298.257222101D0
               SHFT(1) = 0.D0
               SHFT(2) = 0.D0
               SHFT(3) = 0.D0
            endif
          semi = semi*1.d3
          tx = shft(1)*1.d3
          ty = shft(2)*1.d3
          tz = shft(3)*1.d3
c        end of I=1 check
         endif

         IFLAG = 2
         CALL GEOXYZ( SEMI,FINV,TX,TY,TZ,BLAT,BLONG,BH
     1              , X4(1,I),X4(2,I),X4(3,I),IFLAG)

       else
c      assume type is spherical
           CALL CARSPH(X4(1,I),BLAT,BLONG,BH)
       endif
C
C  Rotate Cartesian to NEU
C
       IF( ILOC.EQ.1 ) THEN
            SLAT = DSIN(BLAT*RADCON)
            SLON = DSIN(BLONG*RADCON)
            CLAT = DCOS(BLAT*RADCON)
            CLON = DCOS(BLONG*RADCON)
            NORTH = -SLAT*CLON*X(1,I) -SLAT*SLON*X(2,I) +CLAT*X(3,I)
            EAST  =      -SLON*X(1,I)      +CLON*X(2,I)
            UP    =  CLAT*CLON*X(1,I) +CLAT*SLON*X(2,I) +SLAT*X(3,I)
            IF( OFILE.GT.0 ) WRITE(OFILE,'(a8,3f15.4)') 
     .          SITNAM(I)(1:8),NORTH,EAST,UP
            IF( IPRNT.GT.0 ) WRITE(IPRNT,'(1x,a16,3f15.4)') 
     .           SITNAM(I),NORTH,EAST,UP
         ELSE
            PRINT*,' OPTION NOT AVAILABLE'
         ENDIF

      ENDDO

C
      RETURN
      END
