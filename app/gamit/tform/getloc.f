      SUBROUTINE GETLOC( ILOC,NUMSIT,X,SITNAM,IFILE)
C
C     Read in a set or file of local coordinates.
C     ILOC 1=North,East,Up (meters)
C          2=Delta Lat,Long (decimal degrees), Up (meters)

      implicit none

      include '../includes/tform.h'
         
      logical fcheck
      real*8 x(3,nsdim), x4(3,nsdim)
     .      , north,east,up,shft(3),radcon,slat,slon,alat,alon
     .      , semi,finv,tx,ty,tz,clat,clon,blat,blon,bh
      integer*4 iloc,numsit,ifile,ifile4,ityp4
     .        , idms4,iloc4,ioerr,iflag,i,j 
                   
      character*10 tmpnam
      character*16 sitnam(nsdim)
      character*256 line

      radcon = atan(1.d0)/45.d0

C
C  Obtain local coordinates

c --  read through all records - eof exit to 90   
      numsit = 0 
      do i=1,nsdim
         if(ifile.gt.0) then
           IF( ILOC.EQ.1 ) THEN 
              read(ifile,'(a)',iostat=ioerr,end=90) line
              if( line(1:1).eq.' ') then
                READ(line(2:9),'(a8)') sitnam(i)(1:8)
                read(line(10:256),*) north,east,up
                numsit = numsit + 1                   
                if( iprnt.gt.0 )
     .           write(iprnt,'(1x,a8,3f14.3)') 
     .                      sitnam(numsit)(1:8),north,east,up
              endif
           ELSEIF( ILOC.EQ.2 ) THEN 
              read(ifile,'(a)',iostat=ioerr,end=90) line
              if( line(1:1).eq.' ') then
                READ(line(2:9),'(a8)')  sitnam(i)(1:8) 
                read(line(10:256),*) alat,alon,up
                numsit = numsit + 1            
                if( iprnt.gt.0 ) write(iprnt,'(1x,a8,3f14.3)') 
     .                         sitnam(i)(1:8),alat,alon,up
              endif
           ENDIF
         else
           IF( ILOC.EQ.1 ) THEN
              WRITE(ISCRN,110)
  110         FORMAT(/,1X,'Enter North, East, Up in m')
              READ(ITERM,*) NORTH,EAST,UP
           ELSEIF( ILOC.EQ.2 ) THEN
              WRITE(ISCRN,120)
  120         FORMAT(/,1X,'Enter delta lat, long, up in deg')
              READ(ITERM,*) ALAT,ALON,UP
           ENDIF
        endif

C   Determine Location of Local vector

      IF( I.GT.1 ) GOTO 300

      WRITE(ISCRN,200)
  200 FORMAT(/,1X,
     1  'Specify coordinate input for the local vector location:',/,3x,
     2  'Note:  All system types, except Geodetic (3), assume ',/,3x,
     2  '       the Up vector is parallel to the geocentric ',/,3x,
     3  '       location vector.  Geodetic assumes the up vector',/,3x,
     4  '       is parallel to the local ellipsoidal normal vector.'/)
C
  205 IFILE4 = 4
      CALL FORSYS(ITYP4,IDMS4,ILOC4,IFILE4)

C
      IF( ITYP4.EQ.4 ) THEN
         WRITE(ISCRN,210)
  210    FORMAT(/,1X,'Do not be redundant, try another type.'/)
         GOTO 205

      ENDIF
C

      CALL GETCOR(ITYP4,IDMS4,ILOC4,NUMSIT,X4,SITNAM,IFILE4)

C
  300 CONTINUE
      IF( ITYP4.EQ.3) THEN
        IF( I.EQ.1 ) then
c**          old routine: reinstate eventually with standard tables
c**   1      CALL GEOTAB(ISYS4,SEMI,FINV,TX,TY,TZ)
        tmpnam = 'gdetic.dat'
      if (fcheck(tmpnam)) then
        open(unit=idatum,file=tmpnam,iostat=ioerr,status= 'OLD',err=5)
    5     if (ioerr .ne. 0) then   
            call report_stat('FATAL','TFORM','getloc',' '
     .         ,'Missing datum file (gdetic.dat)',ioerr)
           endif
c*        call gdatum(idatum,iterm,iscrn,iprnt,semi,finv,shft)
c          gdatum returns kilometers
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
       endif
         IFLAG = 2
         CALL GEOXYZ(SEMI,FINV,TX,TY,TZ,BLAT,BLON,BH,
     1     X4(1,I),X4(2,I),X4(3,I),IFLAG)
      ELSE
         CALL CARSPH(X4(1,I),BLAT,BLON,BH)
      ENDIF
C
C  For NEU, rotate into Cartesian
C
      IF( ILOC.EQ.1 ) THEN
         SLAT = DSIN(BLAT*RADCON)
         SLON = DSIN(BLON*RADCON)
         CLAT = DCOS(BLAT*RADCON)
         CLON = DCOS(BLON*RADCON)
C
         X(1,I) = -SLAT*CLON*NORTH -SLON*EAST +CLAT*CLON*UP
         X(2,I) = -SLAT*SLON*NORTH +CLON*EAST +CLAT*SLON*UP
         X(3,I) =       CLAT*NORTH                 +SLAT*UP
C
C  For other formats, add local to degrees, convert to Cartesian
C  and subtract to find local Cartesian vector
C
      ELSE
         ALAT  = ALAT  + BLAT
         ALON  = ALON + BLON
         UP    = UP    + BH
C
         IF( ITYP4.EQ.3 ) THEN
            IFLAG = 1
            CALL GEOXYZ(SEMI,FINV,TX,TY,TZ,ALAT,ALON,UP,
     1        X(1,I),X(2,I),X(3,I),IFLAG)
         ELSE
            CALL SPHXYZ(ALAT,ALON,UP,X(1,I))
         ENDIF
C
         DO J = 1, 3
            X(J,I) = X(J,I) - X4(J,I)
         ENDDO
      ENDIF

C If reading from screen, skip out of loop
      if(ifile4.eq.0) goto 90

c---- end loop on records of table
      enddo
C
   90 CONTINUE

        if( numsit.eq.nsdim ) write(iscrn,96) NSDIM
   96    FORMAT(1X,'Number of sites in input file equals dimension '
     .         ,i2,' in TFORM/getloc.  Have some been missed?')

C
      RETURN
      END
