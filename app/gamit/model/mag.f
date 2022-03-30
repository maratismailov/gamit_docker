      subroutine mag (FNM,ITYPE,IOPT,IDM,DATE,ALT,XLT,XLN,NAME,X,Y,Z,F
     . ,debug)
C IN:	CHARACTER*30 	FNM - output filename
C			ITYPE - 1:geodetic coords (spheroidal earth)
C				2:geocentric coords (spherical earth)
C			IOPT - 	1: want value(s) at locations and dates
C				2: want value yearly intervals at one loc
C				3: want values on a lat/long grid at 1 date
C			IDM - lat long format 	1: degrees and minutes
C						2: decimal degrees
C			DATE - time in decimal years AD
C			ALT - height/km if ITYPE=1
C			    - radial dist from earth centre/km if ITYPE=2
C			XLT - Latitude/decimal degrees
C			XLN - Longitude/decimal degrees
C			NAME - name for location
C
C     This is a program for synthesising geomagnetic field values from the 
C     International Geomagnetic Reference Field series of models as agreed
c     in Decmember 2004 by IAGA Working Group V-MOD. 
C     It is the 10th generation IGRF, ie the 9th revision. 
C     The main-field models for 1900.0, 1905.0,..1940.0 and 2005.0 are 
C     non-definitive, those for 1945.0, 1950.0,...2000.0 are definitive and
C     the secular-variation model for 2005.0 to 2010.0 is non-definitive.
C
C     Main-field models are to degree and order 10 (ie 120 coefficients)
C     for 1900.0-1995.0 and to 13 (ie 195 coefficients) for 2000.0 onwards. 
C     The predictive secular-variation model is to degree and order 8 (ie 80
C     coefficients).
C
C     Options include values at different locations at different
C     times (spot), values at same location at one year intervals
C     (time series), grid of values at one time (grid); geodetic or
C     geocentric coordinates, latitude & longitude entered as decimal
C     degrees or degrees & minutes (not in grid), choice of main field 
C     or secular variation or both (grid only).
C
c     Adapted from 8th generation version to include new maximum degree for
c     main-field models for 2000.0 and onwards and use WGS84 spheroid instead
c     of International Astronomical Union 1966 spheroid as recommended by IAGA
c     in July 2003. Reference radius remains as 6371.2 km - it is NOT the mean
c     radius (= 6371.0 km) but 6371.2 km is what is used in determining the
c     coefficients. Adaptation by Susan Macmillan, August 2003 (for 
c     9th generation) and December 2004.
c     1995.0 coefficients as published in igrf9coeffs.xls and igrf10coeffs.xls
c     used - (Kimmo Korhonen spotted 1 nT difference in 11 coefficients)
c     Susan Macmillan July 2005
C
C Edited by EJP July 2007 to call as subroutine and produce a value for a location.
C 
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      real*8 DTMN, DTMX, FACT, DATE, ALT, XLT, XLN, CLT, X, Y, Z, F, D,
     .       H, S, DX, DY, DZ, F1, DD, DH, DS, DF, XLTI, XLTF, XLTD, 
     .       XLNI, XLNF, XLND
      integer*4 debug
      integer*4 ITYPE , IOPT , IDM , IU , NCOUNT , IDEC , IDECM , INC , 
     .          INCM , IDD , IDS , IH , IDH , IX , IDX , IY , IDY ,
     .          IZ , IDZ , NF , IDF , LTD , LTM , LND , LNM , IMX ,
     .          I , IFL , LTI , LTF , LNI , LNF , LT , LN 

      CHARACTER*1 IA
      CHARACTER*11 TYPE
      CHARACTER*20 NAME
      CHARACTER*30 FNM
      DATA DTMN,DTMX/1900.0,2015.0/
C
C
C      WRITE(6,*)
C      WRITE(6,*)'******************************************************'
        if(debug.ge.3) then
      WRITE(6,*)'MODEL\mag *     IGRF SYNTHESIS PROGRAM        *'
        endif
C      WRITE(6,*)'*                                                    *'
C      WRITE(6,*)'* A program for the computation of geomagnetic       *'
C      WRITE(6,*)'* field elements from the International Geomagnetic  *'
C      WRITE(6,*)'* Reference Field (9th generation) as revised in     *'
C      WRITE(6,*)'* December 2004 by the IAGA Working Group V-MOD.     *'
C      WRITE(6,*)'*                                                    *'
C      WRITE(6,*)'* It is valid for dates from 1900.0 to 2010.0,       *'
C      WRITE(6,*)'* values up to 2015.0 will be computed but with      *'
C      WRITE(6,*)'* reduced accuracy. Values for dates before 1945.0   *'
C      WRITE(6,*)'* and after 2000.0 are non-definitive, otherwise the *'
C      WRITE(6,*)'* values are definitive.                             *'
C     WRITE(6,*)'*                                                    *'
C      WRITE(6,*)'* Susan Macmillan          British Geological Survey *'
C      WRITE(6,*)'*                     Chair IAGA Working Group V-MOD *'
C      WRITE(6,*)'******************************************************'
C      WRITE(6,*)
C      WRITE(6,*)'Enter name of output file (30 characters maximum)'
C      WRITE(6,*)'or press "Return" for output to screen'
C      READ (5,991) FNM
C  991 FORMAT (A30)
      IF (ICHAR(FNM(1:1)).EQ.32) THEN
       IU = 6
      ELSE
       IU = 2
       OPEN (UNIT = IU,FILE = FNM,STATUS = 'UNKNOWN')
      END IF
C EJP - FACT is conversion factor from radians to degrees
      FACT = 180.0/3.141592654
      NCOUNT = 0
C 
C   10 WRITE(6,*)'Enter value for coordinate system:'
C      WRITE(6,*)
C     1'1 - geodetic (shape of Earth is approximated by a spheroid)'
C      WRITE(6,*)
C     1'2 - geocentric (shape of Earth is approximated by a sphere)'
C      READ (5,*) ITYPE
      IF (ITYPE.LT.1.OR.ITYPE.GT.2) Print*,'ITYPE not 1 or 2'
      IF (ITYPE.EQ.1) TYPE = ' geodetic  '
      IF (ITYPE.EQ.2) TYPE = ' geocentric'
C 
C   20 WRITE(6,*) 'Choose an option:'
C      WRITE(6,*) '1 - values at one or more locations & dates'
C      WRITE(6,*) '2 - values at yearly intervals at one location'
C      WRITE(6,*) '3 - values on a latitude/longitude grid at one date'
C      READ (5,*) IOPT
C      IF(IOPT.LT.1.OR.IOPT.GT.3) GO TO 20
      IF(IOPT.LT.1.OR.IOPT.GT.3) Print*,'IOPT is not 1-3'
      IF (IOPT.EQ.3) GO TO 150
C 
C   30 WRITE(6,*)'Enter value for format of latitudes and longitudes:'
C      WRITE(6,*)'1 - in degrees & minutes'
C      WRITE(6,*)'2 - in decimal degrees'
C      READ (5,*) IDM
C      IF (IDM.LT.1.OR.IDM.GT.2) GO TO 30
       IF (IDM.LT.1.OR.IDM.GT.2) Print*,'IDM is not 1 or 2'
      IF (NCOUNT.EQ.0) GOTO 50
C 
C   40 WRITE(6,*) 
C     1'Do you want values for another date & position? (y/n)'
C      READ (5,'(A1)') IA 
C EJP added
   40  IA='N'   
C      IF (IA.NE.'Y'.AND.IA.NE.'y'.AND.IA.NE.'N'.AND.IA.NE.'n')
C     1     GO TO 40
      IF(IA.EQ.'N'.OR.IA.EQ.'n') THEN
C       WRITE(IU,928)
  928  FORMAT (' D is declination (+ve east)'/
     1          ' I is inclination (+ve down)'/
     2          ' H is horizontal intensity'/
     3          ' X is north component'/
     4          ' Y is east component'/
     5          ' Z is vertical component (+ve down)'/
     6          ' F is total intensity')
C       WRITE(IU,929)
  929  FORMAT (/' SV is secular variation (annual rate of change)')
C       IF (ITYPE.EQ.2) THEN
C        WRITE(IU,*)
C     1'These elements are relative to the geocentric coordinate system'
C       ELSE
C        WRITE(IU,*)
C       ENDIF
C       Print*,'Debug mag.f : just before stop'
       Return
      ENDIF
C 
   50 NCOUNT = 1
      IF (IOPT.NE.2) THEN
C       WRITE(6,*) 'Enter date in years A.D.'
C       READ (5,*) DATE
       IF (DATE.LT.DTMN.OR.DATE.GT.DTMX) GO TO 209
      ENDIF

C      IF(ITYPE.EQ.1) THEN
C       WRITE(6,*) 'Enter altitude in km'
C      ELSE  
C       WRITE(6,*) 'Enter radial distance in km (>3485 km)'
C      END IF
C      READ (5,*) ALT
      IF (ITYPE.EQ.2.AND.ALT.LE.3485.0) GO TO 210
C 
C      IF (IDM.EQ.1) THEN
C       WRITE(6,*) 'Enter latitude & longitude in degrees & minutes'
C       WRITE(6,*) '(if either latitude or longitude is between -1'
C       WRITE(6,*) 'and 0 degrees, enter the minutes as negative).'
C       WRITE(6,*) 'Enter 4 integers' 
c       READ (5,*) LTD,LTM,LND,LNM
C       IF (LTD.LT.-90.OR.LTD.GT.90.OR.LTM.LE.-60.OR.LTM.GE.60) GO TO 204
C       IF (LND.LT.-360.OR.LND.GT.360.OR.LNM.LE.-60.OR.LNM.GE.60)
C     1    GO TO 205
C       IF (LTM.LT.0.AND.LTD.NE.0) GO TO 204
C       IF (LNM.LT.0.AND.LND.NE.0) GO TO 205
C       CALL DMDDEC (LTD,LTM,XLT)
C       CALL DMDDEC (LND,LNM,XLN)
C      ELSE
C       WRITE(6,*) 'Enter latitude & longitude in decimal degrees'
C       READ (5,*) XLT,XLN
       IF (XLT.LT.-90.0.OR.XLT.GT.90.0) GO TO 202
       IF (XLN.LT.-360.0.OR.XLN.GT.360.0) GO TO 203
C      ENDIF
C 
C      WRITE(*,*) 'Enter place name (20 characters maximum)'
C      READ (*,'(A)') NAME
C EJP - this forms the co-latitude
      CLT = 90.0 - XLT
      IF (CLT.LT.0.0.OR.CLT.GT.180.0) GO TO 204
      IF (XLN.LE.-360.0.OR.XLN.GE.360.0) GO TO 205
      IF (IOPT.EQ.2) GOTO 60
C 
      CALL IGRF10SYN (0,DATE,ITYPE,ALT,CLT,XLN,X,Y,Z,F)
      D = FACT*ATAN2(Y,X)
      H = SQRT(X*X + Y*Y)
      S = FACT*ATAN2(Z,H)
      CALL DDECDM (D,IDEC,IDECM)
      CALL DDECDM (S,INC,INCM)
C 
      CALL IGRF10SYN (1,DATE,ITYPE,ALT,CLT,XLN,DX,DY,DZ,F1)
      DD = (60.0*FACT*(X*DY - Y*DX))/(H*H)
      DH = (X*DX + Y*DY)/H
      DS = (60.0*FACT*(H*DZ - Z*DH))/(F*F)
      DF = (H*DH + Z*DZ)/F
C 
C EJP - write answer with deg & min or dec deg depending on input
C Have commented out so just get X,Y,Z,F passed back to calling program
C      IF (IDM.EQ.1) THEN
C       WRITE(IU,930) DATE,LTD,LTM,TYPE,LND,LNM,ALT,NAME
  930  FORMAT (1X,F8.3,' Lat',2I4,A11,' Long ',2I4,F10.3,' km ',A20)
C      ELSE
C       WRITE(IU,931) DATE,XLT,TYPE,XLN,ALT,NAME
  931  FORMAT (1X,F8.3,' Lat',F8.3,A11,' Long ',F8.3,F10.3,' km ',A20)
C      ENDIF
C 
      IDD = NINT(DD)
C      WRITE(IU,937) IDEC,IDECM,IDD
  937 FORMAT (15X,'D =',I5,' deg',I4,' min',4X,'SV =',I8,' min/yr')
C 
      IDS = NINT(DS)
C      WRITE(IU,939) INC,INCM,IDS
  939 FORMAT (15X,'I =',I5,' deg',I4,' min',4X,'SV =',I8,' min/yr')
C 
      IH = NINT(H)
      IDH = NINT(DH)
C      WRITE(IU,941) IH,IDH
  941 FORMAT (15X,'H =',I8,' nT     ',5X,'SV =',I8,' nT/yr')
C 
      IX = NINT(X)
      IDX = NINT(DX)
C      WRITE(IU,943) IX,IDX
  943 FORMAT (15X,'X =',I8,' nT     ',5X,'SV =',I8,' nT/yr')
C 
      IY = NINT(Y)
      IDY = NINT(DY)
C      WRITE(IU,945) IY,IDY
  945 FORMAT (15X,'Y =',I8,' nT     ',5X,'SV =',I8,' nT/yr')
C 
      IZ = NINT(Z)
      IDZ = NINT(DZ)
C      WRITE(IU,947) IZ,IDZ
  947 FORMAT (15X,'Z =',I8,' nT     ',5X,'SV =',I8,' nT/yr')
C 
      NF = NINT(F)
      IDF = NINT(DF)
C      WRITE(IU,949) NF,IDF
  949 FORMAT (15X,'F =',I8,' nT     ',5X,'SV =',I8,' nT/yr'/)
C
C EJP: 40 Asks if you want to ask it for another value/location - have set to exit
      GO TO 40
C 
   60 CONTINUE
C
C     SERIES OF VALUES AT ONE LOCATION...
C 
      IF (IDM.EQ.1) THEN
       WRITE(IU,932) LTD,LTM,TYPE,LND,LNM,ALT,NAME
  932  FORMAT ('Lat',2I4,A11,'  Long ',2I4,F10.3,' km ',A20)
      ELSE
       WRITE(IU,933) XLT,TYPE,XLN,ALT,NAME
  933  FORMAT ('Lat',F8.3,A11,'  Long ',F8.3,F10.3,' km ',A20)
      ENDIF
      WRITE (IU,934)
  934 FORMAT (3X,'DATE',7X,'D',3X,'SV',6X,'I',2X,'SV',6X,'H',4X,'SV',
     17X,'X',4X,'SV',7X,'Y',4X,'SV',7X,'Z',4X,'SV',6X,'F',4X,'SV')
      IMX = DTMX - DTMN - 5
      DO 70 I = 1,IMX
      DATE = DTMN - 0.5 + I
      CALL IGRF10SYN (0,DATE,ITYPE,ALT,CLT,XLN,X,Y,Z,F)
      D = FACT*ATAN2(Y,X)
      H = SQRT(X*X + Y*Y)
      S = FACT*ATAN2(Z,H)
      IH = NINT(H)
      IX = NINT(X)
      IY = NINT(Y)
      IZ = NINT(Z)
      NF = NINT(F)
C 
      CALL IGRF10SYN (1,DATE,ITYPE,ALT,CLT,XLN,DX,DY,DZ,F1)
      DD = (60.0*FACT*(X*DY - Y*DX))/(H*H)
      DH = (X*DX + Y*DY)/H
      DS = (60.0*FACT*(H*DZ - Z*DH))/(F*F)
      DF = (H*DH + Z*DZ)/F
      IDD = NINT(DD)
      IDH = NINT(DH)
      IDS = NINT(DS)
      IDX = NINT(DX)
      IDY = NINT(DY)
      IDZ = NINT(DZ)
      IDF = NINT(DF)
C
      WRITE(IU,935)
     1   DATE,D,IDD,S,IDS,IH,IDH,IX,IDX,IY,IDY,IZ,IDZ,NF,IDF
  935 FORMAT(1X,F6.1,F8.2,I5,F7.2,I4,I7,I6,3(I8,I6),I7,I6)
   70 CONTINUE
      IFL = 2
      GOTO 158
C
C     GRID OF VALUES...
C 
  150 WRITE(6,*)'Enter value for MF/SV flag:'
      WRITE(6,*)'0 for main field (MF)'
      WRITE(6,*)'1 for secular variation (SV)'
      WRITE(6,*)'2 for both'
      WRITE(6,*)'9 to quit'
      READ (5,*) IFL
      IF (IFL.EQ.9) STOP
      IF (IFL.NE.0.AND.IFL.NE.1.AND.IFL.NE.2) GOTO 150
C
      WRITE(6,*) 'Enter initial value, final value & increment or'
      WRITE(6,*) 'decrement of latitude, in degrees & decimals'
      READ (5,*) XLTI,XLTF,XLTD
      LTI = NINT(1000.0*XLTI)
      LTF = NINT(1000.0*XLTF)
      LTD = NINT(1000.0*XLTD)
      WRITE(6,*) 'Enter initial value, final value & increment or'
      WRITE(6,*) 'decrement of longitude, in degrees & decimals'
      READ (5,*) XLNI,XLNF,XLND
      LNI = NINT(1000.0*XLNI)
      LNF = NINT(1000.0*XLNF)
      LND = NINT(1000.0*XLND)
      IF (LTI.LT.-90000.OR.LTI.GT.90000) GO TO 206
      IF (LTF.LT.-90000.OR.LTF.GT.90000) GO TO 206
      IF (LNI.LT.-360000.OR.LNI.GT.360000) GO TO 207
      IF (LNF.LT.-360000.OR.LNF.GT.360000) GO TO 207
   98 WRITE(6,*) 'Enter date in years A.D.'
      READ (5,*) DATE
      IF (DATE.LT.DTMN.OR.DATE.GT.DTMX) GO TO 209
      IF (ITYPE.EQ.1) THEN
       WRITE(6,*) 'Enter altitude in km'
      ELSE
       WRITE(6,*) 'Enter radial distance in km (>3485 km)'
      END IF
      READ (5,*) ALT
      IF (ITYPE.EQ.2.AND.ALT.LE.3485.0) GO TO 210
      WRITE(IU,958) DATE,ALT,TYPE
  958 FORMAT (' Date =',F9.3,5X,'Altitude =',F10.3,' km',5X,A11//
     1        '      Lat     Long',7X,'D',7X,'I',7X,'H',7X,'X',7X,'Y',
     2        7X,'Z',7X,'F')
C
      LT = LTI
  151 XLT = LT
      XLT = 0.001*XLT
      CLT = 90.0 - XLT
      IF (CLT.LT.-0.001.OR.CLT.GT.180.001) GO TO 202
      LN = LNI
  152 XLN = LN
      XLN = 0.001*XLN
      IF (XLN.LE.-360.0) XLN = XLN + 360.0
      IF (XLN.GE.360.0) XLN = XLN - 360.0
      CALL IGRF10SYN (0,DATE,ITYPE,ALT,CLT,XLN,X,Y,Z,F)
      D = FACT*ATAN2(Y,X)
      H = SQRT(X*X + Y*Y)
      S = FACT*ATAN2(Z,H)
      IH = NINT(H)
      IX = NINT(X)
      IY = NINT(Y)
      IZ = NINT(Z)
      NF = NINT(F)
      IF (IFL.EQ.0) GOTO 153
      CALL IGRF10SYN (1,DATE,ITYPE,ALT,CLT,XLN,DX,DY,DZ,F1)
      IDX = NINT(DX)
      IDY = NINT(DY)
      IDZ = NINT(DZ)
      DD = (60.0*FACT*(X*DY - Y*DX))/(H*H)
      IDD = NINT(DD)
      DH = (X*DX + Y*DY)/H
      IDH = NINT(DH)
      DS = (60.0*FACT*(H*DZ - Z*DH))/(F*F)
      IDS = NINT(DS)
      DF = (H*DH + Z*DZ)/F
      IDF = NINT(DF)
C
  153 CONTINUE
      IF (IFL.EQ.0) WRITE(IU,959) XLT,XLN,D,S,IH,IX,IY,IZ,NF
      IF (IFL.EQ.1) WRITE(IU,960) XLT,XLN,IDD,IDS,IDH,IDX,IDY,IDZ,IDF
      IF (IFL.EQ.2) THEN
       WRITE(IU,959) XLT,XLN,D,S,IH,IX,IY,IZ,NF
       WRITE(IU,961) IDD,IDS,IDH,IDX,IDY,IDZ,IDF
      ENDIF      
  959 FORMAT (2F9.3,2F8.2,5I8)
  960 FORMAT (2F9.3,7I8)
  961 FORMAT (14X,'SV: ',7I8)
C
  154 LN = LN + LND
      IF (LND.LT.0) GO TO 156
      IF (LN.LE.LNF) GO TO 152
  155 LT = LT + LTD
      IF (LTD.LT.0) GO TO 157
      IF (LT - LTF) 151,151,158
  156 IF (LN - LNF) 155,152,152
  157 IF (LT.GE.LTF) GO TO 151
  158 CONTINUE
      IF (IFL.EQ.0.OR.IFL.EQ.2) THEN
       WRITE(IU,962)
  962  FORMAT (/' D is declination in degrees (+ve east)'/
     1          ' I is inclination in degrees (+ve down)'/
     2          ' H is horizontal intensity in nT'/
     3          ' X is north component in nT'/
     4          ' Y is east component in nT'/
     5          ' Z is vertical component in nT (+ve down)'/
     6          ' F is total intensity in nT')
      IF (IFL.NE.0) WRITE(IU,963)
  963  FORMAT (' SV is secular variation (annual rate of change)'/
     1' Units for SV: minutes/yr (D & I); nT/yr (H,X,Y,Z & F)')
      IF (ITYPE.EQ.2) WRITE(IU,*)
     1'These elements are relative to the geocentric coordinate system'
      ELSE
       WRITE(IU,964)
  964  FORMAT (/' D is SV in declination in minutes/yr (+ve east)'/
     1          ' I is SV in inclination in minutes/yr (+ve down)'/
     2          ' H is SV in horizontal intensity in nT/yr'/
     3          ' X is SV in north component in nT/yr'/
     4          ' Y is SV in east component in nT/yr'/
     5          ' Z is SV in vertical component in nT/yr (+ve down)'/
     6          ' F is SV in total intensity in nT/yr')
      IF (ITYPE.EQ.2) WRITE(IU,*)
     1'These elements are relative to the geocentric coordinate system'
      ENDIF
  159 STOP
C 
  209 WRITE(6,972) DATE
  972 FORMAT (' ***** Error *****'/' DATE =',F9.3,
     1        ' - out of range')
      STOP
C 
  210 WRITE(6,973) ALT,ITYPE
  973 FORMAT (' ***** Error *****'/' A value of ALT =',F10.3,
     1        ' is not allowed when ITYPE =',I2)
      STOP
C 
  202 WRITE(6,966) XLT
  966 FORMAT (' ***** Error *****'/' XLT =',F9.3,
     1        ' - out of range')
      STOP
C 
  203 WRITE(6,967) XLN
  967 FORMAT (' ***** Error *****'/' XLN =',F10.3,
     1        ' - out of range')
      STOP
C 
  204 WRITE(6,968) LTD,LTM
  968 FORMAT (' ***** Error *****'/' Latitude out of range',
     1        ' - LTD =',I6,5X,'LTM =',I4)
      STOP
C 
  205 WRITE(6,969) LND,LNM
  969 FORMAT (' ***** Error *****'/' Longitude out of range',
     1        ' - LND =',I8,5X,'LNM =',I4)
      STOP
C
  206 WRITE(6,970) LTI,LTF
  970 FORMAT (' ***** Error *****'/
     1        ' Latitude limits of table out of range - LTI =',
     2        I6,5X,' LTF =',I6)
      STOP
C
  207 WRITE(6,971) LNI,LNF
  971 FORMAT (' ***** Error *****'/
     1        ' Longitude limits of table out of range - LNI =',
     2        I8,5X,' LNF =',I8)
      STOP
C
      END
* MOD TAH 200124: Moved to igrf_util.f so all mgaxx.f routines can use.
C 
*     SUBROUTINE DMDDEC (I,M,X)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     IMPLICIT NONE
*     real*8 DE, EM, X
*     integer*4 i,m
*     
*     DE = I
*     EM = M
*     IF (I.LT.0) EM = -EM
*     X = DE + EM/60.0
*     RETURN
*     END
C 
*     SUBROUTINE DDECDM (X,I,M)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     IMPLICIT NONE
*     real*8 X, SIG, DR, T
*     integer*4 i,m, isig

*     SIG = SIGN(1.1d0,X)
*     DR = ABS(X)
*     I = INT(DR)
*     T = I
*     M = NINT(60.*(DR - T))
*     IF (M.EQ.60) THEN
*      M = 0
*      I = I + 1
*     ENDIF
*     ISIG = INT(SIG)
*     IF (I.NE.0) THEN
*      I = I * ISIG
*     ELSE
*      IF (M.NE.0) M = M * ISIG
*     ENDIF
*     RETURN
*     END
  
      subroutine igrf10syn (isv,date,itype,alt,colat,elong,x,y,z,f)
c
c     This is a synthesis routine for the 10th generation IGRF as agreed 
c     in December 2004 by IAGA Working Group V-MOD. It is valid 1900.0 to
c     2010.0 inclusive. Values for dates from 1945.0 to 2000.0 inclusive are 
c     definitve, otherwise they are non-definitive.
c   INPUT
c     isv   = 0 if main-field values are required
c     isv   = 1 if secular variation values are required
c     date  = year A.D. Must be greater than or equal to 1900.0 and 
c             less than or equal to 2015.0. Warning message is given 
c             for dates greater than 2010.0. Must be double precision.
c     itype = 1 if geodetic (spheroid)
c     itype = 2 if geocentric (sphere)
c     alt   = height in km above sea level if itype = 1
c           = distance from centre of Earth in km if itype = 2 (>3485 km)
c     colat = colatitude (0-180)
c     elong = east-longitude (0-360)
c     alt, colat and elong must be double precision.
c   OUTPUT
c     x     = north component (nT) if isv = 0, nT/year if isv = 1
c     y     = east component (nT) if isv = 0, nT/year if isv = 1
c     z     = vertical component (nT) if isv = 0, nT/year if isv = 1
c     f     = total intensity (nT) if isv = 0, rubbish if isv = 1
c
c     To get the other geomagnetic elements (D, I, H and secular
c     variations dD, dH, dI and dF) use routines ptoc and ptocsv.
c
c     Adapted from 8th generation version to include new maximum degree for
c     main-field models for 2000.0 and onwards and use WGS84 spheroid instead
c     of International Astronomical Union 1966 spheroid as recommended by IAGA
c     in July 2003. Reference radius remains as 6371.2 km - it is NOT the mean
c     radius (= 6371.0 km) but 6371.2 km is what is used in determining the
c     coefficients. Adaptation by Susan Macmillan, August 2003 (for 
c     9th generation) and December 2004.
c     1995.0 coefficients as published in igrf9coeffs.xls and igrf10coeffs.xls
c     used - (Kimmo Korhonen spotted 1 nT difference in 11 coefficients)
c     Susan Macmillan July 2005
c
c      implicit double precision (a-h,o-z)
      implicit none
      
      real*8 GH, G0,G1, G2, G3, G4, G5, G6, G7, G8, G9, GA, GB, GC, GD,
     .       GE, GF, GG, GI, GJ, GK, GL, GM, GP, X, Y, Z, DATE, T, ONE,
     .       TC, R, ALT, COLAT, CT, ST, ELONG, CL, SL, CD, SD, A2, B2,
     .       TWO, THREE, RHO, RATIO, RR, P, Q, FN, GN, FM, GMM, F
     
      integer*4 ISV ,ITYPE , LL , NMX , NC , KMX , L , M , N ,
     .          K , J , I , LM 

      dimension gh(3060),g0(120),g1(120),g2(120),g3(120),g4(120),
     1          g5(120),g6(120),g7(120),g8(120),g9(120),ga(120),
     2          gb(120),gc(120),gd(120),ge(120),gf(120),gg(120),
     3          gi(120),gj(120),gk(195),gl(195),gm(195),gp(195),
     4          p(105),q(105),cl(13),sl(13)
      equivalence (g0,gh(  1)),(g1,gh(121)),(g2,gh(241)),(g3,gh(361)),
     1            (g4,gh(481)),(g5,gh(601)),(g6,gh(721)),(g7,gh(841)),
     2            (g8,gh(961)),(g9,gh(1081)),(ga,gh(1201)),
     3            (gb,gh(1321)),(gc,gh(1441)),(gd,gh(1561)),
     4            (ge,gh(1681)),(gf,gh(1801)),(gg,gh(1921)),
     5            (gi,gh(2041)),(gj,gh(2161)),(gk,gh(2281)),
     6            (gl,gh(2476)),(gm,gh(2671)),(gp,gh(2866))
c
      data g0/ -31543.,-2298., 5922., -677., 2905.,-1061.,  924., 1121., 1900
     1           1022.,-1469., -330., 1256.,    3.,  572.,  523.,  876., 1900
     2            628.,  195.,  660.,  -69., -361., -210.,  134.,  -75., 1900
     3           -184.,  328., -210.,  264.,   53.,    5.,  -33.,  -86., 1900
     4           -124.,  -16.,    3.,   63.,   61.,   -9.,  -11.,   83., 1900
     5           -217.,    2.,  -58.,  -35.,   59.,   36.,  -90.,  -69., 1900
     6             70.,  -55.,  -45.,    0.,  -13.,   34.,  -10.,  -41., 1900
     7             -1.,  -21.,   28.,   18.,  -12.,    6.,  -22.,   11., 1900
     8              8.,    8.,   -4.,  -14.,   -9.,    7.,    1.,  -13., 1900
     9              2.,    5.,   -9.,   16.,    5.,   -5.,    8.,  -18., 1900
     a              8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., 1900
     b             -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,   -1., 1900
     c             -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., 1900
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1900
     e              0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6./ 1900
      data g1/ -31464.,-2298., 5909., -728., 2928.,-1086., 1041., 1065., 1905
     1           1037.,-1494., -357., 1239.,   34.,  635.,  480.,  880., 1905
     2            643.,  203.,  653.,  -77., -380., -201.,  146.,  -65., 1905
     3           -192.,  328., -193.,  259.,   56.,   -1.,  -32.,  -93., 1905
     4           -125.,  -26.,   11.,   62.,   60.,   -7.,  -11.,   86., 1905
     5           -221.,    4.,  -57.,  -32.,   57.,   32.,  -92.,  -67., 1905
     6             70.,  -54.,  -46.,    0.,  -14.,   33.,  -11.,  -41., 1905
     7              0.,  -20.,   28.,   18.,  -12.,    6.,  -22.,   11., 1905
     8              8.,    8.,   -4.,  -15.,   -9.,    7.,    1.,  -13., 1905
     9              2.,    5.,   -8.,   16.,    5.,   -5.,    8.,  -18., 1905
     a              8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., 1905
     b             -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., 1905
     c             -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., 1905
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1905
     e              0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6./ 1905
      data g2/ -31354.,-2297., 5898., -769., 2948.,-1128., 1176., 1000., 1910
     1           1058.,-1524., -389., 1223.,   62.,  705.,  425.,  884., 1910
     2            660.,  211.,  644.,  -90., -400., -189.,  160.,  -55., 1910
     3           -201.,  327., -172.,  253.,   57.,   -9.,  -33., -102., 1910
     4           -126.,  -38.,   21.,   62.,   58.,   -5.,  -11.,   89., 1910
     5           -224.,    5.,  -54.,  -29.,   54.,   28.,  -95.,  -65., 1910
     6             71.,  -54.,  -47.,    1.,  -14.,   32.,  -12.,  -40., 1910
     7              1.,  -19.,   28.,   18.,  -13.,    6.,  -22.,   11., 1910
     8              8.,    8.,   -4.,  -15.,   -9.,    6.,    1.,  -13., 1910
     9              2.,    5.,   -8.,   16.,    5.,   -5.,    8.,  -18., 1910
     a              8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., 1910
     b             -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., 1910
     c             -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., 1910
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1910
     e              0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6./ 1910
      data g3/ -31212.,-2306., 5875., -802., 2956.,-1191., 1309.,  917., 1915
     1           1084.,-1559., -421., 1212.,   84.,  778.,  360.,  887., 1915
     2            678.,  218.,  631., -109., -416., -173.,  178.,  -51., 1915
     3           -211.,  327., -148.,  245.,   58.,  -16.,  -34., -111., 1915
     4           -126.,  -51.,   32.,   61.,   57.,   -2.,  -10.,   93., 1915
     5           -228.,    8.,  -51.,  -26.,   49.,   23.,  -98.,  -62., 1915
     6             72.,  -54.,  -48.,    2.,  -14.,   31.,  -12.,  -38., 1915
     7              2.,  -18.,   28.,   19.,  -15.,    6.,  -22.,   11., 1915
     8              8.,    8.,   -4.,  -15.,   -9.,    6.,    2.,  -13., 1915
     9              3.,    5.,   -8.,   16.,    6.,   -5.,    8.,  -18., 1915
     a              8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., 1915
     b             -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., 1915
     c             -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., 1915
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1915
     e              0.,   -2.,    1.,    4.,    2.,    0.,    0.,   -6./ 1915
      data g4/ -31060.,-2317., 5845., -839., 2959.,-1259., 1407.,  823., 1920
     1           1111.,-1600., -445., 1205.,  103.,  839.,  293.,  889., 1920
     2            695.,  220.,  616., -134., -424., -153.,  199.,  -57., 1920
     3           -221.,  326., -122.,  236.,   58.,  -23.,  -38., -119., 1920
     4           -125.,  -62.,   43.,   61.,   55.,    0.,  -10.,   96., 1920
     5           -233.,   11.,  -46.,  -22.,   44.,   18., -101.,  -57., 1920
     6             73.,  -54.,  -49.,    2.,  -14.,   29.,  -13.,  -37., 1920
     7              4.,  -16.,   28.,   19.,  -16.,    6.,  -22.,   11., 1920
     8              7.,    8.,   -3.,  -15.,   -9.,    6.,    2.,  -14., 1920
     9              4.,    5.,   -7.,   17.,    6.,   -5.,    8.,  -19., 1920
     a              8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., 1920
     b             -3.,    1.,   -2.,   -2.,    9.,    2.,   10.,    0., 1920
     c             -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., 1920
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1920
     e              0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6./ 1920
      data g5/ -30926.,-2318., 5817., -893., 2969.,-1334., 1471.,  728., 1925
     1           1140.,-1645., -462., 1202.,  119.,  881.,  229.,  891., 1925
     2            711.,  216.,  601., -163., -426., -130.,  217.,  -70., 1925
     3           -230.,  326.,  -96.,  226.,   58.,  -28.,  -44., -125., 1925
     4           -122.,  -69.,   51.,   61.,   54.,    3.,   -9.,   99., 1925
     5           -238.,   14.,  -40.,  -18.,   39.,   13., -103.,  -52., 1925
     6             73.,  -54.,  -50.,    3.,  -14.,   27.,  -14.,  -35., 1925
     7              5.,  -14.,   29.,   19.,  -17.,    6.,  -21.,   11., 1925
     8              7.,    8.,   -3.,  -15.,   -9.,    6.,    2.,  -14., 1925
     9              4.,    5.,   -7.,   17.,    7.,   -5.,    8.,  -19., 1925
     a              8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., 1925
     b             -3.,    1.,   -2.,   -2.,    9.,    2.,   10.,    0., 1925
     c             -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., 1925
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1925
     e              0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6./ 1925
      data g6/ -30805.,-2316., 5808., -951., 2980.,-1424., 1517.,  644., 1930
     1           1172.,-1692., -480., 1205.,  133.,  907.,  166.,  896., 1930
     2            727.,  205.,  584., -195., -422., -109.,  234.,  -90., 1930
     3           -237.,  327.,  -72.,  218.,   60.,  -32.,  -53., -131., 1930
     4           -118.,  -74.,   58.,   60.,   53.,    4.,   -9.,  102., 1930
     5           -242.,   19.,  -32.,  -16.,   32.,    8., -104.,  -46., 1930
     6             74.,  -54.,  -51.,    4.,  -15.,   25.,  -14.,  -34., 1930
     7              6.,  -12.,   29.,   18.,  -18.,    6.,  -20.,   11., 1930
     8              7.,    8.,   -3.,  -15.,   -9.,    5.,    2.,  -14., 1930
     9              5.,    5.,   -6.,   18.,    8.,   -5.,    8.,  -19., 1930
     a              8.,   10.,  -20.,    1.,   14.,  -12.,    5.,   12., 1930
     b             -3.,    1.,   -2.,   -2.,    9.,    3.,   10.,    0., 1930
     c             -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., 1930
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1930
     e              0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6./ 1930
      data g7/ -30715.,-2306., 5812.,-1018., 2984.,-1520., 1550.,  586., 1935
     1           1206.,-1740., -494., 1215.,  146.,  918.,  101.,  903., 1935
     2            744.,  188.,  565., -226., -415.,  -90.,  249., -114., 1935
     3           -241.,  329.,  -51.,  211.,   64.,  -33.,  -64., -136., 1935
     4           -115.,  -76.,   64.,   59.,   53.,    4.,   -8.,  104., 1935
     5           -246.,   25.,  -25.,  -15.,   25.,    4., -106.,  -40., 1935
     6             74.,  -53.,  -52.,    4.,  -17.,   23.,  -14.,  -33., 1935
     7              7.,  -11.,   29.,   18.,  -19.,    6.,  -19.,   11., 1935
     8              7.,    8.,   -3.,  -15.,   -9.,    5.,    1.,  -15., 1935
     9              6.,    5.,   -6.,   18.,    8.,   -5.,    7.,  -19., 1935
     a              8.,   10.,  -20.,    1.,   15.,  -12.,    5.,   11., 1935
     b             -3.,    1.,   -3.,   -2.,    9.,    3.,   11.,    0., 1935
     c             -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., 1935
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1935
     e              0.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6./ 1935
      data g8/ -30654.,-2292., 5821.,-1106., 2981.,-1614., 1566.,  528., 1940
     1           1240.,-1790., -499., 1232.,  163.,  916.,   43.,  914., 1940
     2            762.,  169.,  550., -252., -405.,  -72.,  265., -141., 1940
     3           -241.,  334.,  -33.,  208.,   71.,  -33.,  -75., -141., 1940
     4           -113.,  -76.,   69.,   57.,   54.,    4.,   -7.,  105., 1940
     5           -249.,   33.,  -18.,  -15.,   18.,    0., -107.,  -33., 1940
     6             74.,  -53.,  -52.,    4.,  -18.,   20.,  -14.,  -31., 1940
     7              7.,   -9.,   29.,   17.,  -20.,    5.,  -19.,   11., 1940
     8              7.,    8.,   -3.,  -14.,  -10.,    5.,    1.,  -15., 1940
     9              6.,    5.,   -5.,   19.,    9.,   -5.,    7.,  -19., 1940
     a              8.,   10.,  -21.,    1.,   15.,  -12.,    5.,   11., 1940
     b             -3.,    1.,   -3.,   -2.,    9.,    3.,   11.,    1., 1940
     c             -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., 1940
     d             -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., 1940
     e              0.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6./ 1940
      data g9/ -30594.,-2285., 5810.,-1244., 2990.,-1702., 1578.,  477., 1945
     1           1282.,-1834., -499., 1255.,  186.,  913.,  -11.,  944., 1945
     2            776.,  144.,  544., -276., -421.,  -55.,  304., -178., 1945
     3           -253.,  346.,  -12.,  194.,   95.,  -20.,  -67., -142., 1945
     4           -119.,  -82.,   82.,   59.,   57.,    6.,    6.,  100., 1945
     5           -246.,   16.,  -25.,   -9.,   21.,  -16., -104.,  -39., 1945
     6             70.,  -40.,  -45.,    0.,  -18.,    0.,    2.,  -29., 1945
     7              6.,  -10.,   28.,   15.,  -17.,   29.,  -22.,   13., 1945
     8              7.,   12.,   -8.,  -21.,   -5.,  -12.,    9.,   -7., 1945
     9              7.,    2.,  -10.,   18.,    7.,    3.,    2.,  -11., 1945
     a              5.,  -21.,  -27.,    1.,   17.,  -11.,   29.,    3., 1945
     b             -9.,   16.,    4.,   -3.,    9.,   -4.,    6.,   -3., 1945
     c              1.,   -4.,    8.,   -3.,   11.,    5.,    1.,    1., 1945
     d              2.,  -20.,   -5.,   -1.,   -1.,   -6.,    8.,    6., 1945
     e             -1.,   -4.,   -3.,   -2.,    5.,    0.,   -2.,   -2./ 1945
      data ga/ -30554.,-2250., 5815.,-1341., 2998.,-1810., 1576.,  381., 1950
     1           1297.,-1889., -476., 1274.,  206.,  896.,  -46.,  954., 1950
     2            792.,  136.,  528., -278., -408.,  -37.,  303., -210., 1950
     3           -240.,  349.,    3.,  211.,  103.,  -20.,  -87., -147., 1950
     4           -122.,  -76.,   80.,   54.,   57.,   -1.,    4.,   99., 1950
     5           -247.,   33.,  -16.,  -12.,   12.,  -12., -105.,  -30., 1950
     6             65.,  -55.,  -35.,    2.,  -17.,    1.,    0.,  -40., 1950
     7             10.,   -7.,   36.,    5.,  -18.,   19.,  -16.,   22., 1950
     8             15.,    5.,   -4.,  -22.,   -1.,    0.,   11.,  -21., 1950
     9             15.,   -8.,  -13.,   17.,    5.,   -4.,   -1.,  -17., 1950
     a              3.,   -7.,  -24.,   -1.,   19.,  -25.,   12.,   10., 1950
     b              2.,    5.,    2.,   -5.,    8.,   -2.,    8.,    3., 1950
     c            -11.,    8.,   -7.,   -8.,    4.,   13.,   -1.,   -2., 1950
     d             13.,  -10.,   -4.,    2.,    4.,   -3.,   12.,    6., 1950
     e              3.,   -3.,    2.,    6.,   10.,   11.,    3.,    8./ 1950
      data gb/ -30500.,-2215., 5820.,-1440., 3003.,-1898., 1581.,  291., 1955
     1           1302.,-1944., -462., 1288.,  216.,  882.,  -83.,  958., 1955
     2            796.,  133.,  510., -274., -397.,  -23.,  290., -230., 1955
     3           -229.,  360.,   15.,  230.,  110.,  -23.,  -98., -152., 1955
     4           -121.,  -69.,   78.,   47.,   57.,   -9.,    3.,   96., 1955
     5           -247.,   48.,   -8.,  -16.,    7.,  -12., -107.,  -24., 1955
     6             65.,  -56.,  -50.,    2.,  -24.,   10.,   -4.,  -32., 1955
     7              8.,  -11.,   28.,    9.,  -20.,   18.,  -18.,   11., 1955
     8              9.,   10.,   -6.,  -15.,  -14.,    5.,    6.,  -23., 1955
     9             10.,    3.,   -7.,   23.,    6.,   -4.,    9.,  -13., 1955
     a              4.,    9.,  -11.,   -4.,   12.,   -5.,    7.,    2., 1955
     b              6.,    4.,   -2.,    1.,   10.,    2.,    7.,    2., 1955
     c             -6.,    5.,    5.,   -3.,   -5.,   -4.,   -1.,    0., 1955
     d              2.,   -8.,   -3.,   -2.,    7.,   -4.,    4.,    1., 1955
     e             -2.,   -3.,    6.,    7.,   -2.,   -1.,    0.,   -3./ 1955
      data gc/ -30421.,-2169., 5791.,-1555., 3002.,-1967., 1590.,  206., 1960
     1           1302.,-1992., -414., 1289.,  224.,  878., -130.,  957., 1960
     2            800.,  135.,  504., -278., -394.,    3.,  269., -255., 1960
     3           -222.,  362.,   16.,  242.,  125.,  -26., -117., -156., 1960
     4           -114.,  -63.,   81.,   46.,   58.,  -10.,    1.,   99., 1960
     5           -237.,   60.,   -1.,  -20.,   -2.,  -11., -113.,  -17., 1960
     6             67.,  -56.,  -55.,    5.,  -28.,   15.,   -6.,  -32., 1960
     7              7.,   -7.,   23.,   17.,  -18.,    8.,  -17.,   15., 1960
     8              6.,   11.,   -4.,  -14.,  -11.,    7.,    2.,  -18., 1960
     9             10.,    4.,   -5.,   23.,   10.,    1.,    8.,  -20., 1960
     a              4.,    6.,  -18.,    0.,   12.,   -9.,    2.,    1., 1960
     b              0.,    4.,   -3.,   -1.,    9.,   -2.,    8.,    3., 1960
     c              0.,   -1.,    5.,    1.,   -3.,    4.,    4.,    1., 1960
     d              0.,    0.,   -1.,    2.,    4.,   -5.,    6.,    1., 1960
     e              1.,   -1.,   -1.,    6.,    2.,    0.,    0.,   -7./ 1960
      data gd/ -30334.,-2119., 5776.,-1662., 2997.,-2016., 1594.,  114., 1965
     1           1297.,-2038., -404., 1292.,  240.,  856., -165.,  957., 1965
     2            804.,  148.,  479., -269., -390.,   13.,  252., -269., 1965
     3           -219.,  358.,   19.,  254.,  128.,  -31., -126., -157., 1965
     4            -97.,  -62.,   81.,   45.,   61.,  -11.,    8.,  100., 1965
     5           -228.,   68.,    4.,  -32.,    1.,   -8., -111.,   -7., 1965
     6             75.,  -57.,  -61.,    4.,  -27.,   13.,   -2.,  -26., 1965
     7              6.,   -6.,   26.,   13.,  -23.,    1.,  -12.,   13., 1965
     8              5.,    7.,   -4.,  -12.,  -14.,    9.,    0.,  -16., 1965
     9              8.,    4.,   -1.,   24.,   11.,   -3.,    4.,  -17., 1965
     a              8.,   10.,  -22.,    2.,   15.,  -13.,    7.,   10., 1965
     b             -4.,   -1.,   -5.,   -1.,   10.,    5.,   10.,    1., 1965
     c             -4.,   -2.,    1.,   -2.,   -3.,    2.,    2.,    1., 1965
     d             -5.,    2.,   -2.,    6.,    4.,   -4.,    4.,    0., 1965
     e              0.,   -2.,    2.,    3.,    2.,    0.,    0.,   -6./ 1965
      data ge/ -30220.,-2068., 5737.,-1781., 3000.,-2047., 1611.,   25., 1970
     1           1287.,-2091., -366., 1278.,  251.,  838., -196.,  952., 1970
     2            800.,  167.,  461., -266., -395.,   26.,  234., -279., 1970
     3           -216.,  359.,   26.,  262.,  139.,  -42., -139., -160., 1970
     4            -91.,  -56.,   83.,   43.,   64.,  -12.,   15.,  100., 1970
     5           -212.,   72.,    2.,  -37.,    3.,   -6., -112.,    1., 1970
     6             72.,  -57.,  -70.,    1.,  -27.,   14.,   -4.,  -22., 1970
     7              8.,   -2.,   23.,   13.,  -23.,   -2.,  -11.,   14., 1970
     8              6.,    7.,   -2.,  -15.,  -13.,    6.,   -3.,  -17., 1970
     9              5.,    6.,    0.,   21.,   11.,   -6.,    3.,  -16., 1970
     a              8.,   10.,  -21.,    2.,   16.,  -12.,    6.,   10., 1970
     b             -4.,   -1.,   -5.,    0.,   10.,    3.,   11.,    1., 1970
     c             -2.,   -1.,    1.,   -3.,   -3.,    1.,    2.,    1., 1970
     d             -5.,    3.,   -1.,    4.,    6.,   -4.,    4.,    0., 1970
     e              1.,   -1.,    0.,    3.,    3.,    1.,   -1.,   -4./ 1970
      data gf/ -30100.,-2013., 5675.,-1902., 3010.,-2067., 1632.,  -68., 1975
     1           1276.,-2144., -333., 1260.,  262.,  830., -223.,  946., 1975
     2            791.,  191.,  438., -265., -405.,   39.,  216., -288., 1975
     3           -218.,  356.,   31.,  264.,  148.,  -59., -152., -159., 1975
     4            -83.,  -49.,   88.,   45.,   66.,  -13.,   28.,   99., 1975
     5           -198.,   75.,    1.,  -41.,    6.,   -4., -111.,   11., 1975
     6             71.,  -56.,  -77.,    1.,  -26.,   16.,   -5.,  -14., 1975
     7             10.,    0.,   22.,   12.,  -23.,   -5.,  -12.,   14., 1975
     8              6.,    6.,   -1.,  -16.,  -12.,    4.,   -8.,  -19., 1975
     9              4.,    6.,    0.,   18.,   10.,  -10.,    1.,  -17., 1975
     a              7.,   10.,  -21.,    2.,   16.,  -12.,    7.,   10., 1975
     b             -4.,   -1.,   -5.,   -1.,   10.,    4.,   11.,    1., 1975
     c             -3.,   -2.,    1.,   -3.,   -3.,    1.,    2.,    1., 1975
     d             -5.,    3.,   -2.,    4.,    5.,   -4.,    4.,   -1., 1975
     e              1.,   -1.,    0.,    3.,    3.,    1.,   -1.,   -5./ 1975
      data gg/ -29992.,-1956., 5604.,-1997., 3027.,-2129., 1663., -200., 1980
     1           1281.,-2180., -336., 1251.,  271.,  833., -252.,  938., 1980
     2            782.,  212.,  398., -257., -419.,   53.,  199., -297., 1980
     3           -218.,  357.,   46.,  261.,  150.,  -74., -151., -162., 1980
     4            -78.,  -48.,   92.,   48.,   66.,  -15.,   42.,   93., 1980
     5           -192.,   71.,    4.,  -43.,   14.,   -2., -108.,   17., 1980
     6             72.,  -59.,  -82.,    2.,  -27.,   21.,   -5.,  -12., 1980
     7             16.,    1.,   18.,   11.,  -23.,   -2.,  -10.,   18., 1980
     8              6.,    7.,    0.,  -18.,  -11.,    4.,   -7.,  -22., 1980
     9              4.,    9.,    3.,   16.,    6.,  -13.,   -1.,  -15., 1980
     a              5.,   10.,  -21.,    1.,   16.,  -12.,    9.,    9., 1980
     b             -5.,   -3.,   -6.,   -1.,    9.,    7.,   10.,    2., 1980
     c             -6.,   -5.,    2.,   -4.,   -4.,    1.,    2.,    0., 1980
     d             -5.,    3.,   -2.,    6.,    5.,   -4.,    3.,    0., 1980
     e              1.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6./ 1980
      data gi/ -29873.,-1905., 5500.,-2072., 3044.,-2197., 1687., -306., 1985
     1           1296.,-2208., -310., 1247.,  284.,  829., -297.,  936., 1985
     2            780.,  232.,  361., -249., -424.,   69.,  170., -297., 1985
     3           -214.,  355.,   47.,  253.,  150.,  -93., -154., -164., 1985
     4            -75.,  -46.,   95.,   53.,   65.,  -16.,   51.,   88., 1985
     5           -185.,   69.,    4.,  -48.,   16.,   -1., -102.,   21., 1985
     6             74.,  -62.,  -83.,    3.,  -27.,   24.,   -2.,   -6., 1985
     7             20.,    4.,   17.,   10.,  -23.,    0.,   -7.,   21., 1985
     8              6.,    8.,    0.,  -19.,  -11.,    5.,   -9.,  -23., 1985
     9              4.,   11.,    4.,   14.,    4.,  -15.,   -4.,  -11., 1985
     a              5.,   10.,  -21.,    1.,   15.,  -12.,    9.,    9., 1985
     b             -6.,   -3.,   -6.,   -1.,    9.,    7.,    9.,    1., 1985
     c             -7.,   -5.,    2.,   -4.,   -4.,    1.,    3.,    0., 1985
     d             -5.,    3.,   -2.,    6.,    5.,   -4.,    3.,    0., 1985
     e              1.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6./ 1985
      data gj/ -29775.,-1848., 5406.,-2131., 3059.,-2279., 1686., -373., 1990
     1           1314.,-2239., -284., 1248.,  293.,  802., -352.,  939., 1990
     2            780.,  247.,  325., -240., -423.,   84.,  141., -299., 1990
     3           -214.,  353.,   46.,  245.,  154., -109., -153., -165., 1990
     4            -69.,  -36.,   97.,   61.,   65.,  -16.,   59.,   82., 1990
     5           -178.,   69.,    3.,  -52.,   18.,    1.,  -96.,   24., 1990
     6             77.,  -64.,  -80.,    2.,  -26.,   26.,    0.,   -1., 1990
     7             21.,    5.,   17.,    9.,  -23.,    0.,   -4.,   23., 1990
     8              5.,   10.,   -1.,  -19.,  -10.,    6.,  -12.,  -22., 1990
     9              3.,   12.,    4.,   12.,    2.,  -16.,   -6.,  -10., 1990
     a              4.,    9.,  -20.,    1.,   15.,  -12.,   11.,    9., 1990
     b             -7.,   -4.,   -7.,   -2.,    9.,    7.,    8.,    1., 1990
     c             -7.,   -6.,    2.,   -3.,   -4.,    2.,    2.,    1., 1990
     d             -5.,    3.,   -2.,    6.,    4.,   -4.,    3.,    0., 1990
     e              1.,   -2.,    3.,    3.,    3.,   -1.,    0.,   -6./ 1990
      data gk/ -29692.,-1784., 5306.,-2200., 3070.,-2366., 1681., -413., 1995
     1           1335.,-2267., -262., 1249.,  302.,  759., -427.,  940., 1995
     2            780.,  262.,  290., -236., -418.,   97.,  122., -306., 1995
     3           -214.,  352.,   46.,  235.,  165., -118., -143., -166., 1995
     4            -55.,  -17.,  107.,   68.,   67.,  -17.,   68.,   72., 1995
     5           -170.,   67.,   -1.,  -58.,   19.,    1.,  -93.,   36., 1995
     6             77.,  -72.,  -69.,    1.,  -25.,   28.,    4.,    5., 1995
     7             24.,    4.,   17.,    8.,  -24.,   -2.,   -6.,   25., 1995
     8              6.,   11.,   -6.,  -21.,   -9.,    8.,  -14.,  -23., 1995
     9              9.,   15.,    6.,   11.,   -5.,  -16.,   -7.,   -4., 1995
     a              4.,    9.,  -20.,    3.,   15.,  -10.,   12.,    8., 1995
     b             -6.,   -8.,   -8.,   -1.,    8.,   10.,    5.,   -2., 1995
     c             -8.,   -8.,    3.,   -3.,   -6.,    1.,    2.,    0., 1995
     d             -4.,    4.,   -1.,    5.,    4.,   -5.,    2.,   -1., 1995
     e              2.,   -2.,    5.,    1.,    1.,   -2.,    0.,   -7., 1995
     f           75*0./                                                  1995
      data gl/ -29619.4,-1728.2, 5186.1,-2267.7, 3068.4,-2481.6, 1670.9, 2000
     1           -458.0, 1339.6,-2288.0, -227.6, 1252.1,  293.4,  714.5, 2000
     2           -491.1,  932.3,  786.8,  272.6,  250.0, -231.9, -403.0, 2000
     3            119.8,  111.3, -303.8, -218.8,  351.4,   43.8,  222.3, 2000
     4            171.9, -130.4, -133.1, -168.6,  -39.3,  -12.9,  106.3, 2000
     5             72.3,   68.2,  -17.4,   74.2,   63.7, -160.9,   65.1, 2000
     6             -5.9,  -61.2,   16.9,    0.7,  -90.4,   43.8,   79.0, 2000
     7            -74.0,  -64.6,    0.0,  -24.2,   33.3,    6.2,    9.1, 2000
     8             24.0,    6.9,   14.8,    7.3,  -25.4,   -1.2,   -5.8, 2000
     9             24.4,    6.6,   11.9,   -9.2,  -21.5,   -7.9,    8.5, 2000
     a            -16.6,  -21.5,    9.1,   15.5,    7.0,    8.9,   -7.9, 2000
     b            -14.9,   -7.0,   -2.1,    5.0,    9.4,  -19.7,    3.0, 2000
     c             13.4,   -8.4,   12.5,    6.3,   -6.2,   -8.9,   -8.4, 2000
     d             -1.5,    8.4,    9.3,    3.8,   -4.3,   -8.2,   -8.2, 2000
     e              4.8,   -2.6,   -6.0,    1.7,    1.7,    0.0,   -3.1, 2000
     f              4.0,   -0.5,    4.9,    3.7,   -5.9,    1.0,   -1.2, 2000
     g              2.0,   -2.9,    4.2,    0.2,    0.3,   -2.2,   -1.1, 2000
     h             -7.4,    2.7,   -1.7,    0.1,   -1.9,    1.3,    1.5, 2000
     i             -0.9,   -0.1,   -2.6,    0.1,    0.9,   -0.7,   -0.7, 2000
     j              0.7,   -2.8,    1.7,   -0.9,    0.1,   -1.2,    1.2, 2000
     k             -1.9,    4.0,   -0.9,   -2.2,   -0.3,   -0.4,    0.2, 2000
     l              0.3,    0.9,    2.5,   -0.2,   -2.6,    0.9,    0.7, 2000
     m             -0.5,    0.3,    0.3,    0.0,   -0.3,    0.0,   -0.4, 2000
     n              0.3,   -0.1,   -0.9,   -0.2,   -0.4,   -0.4,    0.8, 2000
     o             -0.2,   -0.9,   -0.9,    0.3,    0.2,    0.1,    1.8, 2000
     p             -0.4,   -0.4,    1.3,   -1.0,   -0.4,   -0.1,    0.7, 2000
     q              0.7,   -0.4,    0.3,    0.3,    0.6,   -0.1,    0.3, 2000
     r              0.4,   -0.2,    0.0,   -0.5,    0.1,   -0.9/         2000
      data gm/ -29556.8,-1671.8, 5080.0,-2340.5, 3047.0,-2594.9, 1656.9, 2005
     1           -516.7, 1335.7,-2305.3, -200.4, 1246.8,  269.3,  674.4, 2005
     2           -524.5,  919.8,  798.2,  281.4,  211.5, -225.8, -379.5, 2005
     3            145.7,  100.2, -304.7, -227.6,  354.4,   42.7,  208.8, 2005
     4            179.8, -136.6, -123.0, -168.3,  -19.5,  -14.1,  103.6, 2005
     5             72.9,   69.6,  -20.2,   76.6,   54.7, -151.1,   63.7, 2005
     6            -15.0,  -63.4,   14.7,    0.0,  -86.4,   50.3,   79.8, 2005
     7            -74.4,  -61.4,   -1.4,  -22.5,   38.6,    6.9,   12.3, 2005
     8             25.4,    9.4,   10.9,    5.5,  -26.4,    2.0,   -4.8, 2005
     9             24.8,    7.7,   11.2,  -11.4,  -21.0,   -6.8,    9.7, 2005
     a            -18.0,  -19.8,   10.0,   16.1,    9.4,    7.7,  -11.4, 2005
     b            -12.8,   -5.0,   -0.1,    5.6,    9.8,  -20.1,    3.6, 2005
     c             12.9,   -7.0,   12.7,    5.0,   -6.7,  -10.8,   -8.1, 2005
     d             -1.3,    8.1,    8.7,    2.9,   -6.7,   -7.9,   -9.2, 2005
     e              5.9,   -2.2,   -6.3,    2.4,    1.6,    0.2,   -2.5, 2005
     f              4.4,   -0.1,    4.7,    3.0,   -6.5,    0.3,   -1.0, 2005
     g              2.1,   -3.4,    3.9,   -0.9,   -0.1,   -2.3,   -2.2, 2005
     h             -8.0,    2.9,   -1.6,    0.3,   -1.7,    1.4,    1.5, 2005
     i             -0.7,   -0.2,   -2.4,    0.2,    0.9,   -0.7,   -0.6, 2005
     j              0.5,   -2.7,    1.8,   -1.0,    0.1,   -1.5,    1.0, 2005
     k             -2.0,    4.1,   -1.4,   -2.2,   -0.3,   -0.5,    0.3, 2005
     l              0.3,    0.9,    2.3,   -0.4,   -2.7,    1.0,    0.6, 2005
     m             -0.4,    0.4,    0.5,    0.0,   -0.3,    0.0,   -0.4, 2005
     n              0.3,    0.0,   -0.8,   -0.4,   -0.4,    0.0,    1.0, 2005
     o             -0.2,   -0.9,   -0.7,    0.3,    0.3,    0.3,    1.7, 2005
     p             -0.4,   -0.5,    1.2,   -1.0,   -0.4,    0.0,    0.7, 2005
     q              0.7,   -0.3,    0.2,    0.4,    0.6,   -0.1,    0.4, 2005
     r              0.4,   -0.2,   -0.1,   -0.5,   -0.3,   -1.0/         2005
      data gp/     8.8,  10.8, -21.3, -15.0,  -6.9, -23.3,  -1.0, -14.0, 2007
     1            -0.3,  -3.1,   5.4,  -0.9,  -6.5,  -6.8,  -2.0,  -2.5, 2007
     2             2.8,   2.0,  -7.1,   1.8,   5.9,   5.6,  -3.2,   0.0, 2007
     3            -2.6,   0.4,   0.1,  -3.0,   1.8,  -1.2,   2.0,   0.2, 2007
     4             4.5,  -0.6,  -1.0,  -0.8,   0.2,  -0.4,  -0.2,  -1.9, 2007
     5             2.1,  -0.4,  -2.1,  -0.4,  -0.4,  -0.2,   1.3,   0.9, 2007
     6            -0.4,   0.0,   0.8,  -0.2,   0.4,   1.1,   0.1,   0.6, 2007
     7             0.2,   0.4,  -0.9,  -0.5,  -0.3,   0.9,   0.3,  -0.2, 2007
     8             0.2,  -0.2,  -0.2,   0.2,   0.2,   0.2,  -0.2,   0.4, 2007
     9             0.2,   0.2,   0.5,  -0.3,  -0.7,   0.5,   0.5,   0.4, 2007
     a         115*0.0/                                                  2007
c
c     set initial values
c
      x     = 0.0
      y     = 0.0
      z     = 0.0
      if (date.lt.1900.0.or.date.gt.2015.0) go to 11
      if (date.gt.2010.0) write (6,960) date
  960 format (/' This version of the IGRF is intended for use up',
     1        ' to 2010.0.'/' values for',f9.3,' will be computed',
     2        ' but may be of reduced accuracy'/)
      if (date.ge.2005.0) go to 1
      t     = 0.2*(date - 1900.0)                                             
      ll    = t
      one   = ll
      t     = t - one
      if (date.lt.1995.0) then
       nmx   = 10
       nc    = nmx*(nmx+2)
       ll    = nc*ll
       kmx   = (nmx+1)*(nmx+2)/2
      else
       nmx   = 13
       nc    = nmx*(nmx+2)
       ll    = 0.2*(date - 1995.0)
       ll    = 120*19 + nc*ll
       kmx   = (nmx+1)*(nmx+2)/2
      endif
      tc    = 1.0 - t
      if (isv.eq.1) then
       tc = -0.2
       t = 0.2
      end if
      go to 2
c
    1 t     = date - 2005.0
      tc    = 1.0
      if (isv.eq.1) then
       t = 1.0
       tc = 0.0
      end if
      ll    = 2670
      nmx   = 13
      nc    = nmx*(nmx+2)
      kmx   = (nmx+1)*(nmx+2)/2
    2 r     = alt
      one   = colat*0.017453292
      ct    = cos(one)
      st    = sin(one)
      one   = elong*0.017453292
      cl(1) = cos(one)
      sl(1) = sin(one)
      cd    = 1.0
      sd    = 0.0
      l     = 1
      m     = 1
      n     = 0
      if (itype.eq.2) go to 3
c
c     conversion from geodetic to geocentric coordinates 
c     (using the WGS84 spheroid)
c 
      a2    = 40680631.6
      b2    = 40408296.0
      one   = a2*st*st
      two   = b2*ct*ct
      three = one + two
      rho   = sqrt(three)
      r     = sqrt(alt*(alt + 2.0*rho) + (a2*one + b2*two)/three)
      cd    = (alt + rho)/r
      sd    = (a2 - b2)/rho*ct*st/r
      one   = ct
      ct    = ct*cd -  st*sd
      st    = st*cd + one*sd
c
    3 ratio = 6371.2/r
      rr    = ratio*ratio
c
c     computation of Schmidt quasi-normal coefficients p and x(=q)
c
      p(1)  = 1.0
      p(3)  = st
      q(1)  = 0.0
      q(3)  =  ct
      do 10 k=2,kmx                                                       
      if (n.ge.m) go to 4
      m     = 0
      n     = n + 1
      rr    = rr*ratio
      fn    = n
      gn    = n - 1
    4 fm    = m
      if (m.ne.n) go to 5
      if (k.eq.3) go to 6
      one   = sqrt(1.0 - 0.5/fm)
      j     = k - n - 1
      p(k)  = one*st*p(j)
      q(k)  = one*(st*q(j) + ct*p(j))
      cl(m) = cl(m-1)*cl(1) - sl(m-1)*sl(1)
      sl(m) = sl(m-1)*cl(1) + cl(m-1)*sl(1)
      go to 6                                                           
    5 gmm    = m*m
      one   = sqrt(fn*fn - gmm)
      two   = sqrt(gn*gn - gmm)/one
      three = (fn + gn)/one
      i     = k - n
      j     = i - n + 1
      p(k)  = three*ct*p(i) - two*p(j)
      q(k)  = three*(ct*q(i) - st*p(i)) - two*q(j)
c
c     synthesis of x, y and z in geocentric coordinates
c 
    6 lm    = ll + l
      one   = (tc*gh(lm) + t*gh(lm+nc))*rr                                     
      if (m.eq.0) go to 9                                                      
      two   = (tc*gh(lm+1) + t*gh(lm+nc+1))*rr
      three = one*cl(m) + two*sl(m)
      x     = x + three*q(k)
      z     = z - (fn + 1.0)*three*p(k)
      if (st.eq.0.0) go to 7
      y     = y + (one*sl(m) - two*cl(m))*fm*p(k)/st
      go to 8
    7 y     = y + (one*sl(m) - two*cl(m))*q(k)*ct
    8 l     = l + 2
      go to 10
    9 x     = x + one*q(k)
      z     = z - (fn + 1.0)*one*p(k)
      l     = l + 1
   10 m     = m + 1
c
c     conversion to coordinate system specified by itype
c
      one   = x
      x     = x*cd +   z*sd
      z     = z*cd - one*sd
      f     = sqrt(x*x + y*y + z*z)
c 
      return
c
c     error return if date out of bounds
c
   11 f     = 1.0d8
      write (6,961) date
  961 format (/' This subroutine will not work with a date of',
     1        f9.3,'.  Date must be in the range 1900.0.ge.date',
     2        '.le.2015.0.  On return f = 1.0d8., x = y = z = 0.')
      return
      end
