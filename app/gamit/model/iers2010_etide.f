       SUBROUTINE IERS2010_etide(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE
     . ,UTCTAI)
*+
*  - - - - - - - - - - - - - - - 
*   D E H A N T T I D E I N E L
*  - - - - - - - - - - - - - - - 
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine computes the tidal corrections of station displacements
*  caused by lunar and solar gravitational attraction (see References). 
*  The computations are calculated by the following steps:
*
*  Step 1): General degree 2 and degree 3 corrections + CALL ST1IDIU 
*  + CALL ST1ISEM + CALL ST1L1.
*  
*  Step 2): CALL STEP2DIU + CALL STEP2LON
*
*  It has been decided that the Step 3 non-correction for permanent tide
*  would not be applied in order to avoid a jump in the reference frame.
*  This Step 3 must be added in order to get the non-tidal station position
*  and to conform with the IAG Resolution.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 1
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
*     XSUN          d(3)   Geocentric position of the Sun (Note 2)
*     XMON          d(3)   Geocentric position of the Moon (Note 2)
*     YR            i      Year (Note 3)
*     MONTH         i      Month (Note 3)
*     DAY           i      Day of Month (Note 3)
*     FHR           d      Hour in the day (Note 4)
* SCM 05/10/2018 modified to work with Model 
*     UTCTAI        d      TAI-UTC correction passed from etide.f
*                          
*
*  Returned:
*     DXTIDE        d(3)   Displacement vector (Note 5)
*
*  Notes:
*
*  1) The IGS station is in ITRF co-rotating frame.  All coordinates,
*     X, Y, and Z, are expressed in meters. 
*  
*  2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
*     coordinates are expressed in meters.
*
*  3) The values are expressed in Coordinated Universal Time (UTC).
*
*  4) The fractional hours in the day is computed as the hour + minutes/60.0
*     + sec/3600.0.  The unit is expressed in Universal Time (UT).
*
*  5) The displacement vector is in the geocentric ITRF.  All components are
*     expressed in meters.
*
*  Called:
*     SPROD             Finds the scalar product and unit vector of two vectors 
*     ZERO_VEC8         Returns the zero vector
*     ST1IDIU           Corrects for the out-of-phase part of Love numbers
*                       for the diurnal band
*     ST1ISEM           Same as above for the semi-diurnal band
*     ST1L1             Corrects for the latitude dependence of Love numbers
*     CAL2JD            Computes Julian Date from Gregorian calendar date
*     DAT               Computes the difference TAI-UTC
*     STEP2DIU          Computes in-phase and out-of-phase corrections in
*                       the diurnal band
*     STEP2LON          Same as above for the long period band
*
*  Test case:
*     given input: XSTA(1) = 4075578.385D0 meters
*                  XSTA(2) =  931852.890D0 meters
*                  XSTA(3) = 4801570.154D0 meters   
*                  XSUN(1) = 137859926952.015D0 meters
*                  XSUN(2) = 54228127881.4350D0 meters
*                  XSUN(3) = 23509422341.6960D0 meters
*                  XMON(1) = -179996231.920342D0 meters
*                  XMON(2) = -312468450.131567D0 meters
*                  XMON(3) = -169288918.592160D0 meters
*                  YR      = 2009
*                  MONTH   = 4
*                  DAY     = 13
*                  FHR     = 0.00D0 seconds 
*                  
*     expected output:  DXTIDE(1) = 0.7700420357108125891D-01 meters
*                       DXTIDE(2) = 0.6304056321824967613D-01 meters
*                       DXTIDE(3) = 0.5516568152597246810D-01 meters
*
*  References:
*
*     Groten, E., 2000, Geodesists Handbook 2000, Part 4,
*     http://www.gfy.ku.dk/~iag/HB2000/part4/groten.htm. See also
*     ''Parameters of Common Relevance of Astronomy, Geodesy, and
*     Geodynamics," J. Geod., 74, pp. 134-140
*
*     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
*     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*     
*     Pitjeva, E. and Standish, E. M., 2009, ''Proposals for the masses
*     of the three largest asteroids, the Moon-Earth mass ratio and the
*     Astronomical Unit," Celest. Mech. Dyn. Astr., 103, pp. 365-372
*
*     Ries, J. C., Eanes, R. J., Shum, C. K. and Watkins, M. M., 1992,
*     ''Progress in the Determination of the Gravitational Coefficient
*     of the Earth," Geophys. Res. Lett., 19(6), pp. 529-531
*
*  Revisions:
*  1996 March    23 V. Dehant      Original code
*                   P. M. Mathews
*                   J. Gipson
*  2000 May      17 V. Dehant      Last modifications
*                   P. M. Mathews
*  2006 February 06 J. Ray         Header comments modified to clarify
*                                  input/output units and systems
*  2006 February 06 J. Ray         Subroutine DUTC modified for leap
*                                  second on 2006.0 and to correct 
*                                  do 5 i=1,87 from 84 to 87
*  2006 August   31 G. Petit       Correct DUTC for dates after 2007
*  2007 June     20 H. Manche      Modified DUTC to correct past mistake
*                                  and corrected DE line in subroutine
*                                  STEP2DIU
*  2007 October  23 H. Manche      Replace subroutines DUTC and FJLDY with
*                   G. Petit       SOFA subroutines iau_CAL2JD and iau_DAT
*                                  and correct time arguments of subroutine
*                                  STEP2DIU
*  2009 February 19 G. Petit       Update routine iau_DAT for 2009.0 leap
*                                  second
*  2009 August   06 B.E. Stetzler  Initial standardization of code 
*  2009 August   07 B.E. Stetzler  Updated MASS_RATIO_SUN, 
*                                  MASS_RATIO_MOON and RE to CBEs and
*                                  provided a test case
*  2009 August  07  B.E. Stetzler  Capitalized all variables for Fortran
*                                  77 compatibility
*  2009 September 01 B.E. Stetzler Removed 'iau_' from redistributed SOFA
*                                  subroutines
*-----------------------------------------------------------------------

      IMPLICIT NONE
      
      INTEGER STATUT,YR,MONTH,DAY,I,J
      DOUBLE PRECISION XSTA(3),XSUN(3),XMON(3),DXTIDE(3),XCORSTA(3),
     .                 FHR,H20,L20,H3,L3,H2,L2,SCS,RSTA,SCM,RSUN,RMON,
     .                 SCSUN,SCMON,COSPHI,P2SUN,P2MON,P3SUN,P3MON,
     .                 X2SUN,X2MON,X3SUN,X3MON,MASS_RATIO_SUN,
     .                 MASS_RATIO_MOON,RE,FAC2SUN,FAC2MON,FAC3SUN,
     .                 FAC3MON,JJM0,JJM1,DTT,T,PI,SINPHI,COSLA,SINLA,
     .                 DR,DN,UTCTAI
      PARAMETER ( PI = 3.1415926535897932384626433D0 ) 
*----------------------------------------------------------------------  
* NOMINAL SECOND DEGREE AND THIRD DEGREE LOVE NUMBERS AND SHIDA NUMBERS  
*----------------------------------------------------------------------  
      DATA H20/0.6078D0/,L20/0.0847D0/,H3/0.292D0/,L3/0.015D0/
*----------------------------------------------------------------------  
* SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR  
*----------------------------------------------------------------------  
      CALL SPROD_iers2010(XSTA,XSUN,SCS,RSTA,RSUN)  
      CALL SPROD_iers2010(XSTA,XMON,SCM,RSTA,RMON)  
      SCSUN=SCS/RSTA/RSUN  
      SCMON=SCM/RSTA/RMON
*----------------------------------------------------------------------   
* COMPUTATION OF NEW H2 AND L2  
*----------------------------------------------------------------------  
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      H2=H20-0.0006D0*(1D0-3D0/2D0*COSPHI**2)
      L2=L20+0.0002D0*(1D0-3D0/2D0*COSPHI**2)  

* P2 term  
      P2SUN=3D0*(H2/2D0-L2)*SCSUN**2-H2/2D0  
      P2MON=3D0*(H2/2D0-L2)*SCMON**2-H2/2D0  

* P3 term  
      P3SUN=5D0/2D0*(H3-3D0*L3)*SCSUN**3+3D0/2D0*(L3-H3)*SCSUN
      P3MON=5D0/2D0*(H3-3D0*L3)*SCMON**3+3D0/2D0*(L3-H3)*SCMON

*----------------------------------------------------------------------  
* TERM IN DIRECTION OF SUN/MOON VECTOR  
*----------------------------------------------------------------------  
      X2SUN=3D0*L2*SCSUN  
      X2MON=3D0*L2*SCMON  
      X3SUN=3D0*L3/2D0*(5D0*SCSUN**2-1D0)  
      X3MON=3D0*L3/2D0*(5D0*SCMON**2-1D0)
*----------------------------------------------------------------------  
* FACTORS FOR SUN/MOON USING IAU CURRENT BEST ESTIMATES (SEE REFERENCES) 
*----------------------------------------------------------------------  
      MASS_RATIO_SUN=332946.0482D0
      MASS_RATIO_MOON=0.0123000371D0
      RE=6378136.6D0
      FAC2SUN=MASS_RATIO_SUN*RE*(RE/RSUN)**3
      FAC2MON=MASS_RATIO_MOON*RE*(RE/RMON)**3
      FAC3SUN=FAC2SUN*(RE/RSUN)
      FAC3MON=FAC2MON*(RE/RMON)
  
* TOTAL DISPLACEMENT  
      DO 10 I=1,3  
      DXTIDE(I)=FAC2SUN*( X2SUN*XSUN(I)/RSUN + P2SUN*XSTA(I)/RSTA ) +
     .          FAC2MON*( X2MON*XMON(I)/RMON + P2MON*XSTA(I)/RSTA ) +  
     .          FAC3SUN*( X3SUN*XSUN(I)/RSUN + P3SUN*XSTA(I)/RSTA ) +   
     .          FAC3MON*( X3MON*XMON(I)/RMON + P3MON*XSTA(I)/RSTA )  
10    CONTINUE  

      CALL ZERO_VEC8_iers2010(XCORSTA)
*+---------------------------------------------------------------------  
* CORRECTIONS FOR THE OUT-OF-PHASE PART OF LOVE NUMBERS (PART H_2^(0)I  
* AND L_2^(0)I )  
*----------------------------------------------------------------------

* FIRST, FOR THE DIURNAL BAND       

      CALL ST1IDIU_iers2010(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
      DO 11 I=1,3
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)  
11    CONTINUE
  
* SECOND, FOR THE SEMI-DIURNAL BAND       
 
      CALL ST1ISEM_iers2010(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
      DO 12 I=1,3
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
12    CONTINUE  
  
*+---------------------------------------------------------------------
* CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )  
*----------------------------------------------------------------------   
      CALL ST1L1_iers2010(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
      DO 13 I=1,3  
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
13    CONTINUE    
  
* CONSIDER CORRECTIONS FOR STEP 2  

*+---------------------------------------------------------------------  
* CORRECTIONS FOR THE DIURNAL BAND:  
* 
*  FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES 
*        
*   1) CALL THE SUBROUTINE COMPUTING THE JULIAN DATE 
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      CALL CAL2JD_iers2010( YR, MONTH, DAY, JJM0, JJM1, STATUT )
      T=((JJM0-2451545.0D0)+JJM1+FHR/24.0D0)/36525D0
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
*   2) CALL THE SUBROUTINE COMPUTING THE CORRECTION OF UTC TIME  
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
* Use UTCTAI passed from etide.f here
*      CALL DAT_iers2010( YR, MONTH, DAY, FHR, DTT, STATUT )
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DTT = UTCTAI + 32.184D0
*     CONVERSION OF T IN TT TIME
      T=T+DTT/(3600.0D0*24.0D0*36525D0)

*  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND
*  CORRECTIONS, (in-phase and out-of-phase frequency dependence):
  
      CALL STEP2DIU_iers2010(XSTA,FHR,T,XCORSTA)  
      DO 14 I=1,3  
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
14    CONTINUE  
  
*  CORRECTIONS FOR THE LONG-PERIOD BAND,
*  (in-phase and out-of-phase frequency dependence):  

      CALL STEP2LON_iers2010(XSTA,T,XCORSTA)
      DO 15 I=1,3  
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
15    CONTINUE  
          
CD      print*,'#I 10: FHR,UTCTAI,DTT,T,DXTIDE:',FHR,UTCTAI,DTT,T,
CD     . DXTIDE*1000.d0

* CONSIDER CORRECTIONS FOR STEP 3  
  
*----------------------------------------------------------------------
* UNCORRECT FOR THE PERMANENT TIDE  
*  
*      SINPHI=XSTA(3)/RSTA  
*      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
*      COSLA=XSTA(1)/COSPHI/RSTA  
*      SINLA=XSTA(2)/COSPHI/RSTA  
*      DR=-DSQRT(5D0/4D0/PI)*H2*0.31460D0*(3D0/2D0*SINPHI**2-0.5D0)
*      DN=-DSQRT(5D0/4D0/PI)*L2*0.31460D0*3D0*COSPHI*SINPHI
*      DXTIDE(1)=DXTIDE(1)-DR*COSLA*COSPHI+DN*COSLA*SINPHI
*      DXTIDE(2)=DXTIDE(2)-DR*SINLA*COSPHI+DN*SINLA*SINPHI  
*      DXTIDE(3)=DXTIDE(3)-DR*SINPHI      -DN*COSPHI
*-----------------------------------------------------------------------
      RETURN  

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
* 
      SUBROUTINE CAL2JD_iers2010 ( IY, IM, ID, DJM0, DJM, J )
*+
*  - - - - - - - - - - -
*   C A L 2 J D
*  - - - - - - - - - - -
*
*  Gregorian Calendar to Julian Date.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     IY,IM,ID    i     year, month, day in Gregorian calendar (Note 1)
*
*  Returned:
*     DJM0        d     MJD zero-point: always 2400000.5
*     DJM         d     Modified Julian Date for 0 hrs
*     J           i     status:
*                           0 = OK
*                          -1 = bad year   (Note 3: JD not computed)
*                          -2 = bad month  (JD not computed)
*                          -3 = bad day    (JD computed)
*
*  Notes:
*
*  1) The algorithm used is valid from -4800 March 1, but this
*     implementation rejects dates before -4799 January 1.
*
*  2) The Julian Date is returned in two pieces, in the usual SOFA
*     manner, which is designed to preserve time resolution.  The
*     Julian Date is available as a single number by adding DJM0 and
*     DJM.
*
*  3) In early eras the conversion is from the "Proleptic Gregorian
*     Calendar";  no account is taken of the date(s) of adoption of
*     the Gregorian Calendar, nor is the AD/BC numbering convention
*     observed.
*
*  Reference:
*
*     Explanatory Supplement to the Astronomical Almanac,
*     P. Kenneth Seidelmann (ed), University Science Books (1992),
*     Section 12.92 (p604).
*
*  This revision:  2001 September 16
*
*  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER IY, IM, ID
      DOUBLE PRECISION DJM0, DJM
      INTEGER J, MY, IYPMY

*  Earliest year allowed (4800BC)
      INTEGER IYMIN
      PARAMETER ( IYMIN = -4799 )

*  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Preset status.
      J = 0

*  Validate year.
      IF ( IY.LT.IYMIN ) THEN
         J = -1
      ELSE

*     Validate month.
         IF ( IM.GE.1 .AND. IM.LE.12 ) THEN

*        Allow for leap year.
            IF ( MOD(IY,4) .EQ. 0 ) THEN
               MTAB(2) = 29
            ELSE
               MTAB(2) = 28
            END IF
            IF ( MOD(IY,100).EQ.0 .AND. MOD(IY,400).NE.0 ) MTAB(2) = 28

*        Validate day.
            IF ( ID.LT.1 .OR. ID.GT.MTAB(IM) ) J = -3

*        Result.
            MY = ( IM - 14 ) / 12
            IYPMY = IY + MY
            DJM0 = 2400000.5D0
            DJM = DBLE( ( 1461 * ( IYPMY + 4800 ) ) / 4
     :                + (  367 * ( IM-2 - 12*MY ) ) / 12
     :                - (    3 * ( ( IYPMY + 4900 ) / 100 ) ) / 4
     :                + ID - 2432076)

*        Bad month
         ELSE
            J = -2
         END IF
      END IF

*  Finished.

*+-----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and 
*     restrictions listed below.
*
*  3. You (the user) may copy and adapt the SOFA software and its 
*     algorithms for your own purposes and you may copy and distribute
*     a resulting "derived work" to others on a world-wide, royalty-free 
*     basis, provided that the derived work complies with the following
*     requirements: 
*
*     a) Your work shall be marked or carry a statement that it (i) uses
*        routines and computations derived by you from software provided 
*        by SOFA under license to you; and (ii) does not contain
*        software provided by SOFA or software that has been distributed
*        by or endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon and/or differs from the
*        original SOFA software.
*
*     c) The name(s) of all routine(s) that you distribute shall differ
*        from the SOFA names, even when the SOFA content has not been
*        otherwise changed.
*
*     d) The routine-naming prefix "iau" shall not be used.
*
*     e) The origin of the SOFA components of your derived work must not
*        be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     f) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have granted 
*        a further right to modify the source code of your derived work.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall acknowledge
*     that the SOFA software was used in obtaining those results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  6. The SOFA software is provided "as is" and the Board makes no 
*     warranty as to its use or performance.   The Board does not and 
*     cannot warrant the performance or results which the user may obtain 
*     by using the SOFA software.  The Board makes no warranties, express 
*     or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms 
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*-----------------------------------------------------------------------

      END
*
      SUBROUTINE DAT_iers2010 ( IY, IM, ID, FD, DELTAT, J )
*+
*  - - - - - - - -
*   D A T
*  - - - - - - - -
*
*  For a given UTC date, calculate delta(AT) = TAI-UTC.
*
*     :------------------------------------------:
*     :                                          :
*     :                 IMPORTANT                :
*     :                                          :
*     :  A new version of this routine must be   :
*     :  produced whenever a new leap second is  :
*     :  announced.  There are five items to     :
*     :  change on each such occasion:           :
*     :                                          :
*     :  1) The parameter NDAT must be           :
*     :     increased by 1.                      :
*     :                                          :
*     :  2) A new line must be added to the set  :
*     :     of DATA statements that initialize   :
*     :     the arrays IDATE and DATS.           :
*     :                                          :
*     :  3) The parameter IYV must be set to     :
*     :     the current year.                    :
*     :                                          :
*     :  4) The "Latest leap second" comment     :
*     :     below must be set to the new leap    :
*     :     second date.                         :
*     :                                          :
*     :  5) The "This revision" comment, later,  :
*     :     must be set to the current date.     :
*     :                                          :
*     :  Change (3) must also be carried out     :
*     :  whenever the routine is re-issued,      :
*     :  even if no leap seconds have been       :
*     :  added.                                  :
*     :                                          :
*     :  Latest leap second:  2008 December 31   :
*     :                                          :
*     :__________________________________________:
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     IY       i     UTC:  year (Notes 1 and 2)
*     IM       i           month (Note 2)
*     ID       i           day (Notes 2 and 3)
*     FD       d           fraction of day (Note 4)
*
*  Returned:
*     DELTAT   d     TAI minus UTC, seconds
*     J        i     status (Note 5):
*                       1 = dubious year (Note 1)
*                       0 = OK
*                      -1 = bad year
*                      -2 = bad month
*                      -3 = bad day (Note 3)
*                      -4 = bad fraction (Note 4)
*
*  Notes:
*
*  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
*     to call the routine with an earlier date.  If this is attempted,
*     zero is returned together with a warning status.
*
*     Because leap seconds cannot, in principle, be predicted in
*     advance, a reliable check for dates beyond the valid range is
*     impossible.  To guard against gross errors, a year five or more
*     after the release year of the present routine (see parameter IYV)
*     is considered dubious.  In this case a warning status is returned
*     but the result is computed in the normal way.
*
*     For both too-early and too-late years, the warning status is J=+1.
*     This is distinct from the error status J=-1, which signifies a
*     year so early that JD could not be computed.
*
*  2) If the specified date is for a day which ends with a leap second,
*     the UTC-TAI value returned is for the period leading up to the
*     leap second.  If the date is for a day which begins as a leap
*     second ends, the UTC-TAI returned is for the period following the
*     leap second.
*
*  3) The day number must be in the normal calendar range, for example
*     1 through 30 for April.  The "almanac" convention of allowing
*     such dates as January 0 and December 32 is not supported in this
*     routine, in order to avoid confusion near leap seconds.
*
*  4) The fraction of day is used only for dates before the introduction
*     of leap seconds, the first of which occurred at the end of 1971.
*     It is tested for validity (zero to less than 1 is the valid range)
*     even if not used;  if invalid, zero is used and status J=-4 is
*     returned.  For many applications, setting FD to zero is
*     acceptable;  the resulting error is always less than 3 ms (and
*     occurs only pre-1972).
*
*  5) The status value returned in the case where there are multiple
*     errors refers to the first error detected.  For example, if the
*     month and day are 13 and 32 respectively, J=-2 (bad month) will be
*     returned.
*
*  6) In cases where a valid result is not available, zero is returned.
*
*  References:
*
*  1) For dates from 1961 January 1 onwards, the expressions from the
*     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
*
*  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
*     the 1992 Explanatory Supplement.
*
*  Called:
*     CAL2JD       Gregorian calendar to Julian Day number
*
*  This revision:  2008 July 5
*
*  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER IY, IM, ID
      DOUBLE PRECISION FD, DELTAT
      INTEGER J

*  Release year for this version of DAT
      INTEGER IYV
      PARAMETER ( IYV = 2009 )

*  Number of Delta(AT) changes (increase by 1 for each new leap second)
      INTEGER NDAT
      PARAMETER ( NDAT = 39 )

*  Number of Delta(AT) expressions before leap seconds were introduced
      INTEGER NERA1
      PARAMETER ( NERA1 = 14 )

*  Dates (year, month) on which new Delta(AT) came into force
      INTEGER IDATE(2,NDAT)

*  New Delta(AT) which came into force on the given dates
      DOUBLE PRECISION DATS(NDAT)

*  Reference dates (MJD) and drift rates (s/day), pre leap seconds
      DOUBLE PRECISION DRIFT(2,NERA1)

*  Miscellaneous local variables
      LOGICAL MORE
      INTEGER JS, I, M, IS
      DOUBLE PRECISION DA, DJM0, DJM

*  Dates and Delta(AT)s
      DATA (IDATE(I, 1),I=1,2),DATS(1)  / 1960,  1,  1.4178180D0 /
      DATA (IDATE(I, 2),I=1,2),DATS(2)  / 1961,  1,  1.4228180D0 /
      DATA (IDATE(I, 3),I=1,2),DATS(3)  / 1961,  8,  1.3728180D0 /
      DATA (IDATE(I, 4),I=1,2),DATS(4)  / 1962,  1,  1.8458580D0 /
      DATA (IDATE(I, 5),I=1,2),DATS(5)  / 1963, 11,  1.9458580D0 /
      DATA (IDATE(I, 6),I=1,2),DATS(6)  / 1964,  1,  3.2401300D0 /
      DATA (IDATE(I, 7),I=1,2),DATS(7)  / 1964,  4,  3.3401300D0 /
      DATA (IDATE(I, 8),I=1,2),DATS(8)  / 1964,  9,  3.4401300D0 /
      DATA (IDATE(I, 9),I=1,2),DATS(9)  / 1965,  1,  3.5401300D0 /
      DATA (IDATE(I,10),I=1,2),DATS(10) / 1965,  3,  3.6401300D0 /
      DATA (IDATE(I,11),I=1,2),DATS(11) / 1965,  7,  3.7401300D0 /
      DATA (IDATE(I,12),I=1,2),DATS(12) / 1965,  9,  3.8401300D0 /
      DATA (IDATE(I,13),I=1,2),DATS(13) / 1966,  1,  4.3131700D0 /
      DATA (IDATE(I,14),I=1,2),DATS(14) / 1968,  2,  4.2131700D0 /
      DATA (IDATE(I,15),I=1,2),DATS(15) / 1972,  1, 10D0 /
      DATA (IDATE(I,16),I=1,2),DATS(16) / 1972,  7, 11D0 /
      DATA (IDATE(I,17),I=1,2),DATS(17) / 1973,  1, 12D0 /
      DATA (IDATE(I,18),I=1,2),DATS(18) / 1974,  1, 13D0 /
      DATA (IDATE(I,19),I=1,2),DATS(19) / 1975,  1, 14D0 /
      DATA (IDATE(I,20),I=1,2),DATS(20) / 1976,  1, 15D0 /
      DATA (IDATE(I,21),I=1,2),DATS(21) / 1977,  1, 16D0 /
      DATA (IDATE(I,22),I=1,2),DATS(22) / 1978,  1, 17D0 /
      DATA (IDATE(I,23),I=1,2),DATS(23) / 1979,  1, 18D0 /
      DATA (IDATE(I,24),I=1,2),DATS(24) / 1980,  1, 19D0 /
      DATA (IDATE(I,25),I=1,2),DATS(25) / 1981,  7, 20D0 /
      DATA (IDATE(I,26),I=1,2),DATS(26) / 1982,  7, 21D0 /
      DATA (IDATE(I,27),I=1,2),DATS(27) / 1983,  7, 22D0 /
      DATA (IDATE(I,28),I=1,2),DATS(28) / 1985,  7, 23D0 /
      DATA (IDATE(I,29),I=1,2),DATS(29) / 1988,  1, 24D0 /
      DATA (IDATE(I,30),I=1,2),DATS(30) / 1990,  1, 25D0 /
      DATA (IDATE(I,31),I=1,2),DATS(31) / 1991,  1, 26D0 /
      DATA (IDATE(I,32),I=1,2),DATS(32) / 1992,  7, 27D0 /
      DATA (IDATE(I,33),I=1,2),DATS(33) / 1993,  7, 28D0 /
      DATA (IDATE(I,34),I=1,2),DATS(34) / 1994,  7, 29D0 /
      DATA (IDATE(I,35),I=1,2),DATS(35) / 1996,  1, 30D0 /
      DATA (IDATE(I,36),I=1,2),DATS(36) / 1997,  7, 31D0 /
      DATA (IDATE(I,37),I=1,2),DATS(37) / 1999,  1, 32D0 /
      DATA (IDATE(I,38),I=1,2),DATS(38) / 2006,  1, 33D0 /
      DATA (IDATE(I,39),I=1,2),DATS(39) / 2009,  1, 34D0 /

*  Reference dates and drift rates
      DATA (DRIFT(I, 1),I=1,2) / 37300D0, 0.001296D0 /
      DATA (DRIFT(I, 2),I=1,2) / 37300D0, 0.001296D0 /
      DATA (DRIFT(I, 3),I=1,2) / 37300D0, 0.001296D0 /
      DATA (DRIFT(I, 4),I=1,2) / 37665D0, 0.0011232D0 /
      DATA (DRIFT(I, 5),I=1,2) / 37665D0, 0.0011232D0 /
      DATA (DRIFT(I, 6),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I, 7),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I, 8),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I, 9),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I,10),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I,11),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I,12),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I,13),I=1,2) / 39126D0, 0.002592D0 /
      DATA (DRIFT(I,14),I=1,2) / 39126D0, 0.002592D0 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Initialize the result to zero and the status to OK.
      DA = 0D0
      JS = 0

*  If invalid fraction of a day, set error status and give up.
      IF ( FD.LT.0D0 .OR. FD.GE.1D0 ) THEN
         JS = -4
         GO TO 9000
      END IF

*  Convert the date into an MJD.
      CALL CAL2JD_iers2010 ( IY, IM, ID, DJM0, DJM, JS )

*  If invalid year, month, or day, give up.
      IF ( JS .LT. 0 ) GO TO 9000

*  If pre-UTC year, set warning status and give up.
      IF ( IY .LT. IDATE(1,1) ) THEN
         JS = 1
         GO TO 9000
      END IF

*  If suspiciously late year, set warning status but proceed.
      IF ( IY .GT. IYV+5 ) JS = 1

*  Combine year and month.
      M = 12*IY+IM

*  Prepare to search the tables.
      MORE = .TRUE.

*  Find the most recent table entry.
      DO 1 I=NDAT,1,-1
         IF ( MORE ) THEN
            IS = I
            MORE = M .LT. ( 12*IDATE(1,I) + IDATE(2,I) )
         END IF
 1    CONTINUE

*  Get the Delta(AT).
      DA = DATS(IS)

*  If pre-1972, adjust for drift.
      IF ( IS .LE. NERA1 ) DA = DA +
     :                          ( DJM + FD - DRIFT(1,IS) ) * DRIFT(2,IS)

*  Return the Delta(AT) value and the status.
 9000 CONTINUE
      DELTAT = DA
      J = JS

*  Finished.

*+-----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and 
*     restrictions listed below.
*
*  3. You (the user) may copy and adapt the SOFA software and its 
*     algorithms for your own purposes and you may copy and distribute
*     a resulting "derived work" to others on a world-wide, royalty-free 
*     basis, provided that the derived work complies with the following
*     requirements: 
*
*     a) Your work shall be marked or carry a statement that it (i) uses
*        routines and computations derived by you from software provided 
*        by SOFA under license to you; and (ii) does not contain
*        software provided by SOFA or software that has been distributed
*        by or endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon and/or differs from the
*        original SOFA software.
*
*     c) The name(s) of all routine(s) that you distribute shall differ
*        from the SOFA names, even when the SOFA content has not been
*        otherwise changed.
*
*     d) The routine-naming prefix "iau" shall not be used.
*
*     e) The origin of the SOFA components of your derived work must not
*        be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     f) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have granted 
*        a further right to modify the source code of your derived work.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall acknowledge
*     that the SOFA software was used in obtaining those results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  6. The SOFA software is provided "as is" and the Board makes no 
*     warranty as to its use or performance.   The Board does not and 
*     cannot warrant the performance or results which the user may obtain 
*     by using the SOFA software.  The Board makes no warranties, express 
*     or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms 
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*-----------------------------------------------------------------------

      END
*
      SUBROUTINE SPROD_iers2010 (X,Y,SCAL,R1,R2) 
*+
*  - - - - - - - - - - -
*   S P R O D 
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine computes the scalar product of two vectors and 
*  their norms.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     X            d(3)      components of vector x 
*     Y            d(3)      components of vector y
*
*  Returned:
*     SCAL         d      scalar product of vector x and vector y
*     R1           d      length of vector x 
*     R2           d      length of vector y
*
*  Called:
*     None
*
*  Test case:
*     given input: X(1) = 2D0	Y(1) = 1D0
*                  X(2) = 2D0	Y(2) = 3D0
*                  X(3) = 3D0	Y(3) = 4D0
*     
*     expected output: SCAL = 20D0
*                      R1 = 4.123105625617660586D0
*                      R2 = 5.099019513592784492D0
*
*  References:
*
*     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
*     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 July 10 B.E.Stetzler Initial standardization of function,
*                            explicit exponential notation and
*                            provided a test case 
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION X(3), Y(3), R1, R2, SCAL

      R1=DSQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)) 
      R2=DSQRT(Y(1)*Y(1)+Y(2)*Y(2)+Y(3)*Y(3)) 
      SCAL=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3) 

      RETURN 

* Finished.
  
*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
*
      SUBROUTINE ST1IDIU_iers2010 (XSTA,XSUN,XMON,FAC2SUN,FAC2MON
     .,XCORSTA)
*+
*  - - - - - - - - - - -
*   S T 1 I D I U
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine gives the out-of-phase corrections induced by
*  mantle anelasticity in the diurnal band. 
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 1
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
*     XSUN          d(3)   Geocentric position of the Sun (Note 2)
*     XMON          d(3)   Geocentric position of the Moon (Note 2)
*     FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)      
*     FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3) 
*
*  Returned:
*     XCORSTA       d(3)   Out of phase station corrections for diurnal band
*
*  Notes:
*
*  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
*     expressed in meters. 
*  
*  2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
*     coordinates are expressed in meters.
*
*  3) The expressions are computed in the main program.  TGP is the tide
*     generated potential.  The units are inverse meters. 
*
*  Test case:
*     given input: XSTA(1) = 4075578.385D0 meters
*                  XSTA(2) =  931852.890D0 meters
*                  XSTA(3) = 4801570.154D0 meters   
*                  XSUN(1) = 137859926952.015D0 meters
*                  XSUN(2) = 54228127881.4350D0 meters
*                  XSUN(3) = 23509422341.6960D0 meters
*                  XMON(1) = -179996231.920342D0 meters
*                  XMON(2) = -312468450.131567D0 meters
*                  XMON(3) = -169288918.592160D0 meters
*                  FAC2SUN =  0.163271964478954D0 1/meters     
*                  FAC2MON =  0.321989090026845D0 1/meters    
*                  
*     expected output:  XCORSTA(1) = -0.2836337012840008001D-03 meters
*                       XCORSTA(2) =  0.1125342324347507444D-03 meters
*                       XCORSTA(3) = -0.2471186224343683169D-03 meters
*
*  References:
*
*     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
*     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  1996 March    23 V. Dehant      Original code
*  2009 July     30 B.E. Stetzler  Initial standardization of code 
*  2009 July     31 B.E. Stetzler  Provided a test case
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION NORM8,RSTA,SINPHI,COSPHI,COS2PHI,SINLA,
     .                 COSLA, RMON, RSUN, DRSUN, DRMON, DNSUN, DNMON,
     .                 DESUN, DEMON, DR, DN, DE, XSTA, XSUN, XMON,
     .                 XCORSTA, DHI, DLI, FAC2SUN, FAC2MON
      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)  
      DATA DHI/-0.0025D0/,DLI/-0.0007D0/  

* Compute the normalized position vector of the IGS station.
      RSTA = NORM8(XSTA)
      SINPHI = XSTA(3)/RSTA  
      COSPHI = DSQRT(XSTA(1)*XSTA(1)+XSTA(2)*XSTA(2))/RSTA
      COS2PHI = COSPHI*COSPHI-SINPHI*SINPHI
      SINLA = XSTA(2)/COSPHI/RSTA  
      COSLA = XSTA(1)/COSPHI/RSTA  
* Compute the normalized position vector of the Moon.
      RMON=NORM8(XMON)
* Compute the normalized position vector of the Sun.
      RSUN=NORM8(XSUN)

      DRSUN=-3D0*DHI*SINPHI*COSPHI*FAC2SUN*XSUN(3)*(XSUN(1)*
     .            SINLA-XSUN(2)*COSLA)/RSUN**2

      DRMON=-3D0*DHI*SINPHI*COSPHI*FAC2MON*XMON(3)*(XMON(1)*
     .            SINLA-XMON(2)*COSLA)/RMON**2

      DNSUN=-3D0*DLI*COS2PHI*FAC2SUN*XSUN(3)*(XSUN(1)*SINLA-
     .            XSUN(2)*COSLA)/RSUN**2

      DNMON=-3D0*DLI*COS2PHI*FAC2MON*XMON(3)*(XMON(1)*SINLA-
     .            XMON(2)*COSLA)/RMON**2

      DESUN=-3D0*DLI*SINPHI*FAC2SUN*XSUN(3)*
     . (XSUN(1)*COSLA+XSUN(2)*SINLA)/RSUN**2

      DEMON=-3D0*DLI*SINPHI*FAC2MON*XMON(3)*
     . (XMON(1)*COSLA+XMON(2)*SINLA)/RMON**2

      DR = DRSUN+DRMON 
      DN = DNSUN+DNMON  
      DE = DESUN+DEMON 

*  Compute the corrections for the station.
      XCORSTA(1)=DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA  
      XCORSTA(2)=DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA  
      XCORSTA(3)=DR*SINPHI+DN*COSPHI  

      RETURN 

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
*
      SUBROUTINE ST1ISEM_iers2010 (XSTA,XSUN,XMON,FAC2SUN,FAC2MON
     .,XCORSTA)
*+
*  - - - - - - - - - - -
*   S T 1 I S E M
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine gives the out-of-phase corrections induced by
*  mantle anelasticity in the semi-diurnal band. 
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 1
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
*     XSUN          d(3)   Geocentric position of the Sun (Note 2)
*     XMON          d(3)   Geocentric position of the Moon (Note 2)
*     FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)      
*     FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3) 
*
*  Returned:
*     XCORSTA       d(3)   Out of phase station corrections for
*                          semi-diurnal band
*
*  Notes:
*
*  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
*     expressed in meters. 
*  
*  2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
*     coordinates are expressed in meters.
*
*  3) The expressions are computed in the main program.  TGP is the tide
*     generated potential.  The units are inverse meters. 
*
*  Test case:
*     given input: XSTA(1) = 4075578.385D0 meters
*                  XSTA(2) =  931852.890D0 meters
*                  XSTA(3) = 4801570.154D0 meters   
*                  XSUN(1) = 137859926952.015D0 meters
*                  XSUN(2) = 54228127881.4350D0 meters
*                  XSUN(3) = 23509422341.6960D0 meters
*                  XMON(1) = -179996231.920342D0 meters
*                  XMON(2) = -312468450.131567D0 meters
*                  XMON(3) = -169288918.592160D0 meters
*                  FAC2SUN =  0.163271964478954D0 1/meters     
*                  FAC2MON =  0.321989090026845D0 1/meters    
*                  
*     expected output:  XCORSTA(1) = -0.2801334805106874015D-03 meters
*                       XCORSTA(2) =  0.2939522229284325029D-04 meters
*                       XCORSTA(3) = -0.6051677912316721561D-04 meters
*
*  References:
*
*     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
*     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  1996 March    23 V. Dehant      Original code
*  2009 July     31 B.E. Stetzler  Initial standardization of code 
*  2009 July     31 B.E. Stetzler  Provided a test case
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION NORM8,RSTA,SINPHI,COSPHI,COSTWOLA,
     .                 SINLA,COSLA,RMON,RSUN,DRSUN,DRMON,DNSUN,DNMON,
     .                 DESUN, DEMON, DR, DN, DE, XSTA, XSUN, XMON,
     .                 XCORSTA, DHI, DLI, FAC2SUN, FAC2MON, SINTWOLA
      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)  
      DATA DHI/-0.0022D0/,DLI/-0.0007D0/  

* Compute the normalized position vector of the IGS station.
      RSTA = NORM8(XSTA)
      SINPHI = XSTA(3)/RSTA  
      COSPHI = DSQRT(XSTA(1)*XSTA(1)+XSTA(2)*XSTA(2))/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA  
      COSLA=XSTA(1)/COSPHI/RSTA  
      COSTWOLA=COSLA*COSLA-SINLA*SINLA 
      SINTWOLA=2D0*COSLA*SINLA  
* Compute the normalized position vector of the Moon.
      RMON=NORM8(XMON)
* Compute the normalized position vector of the Sun.
      RSUN=NORM8(XSUN)

      DRSUN=-3D0/4D0*DHI*COSPHI**2*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*  
     . SINTWOLA-2D0*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2

      DRMON=-3D0/4D0*DHI*COSPHI**2*FAC2MON*((XMON(1)**2-XMON(2)**2)*  
     . SINTWOLA-2D0*XMON(1)*XMON(2)*COSTWOLA)/RMON**2

      DNSUN=3D0/2D0*DLI*SINPHI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)* 
     . SINTWOLA-2D0*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2

      DNMON=3D0/2D0*DLI*SINPHI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)* 
     . SINTWOLA-2D0*XMON(1)*XMON(2)*COSTWOLA)/RMON**2

      DESUN=-3D0/2D0*DLI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*  
     . COSTWOLA+2D0*XSUN(1)*XSUN(2)*SINTWOLA)/RSUN**2

      DEMON=-3D0/2D0*DLI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*  
     . COSTWOLA+2D0*XMON(1)*XMON(2)*SINTWOLA)/RMON**2

      DR=DRSUN+DRMON 
      DN=DNSUN+DNMON  
      DE=DESUN+DEMON 

      XCORSTA(1)=DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA  
      XCORSTA(2)=DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA  
      XCORSTA(3)=DR*SINPHI+DN*COSPHI  

      RETURN 

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END  
*
      SUBROUTINE ST1L1_iers2010 (XSTA,XSUN,XMON,FAC2SUN,FAC2MON
     .,XCORSTA)
*+
*  - - - - - - - - - - -
*   S T 1 L 1
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine gives the corrections induced by the latitude 
*  dependence given by L^1 in Mathews et al. 1991 (See References).
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 1
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
*     XSUN          d(3)   Geocentric position of the Sun (Note 2)
*     XMON          d(3)   Geocentric position of the Moon (Note 2)
*     FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)      
*     FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3) 
*
*  Returned:
*     XCORSTA       d(3)   Out of phase station corrections for
*                          semi-diurnal band
*
*  Notes:
*
*  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
*     expressed in meters. 
*  
*  2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
*     coordinates are expressed in meters.
*
*  3) The expressions are computed in the main program. TGP is the tide
*     generated potential.  The units are inverse meters. 
*
*  Test case:
*     given input: XSTA(1) = 4075578.385D0 meters
*                  XSTA(2) =  931852.890D0 meters
*                  XSTA(3) = 4801570.154D0 meters   
*                  XSUN(1) = 137859926952.015D0 meters
*                  XSUN(2) = 54228127881.4350D0 meters
*                  XSUN(3) = 23509422341.6960D0 meters
*                  XMON(1) = -179996231.920342D0 meters
*                  XMON(2) = -312468450.131567D0 meters
*                  XMON(3) = -169288918.592160D0 meters
*                  FAC2SUN =  0.163271964478954D0 1/meters     
*                  FAC2MON =  0.321989090026845D0 1/meters    
*                  
*     expected output:  XCORSTA(1) = 0.2367189532359759044D-03 meters
*                       XCORSTA(2) = 0.5181609907284959182D-03 meters
*                       XCORSTA(3) = -0.3014881422940427977D-03 meters
*
*  References:
*
*     Mathews, P. M., Buffett, B. A., Herring, T. A., Shapiro, I. I.,
*     1991b, Forced nutations of the Earth: Influence of inner core
*     Dynamics 2. Numerical results and comparisons, J. Geophys. Res.,
*     96, 8243-8257
*
*     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
*     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  1996 March    23 V. Dehant      Original code
*  2009 July     31 B.E. Stetzler  Initial standardization of code 
*  2009 July     31 B.E. Stetzler  Provided a test case and Mathews
*                                  reference
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION NORM8, RSTA, SINPHI, COSPHI, COSTWOLA, SINLA,
     .                 COSLA, RMON, RSUN, DRSUN, DRMON, DNSUN, DNMON,
     .                 DESUN, DEMON, DR, DN, DE, XSTA, XSUN, XMON,
     .                 XCORSTA, DHI, DLI, FAC2SUN, FAC2MON, SINTWOLA,
     .                 L1, L1D, L1SD

      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)  
      DATA L1D/0.0012D0/,L1SD/0.0024D0/

* Compute the normalized position vector of the IGS station.
      RSTA = NORM8(XSTA)
      SINPHI = XSTA(3)/RSTA  
      COSPHI = DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA  
      SINLA = XSTA(2)/COSPHI/RSTA  
      COSLA = XSTA(1)/COSPHI/RSTA  

* Compute the normalized position vector of the Moon.
      RMON = NORM8(XMON)

* Compute the normalized position vector of the Sun.
      RSUN = NORM8(XSUN)

* Compute the station corrections for the diurnal band.

      L1=L1D  
      DNSUN=-L1*SINPHI**2*FAC2SUN*XSUN(3)*(XSUN(1)*COSLA+XSUN(2)*SINLA)
     .            /RSUN**2
      DNMON=-L1*SINPHI**2*FAC2MON*XMON(3)*(XMON(1)*COSLA+XMON(2)*SINLA)
     .            /RMON**2
      DESUN=L1*SINPHI*(COSPHI**2-SINPHI**2)*FAC2SUN*XSUN(3)*
     . (XSUN(1)*SINLA-XSUN(2)*COSLA)/RSUN**2
      DEMON=L1*SINPHI*(COSPHI**2-SINPHI**2)*FAC2MON*XMON(3)*
     . (XMON(1)*SINLA-XMON(2)*COSLA)/RMON**2

      DE = 3D0*(DESUN+DEMON)  
      DN = 3D0*(DNSUN+DNMON)  

      XCORSTA(1) = -DE*SINLA-DN*SINPHI*COSLA  
      XCORSTA(2) = DE*COSLA-DN*SINPHI*SINLA  
      XCORSTA(3) = DN*COSPHI  
   
* Compute the station corrections for the semi-diurnal band.
  
      L1=L1SD  
      COSTWOLA=COSLA**2-SINLA**2  
      SINTWOLA=2.*COSLA*SINLA  

      DNSUN=-L1/2D0*SINPHI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     . COSTWOLA+2D0*XSUN(1)*XSUN(2)*SINTWOLA)/RSUN**2

      DNMON=-L1/2D0*SINPHI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     . COSTWOLA+2D0*XMON(1)*XMON(2)*SINTWOLA)/RMON**2

      DESUN=-L1/2D0*SINPHI**2*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     . SINTWOLA-2D0*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2

      DEMON=-L1/2D0*SINPHI**2*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     . SINTWOLA-2D0*XMON(1)*XMON(2)*COSTWOLA)/RMON**2

      DE = 3D0*(DESUN+DEMON)  
      DN = 3D0*(DNSUN+DNMON)  

      XCORSTA(1)=XCORSTA(1)-DE*SINLA-DN*SINPHI*COSLA  
      XCORSTA(2)=XCORSTA(2)+DE*COSLA-DN*SINPHI*SINLA  
      XCORSTA(3)=XCORSTA(3)+DN*COSPHI  

      RETURN

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END  
*
      SUBROUTINE STEP2DIU_iers2010 (XSTA,FHR,T,XCORSTA)  
*+
*  - - - - - - - - - - -
*   S T E P 2 D I U
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine gives the in-phase and out-of-phase corrections
*  induced by mantle anelasticity in the diurnal band. 
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 1
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
*     FHR           d      Fractional hours in the day (Note 2)
*     T             d      Centuries since J2000
*
*  Returned:
*     XCORSTA       d(3)   In phase and out of phase station corrections
*                          for diurnal band (Note 4)
*
*  Notes:
*
*  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
*     expressed in meters. 
*  
*  2) The fractional hours in the day is computed as the hour + minutes/60.0
*     + sec/3600.0.  The unit is expressed in Universal Time (UT).
*
*  4) All coordinates are expressed in meters.
*
*  Test case:
*     given input: XSTA(1) = 4075578.385D0 meters
*                  XSTA(2) =  931852.890D0 meters
*                  XSTA(3) = 4801570.154D0 meters 
*                  FHR     = 0.00D0 hours
*                  T       = 0.1059411362080767D0 Julian centuries
*                  
*     expected output:  XCORSTA(1) = 0.4193085327321284701D-02 meters
*                       XCORSTA(2) = 0.1456681241014607395D-02 meters
*                       XCORSTA(3) = 0.5123366597450316508D-02 meters
*
*  References:
*
*     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
*     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  1996 March    23 V. Dehant      Original code
*  2009 July     31 B.E. Stetzler  Initial standardization of code 
*  2009 August   06 B.E. Stetzler  Provided a test case
*  2009 August   06 B.E. Stetzler  Capitalized all variables for 
*                                  Fortran 77 compatibility
*  2010 October  20 B.E. Stetzler  Input T corrected to be number of
*                                  centuries since J2000
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I, J
      DOUBLE PRECISION XSTA(3), XCORSTA(3), DATDI(9,31), DEG2RAD, FHR,
     .                 T, S, TAU, PR, H, P, ZNS, PS, RSTA, SINPHI,
     .                 COSPHI, COSLA, SINLA, ZLA, THETAF, DR, DN, DE
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

      DATA ((DATDI(I,J),I=1,9),J=1,31)/  

     . -3D0, 0D0, 2D0, 0D0, 0D0,-0.01D0, 0D0, 0D0, 0D0,   
     . -3D0, 2D0, 0D0, 0D0, 0D0,-0.01D0, 0D0, 0D0, 0D0,   
     . -2D0, 0D0, 1D0,-1D0, 0D0,-0.02D0, 0D0, 0D0, 0D0,   
     . -2D0, 0D0, 1D0, 0D0, 0D0,-0.08D0, 0D0,-0.01D0, 0.01D0,
     . -2D0, 2D0,-1D0, 0D0, 0D0,-0.02D0, 0D0, 0D0, 0D0,
     . -1D0, 0D0, 0D0,-1D0, 0D0,-0.10D0, 0D0, 0D0, 0D0,
     . -1D0, 0D0, 0D0, 0D0, 0D0,-0.51D0, 0D0,-0.02D0, 0.03D0,
     . -1D0, 2D0, 0D0, 0D0, 0D0, 0.01D0, 0D0, 0D0, 0D0,
     .  0D0,-2D0, 1D0, 0D0, 0D0, 0.01D0, 0D0, 0D0, 0D0,
     .  0D0, 0D0,-1D0, 0D0, 0D0, 0.02D0, 0D0, 0D0, 0D0,
     .  0D0, 0D0, 1D0, 0D0, 0D0, 0.06D0, 0D0, 0D0, 0D0,
     .  0D0, 0D0, 1D0, 1D0, 0D0, 0.01D0, 0D0, 0D0, 0D0,
     .  0D0, 2D0,-1D0, 0D0, 0D0, 0.01D0, 0D0, 0D0, 0D0,
     .  1D0,-3D0, 0D0, 0D0, 1D0,-0.06D0, 0D0, 0D0, 0D0,
     .  1D0,-2D0, 0D0,-1D0, 0D0, 0.01D0, 0D0, 0D0, 0D0,
     .  1D0,-2D0, 0D0, 0D0, 0D0,-1.23D0,-0.07D0, 0.06D0, 0.01D0,
     .  1D0,-1D0, 0D0, 0D0,-1D0, 0.02D0, 0D0, 0D0, 0D0,
     .  1D0,-1D0, 0D0, 0D0, 1D0, 0.04D0, 0D0, 0D0, 0D0,
     .  1D0, 0D0, 0D0,-1D0, 0D0,-0.22D0, 0.01D0, 0.01D0, 0D0,
     .  1D0, 0D0, 0D0, 0D0, 0D0,12.00D0,-0.80D0,-0.67D0,-0.03D0,
     .  1D0, 0D0, 0D0, 1D0, 0D0, 1.73D0,-0.12D0,-0.10D0, 0D0,
     .  1D0, 0D0, 0D0, 2D0, 0D0,-0.04D0, 0D0, 0D0, 0D0, 
     .  1D0, 1D0, 0D0, 0D0,-1D0,-0.50D0,-0.01D0, 0.03D0, 0D0,
     .  1D0, 1D0, 0D0, 0D0, 1D0, 0.01D0, 0D0, 0D0, 0D0,
     .  0D0, 1D0, 0D0, 1D0,-1D0,-0.01D0, 0D0, 0D0, 0D0,
     .  1D0, 2D0,-2D0, 0D0, 0D0,-0.01D0, 0D0, 0D0, 0D0,
     .  1D0, 2D0, 0D0, 0D0, 0D0,-0.11D0, 0.01D0, 0.01D0, 0D0,
     .  2D0,-2D0, 1D0, 0D0, 0D0,-0.01D0, 0D0, 0D0, 0D0,
     .  2D0, 0D0,-1D0, 0D0, 0D0,-0.02D0, 0D0, 0D0, 0D0,
     .  3D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0,
     .  3D0, 0D0, 0D0, 1D0, 0D0, 0D0, 0D0, 0D0, 0D0/

      DEG2RAD = D2PI/360D0

*  Compute the phase angles in degrees.
      S = 218.31664563D0
     .  + (481267.88194D0
     .  + (-0.0014663889D0 
     .  + (0.00000185139D0)*T)*T)*T 

      TAU = FHR*15D0
     .    + 280.4606184D0
     .    + (36000.7700536D0
     .    + (0.00038793D0
     .    + (-0.0000000258D0)*T)*T)*T
     .    + (-S)

      PR = (1.396971278D0
     .   + (0.000308889D0
     .   + (0.000000021D0
     .   + (0.000000007D0)*T)*T)*T)*T 

      S = S + PR 

      H = 280.46645D0
     .  + (36000.7697489D0
     .  + (0.00030322222D0 
     .  + (0.000000020D0
     .  + (-0.00000000654D0)*T)*T)*T)*T  

      P = 83.35324312D0
     .  + (4069.01363525D0
     .  + (-0.01032172222D0
     .  + (-0.0000124991D0
     .  + (0.00000005263D0)*T)*T)*T)*T  

      ZNS = 234.95544499D0
     .    + (1934.13626197D0
     .    + (-0.00207561111D0
     .    + (-0.00000213944D0
     .    + (0.00000001650D0)*T)*T)*T)*T 

      PS = 282.93734098D0
     .   + (1.71945766667D0
     .   + (0.00045688889D0
     .   + (-0.00000001778D0
     .   + (-0.00000000334D0)*T)*T)*T)*T

* Reduce angles to between the range 0 and 360.
      S =  DMOD(S,360D0)
      TAU = DMOD(TAU,360D0)
      H =  DMOD(H,360D0)
      P =  DMOD(P,360D0)
      ZNS = DMOD(ZNS,360D0)
      PS = DMOD(PS,360D0)

      RSTA = DSQRT(XSTA(1)**2+XSTA(2)**2+XSTA(3)**2)  
      SINPHI = XSTA(3)/RSTA  
      COSPHI = DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA  

      COSLA = XSTA(1)/COSPHI/RSTA
      SINLA = XSTA(2)/COSPHI/RSTA
      ZLA = DATAN2(XSTA(2),XSTA(1))

      DO 99 I=1,3  
* Initialize.
      XCORSTA(I) = 0D0
99    CONTINUE
      DO 98 J=1,31
* Convert from degrees to radians.
      THETAF=(TAU+DATDI(1,J)*S+DATDI(2,J)*H+DATDI(3,J)*P+
     . DATDI(4,J)*ZNS+DATDI(5,J)*PS)*DEG2RAD

      DR=DATDI(6,J)*2D0*SINPHI*COSPHI*SIN(THETAF+ZLA)+
     . DATDI(7,J)*2D0*SINPHI*COSPHI*COS(THETAF+ZLA)

      DN=DATDI(8,J)*(COSPHI**2-SINPHI**2)*SIN(THETAF+ZLA)+
     . DATDI(9,J)*(COSPHI**2-SINPHI**2)*COS(THETAF+ZLA)
*      DE=DATDI(8,J)*SINPHI*COS(THETAF+ZLA)+
*     Modified 20 June 2007

      DE=DATDI(8,J)*SINPHI*COS(THETAF+ZLA)-
     . DATDI(9,J)*SINPHI*SIN(THETAF+ZLA)


      XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA  
     . -DN*SINPHI*COSLA  
      XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA  
     . -DN*SINPHI*SINLA  
      XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI  
98    CONTINUE   

      DO 97 I=1,3
      XCORSTA(I)=XCORSTA(I)/1000D0
97    CONTINUE  
      RETURN

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END  
*
      SUBROUTINE STEP2LON_iers2010 (XSTA,T,XCORSTA)  
*+
*  - - - - - - - - - - -
*   S T E P 2 L O N
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine gives the in-phase and out-of-phase corrections
*  induced by mantle anelasticity in the long period band. 
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 1
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
*     T             d      Centuries since J2000
*
*  Returned:
*     XCORSTA       d(3)   In phase and out of phase station corrections
*                          for diurnal band (Note 2)
*
*  Notes:
*
*  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
*     expressed in meters. 
*  
*  2) All coordinates are expressed in meters.
*
*  Test case:
*     given input: XSTA(1) = 4075578.385D0 meters
*                  XSTA(2) =  931852.890D0 meters
*                  XSTA(3) = 4801570.154D0 meters 
*                  T       = 0.1059411362080767D0 Julian centuries
*                  
*     expected output:  XCORSTA(1) = -0.9780962849562107762D-04 meters
*                       XCORSTA(2) = -0.2236349699932734273D-04 meters
*                       XCORSTA(3) =  0.3561945821351565926D-03 meters
*
*  References:
*
*     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
*     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  1996 March    23 V. Dehant      Original code
*  2009 August   07 B.E. Stetzler  Initial standardization of code
*                                  and found unnecessary variables tau
*                                  and fhr 
*  2009 August   07 B.E. Stetzler  Provided a test case
*  2009 August   07 B.E. Stetzler  Capitalized all variables for 
*                                  Fortran 77 compatibility
*  2010 October  20 B.E. Stetzler  Input T corrected to be number of 
*                                  centuries since J2000
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I, J
      DOUBLE PRECISION XSTA(3), XCORSTA(3), DATDI(9,5), DEG2RAD, 
     .                 T, S, PR, H, P, ZNS, PS, RSTA, SINPHI,
     .                 COSPHI, COSLA, SINLA, THETAF, DR, DN, DE,
     .                 DR_TOT, DN_TOT
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

      DATA ((DATDI(I,J),I=1,9),J=1,5)/  
     .   0, 0, 0, 1, 0,   0.47D0, 0.23D0, 0.16D0, 0.07D0,
     .   0, 2, 0, 0, 0,  -0.20D0,-0.12D0,-0.11D0,-0.05D0,
     .   1, 0,-1, 0, 0,  -0.11D0,-0.08D0,-0.09D0,-0.04D0,
     .   2, 0, 0, 0, 0,  -0.13D0,-0.11D0,-0.15D0,-0.07D0,
     .   2, 0, 0, 1, 0,  -0.05D0,-0.05D0,-0.06D0,-0.03D0/


      DEG2RAD = D2PI/360D0

*  Compute the phase angles in degrees.
      S = 218.31664563D0
     .  + (481267.88194D0
     .  + (-0.0014663889D0 
     .  + (0.00000185139D0)*T)*T)*T 

      PR = (1.396971278D0
     .   + (0.000308889D0
     .   + (0.000000021D0
     .   + (0.000000007D0)*T)*T)*T)*T 

      S = S + PR 

      H = 280.46645D0
     .  + (36000.7697489D0
     .  + (0.00030322222D0 
     .  + (0.000000020D0
     .  + (-0.00000000654D0)*T)*T)*T)*T  

      P = 83.35324312D0
     .  + (4069.01363525D0
     .  + (-0.01032172222D0
     .  + (-0.0000124991D0
     .  + (0.00000005263D0)*T)*T)*T)*T  

      ZNS = 234.95544499D0
     .    + (1934.13626197D0
     .    + (-0.00207561111D0
     .    + (-0.00000213944D0
     .    + (0.00000001650D0)*T)*T)*T)*T 

      PS = 282.93734098D0
     .   + (1.71945766667D0
     .   + (0.00045688889D0
     .   + (-0.00000001778D0
     .   + (-0.00000000334D0)*T)*T)*T)*T

      RSTA=DSQRT(XSTA(1)**2+XSTA(2)**2+XSTA(3)**2)  
      SINPHI=XSTA(3)/RSTA  
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA  

      COSLA=XSTA(1)/COSPHI/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA

* Reduce angles to between the range 0 and 360.
      S =  DMOD(S,360D0)
*      TAU = DMOD(TAU,360D0)
      H =  DMOD(H,360D0)
      P =  DMOD(P,360D0)
      ZNS = DMOD(ZNS,360D0)
      PS = DMOD(PS,360D0)

      DR_TOT = 0D0
      DN_TOT = 0D0

      DO 99 I=1,3  
      XCORSTA(I)=0D0
99    CONTINUE
      DO 98 J=1,5
      THETAF=(DATDI(1,J)*S+DATDI(2,J)*H+DATDI(3,J)*P+
     . DATDI(4,J)*ZNS+DATDI(5,J)*PS)*DEG2RAD

      DR=DATDI(6,J)*(3D0*SINPHI**2-1D0)/2D0*COS(THETAF)+
     .     DATDI(8,J)*(3D0*SINPHI**2-1D0)/2D0*SIN(THETAF)

      DN=DATDI(7,J)*(COSPHI*SINPHI*2D0)*COS(THETAF)+
     .     DATDI(9,J)*(COSPHI*SINPHI*2D0)*SIN(THETAF)

      DE = 0D0 
      DR_TOT = DR_TOT+DR
      DN_TOT = DN_TOT+DN

      XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA  
     . -DN*SINPHI*COSLA  
      XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA  
     . -DN*SINPHI*SINLA  
      XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI  
98    CONTINUE   

      DO 97 I=1,3
      XCORSTA(I)=XCORSTA(I)/1000D0
97    CONTINUE  
      RETURN

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
* 
      SUBROUTINE ZERO_VEC8_iers2010(V)
*+
*  - - - - - - - - - - -
*   Z E R O _ V E C 8
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine zeroes a vector.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Returned:
*     V            d(3)      vector V 
*
*  Called:
*     None
*
*  Test case: This is a support routine of the main program DEHANTTIDEINEL.F.
*     given input: V(1) = 2D0
*                  V(2) = 2D0
*                  V(3) = 3D0
*     
*     expected output: V = 0D0
*
*  References:
*
*     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
*     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 July 29 B.E.Stetzler Initial standardization of subroutine,
*                            provided a test case, and capitalized
*                            all variables for FORTRAN 77 compatibility  
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION V(3)

*  Get the zeros of the components of the vector.
      V(1)= 0D0
      V(2)= 0D0
      V(3)= 0D0

      RETURN

* Finished.
  
*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
