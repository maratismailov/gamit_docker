SUBROUTINE crs_trs ( iau_model, mjd,xpole, ypole, ut1utc, dx_eop, dy_eop, &
                    CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS, &
                    RBP, RN, RPOM, GST, R_GST, ERA, R_ERA )

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine:  crs_trs.f90
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Purpose:
!  ICRF-ITRF transformation matrix (direct/inverse) for position/velocity vectors
!  Earth Orientation Matrix based on the EOP data corrected for short period tidal 
!  variations
! 
! Requires:
! IAU SOFA Earth Attitude software library http://www.iausofa.org
! Version tested - SOFA Library Issue 2018-01-30 
!
! Remark:
!  The matrices required for the ICRF-ITRF direct and inverse transformation
!  are computed for the position and velocity vectors individually.
!
! Input arguments:
! - mjd:                Modified Julian Day number of the epoch (in TT scale)
! - iau_model:		Precession-Nutation model by International Astronomical Union (IAU)
!			   IAU_model = 1976 refers to the IAU 1976/1980A model
!			   IAU_model = 2000a refers to the IAU 2000A CIO based model
!			   IAU_model = 2000e refers to the IAU 2000A Equinox based model
!			   IAU_model = 2006a refers to the IAU 2006/2000A CIO based model
!			   IAU_model = 2006ab refers to the IAU 2006 x,y series based model
!  
! - xpole, ypole:       Polar motion coordinates (rad) 
! - ut1utc:             Difference between UT1 and UTC (sec)
! - dx_eop, dy_eop:     EOP corrections to nutation model (arsec) (usually 0)
!   		
! Output arguments:
! - CRS2TRS:		GCRS to ITRS transformation matrix
! - TRS2CRS:		ITRS to GCRS transformation matrix 
! - d_CRS2TRS:		Derivative of GCRS to ITRS transformation matrix
! - d_TRS2CRS:		Derivative of ITRS to GCRS transformation matrix
! - RBP:                CRS to TRS Frame Bias + Precession rotation matrixc (P)
! - RN:                 CRS to TRS Nutation rotation matrix (N)
! - RPOM                CRS to TRS Polar Motion (W) rotation matrix
! - GST, ERA            CRS to TRS Greenwitch Aparant Siderial Time (GAST), CIO based Earth Rotation Angle (ERA)
! - R_GST, R_ERA        CRS to TRS GAST and ERA rotation (S) matrix
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! S Mcclusky             March 2019
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IMPLICIT NONE

! ----------------------------------------------------------------------
! Argument Declarations
! ----------------------------------------------------------------------
! IN
      REAL (KIND=8), INTENT(IN) :: mjd, xpole, ypole, dx_eop, dy_eop, ut1utc
      CHARACTER (LEN=10), INTENT(IN) :: iau_model
! OUT
      REAL (KIND=8), INTENT(OUT) :: CRS2TRS(3,3), TRS2CRS(3,3)
      REAL (KIND=8), INTENT(OUT) :: d_CRS2TRS(3,3), d_TRS2CRS(3,3)

! SOFA Called Functions (F77)
      DOUBLE PRECISION iau_S06, iau_ERA00, iau_ANP, iau_SP00, iau_OBL80
      DOUBLE PRECISION iau_EE00, iau_GMST00,iau_EQEQ94, iau_GMST82
      DOUBLE PRECISION iau_GST00A,iau_GST06A

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = 8) :: mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC, mjd_UT1
      DOUBLE PRECISION TT1, TT2, TT1_UT1, TT2_UT1, TAI1, TAI2
      INTEGER (KIND = 8) :: mjd_UTC_day
      INTEGER (KIND = 4) :: i,j

      DOUBLE PRECISION arcsec2rad, sp_iau00
      REAL (KIND = 8) :: X_iau00A,Y_iau00A,s_iau00A , X_iau06,Y_iau06,s_iau06
      REAL (KIND = 8) :: X_iau,Y_iau,s_iau , X_pn,Y_pn
      REAL (KIND = 8) :: pi
      DOUBLE PRECISION RPOM(3,3), RPOM_T(3,3)
      DOUBLE PRECISION RC2I(3,3), RC2I_T(3,3)
      DOUBLE PRECISION R_era(3,3), era
      DOUBLE PRECISION R_gst(3,3), gst
      DOUBLE PRECISION Qt(3,3), Qt_inv(3,3), Rt(3,3), Rt_inv(3,3), Wt(3,3), Wt_inv(3,3)
      DOUBLE PRECISION QR(3,3), Wi_Ri(3,3)
      DOUBLE PRECISION dtheta, P_dR (3,3), P_dR_T(3,3)
      DOUBLE PRECISION QP(3,3), QPR(3,3), QPRW(3,3), Ri_Qi(3,3), PT_Ri_Qi(3,3), Wi_PT_RiQi(3,3)
      DOUBLE PRECISION DX00, DY00, DX06, DY06
      
!     IAU1976/1980 local variables
      DOUBLE PRECISION RP(3,3),RB(3,3),RBP(3,3),RN(3,3),RBPN(3,3)
      DOUBLE PRECISION EE,GMST,DP80,DDP80,DE80,DDE80,DPSI,DEPS,EPSA
      DOUBLE PRECISION RC2TI(3,3),RC2IT(3,3),RC2IT_T(3,3)
      
      DOUBLE PRECISION V1(3),V2(3),X, Y, S,DDP00,DDE00,DE00,DP00
      
      DOUBLE PRECISION jd,tu2000,sidvel,sidtim

! ----------------------------------------------------------------------
      PARAMETER ( pi = 3.1415926535897932D0 )
! ----------------------------------------------------------------------
      arcsec2rad = pi / (3600.0D0 * 180.0D0)
! ----------------------------------------------------------------------
!      print*,'CRS_TRS: iau_model, mjd_tt, xpole, ypole, ut1utc: ',iau_model, mjd, xpole, ypole, ut1utc  
! ----------------------------------------------------------------------
! Time Systems transformation										 
! ----------------------------------------------------------------------
      CALL time_TT (mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC)
! ----------------------------------------------------------------------
! TT
      TT1 = 2400000.5D0
      TT2 = mjd_TT
! ----------------------------------------------------------------------
! TAI
      TAI1 = 2400000.5D0
      TAI2 = mjd_TAI
! ----------------------------------------------------------------------
! UTC
      mjd_UTC_day = INT (mjd_UTC)
! ----------------------------------------------------------------------
! UT1
      mjd_UT1 = mjd_UTC + ut1utc / 86400D0
      TT1_UT1 = 2400000.5D0
      TT2_UT1 = mjd_UT1  
! ----------------------------------------------------------------------
! Precession-Nutation model:  X, Y (in radians)
! ----------------------------------------------------------------------
! 1.  IAU 1976
      if ( index(iau_model,'1976') .ne. 0 ) then
      
         DDP80 = -55.0655D0 * arcsec2rad/1000.0D0
         DDE80 = -6.3580D0 * arcsec2rad/1000.0D0
! ----------------------------------------------------------------------	       
! 2.  IAU 2000A 
      else if ( index(iau_model,'2000') .ne. 0 ) then
      
!         DX00 =  0.1725D0 * arcsec2rad/1000D0
!         DY00 = -0.2650D0 * arcsec2rad/1000D0
          DX00 =  0.D0 * arcsec2rad/1000D0
          DY00 = -0.D0 * arcsec2rad/1000D0	 
! ----------------------------------------------------------------------
! 3.  IAU 2006/2000A     
      else if ( index(iau_model,'2006') .ne. 0 ) then  
           
!         DX06 =  0.1750D0 * arcsec2rad/1000D0
!         DY06 = -0.2259D0 * arcsec2rad/1000D0
          DX06 =  0.D0 * arcsec2rad/1000D0
          DY06 =  0.D0 * arcsec2rad/1000D0
! ----------------------------------------------------------------------
! 4. Unknown IAU model type	  
      else
  
         print*,'IAU Precession-Nutation model: ',iau_model,' Not enabled'
	 stop     
	 	 
      endif
      
! Initialise GST and ERA Rotation Matrices : Set them to identity I3x3
! R3(GST)
      r_gst = 0.0d0
      forall (i = 1:3) r_gst(i,i) = 1.0d0
      r_era(1:3,1:3) = 0.0d0
      forall (i = 1:3) r_era(i,i) = 1.0d0

! ----------------------------------------------------------------------
! Compute TRF <--> CRF transformation depending on IAU model eslected
! ----------------------------------------------------------------------

      if ( trim(iau_model) == '1976') then
            
! ===============================================
! ======= SOFA COOKBOOK - IAU 1976/1980 =========
! ===============================================
!
! IAU 1976 precession matrix, J2000.0 to date.
         CALL iau_PMAT76 ( TT1, TT2, RP )
! IAU 1980 nutation.
         CALL iau_NUT80 ( TT1, TT2, DP80, DE80 )
! Add adjustments: frame bias, precession-rates, geophysical.
         DPSI = DP80 + DDP80
         DEPS = DE80 + DDE80
! Mean obliquity.
         EPSA = iau_OBL80 ( TT1, TT2 )
! Build the rotation matrix.
         CALL iau_NUMAT ( EPSA, DPSI, DEPS, RN )
! Combine the matrices: PN = N x P.
         CALL iau_RXR ( RN, RP, RBPN )
! Equation of the equinoxes, including nutation correction.
         EE = iau_EQEQ94 ( TT1, TT2 ) + DDP80 * COS ( EPSA )
! Greenwich mean sidereal time (IAU 1982/1994).
         GMST = iau_GMST82(TT1_UT1,TT2_UT1) 
! Greenwich apparent sidereal time (IAU 1982/1994).
         GST = iau_ANP ( GMST + EE )
         CALL iau_RZ ( GST, R_GST )
! Form celestial-terrestrial matrix (no polar motion yet).
	 RC2I = RBPN
	 CALL iau_CR ( RBPN, RC2TI )
         CALL iau_RZ ( GST, RC2TI )
! Polar motion matrix (TIRS->ITRS, IERS 1996).
         CALL iau_IR ( RPOM )
         CALL iau_RX ( -ypole, RPOM )
         CALL iau_RY ( -xpole, RPOM )
! Form celestial-terrestrial matrix (including polar motion).
         CALL iau_RXR ( RPOM, RC2TI, RC2IT )      
! Form terrestrial-celestial matrix (including polar motion).
         CALL iau_TR ( RC2IT, RC2IT_T )
! Debug
!	 print*,'IAU2000A CIO CtoT: ',iau_model,RC2IT
!	 print*,'IAU2000A CIO TtoC: ',iau_model,RC2IT_T
	 	 	 
      else if ( trim(iau_model) == '2000a' ) then	 
	 
! ====================================================
! ======= SOFA COOKBOOK - IAU 2000A, CIO based =======
! ====================================================
! CIP and CIO, IAU 2000A.
         CALL iau_XYS00A ( TT1, TT2, X, Y, S )
! Add CIP corrections.
         X = X + DX00 + dX_eop * arcsec2rad
         Y = Y + DY00 + dY_eop * arcsec2rad
! Precession + Nutation 
         CALL iau_PN00A (TT1, TT2, DPSI, DEPS, EPSA, RB, RP, RBP, RN, RBPN)
! Greenwich Apparent Sidereal Time  (GST) 
         GST = iau_GST00A (TT1_UT1,TT2_UT1,TT1,TT2)
         CALL iau_RZ ( GST, R_GST )
! GCRS to CIRS matrix.
         CALL iau_C2IXYS ( X, Y, S, RC2I )
! Earth rotation angle.
         ERA = iau_ERA00 ( TT1_UT1,TT2_UT1 )
         CALL iau_RZ ( ERA, R_ERA )	 
! Form celestial-terrestrial matrix (no polar motion yet).
         CALL iau_CR ( RC2I, RC2TI )
         CALL iau_RZ ( ERA, RC2TI )
! Polar motion matrix (TIRS->ITRS, IERS 2003).
         CALL iau_POM00 ( xpole, ypole, iau_SP00(TT1,TT2), RPOM )
! Form celestial-terrestrial matrix (including polar motion).
         CALL iau_RXR ( RPOM, RC2TI, RC2IT )
         CALL iau_TR ( RC2IT, RC2IT_T )	 
! Debug
!	 write(*,*)   'CRS_TRS - CIO X,Y,S: ',iau_model,X,Y,S
!         write(*,100) 'CRS_TRS - RB:    ',iau_model,((RB(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RP:    ',iau_model,((RP(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RBP:   ',iau_model,((RBP(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RN:    ',iau_model,((RN(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RBPN:  ',iau_model,((RBPN(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RC2I:  ',iau_model,((RC2I(i,j),j=1,3),i=1,3)
!	 write(*,*)   'CRS_TRS - ERA/GST: ',ERA,GST
!         write(*,100) 'CRS_TRS - RC2TI: ',iau_model,((RC2TI(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RPOM:  ',iau_model,((RPOM(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - CtoT:  ',iau_model,((RC2IT(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - TtoC:  ',iau_model,((RC2IT_T(i,j),j=1,3),i=1,3)
100     format(a,1x,a,1x,/,3(1x,3D22.14,/))	 
	 	 
      else if ( trim(iau_model) == '2000e' ) then	 
      
! ========================================================
! ======= SOFA COOKBOOK - IAU 2000A, equinox based =======
! ========================================================
! Nutation, IAU 2000A.
         CALL iau_NUT00A ( TT1, TT2, DP00, DE00 )
! Precession-nutation quantities, IAU 2000.
         CALL iau_PN00 ( TT1, TT2, DP00, DE00, EPSA, RB, RP, RBP, RN, RBPN )
! Transform dX,dY corrections frqom GCRS to mean of date.
         V1(1) = DX00 + dX_eop * arcsec2rad
         V1(2) = DY00 + dY_eop * arcsec2rad
         V1(3) = 0D0
         CALL iau_RXP ( RBPN, V1, V2 )
         DDP00 = V2(1) / SIN ( EPSA )
         DDE00 = V2(2)
! Corrected nutation.
         DPSI = DP00 + DDP00
         DEPS = DE00 + DDE00
! Build the rotation matrix.
         CALL iau_NUMAT ( EPSA, DPSI, DEPS, RN )
! Combine the matrices: PN = N x P.
         CALL iau_RXR ( RN, RBP, RBPN )
! Greenwich apparent sidereal time (IAU 2000).
         GMST = iau_GMST00(TT1_UT1,TT2_UT1,TT1,TT2)
	 EE   = iau_EE00(TT1,TT2,EPSA,DPSI)
	 GST  = iau_ANP (GMST + EE)
         CALL iau_RZ ( GST, R_GST )
! Form celestial-terrestrial matrix (no polar motion yet).
 	 RC2I = RBPN
	 CALL iau_CR ( RBPN, RC2TI )
         CALL iau_RZ ( GST, RC2TI )
! Polar motion matrix (TIRS->ITRS, IERS 2003).
         CALL iau_POM00 ( xpole, ypole, iau_SP00(TT1,TT2), RPOM )
! Form celestial-terrestrial matrix (including polar motion).
         CALL iau_RXR ( RPOM, RC2TI, RC2IT )   
         CALL iau_TR ( RC2IT, RC2IT_T )	
! Debug
!	 write(*,*)   'CRS_TRS - EQE DPSI,DEPS: ',iau_model,DPSI,DEPS
!         write(*,100) 'CRS_TRS - RN:     ',iau_model,((RN(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RBP:    ',iau_model,((RBP(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RBPN:   ',iau_model,((RBPN(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RC2TI:  ',iau_model,((RC2TI(i,j),j=1,3),i=1,3)
!	 write(*,*)   'CRS_TRS - EQE GMST,EE,GST: ',iau_model,GMST,EE,GST
!         write(*,100) 'CRS_TRS - RPOM:   ',iau_model,((RPOM(i,j),j=1,3),i=1,3) 
!         write(*,100) 'CRS_TRS - CtoT:   ',iau_model,((RC2IT(i,j),j=1,3),i=1,3) 
!         write(*,100) 'CRS_TRS - TtoC:   ',iau_model,((RC2IT_t(i,j),j=1,3),i=1,3) 

      else if (trim(iau_model) == '2006a') then	 	

! =========================================================
! ======= SOFA COOKBOOK - IAU 2006/2000A, CIO based =======
! =========================================================
! CIP and CIO, IAU 2006/2000A.
         CALL iau_XYS06A ( TT1, TT2, X, Y, S )
! Add CIP corrections.
         X = X + DX06 + dX_eop * arcsec2rad
         Y = Y + DY06 + dY_eop * arcsec2rad
! Precession + Nutation 
         CALL iau_PN06A (TT1, TT2, DPSI, DEPS, EPSA, RB, RP, RBP, RN, RBPN)
! Greenwich Apparent Sidereal Time  (GST) 
         GST = iau_GST06A (TT1_UT1,TT2_UT1,TT1,TT2)
         CALL iau_RZ ( GST, R_GST )
! GCRS to CIRS matrix.
         CALL iau_C2IXYS ( X, Y, S, RC2I )
! Earth rotation angle.
         ERA = iau_ERA00 ( TT1_UT1,TT2_UT1 )
         CALL iau_RZ ( ERA, R_ERA )
! Form celestial-terrestrial matrix (no polar motion yet).
         CALL iau_CR ( RC2I, RC2TI )
         CALL iau_RZ ( ERA, RC2TI )
! Polar motion matrix (TIRS->ITRS, IERS 2003).
         CALL iau_POM00 ( xpole, ypole, iau_SP00(TT1,TT2), RPOM )
! Form celestial-terrestrial matrix (including polar motion).
         CALL iau_RXR ( RPOM, RC2TI, RC2IT )   
         CALL iau_TR ( RC2IT, RC2IT_T )	 	 
! Debug
!         write(*,*)   'CRS_TRS - CIO X,Y,S: ',iau_model,X,Y,S
!         write(*,100) 'CRS_TRS - RB:    ',iau_model,((RB(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RP:    ',iau_model,((RP(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RBP:   ',iau_model,((RBP(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RN:    ',iau_model,((RN(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RBPN:  ',iau_model,((RBPN(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RC2I:  ',iau_model,((RC2I(i,j),j=1,3),i=1,3)
!         write(*,*)   'CRS_TRS - ERA/GST: ',ERA,GST
!         write(*,100) 'CRS_TRS - RC2TI: ',iau_model,((RC2TI(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RPOM:  ',iau_model,((RPOM(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - CtoT:  ',iau_model,((RC2IT(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - TtoC:  ',iau_model,((RC2IT_T(i,j),j=1,3),i=1,3)

      else if (trim(iau_model) == '2006ab') then	 
      	      
! ==========================================================================
! ======= SOFA COOKBOOK -IAU 2006/2000A, CIO based, using X,Y series =======
! ==========================================================================
! CIP and CIO, IAU 2006/2000A.
         CALL iau_XY06 ( TT1, TT2, X, Y )
         S = iau_S06 ( TT1, TT2, X, Y )
! Add CIP corrections.
         X = X + DX06 + dX_eop * arcsec2rad
         Y = Y + DY06 + dY_eop * arcsec2rad
! Precession + Nutation 
         CALL iau_PN06A (TT1, TT2, DPSI, DEPS, EPSA, RB, RP, RBP, RN, RBPN)
! Greenwich Apparent Sidereal Time  (GST) 
         GST = iau_GST06A (TT1_UT1,TT2_UT1,TT1,TT2)
         CALL iau_RZ ( GST, R_GST )
! GCRS to CIRS matrix.
         CALL iau_C2IXYS ( X, Y, S, RC2I )
! Earth rotation angle.
         ERA = iau_ERA00 ( TT1_UT1,TT2_UT1 )
         CALL iau_RZ ( ERA, R_ERA )
! Form celestial-terrestrial matrix (no polar motion yet).
         CALL iau_CR ( RC2I, RC2TI )
         CALL iau_RZ ( ERA, RC2TI )
! Polar motion matrix (TIRS->ITRS, IERS 2003).
         CALL iau_POM00 ( xpole, ypole, iau_SP00(TT1,TT2), RPOM )
! Form celestial-terrestrial matrix (including polar motion).
         CALL iau_RXR ( RPOM, RC2TI, RC2IT )      
         CALL iau_TR ( RC2IT, RC2IT_T )	 	 
! Debug
!         write(*,*)   'CRS_TRS - CIO X,Y,S: ',iau_model,X,Y,S
!         write(*,100) 'CRS_TRS - RB:    ',iau_model,((RB(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RP:    ',iau_model,((RP(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RBP:   ',iau_model,((RBP(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RN:    ',iau_model,((RN(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RBPN:  ',iau_model,((RBPN(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RC2I:  ',iau_model,((RC2I(i,j),j=1,3),i=1,3)
!         write(*,*)   'CRS_TRS - ERA:   ',ERA
!         write(*,100) 'CRS_TRS - RC2TI: ',iau_model,((RC2TI(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - RPOM:  ',iau_model,((RPOM(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - CtoT:  ',iau_model,((RC2IT(i,j),j=1,3),i=1,3)
!         write(*,100) 'CRS_TRS - TtoC:  ',iau_model,((RC2IT_T(i,j),j=1,3),i=1,3)
 
      else
                 
         print*,'IAU Precession-Nutation model: ',iau_model,' Not enabled'
         stop
      	                      
      end if  
! ----------------------------------------------------------------------
! Save TRS - CRS and CRS - TRS matricies
! ----------------------------------------------------------------------

      TRS2CRS = RC2IT_T
      CRS2TRS = RC2IT 

! ----------------------------------------------------------------------
! Compute Derivatives of CRS-TRS transformation matrix (direct/inverse)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Earth angular velocity
! ----------------------------------------------------------------------
! IERS 2010;
!      dtheta = 0.7292115D0 * 1.0D-04
!      write(*,*)   'CRS_TRS - dtheta IERS2010: ',iau_model,dtheta
! ----------------------------------------------------------------------
! Dtheta is the derivative of Earth Rotation Angle (ERA), dERA/dt in rad/s
      dtheta = 2.0D0 * pi * 1.00273781191135448D0 * (1.0D0 / 86400D0)
!      write(*,*)   'CRS_TRS - dtheta dERA/dt: ',iau_model,dtheta
! ----------------------------------------------------------------------
! Old GAMIT sidmat/sidtim Definitions - for checking only
!      jd = mjd_TT+2400000.5d0       
!      tu2000 = (jd - 2451545.d0 )/36525.d0
!      sidvel = (1.002737909350795d0 +5.9006d-11*tu2000 -5.9d-15*tu2000**2)*2.0D0*pi/86400.d0
!      sidtim = (24110.54841d0+ 8640184.812866d0*tu2000 + 0.093104d0*tu2000**2 - 6.2d-6*tu2000**3)*2.0D0*pi/86400.d0
!      write(*,*)   'CRS_TRS - sidvel,sidtime: ',iau_model,sidvel,sidtim     
!      print*,'jd,mjd_tt,tu2000,sidtim,sidvel',jd,mjd_tt,tu2000,sidtim,sidvel     
! ----------------------------------------------------------------------
! d_CRS2TRS 
! ----------------------------------------------------------------------
! d_CRS2TRS = d(CRS2TRS)/dt = dtheta * W * P_dr_T * R * Q 
! ----------------------------------------------------------------------
!P = [ 0  -1   0
!      1   0   0
!      0   0   0 ] ;
! ----------------------------------------------------------------------
      P_dR  =  0.0D0
      P_dR (1,2) = -1.0D0
      P_dR (2,1) =  1.0D0

! ----------------------------------------------------------------------
! d_CRS2TRS = d(CRS2TRS)/dt = dtheta * W * P_dr_T * R * Q 
! ----------------------------------------------------------------------
      CALL iau_TR ( P_dR, P_dR_T )
      if ( trim(iau_model) .ne. '2000e' .and. trim(iau_model) .ne. '1976' ) then
         CALL iau_RXR ( R_ERA, RC2I, Ri_Qi )
      else
         CALL iau_RXR ( R_GST, RC2I, Ri_Qi )
      endif
      CALL iau_RXR ( P_dR_T, Ri_Qi, PT_Ri_Qi ) 
      CALL iau_RXR ( RPOM, PT_Ri_Qi, Wi_PT_RiQi )
      d_CRS2TRS = dtheta * Wi_PT_RiQi
!      print*,'d_CRS2TRS: ',iau_model,d_CRS2TRS
! d_TRS2CRS is transpose of d_CRS2TRS
      CALL iau_TR (d_CRS2TRS, d_TRS2CRS)
!      print*,'d_TRS2CRS: ',iau_model,d_CRS2TRS

RETURN
END
