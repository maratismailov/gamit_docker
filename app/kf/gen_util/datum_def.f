CTITLE DATUM_DEF

      subroutine datum_def( in_datum )

      implicit none 

*     Routine to set the parameters for datum definitions.  Values
*     are made available through the utmut.h common

      include '../includes/utmut.h'

* PASSED VARIABLES
      character*(*) in_datum   ! Input datum

* LOCAL VARIABLES
      integer*4 i  ! Loop counter
      real*8 C2, C4, C6, C8, R, COSOR, OMO, SO

*
*     Set the parameters based on the datm choices
      do i = 1,3
         dXYZ_ut(i) = 0.d0
      end do 
      if ( in_datum(1:5).eq.'WGS84' ) then
         datum_ut = 'WGS84'
         ER_ut = 6378137.D0
	 RF_ut = 298.257222101D0
      elseif ( in_datum(1:5).eq.'WGSFL' ) then
         datum_ut = 'WGSFL'
         ER_ut = 6378137.D0
	 RF_ut = 298.257223563d0

      elseif ( in_datum(1:5).eq.'NAD27' ) then
         datum_ut = 'NAD27'
         ER_ut = 6378206.4D0  ! Clarke 1866 Ellipdsiod
	 RF_ut = 294.978698D0
*        Datum shift X : -8.00000  ! Values for Western US
*        shift Y : 159.00000
*        shift Z : 175.00000
C         dXYZ_ut(1) =   -8.0d0 
C         dXYZ_ut(2) = 159.0d0
C         dXYZ_ut(3) = 294.0d0
         dXYZ_ut(1) =  -12.01d0  ! These are the GAMIT values
         dXYZ_ut(2) = -162.97d0
         dXYZ_ut(3) = 189.740d0

      elseif ( in_datum(1:4).eq.'OMAN' ) then
         datum_ut = 'OMAN'    ! Clarke 1880 Ellipsoid
         ER_ut = 6378249.145d0
	 RF_ut = 293.465d0
*        Datum shift X : -346.00000
*        shift Y : -1.00000
*        shift Z : 224.00000
         dXYZ_ut(1) = -346.0d0
         dXYZ_ut(2) =   -1.0d0
         dXYZ_ut(3) = 224.0d0
      else    ! Default to WGS84
         datum_ut = 'WGS84'
         ER_ut = 6378137.D0
	 RF_ut = 298.257222101D0
      end if

*     Compute the derived properties
      F_ut = 1.D0/RF_ut
      ESQ_ut =( F_ut+F_ut-F_ut*F_ut)

      EPS_ut=ESQ_ut/(1.-ESQ_ut)
      PR_ut = (1.-F_ut)*ER_ut
      EN_ut = (ER_ut-PR_ut)/(ER_ut+PR_ut)
      A_ut  = -1.5D0*EN_ut + (9.d0/16.d0)*EN_ut**3
      B_ut  = 0.9375D0*EN_ut**2	- (15.d0/32.d0)*EN_ut**4
      C_ut  = -(35./48.)*EN_ut**3
      U_ut  = 1.5D0*EN_ut - (27.d0/32.d0)*EN_ut**3
      V_ut  = 1.3125D0*EN_ut**2 - (55.d0/32.d0)*EN_ut**4
      W_ut  = (151.d0/96.d0)*EN_ut**3
      R_ut  =  ER_ut*(1.-EN_ut)*(1.-EN_ut**2)*
     .         (1.+2.25D0*EN_ut**2+(225.d0/64.d0)*EN_ut**4)
      OMO_ut = OR_ut + A_ut*SIN(2*OR_ut) + 
     .          B_ut*SIN(4*OR_ut) + C_ut*SIN(6*OR_ut)
      SO_ut  = SF_ut*R_ut*OMO_ut

      FE_ut = 500000.0D0
      SF_ut = 0.9996D0
      OR_ut = 0.0D0

****  Factors for conversion from UTM to lat and long
      EN2_ut=EN_ut*EN_ut
      EN3_ut=EN_ut*EN_ut*EN_ut
      EN4_ut=EN2_ut*EN2_ut

      C2 = -3.D0*EN_ut/2.D0+9.D0*EN3_ut/16.D0
      C4 = 15.D0*EN2_ut/16.D0-15.D0*EN4_ut/32.D0
      C6 = -35.D0*EN3_ut/48.D0
      C8 = 315.D0*EN4_ut/512.D0
      U0_ut=2.D0*(C2-2.D0*C4+3.D0*C6-4.D0*C8)
      U2_ut=8.D0*(C4-4.D0*C6+10.D0*C8)
      U4_ut=32.D0*(C6-6.D0*C8)
      U6_ut=128.D0*C8

      C2=3.D0*EN_ut/2.D0-27.D0*EN3_ut/32.D0
      C4=21.D0*EN2_ut/16.D0-55.D0*EN4_ut/32.D0
      C6=151.D0*EN3_ut/96.D0
      C8=1097.D0*EN4_ut/512.D0
      V0_ut=2.D0*(C2-2.D0*C4+3.D0*C6-4.D0*C8)
      V2_ut=8.D0*(C4-4.D0*C6+10.D0*C8)
      V4_ut=32.D0*(C6-6.D0*C8)
      V6_ut=128.D0*C8

      R=ER_ut*(1.D0-EN_ut)*(1.D0-EN_ut*EN_ut)*(1.D0+2.25D0*EN_ut*EN_ut+
     .	   (225.D0/64.D0)*EN4_ut)
      COSOR=DCOS(OR_ut)
      OMO=OR_ut+DSIN(OR_ut)*COSOR*(U0_ut+U2_ut*COSOR*COSOR+
     .    U4_ut*COSOR**4+U6_ut*COSOR**6)
      SO=SF_ut*R_ut*OMO


****  That all
      return
      end

CTITLE CM_DEF

      subroutine cm_def( lng, zone, cm )

*     Routine to get the Central Meridian and Zone from the longitude
*     or given the zone, return the central meridian

      include '../includes/const_param.h'

      real*8 lng   ! East longtiude in radians
      real*8 cm    ! Central Meridian (radians)
      integer*4 zone   ! Zone number for UTM

      real*8 lng_mod   ! Longitude between 0 and 2pi and shifted so that 
                       ! -pi is 0 (start of UTM Zones)

*

****  OK: Convert the longitude to degrees and get the zone of the 
*     input zone is <= 0.
      if( zone.le.0 ) then
*        Make sure that longtide is between 0 and 2*pi radians
         lng_mod = lng
         if( lng.lt. 0 ) lng_mod = lng + 2*pi
         lng_mod = mod(lng_mod,2*pi)
         if( lng_mod.gt.pi ) lng_mod = lng_mod - 2*pi
*        
         lng_mod = lng_mod + pi
         zone = int(lng_mod*30/pi) + 1 ! 6-degree zones
      end if
      cm = (183 - zone*6)*pi/180.d0

****  Thats all
      return
      end
