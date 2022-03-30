CTITLE GEOD_to_UTM
      subroutine GEOD_to_UTM(geod, utm, zone, hemi, datum)

      implicit none 

*     Routine to convert Geodetic Lat, lomg and height to UTM coordinates
*     Based on the NGS tmgrid.f routine. 

      include '../includes/utmut.h' 
      include '../includes/const_param.h'

* PASSED VARIABLES
      real*8 geod(3)   ! Geodetic co-latitude, longitude (+ve East) 
                       ! and height (deg, deg, m)
      real*8 utm(3)    ! UTM Northing, Easting and Ellipoidal height

      integer*4 zone   ! UTM Zone number (derived from longitude)

      character*(*) datum  ! Name of Datum to be used (default is WGS84)
                       ! (See Datum_def.f for current systems allowed)
      character*(*) hemi   ! Hemisphere for UTM Values


* LOCAL VARIABLES
      real*8 omega, S, sinlat, coslat, TN, TS, ETS, L, LS, RN,
     .       A1, A2, A3, A4, A5, A6, A7    ! Constants used in calcutions
      real*8 north, east  ! Northing and easting
      real*8 cm
      real*8 latr, lngw   ! Latiude in radians and west longitude in rads
      real*8 dlng  ! Longitude difference central meridian and point
      integer*4 ncyc ! Number of 2pi cycles in dlng
 


****  Get the datum parameters for the input GEOD coordinates
      call datum_def( datum)
*     Get the longitude of the central meridian and the zone number
      zone = 0
      call cm_def( geod(2), zone, cm)

***** Now run through the calculations in tmgrid.f
      latr = (pi/2-geod(1))
      if( latr.gt.0 ) then
          hemi = 'N'
      else
          hemi = 'S'
      endif

      lngw = -geod(2)

      omega = latr + A_ut*sin(2*latr) + 
     .               B_ut*sin(4*latr) + C_ut*sin(6*latr)
      S = R_ut*omega*SF_ut
      sinlat = sin(latr)
      coslat = cos(latr)
      TN = sinlat/coslat
      TS = TN**2
      ETS = EPS_ut*coslat**2
      dlng = lngw-cm
      ncyc = nint(dlng/(2*pi))
      dlng = dlng - ncyc*2*pi
      L = dlng*coslat

      LS = L**2
      RN = SF_ut*ER_ut/sqrt(1-esq_ut*sinlat**2)
*
      A2 = RN*TN/2.d0
      A4 = (5-TS+ETS*(9.+4.*ETS))/12.
      A6 = (61+TS*(TS-58)+ETS*(270-330*TS))/360
      A1 = -RN
      A3 = (1-TS+ETS)/6
      A5 = (5+TS*(TS-18)+ETS*(14-58*TS))/120
      A7 = (61-479*TS+179*TS**2-TS**3)/5040
      if( latr.lt.0 ) then
         FN_ut = 10000000.D0
      else
         FN_ut = 0.0d0
      end if

      NORTH = S - SO_ut + A2*LS*(1+LS*(A4+A6*LS)) + FN_ut
      EAST = FE_ut +      A1*L*(1+ LS*(A3+LS*(A5+A7*LS)))

      IF(NORTH.LT.0.0) THEN
        NORTH = -NORTH
      ENDIF

****  OK Save the values
      utm(1) = north
      utm(2) = east
      utm(3) = geod(3)

***** Thats all 
      return 
      end

CTITLE UTM_to_GEOD

      subroutine UTM_to_GEOD(utm, geod, zone, hemi, datum)

*     Routine to convert utm coordinates to geodetic latitude, longitude
*     and height.  Based on NGS routine TMGEOD.

      include '../includes/utmut.h' 
      include '../includes/const_param.h'

* PASSED VARIABLES
      real*8 geod(3)   ! Geodetic co-latitude, longitude (+ve East) 
                       ! and height (deg, deg, m)
      real*8 utm(3)    ! UTM Northing, Easting and Ellipoidal height

      integer*4 zone   ! UTM Zone number (derived from longitude)

      character*(*) datum  ! Name of Datum to be used (default is WGS84)
                       ! (See Datum_def.f for current systems allowed)
      character*(*) hemi   ! Either N or S for Northern Southern hemispahere

* LOCAL VARIABLES
      real*8 N, E   ! Northing and Easting
      real*8 cm     ! Central Meridian longitude (rads
      real*8 OM
      real*8 COSOM
      real*8 FOOT
      real*8 SINF
      real*8 COSF
      real*8 TN
      real*8 TS
      real*8 ETS
      real*8 RN
      real*8 Q
      real*8 QS
      real*8 B2
      real*8 B4
      real*8 B6
      real*8 B1
      real*8 B3
      real*8 B5
      real*8 B7
      real*8 LAT
      real*8 L
      real*8 LON

****  Get the datum parameters for the input GEOD coordinates
      call datum_def( datum)
*     Get the longitude of the central meridian from the zone number
      call cm_def( geod(2), zone, cm)
*     See which hemisphere
      if( hemi(1:1).eq.'S' .or. hemi(1:1).eq.'s' ) then
         FN_ut = 10000000.D0
      else
         FN_ut = 0.d0
      end if
*
      N = utm(1)
      E = utm(2)

      OM=(N-FN_ut+SO_ut)/(R_ut*SF_ut)
      COSOM=DCOS(OM)
      FOOT=OM+DSIN(OM)*COSOM*(V0_ut+V2_ut*COSOM*COSOM+V4_ut*COSOM**4+
     &	   V6_ut*COSOM**6)
      SINF=DSIN(FOOT)
      COSF=DCOS(FOOT)
      TN=SINF/COSF
      TS=TN*TN
      ETS=EPS_ut*COSF*COSF
      RN=ER_ut*SF_ut/SQRT(1.D0-ESQ_ut*SINF*SINF)
      Q=(E-FE_ut)/RN
      QS=Q*Q
      B2=-TN*(1.D0+ETS)/2.D0
      B4=-(5.D0+3.D0*TS+ETS*(1.D0-9.D0*TS)-4.D0*ETS*ETS)/12.D0
      B6=(61.D0+45.D0*TS*(2.D0+TS)+ETS*(46.D0-252.D0*TS-
     &	  60.D0*TS*TS))/360.D0
      B1=1.D0
      B3=-(1.D0+TS+TS+ETS)/6.D0
      B5=(5.D0+TS*(28.D0+24.D0*TS)+ETS*(6.D0+8.D0*TS))/120.D0
      B7=-(61.D0+662.D0*TS+1320.D0*TS*TS+720.D0*TS**3)/5040.D0
      LAT=FOOT+B2*QS*(1.D0+QS*(B4+B6*QS))
      L=B1*Q*(1.D0+QS*(B3+QS*(B5+B7*QS)))
      LON=-L/COSF+CM

****  Now convert values to output array
      geod(1) = (pi/2-LAT)
      geod(2) = -LON
      if( geod(2).lt.0 ) then 
         geod(2) = geod(2) + 2*pi
      endif
      geod(3) = utm(3)

***** Thats all
      return
      end






      

