CTITLE ATM_AZ_PARTIAL
 
      subroutine atm_az_partial( type, azimuth, elev, in, part)
 

      implicit none 
 
*     Routine to compute the azimuthal dependent partials for delay
*     rate.   These partials are assumed to 1/(tan e sin e)
*     dependent with cos az and sin az defining the NS and EW
*     components of the assymetry.
*
*     NOTE: Since we only allow random walks for these parameters
*     we only need to do the 'offset' partials.  Therefore we can
*     generate the rate partials independently of the delay partials.
*     If a rate process we allowed, then the delay and rate partials
*     would need to be generated at the same time for the delay-
*     rate observable.
*
* MOD TAH 910918: Added 0.0032 to sin tan in DELAY partial to
*     stpo function going to inifinity.  (This value is mathed
*     for a tilted atmosphere down to 5 deg elevation).
*
*                                      09:48 PM SUN., 15 Feb., 1987
 
*   in      - Site number in baseline (i.e., 1 or 2)
 
      integer*4 in
 
*   azimuth(2,2)    - Azimith and rate (first index is site,
*           - second index is value and rate [rad, rad/sec])
*   cosa    - cos of aziumth
 
*   elev(2,2)   - Elevation angles and rate (see above)
 
*   part(2) - NS and EW partial derivatives, returned either
*           - as delay or rate depending on type.  Units are
*           - ps/ps for delay, and (fs/sec)/ps for rate.
*   sign    - Sign of partial (plus for in=2, minus for in=1)
*   sina    - sin of aziumuth
*   tanesine    - Value of tan e sin e
 
      real*4 azimuth(2,2), cosa, elev(2,2), part(2), sign, sina,
     .    tanesine
 
*   type    - type of data (either 'delay' or 'rate')
 
 
      character*(*) type
 
***** Do some start calculations
 
      if( in.eq.1 ) then
          sign = -1
      else
          sign = +1
      end if
 
      tanesine = tan(elev(in,1)) * sin(elev(in,1))
      cosa     = cos(azimuth(in,1))
      sina     = sin(azimuth(in,1))
 
***** Do NS and EW derivatives for delay or rate depending on 'type'
 
*                                     ! Delay partial
      if( type(1:2).eq.'de' ) then
 
          part(1) = sign*cosa / (tanesine+0.0032d0)
          part(2) = sign*sina / (tanesine+0.0032d0)
 
*                                     ! Rate partial
      else
 
          part(1) = sign*(-sina*azimuth(in,2)    /tanesine -
     .            elev(in,2)*cosa*sin( elev(in,1) )/tanesine**2*
     .                     (1.d0+1.d0/cos( elev(in,1) )**2) )*1000.d0
          part(2) = sign*( cosa*azimuth(in,2)/tanesine -
     .            elev(in,2)*sina*sin( elev(in,1) )/tanesine**2*
     .                     (1.d0+1.d0/cos( elev(in,1) )**2) )*1000.d0
      end if

****  Thats all
      return
      end
 
