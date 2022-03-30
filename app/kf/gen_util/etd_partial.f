CTITLE ETD_PARTIAL
 
      subroutine etd_partial
 

      implicit none 
 
c     routine to compute the extended earth tide partials
c     K1 frequency used for the diurnal tides and M2 is used
c     for the semi-diurnal tides
c MOD TAH 870616:  Converted units of partials to ps/m for
c     delays and (fs/s)/m for rates.
c     The ordering of the partials is
c     1 -- diurnal Radial in phase
c     2 -- diurnal Radial out of phase
c     3 -- diurnal Horizonal East in phase
c     4 -- diurnal Horizonal East out of phase
c     5 -- diurnal Horizonal South in phase
c     6 -- diurnal Horizonal South out of phase
c     7 -- semi-diurnal Radial in phase
c     8 -- semi-diurnal Radial out of phase
c     9 -- semi-diurnal Horizontal East in phase
c     10-- semi-diurnal Horizontal east out of phase
c     11-- semi-diurnal Horizontal South in phase
c     12-- semi-diurnal Horizontal South out of phase
c
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
      include '../includes/const_param.h'
 
 
 
*          i          - local site number
*          idex       - index for storing tidal partials
*          itide      - diurnal (1) and semi-diurnal (2) tides
 
      integer*4 i, idex, itide
 
*       du(6)         - partials of local disp. wrt coefficients
 
 
      real*4 du(6)
 
*       colat         - co-latitude of site
*       cos_colat     - cosine of the co-latitude
*       phase         - phase of the normal mode
*       plm           - Legendre function (function)
*       plmdz         - derivative of Legendre function (function)
*       sgn           - signum function
*       tid_freq(2)   - tidal frequencies (K1, M2 in cycles/solar day) (not used)
*       xplm          - Legendre function (variable)
*       xplmdz        - derivative of Legendre function (variable)
*       arg           - GST of observation (rad)
 
 
      real*8 colat, cos_colat, phase, plm, plmdz,
     .    xplm, xplmdz, arg  
c     real*8 tid_freq(2)   -not used
 
c      save tid_freq   - not used
 
c.... Set tidal frequencies
c      data tid_freq / 1.0027379093, 1.9322736136 /  -- not used
c
c.... Loop over two sites
      do i=1,2
 
c....   Loop over diurnal and semi diurnal tides
        do itide=1,2
 
 
c         calculate phase of mode (reference epoch JD 2451545.0)
          call gst_jd(epoch, arg)
          phase = (arg+longitudes(site(i)))*itide

*         Adjust phase of diurnals by -pi/2 to keep in-phase ie
*         arguments are for:
*         semidirnal    cos(arg)  -i sin(arg)
*         diurnal       sin(arg)  -i cos(arg)
          if( itide.eq.1 ) then
              phase = phase - pi/2.d0
          end if
 
          phase = mod(phase,2.d0*pi)
 
c         calculate Legendre functions
          colat =   pi/2.d0 - latitudes(site(i))
          cos_colat = cos(colat)
          xplm = plm(2,itide,cos_colat,0)
          xplmdz = plmdz(2,itide,cos_colat,0)
 
c....     Partials of local displacements wrt tidal coefficients
c         radial displacements
          du(1) =  xplm*cos(phase)
          du(2) = -xplm*sin(phase)
 
c         eastward displacements.  NOTE the phasing of these terms
*         is such that the tidal potential driving the horizontal
*         displacement is cos(phase) and sin(phase)

          du(3) = -itide*xplm*sin(phase)/sin(colat)
          du(4) = -itide*xplm*cos(phase)/sin(colat)
 
c         southward displacements
          du(5) = -xplmdz*sin(colat)*cos(phase)
          du(6) =  xplmdz*sin(colat)*sin(phase)
 
 
c....     Partials of group delay  wrt displacement coefficients
          idex = (itide - 1) * 6
          etd_ext_part(1+idex,i,1) = -(2*i-3) * sin(elev(i,1))
     .                        * du(1) / vel_light * 1.d12
 
          etd_ext_part(2+idex,i,1) = -(2*i-3) * sin(elev(i,1))
     .                        * du(2) / vel_light * 1.d12
 
          etd_ext_part(3+idex,i,1) = -(2*i-3) * ( cos(elev(i,1))
     .            * sin(azimuth(i,1)) * du(3) ) / vel_light*1.d12
          etd_ext_part(4+idex,i,1) = -(2*i-3) * ( cos(elev(i,1))
     .            * sin(azimuth(i,1)) * du(4) ) / vel_light*1.d12

          etd_ext_part(5+idex,i,1) =  (2*i-3) * ( cos(elev(i,1))
     .            * cos(azimuth(i,1)) * du(5) ) / vel_light *1.d12
          etd_ext_part(6+idex,i,1) =  (2*i-3) * ( cos(elev(i,1))
     .            * cos(azimuth(i,1)) * du(6) ) / vel_light *1.d12

c....     Partials of delay rate wrt displacement coefficients
 
          etd_ext_part(1+idex,i,2) = -(2*i-3) * ( elev(i,2) *
     .              cos(elev(i,1)) * du(1) ) / vel_light *1.d15
 
          etd_ext_part(2+idex,i,2) = -(2*i-3) * ( elev(i,2) *
     .              cos(elev(i,1)) * du(2) ) / vel_light *1.d15
 
          etd_ext_part(3+idex,i,2) = -(2*i-3) * ( (-elev(i,2) *
     .              sin(elev(i,1)) * sin(azimuth(i,1)) +
     .              azimuth(i,2) * cos(elev(i,1)) * cos(azimuth(i,1)))
     .            * du(3) ) / vel_light * 1.d15 

          etd_ext_part(4+idex,i,2) = -(2*i-3) * ( (-elev(i,2) *
     .              sin(elev(i,1)) * sin(azimuth(i,1)) +
     .              azimuth(i,2) * cos(elev(i,1)) * cos(azimuth(i,1)))
     .            * du(4) ) / vel_light * 1.d15 

          etd_ext_part(5+idex,i,2) = -(2*i-3) * ( (elev(i,2) *
     .              sin(elev(i,1)) * cos(azimuth(i,1)) +
     .              azimuth(i,2) * cos(elev(i,1)) * sin(azimuth(i,1)))
     .            * du(5) ) / vel_light * 1.d15 

          etd_ext_part(6+idex,i,2) = -(2*i-3) * ( (elev(i,2) *
     .              sin(elev(i,1)) * cos(azimuth(i,1)) +
     .              azimuth(i,2) * cos(elev(i,1)) * sin(azimuth(i,1)))
     .            * du(6) ) / vel_light * 1.d15 
 
         end do
 
      end do
 
 
      return
      end
 
 
 
