      Subroutine IMFW2( a, b, c, elev, wmf )

c     Received from Arthur Niell 12 January 2004.
c     Modified to accept the mapping function coefficients as arguments
c     and also to pass the height in.
c     -- PT 040113                  
c     RWK 060722: Don't pass the height in, but correct the 'a' coefficient before calling



*   a,b,c       - the a,b,and c coefficients in the continued fraction
*                 form of Marini
*   beta        - intermediate term in calculation
*   gamma       - intermediate term in calculation
*   sine        - Sine of elevation angle
*   cose        - Cos of elevation angle
*   wmf(1)      - wet delay mapping function
*   wmf(2)      - d_wet_mapping_function/d_elevation
*   topcon      - Constant of top of mapping fuinction to ensure
*                 that value is 1.0000 at zenith

      real*8 a,b,c, beta, cose, wmf(2), gamma, sine, topcon

*   elev       - elevation (degrees)
c   height     - height of the station in metres (ellipsoidal I think!)
      real*8 elev, deg2rad
      parameter (deg2rad = 3.14159265358979/180.d0)

c debug
c      print*,'input arguments:'
c      print*,'mapcof',mapcof
c      print*,'elev ',elev

c  the IMF mapping function has a station height dependency
c  on the A coefficient.
c ***  call changed to pass a,b,c with 'a' corrected in ATMDEL   rwk 060722
c      a=mapcof(4) + height * 1.6580e-7 
c      b=mapcof(5)
c      c=mapcof(6)
c      print*,'a,b and c are',a,b,c

*   Now the coefficients exist; calculate the mapping function, wmf(1),
*       and the change of mapping function with elevation,
*       dwmf/d_el =wmf(2).
*   To calculate the delay-rate correction, d_tau/dt:
*       d_tau/dt = d_tau_zen/dt * wmf(1) + tau_zen * dwmf/d_el * d_el/dt
c      print*,'elev and deg2rad',elev,deg2rad,elev*deg2rad
      sine  = sin( elev * deg2rad)
      cose  = cos( elev * deg2rad)
      beta  = b/( sine + c )
      gamma = a/( sine + beta)
      topcon = (1.d0 + a/(1.d0 + b/(1.d0 + c)))

c      write(*,'("sine, cose, beta, gamma, topcon = ", 5f10.5)')
c     .           sine, cose, beta, gamma, topcon

      wmf(1) = topcon / ( sine + gamma )

      wmf(2) = -topcon / ( sine + gamma )**2 *
     .         ( cose - a/( sine + beta)**2 * cose *
     .         ( 1.d0 - b/( sine + c )**2 ) )

c      write(*,'("wmf(1), wmf(2) = ", 2f10.4)') wmf(1), wmf(2)
c      write(*,'("wmf(1), wmf(2) = ", 2f10.4)') wmf

caen   write out diagnostic info.
c     write(*,'("  elev, latitude, wmf2.0, dwmf/del = ",4f15.6)')
c    .             elev, latitude, wmf(1), wmf(2)

      return
      end

