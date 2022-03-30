C Calculate the earth radiation fluxes according to the analytical models of Rodrigues Solano (2009)
C with the addition of an adjustment for how far the earth is from the sun.
      Subroutine eflux(angphi,twopi,aunit,rc2,rsb,fluxvis,fluxir)

C Input :   	angphi: satellite - earth -sun angle, radians
C		twopi: two pi
C		aunit:  Astronomical unit in km (mean distance from Sun to Earth)
C		rc2:	magnitude squared sun -> earth   (km^2)
C		rsb:	magnitude earth -> earth satellite (km)
C
C
C Output:       fluxvis (visible/shortwave reflected radiation at the satellite W/m^2	
C	        fluxir (infrared/longwave emitted radiation at the satellite W/m^2 )
C               both in the direction of rvec.

      implicit none

      Real*8 angphi,twopi,aunit,rc2,rsb,fluxvis,fluxir

C      Mean earth radius(km),albedo,solar flux at 1 AU (W/m^2)
      Real*8 meanertrad,albedo,solcon
C      Adjusted incident solar radiation
       Real*8 solflux

C      Declare needed constants
      meanertrad=6371.d0
      albedo = 0.3d0
      solcon = 1367.d0

C   Adjusted incident radiation to Earth
      solflux=solcon*aunit*aunit/rc2

C   Reflected radiation
      fluxvis=(2.d0*albedo*meanertrad*meanertrad*solflux)*
     .          ( ((twopi/2.d0)-angphi)*dcos(angphi)+dsin(angphi) ) /
     .          (3.d0*(twopi/2.d0)*rsb*rsb)


C   Emitted radiation
      fluxir=(1-albedo)*meanertrad*meanertrad*solflux /
     .         (4.d0*rsb*rsb)

cd      Print*,'meanertrad',meanertrad,'albedo',albedo,'aunit',aunit
cd      Print*,'rsb',rsb,'rc2',rc2, 'twopi',twopi


      Return
      end
