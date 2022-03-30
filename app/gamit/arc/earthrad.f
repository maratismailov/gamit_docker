C Calculate forces from earth radiation using a box wing model and an earth radiation model.
Called by ertorb.f
C Calls satprop.f and eflux.f
      Subroutine earthrad(angb,angphi,twopi,sbmass,antbody,ltvel,
     .           aunit,rc2r,rsb,rvecxvec,idbprt,idebug)
C Inputs: 
C		B-angle (angb) between SVEC and RVEC (radians)
C       Sat-Earth-Sun angle (angphi) between ESVEC and RVEC (radians)
C		2pi (twopi)
C		Satellite mass (sbmass) kg
C		satellite body type (works for all GNSS)
C 		Speed of light (ltvel) km s-1
C		aunit (mean earth sun dist) km
C		rc2r (earth sun dist squared) km^2
C       rsb        dist earth - satellite  km
C       idbprt - debug output specification (char 14-16 of arcout filename)
C       idebug - debug file unit number
C
C Outputs:
C		rvecxvec(2):  accelerations in kms-2 
C                     (1) in direction of rvec vector, 
C                     (2) in direction of r perpendicular (xvec) vector

      implicit none

      Real*8 angb,angphi,twopi,sbmass,ltvel,aunit,rc2r,rsb
      integer*4  i,idbprt,idebug                            
      character*20 antbody

      logical kbit
      Real*8 theta
C     Coefficients for bus+masts and solar panels (front or back) for visible and IR wavelengths 
      Real*8 bmradvis(2),bmradir(2),panelFradvis(2),panelFradir(2),
     .       panelBradvis(2),panelBradir(2)
C     Calculated earth radiation at the satellite (visible,IR), Combined effect in 
C     rvec direction, and in perpendicular to rvec (xvec) direction.
      Real*8 fluxvis,fluxir,rvecvis,rvecir,rvecperpvis,rvecperpir
 
      Real*8 rvecxvec(2)

C     Initialise
            do i = 1,2
               bmradvis(i)=0.d0
               bmradir(i)=0.d0
               panelFradvis(i)=0.d0
               panelFradir(i)=0.d0
               panelBradvis(i)=0.d0
               panelBradir(i)=0.d0
            enddo 
C testing
C      twopi = 2.d0*3.14159265358979d0
C      angb=3.d0*(twopi/8.d0)
C      iblock = 3  ! Block IIA
C      iblock = 4  ! Block IIR-A
C      sbmass = 1100.d0
C      ltvel= 299792.458D0  ! km/sec

C	num 	(1 = bus +  masts shortwave)
C			(2 = bus + masts longwave)
C			(3 = solar panels front shortwave)
C			(4 = solar panels front longwave)
C			(5 = solar panels back shortwave)
C			(6 = solar panels back longwave)


C	Calculate theta (angle between the earth radiation vector and 
C   the normal to the solar panel surface)
       if (angb.le.twopi/4.d0) then
            theta=angb
       else 
            theta=(twopi/2.d0)-angb
       endif

cd      Print*, 'theta', theta
cd      Print*, 'angb', angb
cd      Print*, 'twopi', twopi
cd      print *,'EARTHRAD antbody ',antbody
C Get coefficient for rvec for bus / masts
C	Shortwave
      call satprop (theta,angb,twopi,antbody,1,bmradvis)
C       Longwave
      call satprop (theta,angb,twopi,antbody,2,bmradir)
C Get coefficient for rvec for solar panels
      if (angb.lt.twopi/4.d0) then
C       Get coefficient for rvec for front of solar panels
C	Shortwave
        call satprop (theta,angb,twopi,antbody,3,panelFradvis)
C       Longwave
        call satprop (theta,angb,twopi,antbody,4,panelFradir)
      else
C       Get coefficient for rvec for back of solar panels
C	Shortwave
        call satprop (theta,angb,twopi,antbody,5,panelBradvis)
C       Longwave
        call satprop (theta,angb,twopi,antbody,6,panelBradir)
      endif
C     Print debug (all epochs/eclipses only) if want earth radiation debug
      if(((kbit(idbprt,1)).or.(kbit(idbprt,2))).and.kbit(idbprt,9)) then
cd         idebug = 6 
         write(idebug,*) 'theta', theta
         write(idebug,*) 'bmradvis', bmradvis
         write(idebug,*) 'bmradir', bmradir
         write(idebug,*) 'panelFradvis',panelFradvis
         write(idebug,*) 'panelFradir',panelFradir
         write(idebug,*) 'panelBradvis',panelBradvis
         write(idebug,*) 'panelBradir',panelBradir
      endif 

C     Earth radiation :
C     Visible, IR, W/m2
C     Initialise as zero
      fluxvis=0.d0
      fluxir=0.d0
C     calculate earth radiation according to simple analytical approximations from Rodrigues 2009
      call eflux(angphi,twopi,aunit,rc2r,rsb,fluxvis,fluxir)
cd       Print*,'fluxvis',fluxvis
cd       Print*,'fluxir',fluxir
cd       print *
cd     .   ,'bmradvis panelfradvis panelfradir panelbradvis panelbradir'
cd     .   ,bmradvis,panelfradvis,panelfradir,panelbradvis,panelbradir
C     Calculate total effects of Earth radiation
C     If earth radiation will hit the front of the solar panels
      if (angb.lt.twopi/4.d0) then
C       // rvec
        rvecvis=fluxvis*(bmradvis(1)+panelFradvis(1))
        rvecir=fluxir*(bmradir(1)+panelFradir(1))
C       // xvec
        rvecperpvis=fluxvis*panelFradvis(2)
        rvecperpir=fluxir*panelFradir(2)
      else
C       If earth radiation will hit the back of the solar panels
C       // rvec
        rvecvis=fluxvis*(bmradvis(1)+panelBradvis(1))
        rvecir=fluxir*(bmradir(1)+panelBradir(1))
C       // xvec
        rvecperpvis=fluxvis*panelBradvis(2)
        rvecperpir=fluxir*panelBradir(2)
      endif
C      complete calculation of accelerations and convert to km to be consistent with arc
       rvecxvec(1)=(rvecvis+rvecir)*1.d-6/(sbmass*ltvel)
       rvecxvec(2)=(rvecperpvis+rvecperpir)*1.d-6/(sbmass*ltvel)
C     output debug data
      if(((kbit(idbprt,1)).or.(kbit(idbprt,2))).and.kbit(idbprt,9)) then
cd          idebug = 6 
         write(idebug,*) 'fluxvis W/m^2',fluxvis
         write(idebug,*) 'fluxir W/m^2',fluxir
         write(idebug,*) 'rvecvis',rvecvis
         write(idebug,*) 'rvecir',rvecir
         write(idebug,*) 'rvecperpvis',rvecperpvis
         write(idebug,*) 'rvecperpir',rvecperpir
         write(idebug,*) 'rvecxvec ms-2',rvecxvec(1)*1000.d0,
     .             rvecxvec(2)*1000.d0
      endif

      return
      end
