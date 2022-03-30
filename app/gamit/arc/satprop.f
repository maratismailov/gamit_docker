c Satellite properties - for calculating earth radiation effects.
      Subroutine satprop (theta,angb,twopi,antbody,num,out)

      implicit none
C Inputs 
C		theta (smallest angle between rvec and angle to solar panels)
C		angb (angle between svec and rvec)
C		antbody - SV body type (works for all GNSS)
C		2pi (twopi)
C		num 	(1 = bus +  masts shortwave)
C			(2 = bus + masts longwave)
C			(3 = solar panels front shortwave)
C			(4 = solar panels front longwave)
C			(5 = solar panels back shortwave)
C			(6 = solar panels back longwave)
C Output
C		out ( parallel to rvec, perpendicular to rvec (//to xvec)
C       
      Real*8 theta,angb,twopi,out(2)
      integer*4  num
      character*20 antbody 
      character*256 message

C Area of earth facing side of satelllite bus, masts(*2)  and solar panels (*2?) in m^2
      Real*8 azbus, a2mast, a2panel
c Specularity(1) and reflectivity(2) of bus, masts, solar panel fronts, solar panel rears to short wave light
      Real*8 srzbus(2), sr2mast(2), srpanelF(2), srpanelB(2)
C Estimate of long wave properties from Rodrigues for all surfaces (see p46)
      Real*8 srlong(2)
      Real*8 area,theta2,coeff,coeff1,coeff2,coeffrhat,coeffrperp

      logical BIIIA_norep / .true. /  ! Warning Block IIIA not reported yet.  Set false
                                      ! once reported, 

      logical firstcall/.true./

C Initialise output
        out(1)=0.d0
        out(2)=0.d0

CValues from Rodrigues-Solano (2009) Masters thesis, mainly from Fliegel et al. (1992) and Fliegel and Gallini (1996)
        srlong(1)=0.5d0
        srlong(2)=0.2d0
C Area of solar panels (m^2)
C       Block I
cc        if( iblock.eq.1 ) then  
        if( antbody(1:9).eq.'BLOCK I  ' ) then
C          Satellite bus
           azbus=1.510d0
           srzbus(1)=0.75d0
           srzbus(2)=0.86d0
C 	   Solar panel masts
           a2mast=0.470d0
           sr2mast(1)=0.85d0
           sr2mast(2)=0.85d0
C          Solar panels
           a2panel=5.583d0
           srpanelF(1)=0.85d0
           srpanelF(2)=0.23d0
           srpanelB(1)=0.5d0
           srpanelB(2)=0.11d0
C       Block II, IIA
cc        elseif ( iblock.eq.2.or.iblock.eq.3 ) then
        elseif ( antbody(1:9).eq.'BLOCK II ' .or.
     .           antbody(1:9).EQ.'BLOCK IIA' ) then
C          Satellite bus
           azbus=2.881d0
           srzbus(1)=0.20d0
           srzbus(2)=0.56d0
C 	   Solar panel masts
           a2mast=0.985d0
           sr2mast(1)=0.41d0
           sr2mast(2)=0.52d0
C          Solar panels
           a2panel=10.866d0
           srpanelF(1)=0.85d0
           srpanelF(2)=0.23d0
           srpanelB(1)=0.5d0
           srpanelB(2)=0.11d0
C       Block IIR_A, IIR_B, IIR-M
cc        elseif ( iblock.eq.4.or.iblock.eq.5.or.iblock.eq.6 ) then
        elseif ( antbody(1:11).eq.'BLOCK IIR-A' .or.
     .           antbody(1:11).eq.'BLOCK IIR-B' .or.
     .           antbody(1:11).eq.'BLOCK IIR-M' ) then 
C          Satellite bus
           azbus=3.750d0
           srzbus(1)=0.0d0
           srzbus(2)=0.06d0
C 	   Solar panel masts
           a2mast=0.320d0
           sr2mast(1)=0.85d0
           sr2mast(2)=0.85d0
C          Solar panels
           a2panel=13.6d0
           srpanelF(1)=0.85d0
           srpanelF(2)=0.28d0
           srpanelB(1)=0.5d0
           srpanelB(2)=0.11d0
c      For IIF, no solid information as yet (Aug 2012) 
C      putting the values from Rodrigues-Solano IGS recommended code, see 'PROPBOXW.f'
C      N.B that has 'SEE IIF PRESENTATION' and 'OPTICAL PROPERTIES NOT KNOWN, SAME ASSUMPTIONS
c      AS ZIEBART(2001)' Think this is Marek Ziebart's PhD thesis, but not sure which 
C      assumptions are referred to.
cc        elseif( iblock.eq.7 ) then 
cc      elseif( antbody(1:9).eq.'BLOCK IIF' ) then
        elseif( antbody(1:9).eq.'BLOCK IIF'.or.
     .          antbody(1:10).eq.'BLOCK IIIA' )then
* MOD TAH 200504: Added test to print message just once.
          if(antbody(1:10).eq.'BLOCK IIIA' .and. BIIIA_norep) then
             call report_stat('WARNING'  ,'ARC','satprop', ' ', 
     .      'Block IIIA not coded for Earth-radiation, use IIF values',
     .       0)
             BIIIA_norep = .false.
           endif

CC           call report_stat('FATAL','ARC','satprop',' ',
CC     .     'Satellite properties for IIF earth radiation unknown',0)
C          Satellite bus
           azbus=5.40d0
           srzbus(1)=0.d0
           srzbus(2)=0.d0
C 	   Solar panel masts
           a2mast=0.d0
           sr2mast(1)=0.d0
           sr2mast(2)=0.d0
C          Solar panels
           a2panel=22.25d0
           srpanelF(1)=0.85d0
           srpanelF(2)=0.23d0
           srpanelB(1)=0.5d0
           srpanelB(2)=0.11d0
        else
          if(firstcall) then   
            write(message,'(2a)') 'Earth radiation not coded for '
     .                        ,antbody
            call report_stat('FATAL','ARC','satprop',' ',message,0)
            firstcall =.false.
          endif
        endif

C Calculate the variable part of the radiation (Area and angle of panels to Earth radiation)

C       Calculate for bus and masts
      if ((num .eq. 1).or.(num.eq.2)) then 
	    area = azbus + a2mast
	    if (num.eq.1) then
	      coeff1=(1.d0+srzbus(1)*srzbus(2)+(2.d0/3.d0)*
     .         srzbus(2)*(1.d0-srzbus(1)))
	      coeff2=(1.d0+sr2mast(1)*sr2mast(2)+(2.d0/3.d0)*
     .         sr2mast(2)*(1.d0-sr2mast(1)))
          out(1) = azbus*coeff1+a2mast*coeff2
        elseif (num.eq.2) then
          coeff=(1.d0+srlong(1)*srlong(2)+(2.d0/3.d0)*
     .         srlong(2)*(1.d0-srlong(1)))  
          out(1) = area*coeff
        else
CC	           call report_stat('FATAL','ARC','satprop',' ',
CC     .     'Num error 1',0)
        endif 

C       Calculate for solar panels 
      elseif ((num.eq.3).or.(num.eq.4).or.(num.eq.5).or.(num.eq.6)) then
        area = a2panel
        theta2=theta*2.d0
C          Front of panels, shortwave
	if (num.eq.3) then
         coeffrhat=dcos(theta)*(1+srpanelF(1)*srpanelF(2)*dcos(theta2)
     .          +(2.d0/3.d0)*srpanelF(2)*dcos(theta)*(1.d0-srpanelF(1)))

         coeffrperp=dcos(theta)*(srpanelF(1)*srpanelF(2)*dsin(theta2)
     .           +(2.d0/3.d0)*srpanelF(2)*(1-srpanelF(1))*dsin(theta))
C          Back of panels, shortwave
        elseif (num.eq.5) then
         coeffrhat=dcos(theta)*(1+srpanelB(1)*srpanelB(2)*dcos(theta2)
     .         +(2.d0/3.d0)*srpanelB(2)*dcos(theta)*(1.d0-srpanelB(1)))

         coeffrperp=dcos(theta)*(srpanelB(1)*srpanelB(2)*dsin(theta2)
     .           +(2.d0/3.d0)*srpanelB(2)*(1-srpanelB(1))*dsin(theta))
C	   Either side of panels, longwave (specularity, reflectivity estimates are the same)
        elseif ((num.eq.4).or.(num.eq.6)) then

         coeffrhat=dcos(theta)*(1+srlong(1)*srlong(2)*dcos(theta2)
     .           +(2.d0/3.d0)*srlong(2)*dcos(theta)*(1.d0-srlong(1)))

         coeffrperp=dcos(theta)*(srlong(1)*srlong(2)*dsin(theta2)
     .           +(2.d0/3.d0)*srlong(2)*(1-srlong(1))*dsin(theta))
        else
CC	           call report_stat('FATAL','ARC','satprop',' ',
CC     .     'Num error 2',0)
        endif 
C          Make sure correct sign for perpendicular
        if (angb.lt.(twopi/4.d0)) then
		coeffrperp=-coeffrperp
        endif
C       set output
        out(1) = area * coeffrhat
        out(2) = area * coeffrperp
      else
CC	           call report_stat('FATAL','ARC','satprop',' ',
CC     .     'Num value not recognised',0) 
      endif
      return
      end
