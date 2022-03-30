Copyright (c) Massachusetts Institute of Technology and The University of
California at San Diego, 1994. All rights reserved.

      Subroutine rotsnp_sofa( idir,jd,t,tdtgpst,iut1pol
     .                 , frame,precmod,iut1,ipole,inut
     .                 , rot,rotdot,sidtm,xpole,ypole
     .                 , prec,rnut,prot,sidten,era,erarot)
C MOD TAH 200504: F90 module implemenation of the Deasai ans Sibois ocean tide HFEOP model.
      use hfeop_desai      ! module that contains desa model data and calcuation routines.
C MOD TAH 200509: F90 module with Gipson VLBI model with libration apriori.
      use hfeop_gipson     ! module constains Gipson VLBI model.
C MOD TAH 200504: Also updated logic so the PM_GRAVI and UTLIBR are treated together as
C     librarion terms.  All HF-EOP models for oceans should only be oceans (IERS2010 and
C     Desai and Sibois models).
c
C     Based on routin rotsnp by: Y. Bock and R. King 1986-1994; arguments changed RWK 170412 .
C     Modified by S. McClusky to use IAU SOFA library Earth Attitude routines .
C
C     Calculate the rotation matrix between earth-fixed and inertial
c
c     Input:
C       idir = 1   Earth-fixed to inertial
C             -1   Inertial to Earth-fixed
c       jd         Julian day (PEP JD, beginning at midnight, = MJD + 1)
c       t          Seconds-of-day  (GPST)
c       tdtgpst    TDT - GPST  (seconds)
c       iut1pol    Binary coded control for diurnal pole (bit 1) and UT1 (bit 2)
c       iut1, ipole, inut - unit numbers for tables 
c       frame,precmod - Inertial frame and precession model to be used in rotations

c      Output:
c        rot(3,3)     Rotation matix (rad)
c        rotdot(3,3)  Derivative of rotation matrix (rad/s)
c        sidtm        Sidereal time (GAST) (rad) (when EE based frame)
c        sidten       Sidereal time rotation matrix (rad)
c        ERA          Earth Rotation Angle (rad) (when CIO based frame)
c        ERArot       Earth Rotation Angle (rad) rotation matrix
c        xpole        X-pole position (rad)
c        ypole        Y-pole position (rad)
c        prec         Precession rotation matrix (rad)
c        rnut         Nutation rotation matrix (rad)
c        prot         Polar Motion rotation matrix (rad)

      implicit none
C
      character*5 frame, precmod

      integer*4 idir,jd,jds,jdutc,iut1pol,iut1,ipole,inut,i,j
      integer*4 iuttyp
      integer*4 len,rcpar
      
      real*8 t,tdtutc,tdtgpst,rot,rotdot,sidtm,xpole,ypole,xpdot,ypdot
     .     ,tutc,gpstutc,taiutc,taiut1,fract,lod,ut1dot,ut1utc

      dimension rot(3,3),rotdot(3,3)
     
      REAL*8 CRS2TRS(3,3),TRS2CRS(3,3)
      REAL*8 d_CRS2TRS(3,3),d_TRS2CRS(3,3)
      REAL*8 mjd_gpst,mjd_tt,mjd_utc,dx_eop,dy_eop
      REAL*8 rnut(3,3),prec(3,3),prot(3,3)
      REAL*8 sidten(3,3),era,erarot(3,3)
      REAL*8 rnut_t(3,3),prec_t(3,3),srot(3,3)
      REAL*8 sidten_t(3,3),erarot_t(3,3)
      real*8 pi,casr

* MOD TAH 200504: Added aaray for the hfeop_deasi model.
      real*8 eop(4,2)  ! Xp, Yp, UT1, LOD  and their rates. 
                       ! Units are (uas,uas,uts, uas/Day)
      real*8 mjday_TT, Delta_T ! TT time (MJD format) and time
                       ! difference between UTC and TT 
                       ! (Delta_TT = 68.931118 seconds, 2017Nov28)
           
      real*8 cor_x, cor_y, cor_lod,dx,dy,dUT1,dUTlib,dLODlib ! returns from IERS routine

      CHARACTER*10 iau_model
      character*80 p_name
      character*256 message
      
      LOGICAL kbit, debug/.false./

       pi = 4.D0*datan(1.D0)
       casr = pi/180.d0/3600.d0
c
c Get the program name calling rotsnp
      len = rcpar(0,p_name)
      
cd      write(*,'(a,1x,a5,1x,a5)')' ROTSNP frame precmod ',frame,precmod

c Compute jd and t in UTC (needed for interpolating x/y- pole and UT1-UTC values)
       jdutc = jd
       tutc = t
       gpstutc = taiutc(jdutc) - 19.d0
       call timinc(jdutc,tutc,-gpstutc)
       
c Catch possible leap second
       if( jdutc.ne.jd ) then
          gpstutc = taiutc(jdutc) - 19.d0
          jdutc = jd
          tutc = t
          call timinc(jdutc,tutc,-gpstutc)
       endif
       
       tdtutc = taiutc(jdutc) + 32.184d0
       mjd_utc = dfloat(jdutc-2400001) +  tutc/86400.d0 

       fract = tutc/86400.d0

c Get polar motion values from table      
       call polred(ipole,jdutc,fract,xpole,ypole,xpdot,ypdot)    

c Get ut1-utc values from table             
       call ut1red(iut1,jdutc,fract,taiut1,ut1dot,iuttyp) 
       if( taiut1.ge.0.d0 ) ut1utc=taiutc(jdutc)-taiut1
       
       if( kbit(iut1pol,4) ) then   ! IERS model
           call PMUT1_OCEANS (mjd_utc,dX, dY, dUT1,cor_lod)             
*          Solid Earth terms          
           call PM_GRAVI (mjd_utc,cor_x,cor_y)
           dX = (dX + cor_x)*1000  ! Convert to mas
           dY = (dY + cor_y)*1000  ! Convert to mas
           dUT1 = dUT1*1000        ! Convert to ms
c	   print*,'IERS dx,dy,dut1: ',dx,dy,dut1
       elseif( kbit(iut1pol,6) ) then   ! Gipson model 
*          Not clear if this model has solid-earth term included?
C          call gipson(fjdutc-0.5d0,dX, dY, dUT1)
C          dX = dX*1000  ! Convert to mas
C          dY = dY*1000  ! Convert to mas
C          dUT1 = dUT1*1000        ! Convert to ms
* MOD TAH 200509: Replaced with latest Gipson model from IERS HF EOP WG site
*          https://ivscc.gsfc.nasa.gov/hfeop_wg/ 
*          Code above rmjd = fjdutc - 2400001.d0 (MJD_UTC)
           Delta_T = tdtutc
           mjday_TT = mjd_utc + Delta_T/86400.d0 
*          Call routine to compute Ocean tide contribution.  Results
*          are retuned in uas and us.
           call calc_hf_gip_xyu(MJDAY_TT, Delta_T, eop)
*          Solid Earth terms: And apply to ocean contribution.         
           call PM_GRAVI (mjd_utc,cor_x,cor_y)
           dX = eop(1,1)/1000+cor_x*1000        ! Convert to mas
           dY = eop(2,1)/1000+cor_y*1000        ! Convert to mas
           dUT1 = eop(3,1)/1000      ! Convert to ms
       elseif( kbit(iut1pol,7) ) then   ! Desai and Sibois model 
*          Code above rmjd = fjdutc - 2400001.d0 (MJD_UTC)
*          Delta_T = tdtutc = 32.184d0 + taiutc(jd) (seconds)
           Delta_T = tdtutc
           mjday_TT = mjd_utc + Delta_T/86400.d0 
*          Call routine to compute Ocean tide contribution.  Results
*          are retuned in uas and us.
           call calc_hf_eop_xyu(MJDAY_TT, Delta_T, eop)
*          Solid Earth terms: And apply to ocean contribution.         
           call PM_GRAVI (mjd_utc,cor_x,cor_y)
           dX = eop(1,1)/1000+cor_x*1000        ! Convert to mas
           dY = eop(2,1)/1000+cor_y*1000        ! Convert to mas
           dUT1 = eop(3,1)/1000      ! Convert to ms          
       else 
           call ray(mjd_utc,dx,dy,dut1 )  ! Ray Model
c	   print*,'Ray dx,dy,dut1: ',dx,dy,dut1	   
       end if  
           
      if (kbit(iut1pol,2)) then
        ut1utc = ut1utc + (dUT1/1000.d0)
      endif            
      
      if( kbit(iut1pol,5) ) then 
         call UTLIBR(mjd_utc,dUTlib,dLODlib )
	 ut1utc = ut1utc + (dUTlib/1.d6)    
!	 print*,'UTLibr dut1: ',dutlib	   
      endif	
       
      if (kbit(iut1pol,1)) then
        xpole = xpole + (dx/1000.d0)
        ypole = ypole + (dy/1000.d0)
      endif
         
c IAU SOFA Earth Attitude software library http://www.iausofa.org
c Version tested - SOFA Library Issue 2018-01-30 
c iau_model: Precession-Nutation model by International Astronomical Union (IAU)
c IAU_model = 1976 refers to the IAU 1976/1980A model
c IAU_model = 2000a refers to the IAU 2000A CIO based model
c IAU_model = 2000e refers to the IAU 2000A Equinox based model
c IAU_model = 2006a refers to the IAU 2006/2000A CIO based model
c IAU_model = 2006ab refers to the IAU 2006 x,y series based model		   

      if ( precmod .eq. 'IAU0A' ) then
         iau_model = '2000e'
      else if ( precmod .eq. 'IAU0C' ) then
         iau_model = '2000a'
      else if ( precmod .eq. 'IAU06' ) then
         iau_model = '2006a'
      else if ( precmod .eq. 'IAU6A' ) then
         iau_model = '2006ab'
      else
         write(message,'(a,a5)') 
     .   'Error, unknown precession model',precmod
         call report_stat('FATAL',p_name,'orbits/rotsnp','',message,0)
      endif   
          
c mjd_tt: Modified Julian Day number of the epoch (in TT scale) 
c Converting from PEP JD to mjd
       mjd_tt   = dfloat(jd-2400001) +  t/86400.d0 + tdtgpst/86400.d0
       mjd_gpst = dfloat(jd-2400001) +  t/86400.d0 
       
c==================================================
cd SOFA Cookbook Test data
cd       mjd_tt   = 54195.500754444444444D0     
cd       xpole = 0.0349282D0
cd       ypole = 0.4833163D0 
cd       ut1utc = -0.072073685D0
c==================================================

c Set IERS nutation model corrections to 0. Maybe use at later stage?
      dx_eop = 0.d0
      dy_eop = 0.d0
c Convert pole position to radians
      xpole = xpole*casr
      ypole = ypole*casr
       
      call crs_trs(iau_model,mjd_tt,xpole,ypole,ut1utc,dx_eop,dy_eop
     .             , CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS
     .             , prec,rnut,prot,sidtm,sidten,era,erarot)
     
c Compute S-rotation matrix
      call matmpy(prot,sidten,srot,3,3,3)
c      call matmpy(prot,erarot,srot,3,3,3)
c      print*,'srot:   ',srot

c Move to GAMIT naming conventions with correct direction applied
      if(idir.eq.-1 ) then
         rot     =   CRS2TRS     
         rotdot  = d_CRS2TRS
        
      else      
         rot     =   TRS2CRS
         rotdot  = d_TRS2CRS
      end if

c Debug output
      if ( debug ) then
        print*,'ROTSNP_SOFA mjd_gpst,mjd_tt,mjd_utc,xpole,ypole,ut1utc:'
     .                     ,mjd_gpst,mjd_tt,mjd_utc,xpole,ypole,ut1utc
        print*,'ROTSNP_SOFA jd,t,tdtgpst,idir: ',jd,t,tdtgpst,idir
        if (idir .eq. 1 ) then
          write(*,100) 'ROTSNP_SOFA - TRS2CRS:   ',precmod,iau_model,
     .                 ((TRS2CRS(i,j),j=1,3),i=1,3) 
          write(*,100) 'ROTSNP_SOFA - d_TRS2CRS: ',precmod,iau_model,
     .               ((d_TRS2CRS(i,j),j=1,3),i=1,3) 
        else
          write(*,100) 'ROTSNP_SOFA - CRS2TRS:   ',precmod,iau_model,
     .               ((CRS2TRS(i,j),j=1,3),i=1,3) 
          write(*,100) 'ROTSNP_SOFA - d_CRS2TRS: ',precmod,iau_model,
     .             ((d_CRS2TRS(i,j),j=1,3),i=1,3)
        endif
        write(*,100) 'ROTSNP_SOFA - P_rot: ',precmod,iau_model,
     .             ((prec(i,j),j=1,3),i=1,3) 
        write(*,100) 'ROTSNP_SOFA - N_rot: ',precmod,iau_model,
     .             ((rnut(i,j),j=1,3),i=1,3) 
        write(*,100) 'ROTSNP_SOFA - PM_rot: ',precmod,iau_model,
     .             ((prot(i,j),j=1,3),i=1,3) 
        write(*,100) 'ROTSNP_SOFA - SID_rot: ',precmod,iau_model,
     .             ((sidten(i,j),j=1,3),i=1,3) 
        write(*,100) 'ROTSNP_SOFA - ERA_rot: ',precmod,iau_model,
     .             ((erarot(i,j),j=1,3),i=1,3) 
        write(*,100) 'ROTSNP_SOFA - S_rot: ',precmod,iau_model,
     .             ((srot(i,j),j=1,3),i=1,3) 
        print*,'ROTSNP_SOFA - sidtm: ',precmod,iau_model,sidtm
        print*,'ROTSNP_SOFA - era:   ',precmod,iau_model,era
100     format(a,1x,a,1x,a,/,3(1x,3D22.14,/))
      endif

      return
      end
