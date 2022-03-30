Copyright (c) Massachusetts Institute of Technology and University of
California at San Diego, 1994. All rights reserved.

      subroutine srotat( jd,t,tdtutc,eqe,iut1pol,srot,sdrot,sidtm
     .                 , xpole,ypole,sidten,prot,iut1,ipole,precmod )

C MOD TAH 200504: F90 module implemenation of the Deasai ans Sibois ocean tide HFEOP model.
      use hfeop_desai      ! module that contains desa model data and calcuation routines.
C MOD TAH 200509: F90 module with Gipson VLBI model with libration apriori.
      use hfeop_gipson     ! module constains Gipson VLBI model.
C MOD TAH 200504: Also updated logic so the PM_GRAVI and UTLIBR are treated together as
C     librarion terms.  All HF-EOP models for oceans should only be oceans (IERS2010 and
C     Desai and Sibois models).
C
C      Compute S rotation matrix from CRS true-of-date to CTRS

c     Note that this routine (uniquely in GAMIT) uses UTC not GPST as
c     a calling argument since the former is more natural for the Earth's
c     rotation (rwk 940714)

c     Input:
c       jd         Julian day (PEP JD, beginning at midnight, = MJD + 1)
c       t          Seconds-of-day
c       tdtutc     TDT - UTC  (seconds)
c       eqe        Equation of equinoxes
c       iut1pol    Binary coded control for diurnal pole (bit 1) and UT1 (bit 2)
c                   bit 3 (8) is Ray model; bit 4 (16) is libration terms for UT1
c                  (iut1pol is isptide in 'model' and 'solve' 
c       iut1, ipole, inut - unit numbers for tables
c       precmod - precession model to use for rotations

c      Output:
c        srot(3,3)    Rotation matix - wobble and sideral (rad)
c        sdrot(3,3)   Derivative of rotation matrix (rad/s)
c        sidtm        Sidereal time (rad)
c        xpole, ypole Pole position (rad)
c        sidten(3,3)  Sideral time matrix |  These needed in MODEL for
c        prot(3,3)    Wobble matrix       |  EOP partials

      implicit none

      logical kbit

      character*5 precmod
      character*80 prog_name

      integer*4 jd,iuttyp,iut1,ipole,iut1pol,len,rcpar
c      integer*4 i,j

      real*8 t,taiut1,taiutc,ut1utc,ut1dot,fract,tjd,xl,f,d,ascm
     .     , xrot,yrot,prot,srot,sdrot,sidten,sdtmat
     .     , xpole,xpdot,ypole,ypdot,eqe,pi,casr,sidtm
     .     , tdtutc,fjdutc,rmjd,fjdtdt,dx,dy,dUT1,dUTlib,dLODlib

c MOD TAH 110207: Implemented IERS model
      real*8 cor_x, cor_y, cor_lod  ! returns from IERS routine
              ! pmut1_oceans and pm_gravi.  Units arc-sec and secs.
c MOD TAH 200504: Added aaray for the hfeop_deasi model.
      real*8 eop(4,2)  ! Xp, Yp, UT1, LOD  and their rates. 
                       ! Units are (uas,uas,uts, uas/Day)
      real*8 mjday_TT, Delta_T ! TT time (MJD format) and time
                       ! difference between UTC and TT 
                       ! (Delta_TT = 68.931118 seconds, 2017Nov28)

      dimension xrot(3,3),yrot(3,3),prot(3,3),srot(3,3),sdrot(3,3)
      dimension sidten(3,3),sdtmat(3,3)
                        
c      get calling program name and m-file name for report_stat
       len = rcpar(0,prog_name)

      fract = t/86400.d0

c compute diurnal and seimdiurnal contributions to the pole position
c and ut1 using VLBI derived values for tidal variations. dx, dy, are
c output in milliarc seconds, and dut1 in millitime seconds.

c **  this recalculation should be unncessary since tdtutc is passed in
c **  add the following as a trap until we've verified this in all modules:
c **  (tdt-utc was 52.184 in 1981 and has increased to 61.184 by 1994)
      if( dabs(tdtutc-56.d0).ge.20.d0 )
     .  call report_stat('WARNING',prog_name,'lib/srotat',' '
     .                  , 'Input tdtutc undefined ',0)
      tdtutc = 32.184d0 + taiutc(jd)
c     JD here is PEP JD (= true JD + 0.5  and MJD -1 + 2400000 )
      fjdtdt= dble(jd) + fract + tdtutc/86400.d0
      fjdutc= dble(jd) + fract 
      rmjd = fjdutc - 2400001.d0 
cd     print *, 'jd fract fjdutc rmjd ',jd,fract,fjdutc,rmjd 
c      if ((kbit(iut1pol,1)).or.(kbit(iut1pol,2))) then
c       call sd_comp(fjdutc-0.5d0,dx,dy,dut1) 
c       write(*,'(a,f12.3,3f7.3)') 'SD_COMP: ',fjdutc,dx,dy,dut1
c **    force use of new short-period ut1/pole terms -- tah/rwk 990320
c         call ray( fjdutc-0.5d0,dx,dy,dut1 ) 
c       write(*,'(a,f12.3,3f7.3)') 'RAY   P: ',fjdutc,dx,dy,dut1
c      endif

*      See if Ray or IERS model (bit 4 sets IERS)
c      The Gipson model is implementbed (bit 6 sets Gipson model) -- lei 151007
       if( kbit(iut1pol,4) ) then   ! IERS model
           call PMUT1_OCEANS (fjdutc-0.5d0,dX, dY, dUT1,cor_lod)             
*          Solid Earth terms          
           call PM_GRAVI (fjdutc-0.5d0,cor_x,cor_y)
           dX = (dX + cor_x)*1000  ! Convert to mas
           dY = (dY + cor_y)*1000  ! Convert to mas
           dUT1 = dUT1*1000        ! Convert to ms
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
           mjday_TT = rmjd + Delta_T/86400.d0 
*          Call routine to compute Ocean tide contribution.  Results
*          are retuned in uas and us.
           call calc_hf_gip_xyu(MJDAY_TT, Delta_T, eop)
*          Solid Earth terms: And apply to ocean contribution.         
           call PM_GRAVI (fjdutc-0.5d0,cor_x,cor_y)
           dX = eop(1,1)/1000+cor_x*1000        ! Convert to mas
           dY = eop(2,1)/1000+cor_y*1000        ! Convert to mas
           dUT1 = eop(3,1)/1000      ! Convert to ms
       elseif( kbit(iut1pol,7) ) then   ! Desai and Sibois model 
*          Code above rmjd = fjdutc - 2400001.d0 (MJD_UTC)
*          Delta_T = tdtutc = 32.184d0 + taiutc(jd) (seconds)
           Delta_T = tdtutc
           mjday_TT = rmjd + Delta_T/86400.d0 
*          Call routine to compute Ocean tide contribution.  Results
*          are retuned in uas and us.
           call calc_hf_eop_xyu(MJDAY_TT, Delta_T, eop)
*          Solid Earth terms: And apply to ocean contribution.         
           call PM_GRAVI (fjdutc-0.5d0,cor_x,cor_y)
           dX = eop(1,1)/1000+cor_x*1000        ! Convert to mas
           dY = eop(2,1)/1000+cor_y*1000        ! Convert to mas
           dUT1 = eop(3,1)/1000      ! Convert to ms
       else 
           call ray( fjdutc-0.5d0,dx,dy,dut1 )  ! Ray Model
       end if
       

c Read UT1 from the input file
      call ut1red( iut1,jd,fract,taiut1,ut1dot,iuttyp )
cd    print *,'jd fract taiut1 ',jd,fract,taiut1
      if( taiut1.ge.0.d0 ) ut1utc=taiutc(jd)-taiut1

C Compute fundamental arguments for tidal effects on UT1
      tjd= jd + fract - 0.5D0
      call funarg( tjd,xl,f,d,ascm )
C     Add tidal correction if UT1-UTC is regularized (UT1R)
      if( iuttyp.eq.2 ) then
        call ut1tid( ut1utc,xl,f,d,ascm )
      else if (iuttyp.eq.0 ) then
        call report_stat('FATAL',prog_name,'lib/srotat',' '
     .         ,'UT1 type = 0, set = 2 or 4 in UT1. table',0)
      endif

c Add in diurnal and semidiurnal tide contribution to UT1 if requested.
      if (kbit(iut1pol,2)) then
         ut1utc = ut1utc + (dUT1/1000.d0)
      endif

c Add in semidiurnal libration contribution to UT1 if requested
      if( kbit(iut1pol,5) ) then 
         call UTLIBR( rmjd,dUTlib,dLODlib )
         ut1utc = ut1utc + (dUTlib/1.d6) 
      endif
 
c Read pole position from table
      call polred( ipole,jd,fract,xpole,ypole,xpdot,ypdot )
cd    print *,'jd fract xpole ypole ',jd,fract,xpole,ypole
      

c Add in diurnal and semi diurnal contributions to pole position if requested.
      if (kbit(iut1pol,1)) then
        xpole = xpole + (dx/1000.d0)
        ypole = ypole + (dy/1000.d0)
      endif

C     write(*,800) fjdutc-0.5d0, xpole, ypole, ut1utc, 
C    .             dx, dy, dut1, iut1pol
C800  format('EOP ',F14.6,1x,2F11.6,1x,F11.8,1x,2F8.3,1x,F8.5,1x,o4)

c Convert pole position to radians
      pi= 4.D0*datan(1.D0)
      casr= pi/180.d0/3600.d0
      xpole = xpole*casr
      ypole = ypole*casr

c Compute GAST and sidereal rotation matrix
      call sidmat(jd,fract,ut1utc,eqe,sidten,sdtmat,sidtm,precmod)
c     print *,'for SIDMAT iut1pol ut1utc ',iut1pol,ut1utc
c
c      print*,'in SROTAT jd,t,tdtutc,eqe,iut1pol,srot,sdrot,sidtm '
c     .                ,jd,t,tdtutc,eqe,iut1pol,srot,sdrot,sidtm
c
c      print*,'          xpole,ypole,sidten,prot,iut1,ipole'
c     .                ,xpole,ypole,sidten,prot,iut1,ipole

c ******** debug *************
c      print*,' in srotat sidtm = ',sidtm*86400.d0/(2*pi)
c      stop
c ****************************

c      write(*,700) ((SIDTEN(i,j),j=1,3),i=1,3)
c 700  Format(' GAST rotation matrix :',/,3(1x,3D22.14,/))
c      write(*,701) ((SDTMAT(i,j),j=1,3),i=1,3)
c 701  format(' GAST dot rotation matrix :',/,3(1x,3D22.14,/))

C Compute polar motion xp rotation
      Call rotmat(-xpole,2,xrot)
C Compute polar motion yp rotation
      Call rotmat(-ypole,1,yrot)
C Compute total polar motion rotation
      Call matmpy(xrot,yrot,prot,3,3,3)

c      write(*,705) ((XROT(i,j),j=1,3),i=1,3)
c 705  format(' X-Pole rotation matrix :',/,3(1x,3D22.14,/))
c      write(*,706) ((YROT(i,j),j=1,3),i=1,3)
c 706  format(' Y-Pole rotation matrix :',/,3(1x,3D22.14,/))
c      write(*,702) ((PROT(i,j),j=1,3),i=1,3)
c 702  format(' Pole rotation matrix :',/,3(1x,3D22.14,/))

c Compute S-rotation matrix
      Call matmpy(prot,sidten,srot,3,3,3)
c Compute SDOT-rotation matrix
      Call matmpy(prot,sdtmat,sdrot,3,3,3)
c      write(*,703) ((srot(i,j),j=1,3),i=1,3)
c 703  format(' S-rotation matrix :',/,3(1x,3D22.14,/))
c      write(*,704) ((sdrot(i,j),j=1,3),i=1,3)
c 704  format(' SDOT-rotation matrix :',/,3(1x,3D22.14,/))

      return
      end

