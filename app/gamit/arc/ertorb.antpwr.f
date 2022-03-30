Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      Subroutine ERTORB ( ncall,kk,kkk,fjd,fn )
c
c     Computations for forces on earth satellite motion due to
c     earth or moon shadow, radiation`667 pressure (direct and reflected),
c     satellite thrustors and sensors, etc.
c     Ash / Amuchastegui - October 1969
c     Rick Abbot December 1984, modifications for GPS
c     King / Ho / McClusky, 1993-94, new radiation pressure model
c     King, 1994, solid earth tide effects (Rotchacher model--not used)
c     King June 1995 - simplify (and change) constants for radiation pressure models
c     Tregoning September 1995 - remove the rock4 component from all but the SRXYZ
c                                models.  Also, remove the lambda dependance from all
c                                but the direct rad. pressure term and the x and z
c                                terms in the SRXYZ (original ROCK4) model.  
c     King  June 1998 - Add the new (CODE9801) Berne radiation pressure model
c     Petrie Oct 2010 - Jan 2011 - Add new models for transmitter thrust, Earth 
c           radiation, and the direct radiation (UCL).  
c     King Jan 2011 - Use arc.h for all commons
c     Petrie Mar 2012 -  Simplify code and add yaw capability
c     King Oct 2015 - Remove obsolete radiation-pressure models and add the
c       2015 extended Empirical CODE Orbit Model as BERN2 elias ECOM2
c
c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
           
c Class of Equation
c          k   = eqn number ( 1 -> neq )
c          kk  = component number for each parameter (1-3)
c                  ( -2 -> 0 for position ; 1 -> 3 for velocity )
c          kkk = eqn group  ( 0 for motion, 1 -> 9 for partials


      implicit none

      include '../includes/dimpar.h'    
      include '../includes/global.h'
      include '../includes/arc.h'
                
      integer*4 kk,kkk,ncall,i,j

      real*8 anga,angb,dt,bt,xt,yt,zt
     .     , fn,dot
     .     , angvec,u,cosu,sinu,u0,cosuu0,sinuu0,cos2uu0,sin2uu0
     .     , cos4uu0,sin4uu0,lambda
     .     , rvec(3),zvec(3),xvec(3),bvec(3),vvec(3),vmag,hvec(3)
     .     , nvec(3),esvec(3),savec(3),spvec(3)
     .     , satsunvec(3),satsunbfsvec(3),satsunbfsnvec(3) 
     .     , beta,distfct,d0,xt1,xt3
c*    .     , k2,fomeg0,asat,anomin,wper0,eccen,fincl
     .     ,decyrs,fjd 
                       
      character*1 ecltyp    
      character*256 message
      character*7  paramtype
      
      logical firstcall,newsv
      integer*4 iprnsave

      real*8 zaxis(3)/0.d0,0.d0,1.d0/
C     Estimate of emitted power from satellite antenna
      real*8 antpwr

      integer*4 iyr,idoy,imonth,iday,ihr,imin  
      
c     Numerical GPS block code now temporary for UCL routines
        
C     Yaw: yaw origin vector,nominal,cosine of nominal yaw, yaw from yaw table (actual yaw of satellite as calculated by yawtab),
C     difference between nominal yaw and actual yaw, cos and sin of difference. 
      integer*4 ievent,iepoch 
      real*8 wvec(3),yatt_nom,cosyatt_nom,yatt_tab,yaw_error,coserr,
     .       sinerr,xvec_yaw(3),yvec_yaw(3)
C     Acceleration along rvec, xvec vectors from earth radiation
      real*8 rvecxvec(2)  
C     Acceleration from TUM earth radiation models, model choice
      real*8  TUMEvec(3)
      integer*4 erm,grd    
C     U-U0
      real*8 uu0
C     Satellite-earth-sun angle (esvec dot rvec)
      real*8 angphicos, angphi
C     earth sun dist ^2 km
      Real*8 rc2r
C     Direction cosines matrix for converting accels to HCL (Height, cross track, along track)
      Real*8 dircos(3,3),eradeci(3),eradhcl(3)
C     Berne or Marek fourier SRP
      Integer*4 srptype
C     Marek fourier accelerations
      Real*8 srpxacc,srpyacc,srpzacc
      Real*8 radtest(3),radtestf(3),dircos1(3,3),radprsSRP(3)
     .       ,radprsfea(3),radout(3)
           
c**   Function
      logical kbit 

      logical debug/.false./,debugsrp/.false./,debugfn/.false./      

      save distfct,d0,lambda,zvec,bvec,xt,zt,cosu,sinu,ecltyp
c     This flag for the first call overall, reporting parameterization
      data firstcall /.true./ 
      save firstcall                                                  
c     This flag for the first call for each SV, reporting radiation scaling
      data newsv /.true./,iprnsave/0/
      save newsv

*** temp

      if( iprn.ne.iprnsave ) then
         newsv = .true.
         iprnsave = iprn
      else
         newsv = .false.
      endif
                          
      if( ncall ) 201,101,1001
c
c-----------Set up once per step----------------------------------------
c
  101 continue    
           
      go to 9999
c
c-----------Set up once per iteration of a given step for motion------
c
  201 if( ncall.lt.-1 ) go to 601

c      Solid-Earth tides according to M. Rotchacher's thesis (this model is not used
c      --rather the IERS standard model using spherical harmonics is added instead)

c     k2 = 0.3d0
c     tfactm = k2*gm(2)*ertrad**5/rpc3/rsb2**2
c     tfacts = k2*gm(3)*ertrad**5/rc3/rsb2**2
c     coszm =  dot(sbcor,pccor)/(rsb*rpc)
c     coszs =  -dot(sbcor,ccor)/(rsb*rc)
c     p2pm = 3.d0*coszm
c     p2ps = 3.d0*coszs
c     p3pm = 3.d0*(5.d0*coszm**2 - 1.d0)/2.d0
c     p3ps = 3.d0*(5.d0*coszs**2 - 1.d0)/2.d0


c ........................Input vectors
c
c    sbcor(6) : coordinates of satellite wrt Earth
c    bcor(3)  : coordinates of satellite wrt Sun   
c    ccor(6)  : coordinates of Earth wrt Sun
           
c     debug for E, SV, Sun vectors
      if ( kbit(idbprt,1).and. kbit(idbprt,4)) then
         call pjdhms(fJD,iyr,idoy,imonth,iday,ihr,imin)
         write(idebug,*) 'PJD Date ',fjd,iyr,idoy,imonth,iday,ihr,imin
         write(idebug,'(a,3f12.3,3f8.4)') 'SV wrt Earth  ',sbcor
         write(idebug,'(a,3f12.0)') 'SV wrt Sun    ',bcor
         write(idebug,'(a,3f12.0)') 'Earth wrt Sun ',(ccor(i),i=1,3)
      endif

c ........................Axis naming conventions
c      
c    The following are unit vectors used to define the primary axes for 
c    the radiation pressure models:
c
c    SNVEC points points from the sun to the satellite and defines the
c       Sun or direct (D-) axis; the sign is consistent with Fleigel et al. but
c       opposite to Springer et al. [1998]).   
c    RVEC points from the Earth to the satellite; it is the negative of the
c       SV-body Z-axis (ZVEC) which points to the Earth.
c    XVEC is the SV-body X-axis, which points toward the half-plane that
c       contains the Sun. 
C        N.B. Kouba (2009 GPS Solut) says that for IIR satellites X-axis points 
C        away from the sun. However, in GAMIT, values and codes are arranged so that the IIR
C        x-axis is defined consistently with the other satellites.
c    YVEC is the SV-body Y-axis, which completes a right-hand-system with
c       XVEC and ZVEC (viz Z = X x Y); it points along the axis of the solar panels
c       and can be formed from RVEC (or ZVEC) and SNVEC:  YVEC = RVEC x SNVEC or 
c       YVEC = SNVEC x ZVEC.   Note that Springer et al. [1998] define the Y-axis 
c       with opposite sign, deriving it by taking Z as the primary ("X") axis.  
c    BVEC completes an orthogonal right-hand-system with D and Y in the 
c       Berne models.  Springer et al. [1998] derive it from D x Y, which is our
c       -SNVEC x -YVEC, or equivalently, SBVEC x YVEC in our derivation.  However,
c       we define it as YVEC x SNVEC (misunderstanding the sign convention in 
c       Beutler et al. [1994]).  With this scheme, our D/Y/B system is left-handed
c       and all three axes have sign opposite to that given by Beutler et al. [1994]
c       and Springer et al. [1998].
            
c    The following are unit vectors used to determine the angle of the Sun 
c    above the orbital plane (beta) and the elongation of the satellite from 
c    the Sun (u - u0) used in the Berne models:
c    
c    ZAXIS is (0. 0. 1.) in the Earth-fixed inertial system
c    VVEC points in the direction of the satellite's motion
c    HVEC is the orbit normal in the direction the angular momentum (R x V) ("up in a RH sense)
c    NVEC is the ascending node of the satellite in the Earth-fixed inertial system
c    ESVEC points from the center of the Earth to the Sun (= approx. -SNVEC)
c    SAVEC is an auxilliary direction perpendicular to ESVEC and HVEC, in the orbital plane
c    SPVEC is the projection of the Sun in the orbital plane (= HVEC x SAVEC) 
c
c    Beta is the angle between ESVEC and SPVEC, positive up (toward HVEC), and thus
c    the angle between the Sun and the orbital plane, or the Sun's latitude in the
c    orbital frame.
c
c    U is the angle between the satellite and node in the orbital plane, thus it's
c    argument of latitude.  U0 is the angle between the Sun's projection onto the
c    orbital plane (SPVEC) and the satellite's node.  Thus U - U0 is the elongation
c    of the satellite from the Sun's projection in the orbital plane.  

c     
c ...............Compute unit vectors as described above.
c                             
      vmag = dsqrt( dot(sbcor(4),sbcor(4)) )
      do i = 1,3     
        rvec(i) = sbcor(i)/rsb  
        vvec(i) = sbcor(i+3)/vmag
        zvec(i) = -rvec(i)  
        snvec(i) = bcor(i)/rb   
c       get a unit vector pointing from the Earth to the Sun for purposes 
c       of computing beta and u - u0 for the Berne models
        esvec(i) = -ccor(i)/rc
      enddo       
      if(debug) print *,'rvec vvec snvec ',rvec,vvec,snvec
      call cross_unit(vvec,zvec,hvec)
      call cross_unit(rvec,snvec,yvec)
      call cross_unit(yvec,snvec,bvec)
      call cross_unit(rvec,yvec,xvec)
c......Compute the satellite-sun vector in the satellite body frame for
C      UCLR2 (grid model)
      if (modrad.eq.8) then
          do i=1,3
            satsunvec(i)=-snvec(i)
          enddo
          satsunbfsvec(1)=dot(satsunvec,xvec)
          satsunbfsvec(2)=dot(satsunvec,yvec)
          satsunbfsvec(3)=dot(satsunvec,zvec)
          call normalise(satsunbfsvec,satsunbfsnvec)
      endif   

c .....Compute the B-angle between the sun and the SV-BODY +Z axis for 
c      the Earth radiation models

      anga = dot(snvec,rvec)
      angb = dacos(anga)

c .....Compute the beta-angle between the sun and the orbital plane for BERN2
       
c     For debug, use Keplerian elements as a rough check (though u is not exactly the same)
c     For the Berne once-per-rev terms, compute the mean amomoly (M)
c     and argument of perigee (w) ), and argument of latitude (u)
c      call keplr(sbcor,asat,eccen,wper0,fomeg0,fincl,anomin)
cd      print *,'w O i M ',wper0,fomeg0,fincl,anomin
c      u = wper0 + anomin
c      if ( u .gt. twopi ) u = u - twopi
c      cosu = dcos(u)
c      sinu = dsin(u)  
cd      print *,'sinu cosu u 1 ',sinu,cosu,u
          
c      beta is the angle (between -90 and +90) of the Sun above (+) the orbital plane
       beta = twopi/4.d0 - dacos(dot(hvec,esvec))
      if(debug)  print *,'esvec hvec twopi temp beta 1 '
     .       ,esvec,hvec,twopi,u,beta
c     get the node direction 
      call cross_unit(zaxis,hvec,nvec)
c     then get a perpendicular to the orbit normal and Sun (will be in orbit plane)
      call cross_unit(hvec,esvec,savec)
c     finally get the projection of the Sun in the orbit plane
      call cross_unit(savec,hvec,spvec)
c     u is the angle from the node to the SV direction, positive if N X R points toward H
      u = angvec(nvec,rvec,hvec) 
      cosu = dcos(u)
      sinu = dsin(u)  
      if(debug) print *,'sinu cosu u 2 ',sinu,cosu,u
c     u0 is the angle from the SV node to the Sun's projection, postive if N x SP points toward H
      u0 = angvec(nvec,spvec,hvec)
C     for testing and plotting -------- calculate angle u-u0, (elongation
c     of the satellite from the Sun's projection in the orbital plane).
      uu0=u-u0
      if ((u-u0).lt.0.d0) uu0=uu0+twopi
c     precalculate cosine and sine
      cosuu0 = dcos(u-u0)
      sinuu0 = dsin(u-u0)
      cos2uu0 = dcos(2.d0*(u-u0))
      sin2uu0 = dsin(2.d0*(u-u0))
      cos4uu0 = dcos(4.d0*(u-u0))
      sin4uu0 = dsin(4.d0*(u-u0))
                                           
C     Calculate satellite-earth-sun angle (radians)     
      angphicos = dot(esvec,rvec)
      angphi = dacos(angphicos)

C-----------------------------------------------------------------------------------------------
c   Binary-code options for debug prinout (hidden in cols 14-16 of afname line  
c      Bit #  Value
c        1     1   print at all epochs
c        2     2   print only during eclipses
c        3     4   shadowing values 
c        4     8   sun, moon, satellite positions
c        5    16   radiation-pressure quantities 
c        6    32   **switch on earth radiation quantities**
c        7    64   print SRP quantities 
c        8   128   ** switch on antenna thrust **
c        9   256   print earth radiation and antenna thrust quantities
c       10   512   second run of arc, now with yaw table input
C-----------------------------------------------------------------------------------------------
C     Calculate yaw angle if debug bit set, needs to be when 
C     running through arc with yawtable available, i.e. second call
      if ( iyaw.gt.0 .or. kbit(idbprt,10)) then
C     Calculate yaw angle error:
C       Nominal yaw angle (code taken from the last part of lib/yaw_attit.f)

C         Calculate yaw origin vector: the vector (wvec) in the orbital plane which is orthogonal
c         to zvec is then the cross product of zvec and hvec
C         wvec is in the approximate direction of vvec.
          call cross(zvec,hvec,wvec)
c
c         Now compute the (ideal or nominal) yaw angle. (counter-clockwise rotation of W around Z)
          cosyatt_nom=dot(xvec,wvec) 
C         cos, so shouldn't be greater than 1
          if (cosyatt_nom .gt. 1.d0) cosyatt_nom=sign(1.0d0,cosyatt_nom)
C         Convert to angle 
          yatt_nom = dacos(cosyatt_nom)
C         Define yatt_nom as having the opposite sign to beta (give yatt_nom sign of beta, then reverse sign)
          yatt_nom = -sign(yatt_nom,beta)
C         Check the output matches the Bar-Sever 1996 J Geod. output
          if(debug) print*,'yatt, beta, mu'
     .       ,yatt_nom*360/twopi,beta*360/twopi,(uu0*360/twopi)-180
C         Check if yaw table exists


C         yaw angle from tables (think fjd is equivalent to jdobs*1.d0+tobs/86400.d0, as called from model.f )
          yatt_tab=0.d0
          
          if(iyawtab.gt.0) call satatt(iyawtab,fjd,yatt_tab,nsats
     .                   , ievent,iepoch)
C         yaw angle error
CC---------------------------------------------------------------------------------
C code adapted from model/svbody_coords.f to calculate yaw error and adjust the vectors (code for partials is there too, though is not active and have not added here)


c         if no input yaw file, set actual yaw = ideal yaw
          if( iyawtab.eq.0 ) yatt_tab = yatt_nom

c         Now compute the yaw error (being yatt from the yaw table minus ideal or nominal yaw)
          yaw_error=yatt_tab-yatt_nom

C Debug
          yaw_error=5.d0/360.d0*twopi

C         Think this should make sure the yaw error is within plus or minus pi (+-180 degrees): 
ccc       yaw_error=yaw_error-nint(yaw_error/360.0d0)*360.0d0  
          yaw_error=yaw_error-nint(yaw_error/twopi)*twopi
c**       debug to look for floating point exception
c         if( dabs(yatt).lt.0.1d0 .or. dabs(yatt).gt.360.d0 .or.
c     .       dabs(ideal_yaw).lt.0.1d0 .or. dabs(ideal_yaw).gt.360.d0 ) then
c         endif
cd         if( dabs(yaw_error).lt.1.d-12) print *,'yaw error',yaw_error  
          if(debug) print *,'yaw_error ',yaw_error

C     effects on satellite axis unit vectors
c** Now adjust the X and Y coordinate by the yaw error. Here we can assume that the
c   yaw angle from the yaw table (yatt_tab) already contains the correct adjustments 
c   for any yaw manoeuvres etc etc. Therefore, irrespective of what the satellite
c   has done, if there is a difference between yatt and ideal_yaw the sc axes xhat
c   and yhat must be adjusted. 
          if(dabs(yaw_error).gt.1.0d-3)then      
            coserr=dcos(yaw_error)
            sinerr=dsin(yaw_error)
            do i=1,3
              xvec_yaw(i)=xvec(i)*coserr+yvec(i)*sinerr
              yvec_yaw(i)=yvec(i)*coserr-xvec(i)*sinerr 
            enddo 
          else
            do i=1,3
              xvec_yaw(i) = xvec(i)
              yvec_yaw(i) = yvec(i)   
            enddo
          endif
C Debug
cd      print*,'yaw_error',yaw_error*360.d0/twopi
cd      print*,'xvec xvec yaw',xvec,xvec_yaw
cd      print*,'xvec xvec yaw',xvec,xvec_yaw
CCC----------------------------------------------------------------------

      endif
                                            

c.....The following models are available for non-gravitational forces ('radiation 
c.....pressure') invoked according to the value of 'modrad', set in ARC:
       
c  RWK NOTE Oct 2015/Apr 2019: We have removed all of the old models that are no longer
c  used (SPHRC, SRDYZ, SRXYZ, SRDYB (modrad 0-3), BERN1 (modrad 5), leaving only
c  the 1998 BERNE (now modrad 1, alias of ECOM1) and substituting for BERN2 (modrad 6)  
c  the 2015 extended ECOM (ECOM2) (now modrad 2).  Additional mod May 2019 introduces
c  'ECOMC' (MIT designation) which integrates both the ECOM1 once-per-rev DYB 
c  terms plus the ECOM2 D 2-per-rev and 4-per-rev terms. Although there is a
c  distinction at the estimation stage between ECOM2 (9 parameters) and ECOMC
c  (13 parameters), we carry all 13 terms through ARC, MODEL, and SOLVE. 
c
c modrad
c
c  1 ECOM1(formerly BERNE) : Emperical CODE Orbit Model [Beutler et al.,Man. Geod. 19, 367, 1994]
c      radcon(1): [DRAD]: coeff. of (radiation) away from sun (nominally 1.0)
c      radcon(2): [YRAD]: coeff. of (radiation) in sv (+) y-axis direction (nom. 0.)
c      radcon(3): [BRAD]: coeff. of (radiation) in the bvec direction (nom. 0.0)
c      radcon(4): [DCOS]: periodic once per rev cosine part in Sun direction
c      radcon(5): [DSIN]: periodic once per rev sine part in Sun direction
c      radcon(6): [YCOS]: periodic once per rev cosine part in y-axis direction
c      radcon(7): [YSIN]: periodic once per rev sine part in y-axis direction
c      radcon(8): [BCOS]: periodic once per rev cosine part in b-axis directcion
c      radcon(9): [BSIN]: periodic once per rev sine part in b-axis direction
c  2 ECOM2/ECOM2 : Extended Emperical CODE orbit Model [Arnold et al., J. Geod., 89, 775, 2015]
c      radcon(1-9) same as ECOM1
c      radcon(10): [DCS2]: periodic twice per rev cosine in Sun direction
c      radcon(11): [DSN2]: periodic twice per rev sine  in Sun direction
c      radcon(12): [DCS4]: periodic four per rev cosine in Sun direction
c      radcon(13): [DSN4]: periodic four per rev sine in Sun direction
c  modread 3-9 reserved for JPL and other a priori models
c  7 UCLR1 : New fourier series model from University College London; ECOM1 parameterization
c  8 UCLR2 : New grid based model from University College London; ECOM1 parameterization
c   
c....Shadowing :  Tests carried out by Paul Tregoning at Scripps and Simon McClusky 
c   at MIT demonstrated that ARC orbits fit better with IGS orbits when only the
c   direct solar radiation terms are turned off during eclipse.  Hence only these
c   terms and the ROCK4 x- and y-axis terms are multiplied by the shadow factor
c   'lambda' in the accerlations and partials.  These tests may be flawed by the
c   fact that shadowing, even with some tapering, introduces a discontinuity in the
c   accelerations that may itself introduce errors into the integration.  The tests
c   should be repeated if additional tapering is coded.
    

c....See if the Earth or Moon is partially or wholly obstructing the Sun (lamba)

      lambda=1.d0
      call shadow(lambda,ecltyp) 
      if(debug) then 
        print *,'After SHADOW fjd beta ',fjd,beta 
        print *,' ecltyp lambda ',ecltyp,lambda 
      endif 
      call shwprt(fjd,lambda,ecltyp,beta)


c....Compute the distance factor for the radiation force 
      
      distfct = (aunit**2)/rb2 
   
c ....Zero out all correction terms
      d0 = 0.d0
      dt = 0.d0 
      xt = 0.d0
      yt = 0.d0
      zt = 0.d0
      bt = 0.d0
      xt1 = 0.d0
      xt3 = 0.d0
 
c .... Compute the coefficients based on reflectivity and area-to-mass ratio 


c .....Now compute the accelerations according to the model specified

c        The coefficient given by Spinger et al. are accelerations in units of 
c        1.e-9 m/s**2, or 1.e-12 km/s**2.  For the direct term, which in GAMIT provides 
c        the units for adjustments to coefficients representing accelerations about 
c        other axes, we compute the value to match that coded here for BERNE in Dec, 
c        1994. namely, the nominal ROCK4 force divided by the S/C mass.  The ROCK4 force 
c        is given in units of 1.e-5 nt (kg m /s**2); dividing by the S/C mass (now in kg) 
c        and multiplying by 1.e-3 gives acceleration in km/s**2, the units of ARC.
                        
      if(debug) print *,'antbody ',antbody 
      if( antbody(1:9).eq.'BLOCK I  ' ) then
        d0 = 4.54d-8 / sbmass
c         = 103.488e-9 m/s**2  for mass of most recent Blk I (521.8 kg) 
      elseif ( antbody(1:9).eq.'BLOCK II ' .or.   
     .          antbody(1:9).eq.'BLOCK IIA' ) then
         d0 = 8.695d-8 / sbmass
c          = 98.449e-9 m/s**2 for most Blk IIs (mass 883.2 kg; exceptions are 
c                              PRNs 2, 14, 15, 21; cf. CODE9801 mean of 99.428)
c          = 89.372e-9  for Blk IIA (mass 972.9 kg, cf. CODE 9801 mean of 90.3  )
      elseif ( antbody(1:9).eq.'BLOCK IIR' ) then
           d0 = 11.15d-8 / sbmass
c             = 101.364e-9 m/s**2 for Blk IIR (no value yet available for Berne models)    
c          For IIF, multiply IIR constant by 1.5 (from MIT IGS runs) to get the direct 
c          coefficient ~ 1.0.  Since the mass is 1.416 times larger, this implies that 
c          the effective cross-sectional area is 2.12 times larger,
      elseif( antbody(1:9).eq.'BLOCK IIF' ) then
        d0 = 16.7d-8 / sbmass    
      elseif( antbody(1:10).eq.'BLOCK IIIA') then                  
****We don't yet now the effective cross section, so put in a dummy value
****and scale from the estimates. rwk 190109 
        d0 = 1.d-7 / sbmass                  
** The DRAD value for 2019 016 was 1.566, so change the coefficient to 1.5d-7 
        d0 = 1.5d-7/sbmass 
*--------------------------------------------------------------------------------
c          ** Beidou and Galileo cross-sectional areas are derived from a scaling
c             of estimates from rwk sp3 fits to IGSF 5 April 2018. The come out
c             to be approximately the same as the mass ratio, implying that the
c             area/mass ratio for all of the SVs is about the same.
*--------------------------------------------------------------------------------
      elseif( antbody(1:9).eq.'BEIDOU-2M' .or.
     .        antbody(1:9).eq.'BEIDOU-3M' ) then
c             Factor for the MEO is 0.8 times the BLOCK IIF value; since the mass
c             ratio is 0.64, the scaling is 0.64*0.8 = 0.5 
        d0 = 8.35d-8 / sbmass
      elseif( antbody(1:9).eq.'BEIDOU.2I' .or. 
     .        antbody(1:9).eq.'BEIDOU.2G' ) then 
c             Factor for the GEOs is 3 times the BLOCK IIF value; since the mass
c             ratio is 3, the scaling is 3.
        d0 = 50.1d-8 / sbmass
      elseif( antbody(1:9).eq.'GALILEO-1' .or.
     .        antbody(1:9).eq.'GALILEO-2' ) then
c             Factor is 0.5 times the BLOCK IIF value; since the mass ratio 
c             is 0.5, the scaling is 0.5
        d0 = 8.35d-8 / sbmass
      elseif( antbody(1:9).eq.'GLONASS  ' .or.
     .        antbody(1:9).eq.'GLONASS-M' ) then
c       temporary: Glonass mass is nearly the same as IIF
c       d0 = 16.7d-8 / sbmass
c       TAH 190702: Increased value to make DRAD=~1 (25% change) 
        d0 = 20.9d-8 / sbmass 
      elseif( antbody(1:10).eq.'GLONASS-K1' ) then
c       TAH: 190702: This value looks about correct (i.e., DRAD=~1.0)
        d0 = 10.0d-8 / sbmass
      else   
        d0 = 16.7d-8 / 1555.256d0
        if( newsv ) then                                 
            write(message,'(a,a16,a)') 'Radiation scaling unknown for '
     .                   ,antbody,': use GPS Block IIF value' 
            call report_stat('WARNING','ARC','ertorb',' ',message,0)
          endif
      endif
      do i=1,3
         radprs(i)=0.d0
      enddo                           
      if(debug) print *,'sbmass d0 ',sbmass,d0 

c .... Compute the accelerations for the model specified 
                

      if( modrad.eq.1 ) then
c       ECOM1 (formerly BERNE)
        do i= 1,3     
c         For now, apply shadowing only to the direct constant term
          radprs(i) = radprs(i)+(lambda*(radcon(1)*d0+dt)* snvec(i)
     .                 +          ( radcon(2) * d0 + yt ) * yvec(i)
     .                 +          ( radcon(3) * d0 + bt ) * bvec(i) 
     .                  ) * distfct
* MOD TAH 190615: Replaced the u argument with uu0 so that BERNE and ECOMC
*         once per rev terms have the same phase (i.e., cosu -> cosuu0 etc.
          radprs(i) = radprs(i) +
     .                   ( radcon(4) * d0 * cosuu0 * snvec(i)
     .                   + radcon(5) * d0 * sinuu0 * snvec(i)
     .                   + radcon(6) * d0 * cosuu0 * yvec(i)
     .                   + radcon(7) * d0 * sinuu0 * yvec(i)
     .                   + radcon(8) * d0 * cosuu0 * bvec(i)
     .                   + radcon(9) * d0 * sinuu0 * bvec(i) )*distfct 
         enddo                 
         if(debugsrp) 
     .     print *,'modrad radcon cosu sinu snvec yvec bvec radprs '
     .            , modrad,radcon,cosu,sinu,snvec,yvec,bvec,radprs

      elseif ( modrad.eq.2 ) then
c       ECOM2 or ECOMC 
        do i=1,3
c         For now, apply shadowing only to the direct constant term
          radprs(i) = radprs(i)+(lambda*(radcon(1)*d0+dt)* snvec(i)
     .          +          ( radcon(2) * d0 + yt ) * yvec(i)
     .          +          ( radcon(3) * d0 + bt ) * bvec(i) 
     .           ) * distfct
c         ECOM2 uses arg from solar node
          radprs(i) = radprs(i) +
     .          (  radcon(4) * d0 * cosuu0 * snvec(i)
     .          +  radcon(5) * d0 * sinuu0 * snvec(i)
     .          +  radcon(6) * d0 * cosuu0 * yvec(i)
     .          +  radcon(7) * d0 * sinuu0 * yvec(i)
     .          +  radcon(8) * d0 * cosuu0 * bvec(i)
     .          +  radcon(9) * d0 * sinuu0 * bvec(i) 
     .          +  radcon(10)* d0 * cos2uu0 * snvec(i)
     .          +  radcon(11)* d0 * sin2uu0 * snvec(i)
     .          +  radcon(12)* d0 * cos4uu0 * snvec(i)
     .          +  radcon(13)* d0 * sin4uu0 * snvec(i) ) * distfct
         enddo  

         if(debugsrp) 
     .     print *,'modrad radcon cosuu0 sinuu0 snvec yvec bvec radprs '
     .            , modrad,radcon,cosuu0,sinuu0,snvec,yvec,bvec,radprs
      elseif (modrad.eq.7 .or. modrad.eq.8) then
c       UCLR1 and UCLR2  UCLR1 

C         EJP option to test novel SRP setups  28 Mar 2012.
C  could set up this choice to be done by a debug bit
C  Choose SRP model - 1, UCL fourier series model (UCLR1)
C                     2, UCL grid model (UCLR2)
C         srptype = 1
C         srptype = 2

        if (modrad.eq.7) then
C Compute UCL fourier SRP model accelerations
C        This only works for certain blocks at present
          call srpfmod(antbody,sbmass,angb,srpxacc,srpyacc,srpzacc)
        endif
        if (modrad.eq.8) then
C Compute UCL grid model SRP model accelerations
cd       Print*,'ertorb: calling UCL grid model'
           call srpfgrid(satsunbfsnvec,srpxacc,srpyacc,srpzacc)
        endif
C 
        if (modrad.eq.7)  then
C  Calculate SRP accelerations using UCL fourier SRP models
C             add fourier SRP models accelerations (km s-2)
C             multiply by distfct to adjust for distance from sun 
           do i=1,3 
             radprs(i) = ((xvec(i)*srpxacc)
     .                   +(yvec(i)*srpyacc)
Cd     .                   +(zvec(i)*srpzacc))
     .                   +(zvec(i)*srpzacc))*distfct*lambda
             radtestf(i)=radprs(i)
           enddo

        elseif (modrad.eq.8) then
c  Calculate SRP accelerations using UCL grid SRP models (identical code to fourier here at present, may combine )
c   add SRP models accelerations (km s-2)
C             multiply by distfct to adjust for distance from sun 
           do i=1,3
             radprs(i) = ((xvec(i)*srpxacc)
     .                   +(yvec(i)*srpyacc)
Cd     .                   +(zvec(i)*srpzacc))
     .                   +(zvec(i)*srpzacc))*distfct*lambda

             radtestf(i)=radprs(i)
           enddo
cd             Print*,'ertorb: code to add UCL grid model accelerations'
        endif
      else 
         call report_stat('FATAL','ARC','ertorb',' ',
     .         'Unknown solar raditaion pressure (SRP) model',0)
      endif
CC    End of different SRP models section (modrad)
      if(debug) print *,'After srpmod radprs ',srpmod,radprs


C-------------------------------------------------------------------------------------

C      Other non-gravitational forces
CC     Earth radiation
      if (eradmod.ne.'NONE '.or.kbit(idbprt,6)) 
     .     then
          if (eradmod.eq.'NCLE1 ') then 
cd              print*, 'NCLE1'
C             Compute earth radiation accelerations according to model coded by EJP
C             Face and wing, analytical earth radiation - see earthrad.f
              rc2r=rc*rc
C               Calculate satellite-earth-sun angle (radians) (calculated earlier)    
cd                angphicos = dot(esvec,rvec)
cd                angphi = dacos(angphicos)
cd             print *,'ERTORB antbody ',antbody
              call earthrad(angb,angphi,twopi,sbmass,antbody,ltvel,
     .               aunit,rc2r,
     .               rsb,rvecxvec,idbprt,idebug)
              if(debug) then 
                print*,'ARC/ertorb earthrad on'
                print *,'radprs ',radprs
                print *,'Direct solar ',(radprs(i),i=1,3)
              endif 
              do i=1,3
C               Add earth radiation accelerations 
                radprs(i) = radprs(i) + rvec(i)*rvecxvec(1)
                radprs(i) = radprs(i) + xvec(i)*rvecxvec(2)                 
                if(idbprt.gt.0) write(idebug,*) 'add earth accels NCLE1'
              end do           
              if( debug ) print *,'Direct + Albedo ',(radprs(i),i=1,3) 
          elseif (eradmod.eq.'TUME1 '.or.eradmod.eq. 'TUME2 ')  then
C              Compute earth radiation accelerations according to TUM model. TUME1 similar to NCLE1, TUME2 CERES model.
               if (eradmod.eq.'TUME1 ') erm = 1
               if (eradmod.eq.'TUME2 ') erm = 2
               grd=3
C               GRD       : 1 = 2.5 degree grid
C                         : 2 = 5.0 degree grid
C                         : 3 = 10  degree grid
               if(debug)  print*,'erm',erm,'grd',grd
               call earthradTUM(erm,grd,fjd,TUMEvec)     
              if( debug ) print *,'tumevec'
     .          ,gnss,iprn,fjd,tumevec
               do i=1,3
C                 Add earth radiation accelerations 
                  radprs(i) = radprs(i) + TUMEvec(i) 
                  if(debug) print *,'add earth accels TUM'
               end do    
          elseif (eradmod.eq.'NONE ' ) then 
            if(newsv) call report_stat('WARNING','ARC','ertorb',' '
     .          ,'Earth radiation pressure model = NONE ',0)     
          else
            call report_stat('FATAL','ARC','ertorb',eradmod,
     .          'Unknown earth radiation model',0)          
          endif
      endif   
      if(debug) print *,'After eradmod radprs ',eradmod,radprs


C      Antenna thrust
           
      antpwr = 0.d0
      if (antradmod.eq.'ANT77 '.or. kbit(idbprt,8) .or. 
     .                                antradmod.eq.'ANTBK ' ) then
c       Add antenna thrust (antpwr is radiated energy in watts, 
c       1 watt = 1 kg m^2/s^3
        if (antradmod.eq.'ANT77 '.or. kbit(idbprt,8)) then
              antpwr = 77.d0
        elseif (antradmod.eq.'ANTBK ') then
c         Block specific antenna thrust (see http://acc.igs.org/orbits/thrust-power.txt)
c         Compares to Carlos Rodriguez Solano implementation
c         Use Block IIA value for Block I and II
cc         if( iblock.eq.1 ) then  
          if( antbody(1:9).eq.'BLOCK I  ' ) then
            antpwr = 76.d0
cc          elseif ( iblock.eq.2.or.iblock.eq.3 ) then 
          elseif( antbody(1:9).eq.'BLOCK II '  .or.
     .            antbody(1:9).eq.'BLOCK IIA' ) then 
            antpwr = 76.d0
cc            elseif ( iblock.eq.4.or.iblock.eq.5 ) then
          elseif ( antbody(1:11).eq.'BLOCK IIR-A' .or.
     .             antbody(1:11).eq.'BLOCK IIR-B' ) then
            antpwr = 85.d0
cc          elseif( iblock.eq.6 ) then
          elseif( antbody(1:11).eq.'BLOCK IIR-M') then
c           IIR-M: this includes the M-code power which is nominally switched on 
c           (would be 108W without M code)           
            antpwr = 198.d0 
cd            elseif( iblock.eq.7 ) then
          elseif( antbody(1:9).eq.'BLOCK IIF' ) then   
c           M-code was switched off for SVN62/PRN25 on 05Apr2011 JD 2455656.5 PJD 2455657) 
            antpwr = 249.d0       
c           (would be 154W without M code)
cc            if (iprn.eq.25 .and. fjd.ge.2455657.d0 ) then
            if( sname(7:8).eq.'62' .and.fjd.ge.2455657.d0) then
               antpwr = 154.d0
             endif     
          elseif( antbody(1:10).eq.'BLOCK IIIA' ) then
**** This is a preliminary estimate of the power based on an IGS AC email
**** from Peter Steigenberg 190104. ---rwk
             antpwr = 300.d0 
cd             if(firstcall)  print *,'sname antpwr ',sname,antpwr 
          elseif( antbody(1:9).eq.'GLONASS-M' ) then
             antpwr = 100.d0 
          else                 
            if(newsv) then                       
              write(message,'(a,a20)') 
     .            'Antenna thrust model not available for ',antbody
              call report_stat('WARNING','ARC','ertorb',' ',message,0)
             endif 
          endif                                                                 
        else
          call report_stat('WARNING','ARC','ertorb',antradmod,
     .                  'Unknown antenna thrust model',0)
        endif
Cd           print*,'anttr BERNE', (antpwr)
Cd            write(idebug,*) 'anttr BERNE', (antpwr)/(sbmass*ltvel)
        if(debug)  write(*,*) 'antenna force '
     .     ,(antpwr*rvec(i)/(ltvel*1.d6),i=1,3)
        
        do i=1,3
c         Add to accelerations
          radprs(i) = radprs(i) + (antpwr*rvec(i))/(sbmass*
     .                ltvel*1000.d0*1000.d0)
        end do
      endif
      if( newsv ) then                                     
        write(message,'(3a,f6.2,a)') 'Antenna power for ',antbody(1:12)
     .      ,'= ',antpwr,' watts'
         call report_stat('STATUS','ARC','ertorb',' ',message,0) 
      endif                                   
      if(debug) print *,'After antradmod radprs '
     .                 ,antradmod,radprs
                  

C------------------------------------------------------------------------------------
CCC For testing, calculate the earth radiation accelerations in an HCL frame
      if (kbit(idbprt,7).and.eradmod.ne.'NONE ') then 
         do i=1,3
C          set up initial acceleration matrix:

           if (eradmod.eq.'NCLE1 ') then 
           eradeci(i)=rvec(i)*rvecxvec(1)+xvec(i)*rvecxvec(2)
           elseif (eradmod.eq.'TUME1 '.or.eradmod.eq. 'TUME2 ')  then
           eradeci(i)=tumevec(i)
           else
               call report_stat('FATAL','ARC','ertorb',eradmod,
     .           'Unknown earth radiation model hcl',0)          
          endif
C          set up the direction cosine matrix
C          radial
           dircos(1,i)=rvec(i)
C          alongtrack
           dircos(2,i)=vvec(i)
C          crosstrack
           dircos(3,i)=hvec(i)
         enddo
C        Multiply together to get output accelerations
         call matmpy(dircos,eradeci,eradhcl,3,3,1)
      endif
CCC-----------------------------------------------------------------
CCC  Have calculated accelerations for SRP, AT, and ER as specified
CCC  Now set up for parameterisation as desired: 
      if (modrad.eq.7 .or. modrad.eq.8) then
C       set up the direction cosine matrix to convert SRP to DYB axes
         do i=1,3
            dircos1(1,i)=snvec(i)
            dircos1(2,i)=yvec(i)
            dircos1(3,i)=bvec(i)
         enddo
C        Multiply together to get output accelerations
C        Only SRP
         call matmpy(dircos1,radtestf,radprsSRP,3,3,1)
C        SRP, erthrad, antenna thrust, if switched on.
         call matmpy(dircos1,radprs,radout,3,3,1)
C
CCCCCC Control for trying different parameterisations and reporting it to GAMIT.status
C           For now, use DYB parameterisation (nb, in BERNE sinu and cosu are multiplied by d0)
C        paramtype = "UCLRx_2"
         paramtype = "UCLRx_1"
         if (firstcall ) then
           call report_stat('STATUS','ARC','ertorb',paramtype,
     .       'SRP parameterisation',0)
         endif 

         if ( paramtype == "UCLRx_1" ) then
            do i=1,3
             radprs(i) = radcon(1)*radout(1)*snvec(i) +
     .			 radcon(2)*radout(2)*yvec(i) +
     .			 radcon(3)*radout(3)*bvec(i) +
C   D axis (snvec)
     .                       radcon(4) * radout(1) * cosu *snvec(i)
     .                     + radcon(5) * radout(1) * sinu *snvec(i)
C   Y axis (yvec)
     .                     + radcon(6) * radout(2) * cosu *yvec(i)
     .                     + radcon(7) * radout(2) * sinu *yvec(i)
C   b axis (bvec)
     .                     + radcon(8) * radout(3) * cosu *bvec(i)
     .                     + radcon(9) * radout(3) * sinu *bvec(i)
            enddo
          
         elseif ( paramtype == "UCLRx_2" ) then
            do i=1,3
             radprs(i) = radcon(1)*radout(1)*snvec(i) +
     .			 radcon(2)*radout(2)*yvec(i) +
     .			 radcon(3)*radout(3)*bvec(i) +
C   D axis (snvec)
     .                       radcon(4) * radout(1) * cosu *snvec(i)
     .                     + radcon(5) * radout(1) * sinu *snvec(i)
C   Y axis (yvec)
     .                     + radcon(6) * radout(1) * cosu *yvec(i)
     .                     + radcon(7) * radout(1) * sinu *yvec(i)
C   b axis (bvec)
     .                     + radcon(8) * radout(1) * cosu *bvec(i)
     .                     + radcon(9) * radout(1) * sinu *bvec(i)
            enddo
         elseif ( paramtype == "UCLRx_3" ) then
            do i=1,3
             radprs(i) = radcon(1)*radout(1)*snvec(i) +
     .			 radcon(2)*radout(2)*yvec(i) +
     .			 radcon(3)*radout(3)*bvec(i) +
C   D axis (snvec)
     .                       radcon(4) * d0 * cosu *snvec(i)
     .                     + radcon(5) * d0 * sinu *snvec(i)
C   Y axis (yvec)
     .                     + radcon(6) * d0 * cosu *yvec(i)
     .                     + radcon(7) * d0 * sinu *yvec(i)
C   b axis (bvec)
     .                     + radcon(8) * d0 * cosu *bvec(i)
     .                     + radcon(9) * d0 * sinu *bvec(i)
            enddo          
         endif
C      end of parametrisation if modrad = 7 or 8
      endif

C -------------- Test output ---------------------------------------
C      SRP quantities
       if(((kbit(idbprt,1)).or.(kbit(idbprt,2))).and.kbit(idbprt,7))then
          call pjdhms(fJD,iyr,idoy,imonth,iday,ihr,imin)
          call jd_to_decyrs(fJD+0.5,decyrs)
C     List of output: 
C       prn,block number
C       [satellitename (e.g. PRN 28)],
C       year, day of year, hour, minute, 
C       [fractionaljulianday,month, day],
C       nominal satellite mass 
C       beta angle(angle of the Sun above the orbital plane)
C       U(angle from node to SV direction), U0(angle from the SV node to the Sun's projection), U-U0
C       angb(sun-satellite-earth angle), angphi(satellite-earth-sun angle ), 
C       nominal satellite mass
C       Three radiation pressure accels
C       Satellite X, Y, Z unit vectors(inertial)and 
C       4 more unit vectors: bvec,vvec,snvec,esvec(see above). Angles in degrees.
C       earth radiation: height alongtrack crosstrack accelerations in kms-2.
C       antenna thrust / watts
C       Also get  'SHADOW r1 r2 phi thet hgt area1 area2 area3 ari ' from shadow.f in the file - use 'grep -v SH' to remove.
          write(idebug,*) iprn,antbody,
C     .       sname,
C     .       fjd,
c     .       decyrs,iyr,idoy,ihr,imin,
C     .       imonth,iday, 
c     .       sbmass,
     .       beta*360.d0/twopi,
     .       u*360.d0/twopi,u0*360.d0/twopi,uu0*360.d0/twopi,
     .       angb*360.d0/twopi,angphi*360.d0/twopi,
     .       sbmass,
     .       radprs,
     .       xvec,yvec,zvec,
     .       bvec,vvec,snvec,esvec,
     .       eradhcl,
     .       antpwr
       endif
cd     endif on test output
Cdd       Write (idebug,*) 'radprs 1 to 3',radprs
Cd       Print*,'radcon 1-9',radcon(1),radcon(2),radcon(3),radcon(4),
Cd     .      radcon(6),radcon(7),radcon(8),radcon(9)
Cd       Print*,'-----------------------'
C-------------------------------------------------------------------
C-------------------------------------------------------------------
c-----end of condition on choice of model group
                       
                

      if(debug) then    
        print *,'End of SRP ERAD calculations fjd ',fjd 
        print *,'r s x y z b beta ',rvec,snvec,xvec,yvec,zvec,bvec,beta
        print*, 'rvec',rvec
        print*,' beta',beta
        print*, 'radprs ',radprs
        print *,'srpmod antbody lambda distfct sbmass'
     . , srpmod,antbody,lambda, distfct,sbmass
        print *,'dt yt bt xt1 xt3 zt ',dt,yt,bt,xt1,xt3,zt
      endif  


      firstcall = .false.
      go to 9999


c-----------Set up once per iteration of a given step for partials------
c
  601 continue
      go to 9999

c-----------Compute acceleration for motion------------------------

 1001 continue
      if( kkk.gt.0 ) go to 2001

c        Solid-earth tides from Berne model (not used--see SBFN, SBFN1)

cd      print *,'kk fn ',kk, fn
cd       fn = fn + tfactm*(p2pm*pccor(kk)/rpc - p3pm*sbcor(kk)/rsb)
cd      print *,'kk fn ',kk, fn
cd       fn = fn + tfacts*(-p2ps*ccor(kk)/rc - p3ps*sbcor(kk)/rsb)
cd      print *,'kk fn ',kk, fn
cd      if(kk.eq.3) 
c        Radiation pressure
                

      fn = fn + radprs(kk)
      if(debugfn) print *,'fjd kkk kk radprs ',fjd,kk,kk,radprs(kk)
cd      stop 
                         
                   
      if( debug.or.debugsrp ) then
         print *,'srpmod distfct radcon snvec yvec bvec radprs'
     .  ,srpmod,distfct,radcon,snvec,yvec,bvec,radprs
cd          stop
      endif
        
      go to 9999

c-----------compute acceleration for partials------------------------

 2001 continue

c     Radiation pressure

c    --Note that we express all of the parameters as a fraction of the the nominal 
c      direct-radiation acceleration (d0 for BERN2, radforce for all others) (EJP Oct 2010 - actually radforce not used)
                               
      if(debugsrp) print *,'ERTORB rad partials modrad kk kkk '
     .   ,modrad,kk,kkk 
c     Check if partial is wrt Sun direction

      if( kkk.eq.7 ) then
        if ( modrad.eq.7 .or. modrad.eq.8) then
           fn = fn + radout(1)*snvec(kk)
        else 
           fn = fn + lambda * d0 * snvec(kk) * distfct
        endif       
        if(debugsrp) print *,' axis 1 distfct lambda d0 snvec '
     .                        ,distfct,lambda,d0,snvec(kk)
                                     

c     Check if partial is wrt Y-axis 

      elseif( kkk.eq.8 ) then
         if ( modrad.eq.7 .or. modrad.eq.8) then
           fn = fn + radout(2)*yvec(kk)
         else
           fn = fn + d0 * yvec(kk) * distfct  
         endif                       
         if(debugsrp) print *,' axis 2 distfct d0 yvec '
     .     ,distfct,d0,yvec(kk)

c     Check if partial is wrt B-axis

      elseif (kkk.eq.9) then          
        if ( modrad.eq.7 .or. modrad.eq.8) then
          fn = fn + radout(3)*bvec(kk) 
        else 
          fn = fn + d0 * bvec(kk)  * distfct  
        endif
        if(debugsrp) print *,' axis 3 distfct d0 bvec '
     .     ,distfct,d0,bvec(kk)
           
c     Parameters 10-15 for the CODE and UCL models are once-per-rev terms in D, Y, B
* MOD TAH 190615: replaced cosu -> cosuu0 sinu -> sinuu0 to be constistent
*     with the BERNE and ECOMC OPR force integration.
      elseif (kkk.eq.10) then    
         if(modrad.eq.7 .or. modrad.eq.8)  then
           if ( paramtype == "UCLRx_1" ) then
             fn = fn + radout(1) * cosu *snvec(kk)
           elseif ( paramtype == "UCLRx_2" ) then
             fn = fn + radout(1) * cosu *snvec(kk)
           elseif ( paramtype == "UCLRx_3" ) then
             fn = fn + radout(1) * cosu *snvec(kk)
           endif            
         else
* MOD TAH 180915 replace u with uu0
C          fn = fn + cosu * d0 * snvec(kk)  * distfct   
           fn = fn + cosuu0 * d0 * snvec(kk)  * distfct   
         endif 
         if(debugsrp) print *,' 1pr D distfct d0 cosuu0 snvec'
     . ,distfct,d0,cosuu0,snvec(kk)
  
      elseif (kkk.eq.11) then
        if(modrad.eq.7 .or. modrad.eq.8)  then
          if ( paramtype == "UCLRx_1" ) then
            fn = fn + radout(1) * sinu *snvec(kk)
          elseif ( paramtype == "UCLRx_2" ) then
            fn = fn + radout(1) * sinu *snvec(kk)
          elseif ( paramtype == "UCLRx_3" ) then
            fn = fn + radout(1) * sinu *snvec(kk)
          endif
        else
* MOD TAH 180915 replace u with uu0
          fn = fn + sinuu0 * d0 * snvec(kk)  * distfct
        endif    
        if(debugsrp) print *,'1pr D distfct d0 sinuu0 snvec'
     .    ,distfct,d0,sinuu0,snvec(kk)
  

      elseif (kkk.eq.12) then
        if(modrad.eq.7 .or. modrad.eq.8)  then
          if ( paramtype == "UCLRx_1" ) then
             fn = fn + radout(2) * cosu *yvec(kk)
          elseif ( paramtype == "UCLRx_2" ) then
            fn = fn + radout(1) * cosu *yvec(kk)
          elseif ( paramtype == "UCLRx_3" ) then
            fn = fn + d0 * cosu *yvec(kk)
          endif
        else 
* MOD TAH 180915 replace u with uu0
          fn = fn + cosuu0 * d0 * yvec(kk)  * distfct   
      endif   
      if(debugsrp) print *,'1pr Y distfct d0 cosuu0 yvec'
     .  ,distfct,d0,cosuu0,yvec(kk)
  

      elseif (kkk.eq.13) then
        if(modrad.eq.7 .or. modrad.eq.8)  then
          if ( paramtype == "UCLRx_1" ) then
            fn = fn + radout(2) * sinu *yvec(kk)
          elseif ( paramtype == "UCLRx_2" ) then
            fn = fn + radout(1) * sinu *yvec(kk)
          elseif ( paramtype == "UCLRx_3" ) then
            fn = fn + d0 * sinu *yvec(kk)
          endif
        else
* MOD TAH 180915 replace u with uu0
          fn = fn + sinuu0 * d0 * yvec(kk)  * distfct
        endif     
        if(debugsrp) print *,'1pr Y distfct d0 sinuu0 yvec'
     .    ,distfct,d0,sinuu0,yvec(kk)


      elseif (kkk.eq.14) then
        if(modrad.eq.7 .or. modrad.eq.8)  then
          if ( paramtype == "UCLRx_1" ) then
             fn = fn + radout(3) * cosu *bvec(kk)
          elseif ( paramtype == "UCLRx_2" ) then
             fn = fn + radout(1) * cosu *bvec(kk)
          elseif ( paramtype == "UCLRx_3" ) then
             fn = fn + d0 * cosu *bvec(kk)
          endif
        else
* MOD TAH 180915 replace u with uu0
          fn = fn + cosuu0 * d0 * bvec(kk)  * distfct
        endif 
        if(debugsrp)  print *,'1pr B distfct d0 cosuu0 yvec'
     .     ,distfct,d0,cosuu0,bvec(kk)
                        

      elseif (kkk.eq.15) then
        if(modrad.eq.7 .or. modrad.eq.8)  then
          if ( paramtype == "UCLRx_1" ) then
            fn = fn + radout(3) * sinu *bvec(kk)
          elseif ( paramtype == "UCLRx_2" ) then
            fn = fn + radout(1) * sinu *bvec(kk)
          elseif ( paramtype == "UCLRx_3" ) then
            fn = fn + d0 * sinu *bvec(kk)
          endif
        else
* MOD TAH 180915 replace u with uu0
          fn = fn + sinuu0 * d0 * bvec(kk)  * distfct
        endif 
        if(debugsrp) print *,'1pr B distfct d0 sinu bvec'
     .     ,distfct,d0,sinu,bvec(kk)
cd        if(debug.or.debugsrp) stop 1 

* MOD TAH 190615: Added distfct (distance function) scaling to 
*     partials to be consistent with force model (all 4 terms)
      elseif (kkk.eq.16) then
        if( modrad.ne.2 ) call report_stat('FATAL','ARC','ertorb'
     .      ,'SRP parameters 10-13 defined only for ECOM2',0)
        fn = fn + cos2uu0 * d0 * snvec(kk) * distfct
      elseif (kkk.eq.17) then
        fn = fn + sin2uu0 * d0 * snvec(kk) * distfct
      elseif (kkk.eq.18) then
        fn = fn + cos4uu0 * d0 * snvec(kk) * distfct
      elseif (kkk.eq.19) then  
        fn = fn + sin4uu0 * d0 * snvec(kk) * distfct
      endif
c     end of condition on partial number

9999  continue
      return
      end





