      subroutine ptide( ipole,ptd_model,jd,t,geoc,delta)
c     Subroutine to compute the pole tide component
c
c     code based on IERS Technical Note 3
c     written by peterm 1-oct-93
c
c     INPUT PARAMETERS (NONE CHANGED)
c     *******************************  
c
c     ipole      unit number for the wobble file 'pole.' 
* MOD TAH 200219: Removed logical lmptide and replaced with model name.
*     lmptid     logical = .true. if mean tide of 2000 to be removed
c     ptd_model  Name of model to apply (IERS20, IERS20, ZERO, NONE)
c     jd         the julian date of the observation
c     t          the time in sec from the start of the day
c     geoc       an array of the station vector in format
c                latitude, longitude, radius
c                lat and lon in radian measure in the terrestrial frame
c                radius is not used (hence units unimportant).
c                NOTE the distinction between geodetic and geocentric
c                latitude has not been maintained. The error is less
c                than 0.1 mm.
c
c
c     OUTPUT PARAMETERS (NONE DEFINED ON INPUT)
c     *****************************************
c
c     delta      the X,Y,Z Cartesian displacements due to the pole tide.
c                The units are millimeters in the terrestrial frame.

      implicit none
              
      integer*4 jd, ipole

      real *8 geoc(3), delta(3), r(3,3), s(3), xpdot, ypdot, pi,
     .        t, fract, along, colat, jdr, dt
      real*8 xp, yp  ! Pole offsets to be used in pole-tide (seconds)
         

* MOD TAH 200219: Added Mean (secular) pole position
      real*8 mxp, myp  ! Mean (secular) pole position for different models.(mas)
  
C     logical lmptid
      character*8 ptd_model  ! Pole tide model.

c      data CONVD/0.017453292519943296D0/,CONVDS/4.8481368110953599D-6/

      pi = 4.0d0*datan(1.d0)

c     Convert UTC seconds-of-day to fraction-of-day

      fract = t/86400.d0

c     Read the pole-position table and compute pole position
c     in arc-seconds.
      call POLRED( ipole,jd,fract,xp,yp,xpdot,ypdot )

c     Remove the mean pole position  
c       Old code uses IERS 1996 standards with a linear term only
c       New code (110524) uses IERS 2010 standards with a cubic model
c         prior to 2010 and a linear term afterward 
c       Both refer to 2000.0 as a reference epoch
* MOD TAH 200219: Replace code with new more complete model
C     if( lmptid ) then   
C       jdr = dfloat(jd) + fract
C       dt = (jdr-2451544.5d0)/365.25d0   ! time in years
c*         Linear trend from IERS Conventions 1996
c*         xp = xp  - (0.054d0+0.00083d0*dt)
c*         yp = yp  - (0.357d0+0.00395d0*dt) 
c         IERS 2010 model is in two parts
C        if( jdr.lt.2455197.5d0 ) then   ! Data before 2010
c            Cubic trend from IERS Conventions 2010.
C            xp = xp - ( 55.974d0 + 1.8243d0*dt + 0.18413d0*dt**2 
C    .                                      + 0.007024d0*dt**3 )/1000.d0
C            yp = yp - ( 346.346d0 + 1.7896d0*dt - 0.10729d0*dt**2
C    .                                      - 0.000908d0*dt**3 )/1000.d0
C        else     ! Linear trend after 2010.00
C            xp = xp - (  23.513d0 + 7.6141d0*dt )/1000.d0
C            yp = yp - ( 358.891d0 - 0.6287d0*dt )/1000.d0
C         endif
C     endif    
* MOD TAH 200219: See what model is requested
      if( ptd_model(1:4).ne.'NONE' ) then
         jdr = dfloat(jd) + fract + 0.5d0  ! PEPJD to standard JD
         call mean_pole( jdr, ptd_model, mxp, myp )
         xp = xp - mxp/1000.d0
         yp = yp - myp/1000.d0
      else  ! No pole tide; so set xp and yp zero
         xp = 0.d0
         yp = 0.d0
      endif

     
cd      print *,'jdr dt xp yp ',jdr,dt,xp,yp
cd     stop

c     compute the colatitude of the site

      colat = pi/2.0d0 - geoc(1)

c     reassign the longitude
      along = geoc(2)

c     compute the rotation matrix in its transposed position
      r(1,1) = dcos(colat)*dcos(along)
      r(2,1) = dcos(colat)*dsin(along)
      r(3,1) = -dsin(colat)
      r(1,2) = -dsin(along)
      r(2,2) = dcos(along)
      r(3,2) = 0.0d0
      r(1,3) = dsin(colat)*dcos(along)
      r(2,3) = dsin(colat)*sin(along)
      r(3,3) = dcos(colat)


c     compute the s vector. caution parameters are in millimeters and
c     arcseconds.
c** pjm code
c     s(1)  = -32.d0*dsin(2.d0*colat)*(xp*dcos(along) - yp*dsin(along))
c     s(2)  =  -9.d0*dcos(2.d0*colat)*(xp*dcos(along) - yp*dsin(along))
c     s(3)  =   9.d0*dcos(colat)     *(xp*dsin(along) + yp*dcos(along))
c** end pjm code

c** Dong correction -- assumes S, E, U; rotation matrix tranforms correctly to X Y Z 
C     s(1)  =  -9.d0*dcos(2.d0*colat)*(xp*dcos(along) - yp*dsin(along))
C     s(2)  =   9.d0*dcos(colat)     *(xp*dsin(along) + yp*dcos(along))
C     s(3)  = -32.d0*dsin(2.d0*colat)*(xp*dcos(along) - yp*dsin(along)) 
* MOD TAH 20140406: Updated coefficients to be more accurate and consistent
*     with the IERS2010 stanards.  (See kf/gen_util/comp_ptide.f for
*     detailed explanation
* NOTES: s(1) here is South motion (dtheta in IERS conventions),  In
*     GLOBK, dNorth is used and hence opposite sign,
      s(1) =  -8.94d0*dcos(2.d0*colat)*(xp*dcos(along) - yp*dsin(along))
      s(2) =   8.94d0*dcos(colat)     *(xp*dsin(along) + yp*dcos(along))
      s(3) = -33.20d0*dsin(2.d0*colat)*(xp*dcos(along) - yp*dsin(along)) 

      call matmpy (r,s,delta,3,3,1)
c
c     the array delta contains the pole tide corrections
c     values are in the terrestrial system in cartesian coordinates
c
CD     write (*, 10) (delta(i),i=1,3)
CD    format (1x, 'The pole tide components are (mm): ',3f9.3 )

      return
      end
