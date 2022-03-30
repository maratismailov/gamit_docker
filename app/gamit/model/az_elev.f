Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1995.   All rights reserved.

         subroutine az_elev(fstelev, naxis,uaxis,xrot,evec0,r1,twopi
     .                     , ichan,azim,elev,nhat_neu_i,uhat_neu_i)
c
c PURPOSE: subroutine to compute the elevation angle and azimuth
c          from a station to satellite.
c
c PARAMETERS:
c         IN:
c             fstelev : T if this is the first call, F otherwise     LOGICAL
c                       <note: this is to save computation time - it works
c                              just as well if it is true all the time>
c             naxis   : local north unit vector                      R*8(3)
c             uaxis   : local up unit vector                         R*8(3)
c             xrot    : rotation matrix efixed-inertial              R*8(3,3)
c             evec0   : terrestrial station coordinates              R*8(6,2)
c             r1      : difference vector of station and satellite
c                       coordinates                                  R*8(3)
c             twopi   : value of two*pi constant                     R*8
c             ichan   : satellite array index                        I*4
c
c        OUT: azim    : terrestrial azimuth from north - station to satellite    R*8(maxsat)
c             elev    : terrestrial elevation angle of satellite from station    R*8(maxsat)
c         nhat_neu_i  : local north unit vector in inertial space                R*8(3)
c         ehat_neu_i  : local east unit vector in inertial space                 R*8(3)
c         uhat_neu_i  : local up unit vector in inertial space                   R*8(3)
c
c SUBROUTINES CALLED: cross, dot
c
c CREATED: 27th DEC 1993               LAST MODIFIED: 14 JUL 1995
c
c AUTHOR: rwk, put in SR by S McClusky.
c
c  modified 28th June 95 by P Tregoning/S McClusky: The azimuth and elevation angles MUST
c    be computed using local terrestrial (ie geodetic) latitude. To this end, we will rotate
c    the site/sat vector into local NEU space, and then compute the elevation and azimuth
c    from site to satellite. To get the rotation matrix from inertial to NUE space we
c    rotate the local unit vectors into inertial space (this gives us a byproduct which
c    can be used in dipole_comp) and form the 3x3 rotation matrix from these direction
c    cosines
c
         implicit none
c
      include '../includes/dimpar.h'
c
         integer*4 ichan,j
c


      real*8 rxyz(3),rotmat(3,3),xrot(3,3)
     .      ,nhat_neu_i(3),uhat_neu_i(3),evec0(6,2),nhat_rt(3)
     .      ,uhat_rt(3),ehat_neu_i(3),elev(maxsat),twopi
     .      ,r1(3),azim(maxsat),rllh(3),r1hat(3),r1leng,amag3,naxis(3)
     .      ,uaxis(3),r1hat_neu(3),inert_neu(3,3)

      logical fstelev

c bookkeeping: get evec0 in right units and normalise r1
      r1leng = amag3(r1)
      do 10 j = 1,3
        rxyz(j) = evec0(j,1) * 1000.d0
        r1hat(j) = r1(j)/r1leng
10    continue

c  --------- do this part once only -----------------------
cSimon removed this if statement because it doesn't work!!!!!!
c       if(fstelev)then  
c RWK temporary dummy statement to avoid compiler warning
      if( .not.fstelev ) then
         print *,'DUMMY fstelev ',fstelev
      endif


c  only need to rotate N,U unit vectors into inertial space once - the
c  matrix doesn't change from satellite to satellite

c Convert receiver antenna N,U unit vectors from local NEU coordinates
c to geocentric earth fixed cartesian vectors (uhat_rt,nhat_rt)
c
        call rotate_geod(naxis,nhat_rt,'NEU','XYZ',rxyz,rllh,rotmat)
        call matmpy(rotmat,uaxis,uhat_rt,3,3,1)
c
c Convert receiver antenna N,U unit vectors from earth fixed
c rotating frame to an inertial non rotating frame.
c
        call matmpy(xrot,nhat_rt,nhat_neu_i,3,3,1)
        call matmpy(xrot,uhat_rt,uhat_neu_i,3,3,1)

c  compute local E direction cosine in inertial space
        call cross(nhat_neu_i,uhat_neu_i,ehat_neu_i)

c  form up a rotation matrix to go from inertial to local NEU coords.
c  It is just the NEU unit vectors (in inertial) as columns in a 3x3 matrix

c        do 20 j=1,3
c          inert_neu(j,1) = nhat_neu_i(j)
c          inert_neu(j,2) = ehat_neu_i(j)
c          inert_neu(j,3) = uhat_neu_i(j)
c20      continue

        do 20 j=1,3
          inert_neu(1,j) = nhat_neu_i(j)
          inert_neu(2,j) = ehat_neu_i(j)
          inert_neu(3,j) = uhat_neu_i(j)
20      continue


c -------   must always compute from here down for all satellites ----------

c  now rotate the site/sat vector into NEU
        call matmpy(inert_neu,r1hat,r1hat_neu,3,3,1)

c  compute the elevation angle
      elev(ichan) = twopi/4.d0 - dacos(r1hat_neu(3))
c       print*,' elev comp = ',elev(ichan)

c compute it the other way
c      cosz = dot(uhat_neu_i,r1hat)
c       print*,' cosz = ',cosz
c      elev(ichan) = twopi/4.d0 - dacos(cosz)
c      print*,' 2nd elev comp = ',elev(ichan)
c      stop
c  compute the azimuth
       azim(ichan) = datan2(r1hat_neu(2),r1hat_neu(1))
       if( azim(ichan).lt.0.d0 ) azim(ichan) = azim(ichan) + twopi
c       print*,' azimuth = ',azim(ichan)*360.d0/twopi
c       stop

         return
         end

