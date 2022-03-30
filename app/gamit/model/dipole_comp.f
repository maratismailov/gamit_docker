      subroutine dipole_comp( evec0,xrot,r1,r1leng,freql1,freql2,ichan
     .                      , xhat_t,yhat_t,dipoldel,iepoch )
c
c     PURPOSE: compute the transmitter / receiver antenna orientation dependent
c              phase corrections for Right Circularly Polarized electro-magnetic waves
c              Reference Wu J.T., et al., Manuscripta Geogetica (1993) 18, pp91-98...
c
c     Variables:
c               IN:  evec0    - terrestrial frame receiver coordinates                 R*8(3,2)
c                    xrot     - rotation matrix from terrestrial to inertial           R*8(3,3)
c                    r1       - vector of satellite wrt L1 receiver phase centre       R*8(3)
c                    r1leng   - distance from receiver L1 phase centre to satellite CM R*8
c                    freql1,freql2 - frequencies of the l1 and l2 carriers.            R*8
c                    ichan    - satellite slot number on X/C file                      I*4
c
c               OUT: dipoldel - correction to computed delay due to transmitter /
c                               receiver antenna dipole orientations.                  R*8(2)
c
c SUBROUTINES CALLED: rotate_geod, dot, cross, matmpy,
c
c CREATED: 02 DEC 1994               LAST MODIFIED: 02 DEC 1994
c
c AUTHOR: S McClusky.
c
c COPYRIGHT: DEPARTMENT OF EARTH AND PLANETRY SCIENCES
c            M.I.T. 1994
c
      implicit none

      include '../includes/dimpar.h'
c
      character*3 xyz, neu
      integer ichan,j,n,iepoch
      logical first_call,first(maxsat)
      real*8  xhat_neu(3),yhat_neu(3)
      real*8  xhat_r(3),yhat_r(3),khat(3)
      real*8  xhat_rt(3),yhat_rt(3)
      real*8  xhat_t(3),yhat_t(3)
      real*8  rxyz(3)
      real*8  rllh(3),rotmat(3,3),xrot(3,3)
      real*8  dot,pi,one
      real*8  dphi,dphi_prev(maxsat),dotdtdr
      real*8  d_r(3),d_t(3)
      real*8  pos_neg,len_d_r,len_d_t
      real*8  tmp1,tmp2(3),tmp3(3),tmp4(3)
      real*8  dipoldel(2)
      real*8  freql1,freql2
      real*8  evec0(3,2),r1(3),r1leng
c
      parameter ( pi = 3.1415926535897932D0 )
      parameter ( one = 1.D0 )

      save first, dphi_prev
c
      data xyz/'XYZ'/, neu/'NEU'/
      data xhat_neu/ 1.0, 0.0, 0.0/, yhat_neu/0.0,-1.0, 0.0/
      data first_call/.true./
                     
c     Initialize the flags for each SV
      if( first_call ) then
        do j=1,maxsat
          first(j) = .true.
        enddo
        first_call = .false.
      endif

c-----Convert receiver antenna dipole unit vectors from local NEU coordinates
c     to geocentric earth fixed cartesian vectors.
c
      do j = 1,3
        rxyz(j) = evec0(j,1) * 1000.d0
      enddo

      call rotate_geod(xhat_neu,xhat_rt,neu,xyz,rxyz,rllh,rotmat)
      call rotate_geod(yhat_neu,yhat_rt,neu,xyz,rxyz,rllh,rotmat)
c
c-----Convert receiver antenna diopole unit vectors from earth fixed
c     rotating frame to an inertial non rotating frame.
c
      call matmpy(xrot,xhat_rt,xhat_r,3,3,1)
      call matmpy(xrot,yhat_rt,yhat_r,3,3,1)
c
c
c-----Compute the Unit vectors pointing from transmitter to receiver...
c
      do j = 1,3
        khat(j) =  -r1(j)/r1leng
      enddo

c-----Compute effective dipoles for a transmitter and receiver, and
c     compute the phase correction between the 2 effective dipoles.
c
c-----Compute the receiver effective Dipole trem:
c     D_r = xhat_r - khat(khat . xhat_r) + khat x yhat_r
c     Where: xhat_r & yhat_r = the receiver antenna unit vectors in the directions
c                              of the two dipole elements.
c          : khat = the unit vector pointing from the transmitter to the receiver
c
      tmp1 = dot(khat,xhat_r)
      call matmpy(khat,tmp1,tmp2,3,1,1)
      call cross(khat,yhat_r,tmp3)

      do j=1,3
        d_r(j) = xhat_r(j) - tmp2(j) + tmp3(j)
      enddo

c-----Compute the transmitter effective Dipole trem:
c     D_t = xhat_t - khat(khat . xhat_t) - khat x yhat_t
c     Where: xhat_t & yhat_t = the transmitter antenna unit vectors in the directions
c                              of the two dipole elements.
c          : khat = the unit vector pointing from the transmitter to the receiver
c
      tmp1 = dot(khat,xhat_t)
      call matmpy(khat,tmp1,tmp2,3,1,1)
      call cross(khat,yhat_t,tmp3)

      do  j=1,3
        d_t(j) = xhat_t(j) - tmp2(j) - tmp3(j)
      enddo

      len_d_t = dsqrt( dot(d_t,d_t) )
      len_d_r = dsqrt( dot(d_r,d_r) )

c   normalises d_t and d_r
       do  j=1,3
         d_t(j) = d_t(j)/len_d_t
         d_r(j) = d_r(j)/len_d_r
       enddo


c-----Compute the SIGN of the phase correction
c     SIGN = khat . (D_t x D_r)
c
      call cross(d_t,d_r,tmp4)
      pos_neg = dot(khat,tmp4)
c
c
c-----Compute the phase correction dphi
c     dphi = SIGN cos-1(D_t . D_r )
c                          
c** tah/rwk 080122: The combination of three intrinsics is causing a problem
c                   with gfortran at optimization level 3, resulting in dphi=NaN
c                   when the dot product is close to -1.000000.  It works without
c                   optimization, and it works in a test program, but not here.
c                   The test is 2008 day 002, BOR1 channel 6, epoch 1230
c         dphi = dsign(one,pos_neg) * dacos( dot(d_t,d_r))    
cd       if( iepoch.eq.1230.and.ichan.eq.6) then            
cd         print *,'khat tmp4 pos_neg ',khat,tmp4,pos_neg
cd         print *,'d_t len  d_r len dphi ',d_t,len_d_t,d_r,len_d_r,dphi
         dotdtdr =  dot(d_t,d_r)   
         if( dotdtdr.lt.-0.999999d0.and.dotdtdr.gt.-1.000001 ) then
            dotdtdr = -1.d0  
cd            print *,'dotdtdr set to -1.d0'
         elseif( dotdtdr.gt.0.99999d0.and.dotdtdr.lt.1.00001 ) then
            dotdtdr = 1.d0            
cd            print *,'dotdtdr set to +1.d0'
         endif
         dphi = dsign(one,pos_neg) * dacos( dotdtdr)  
cd         print *,'dotdtdr,dphi ',dotdtdr,dphi     
cd       endif
c** end 080122 workaround

c
c-----Compute the integer part of the phase correction
c     N = nint[ (DPHI_previous - dphi) / (2*pi) ]
c
      if (first(ichan)) then
        n = 0
        dphi_prev(ichan) = 0 ! set the previous dphi to 0 for first run
c        jpl_prev(ichan) = 0
        first(ichan) = .false.
      else
        n = nint( ((dphi_prev(ichan) - dphi) / (2*pi)) )
      endif
c
c-----Compute the total phase correction
c     DPHI = 2 * N * pi + dphi
c
      dphi =  2 * n * pi + dphi
c
c debug by PT - adding in a call to a JPL routine to compute the phase_wind_up
c to test Simon's code

c       call phase_wind_up(yhat_r,xhat_r,khat,yhat_t,xhat_t,
c     .                     jpl_prev(ichan), 0,0,correc,0,0)


c-----Save dphi in dphi_previous form next computation
c
      dphi_prev(ichan) = dphi
c      jpl_prev(ichan) = correc

c-----Convert Phase correction (unitless) into L1 and L2 path delay
c     correction is required in seconds.
c
c   ** I think this is reversed sign, changed as trial by simon 92/12/2
c      dipoldel(1) = -1.d0* ((dphi / (2*pi)) * (vlight/freql1) ) / vlight
c      dipoldel(2) = -1.d0* ((dphi / (2*pi)) * (vlight/freql2) ) / vlight
c   ** which means this is the original sign
c      dipoldel(1) =  ((dphi )) * (vlight/freql1)
c      dipoldel(2) =  ((dphi )) * (vlight/freql2)

       dipoldel(1) = dphi/(2*pi*freql1)
       dipoldel(2) = dphi/(2*pi*freql2)


cd      print*, 'PATH DELAY CORRECTION'
cd      print*, 'dipoldel(1),dipoldel(2) ',dipoldel(1),dipoldel(2)
cd      print*, '  '
c

c       write(*,123)ichan,dphi/(2*pi),correc,(correc-dphi/(2*pi))
c123    format('Channel ',i3,' Simon ',f8.5,' JPL ',f8.5,f10.5,/)
c       print*,' Channel ',ichan,' dipoldel = ',dipoldel

c     THE END
c
      return
      end




