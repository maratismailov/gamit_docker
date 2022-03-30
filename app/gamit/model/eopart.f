Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994. All rights reserved.

      Subroutine eopart ( xpole, ypole, sidtm, prec, rnut, prot, sidten
     1                  , diff,evec0, polepart,ut1part )
c
c PURPOSE: subroutine to compute partial derivatives for earth rotation parameters
c          with respect to inertial site coordinates.
c
c PARAMETERS:
c         IN: xpole   : xpole rotation angle (arcs)                R*8
c             ypole   : ypole rotation angle (arcs)                R*8
c             sidtm   : geenwich aparent siderialtime GAST (rad)   R*8
c             prec    : precession rotation matrix                 R*8(3,3)
c             rnut    : nutation rotation matrix                   R*8(3,3)
c             prot    : polar motion rotation matrix               R*8(3,3)
c             sidten  : GAST rotation matrix                       R*8(3,3)
c             diff    : time from mid pt of data span to obs epoch R*8
c             casr    : radians/sec                                R*8
c             pi      : pi                                         R*8
c             evec0   : terrestrial site coordinates (x,y,z)       R*8(3,2)
c
c        OUT: polepart  : partial derivatives of the pole position
c                         and rate wrt inertial site coords.       R*8(3,4)
c             ut1part   : partial derivatives of the ut1-tai
c                         and rate wrt inertial site coords.       R*8(3,4)
c
c       USED: xpolemat  : partial derivative of W wrt xpole        R*8(3,3)
c             ypolemat  : partial derivative of W wrt ypole        R*8(3,3)
c             xpdotmat  : partial derivative of W wrt xpole rate   R*8(3,3)
c             ypdotmat  : partial derivative of W wrt ypole rate   R*8(3,3)
c             ut1taimat : partial derivative of S wrt ut1-tai      R*8(3,3)
c             ut1taidotmat: partial derivative of S wrt ut1-tai rate R*8(3,3)
c
c THEORY: Inertial coordinates X are obtained by rotating terrestrial coordinates x
c         using P; precession matrix, N; nutation matrix, S earth rotation matrix, and
c         W; pole rotation matrix.
c                                        X = P*N*S*W*x
c         We wish to obtain the partial derivatives of X wrt, Xp (xpole), Yp (ypole)
c         Xpdot (xpole rate), Ypdot (ypole rate), ut-utc (ut1-utc), and ut1dot ( ut1-utc rate).
c
c         Knowing:      W = R2 [-Xp] R1 [-Yp] and
c                       S = R3 [sidtm]
c
c         Where:        R1, R2, R3 are standard rotation matricies for rotation about
c                       x,y,z axes respectively.
c
c         We Want: dX/dXp =     P*N*S*dW/dXp*x
c                  dX/dYp =     P*N*S*dW/dYp*x
c                  dX/dXpdot =  P*N*S*dW/dXpdot*x
c                  dX/dYpdot =  P*N*S*dW/dYpdot*x
c                  dX/dut1 =    P*N*dS/dut1*W*x
c                  dX/dut1dot = P*N*dS/dut1dot*W*x
c
c        Where: dW/dXp = ---                    ---
c                        | -sin(-Xp) 0 -cos(-Xp)  |
c                        |      0    0     0      |
c                        |  cos(-Xp) 0 -sin(-Xp)  |
c                        ---                    ---
c               dW/dXpdot is the same matrix multiplied by (T)
c               where T is time to interpolation epoch (days).
c
c               dW/dYp = ---                    ---
c                        |  0       0        0    |
c                        |  0 -sin(-Yp)  cos(-Yp) |
c                        |  0 -cos(-Yp) -sin(-Yp) |
c                        ---                    ---
c               dW/dYpdot is the same matrix multiplied by (T)
c               where T is time to interpolation epoch (days).
c
c               dS/dut1 = ---                     ---
c                         | -sin(gast)  cos(gast) 0 |
c                         | -cos(gast) -sin(gast) 0 |
c                         |     0          0      0 |
c                         ---                    ---
c               dS/dut1dot is the same matrix multiplied by (T)
c               where T is time to interpolation epoch (days).
c
c SUBROUTINES CALLED:
c
c CREATED: 25th APR 1994               LAST MODIFIED:
c
c AUTHOR: S McClusky.
c
            implicit none
c
      include '../includes/dimpar.h'
c

      integer   i,j

      real*8    prec(3,3),rnut(3,3),prot(3,3),sidten(3,3),casr,pi,twopi
     .        , xpolemat(3,3),ypolemat(3,3),xpdotmat(3,3),ypdotmat(3,3)
     .        , ut1taimat(3,3),ut1taidotmat(3,3),diff,xpole,ypole
     .        , radsec,polepart(3,4),ut1part(3,2),evec0(3,2),sidtm
     .        , dws_dxp(3,3),dws_dyp(3,3),dws_dxpdot(3,3)
     .        , dws_dypdot(3,3), dsw_dut1(3,3),dsw_dut1dot(3,3)
     .        , temp(3,3),temp1(3),temp2(3),temp3(3),temp4(3),temp5(3)
     .        , temp6(3)
c
      pi = 4.d0*datan(1.d0)
      twopi = pi * 2.d0
      radsec = twopi/86400.d0
      casr = pi/180.d0/3600.d0


c
c zero out arrays
      do 10 i=1,4
        do 10 j=1,3
          polepart(j,i)=0.0d0
10    continue
      do 20 i=1,2
        do 20 j=1,3
          ut1part(j,i)=0.0d0
20    continue
c
c compute xpole, ypole, xpdot and ypdot pole partial derivatives dW/dXp, dW/dYp
c dW/dXp_dot, dW/dYp_dot saved in matricies xpolemat, ypolemat, xpdotmat, ypdotmat
c
cc
      xpolemat(1,1)= -dsin(-xpole)
      xpolemat(1,2)= 0.0d0
      xpolemat(1,3)= -dcos(-xpole)
      xpolemat(2,1)= 0.0d0
      xpolemat(2,2)= 0.0d0
      xpolemat(2,3)= 0.0d0
      xpolemat(3,1)= dcos(-xpole)
      xpolemat(3,2)= 0.0d0
      xpolemat(3,3)= -dsin(-xpole)

cd      write(*,30) ((xpolemat(i,j),j=1,3),i=1,3)
cd 30    format(' dW/dXp partial matrix :',/,3(1x,3D22.14,/))

      xpdotmat(1,1)= -dsin(-xpole)*diff
      xpdotmat(1,2)= 0.0d0
      xpdotmat(1,3)= -dcos(-xpole)*diff
      xpdotmat(2,1)= 0.0d0
      xpdotmat(2,2)= 0.0d0
      xpdotmat(2,3)= 0.0d0
      xpdotmat(3,1)= dcos(-xpole)*diff
      xpdotmat(3,2)= 0.0d0
      xpdotmat(3,3)= -dsin(-xpole)*diff

cd      write(*,40) ((xpdotmat(i,j),j=1,3),i=1,3)
cd 40    format(' dW/dXpdot partial matrix :',/,3(1x,3D22.14,/))

      ypolemat(1,1)= 0.0d0
      ypolemat(1,2)= 0.0d0
      ypolemat(1,3)= 0.0d0
      ypolemat(2,1)= 0.0d0
      ypolemat(2,2)= -dsin(-ypole)
      ypolemat(2,3)= dcos(-ypole)
      ypolemat(3,1)= 0.0d0
      ypolemat(3,2)= -dcos(-ypole)
      ypolemat(3,3)= -dsin(-ypole)

cd      write(*,50) ((ypolemat(i,j),j=1,3),i=1,3)
cd 50    format(' dW/dYp partial matrix :',/,3(1x,3D22.14,/))

      ypdotmat(1,1)= 0.0d0
      ypdotmat(1,2)= 0.0d0
      ypdotmat(1,3)= 0.0d0
      ypdotmat(2,1)= 0.0d0
      ypdotmat(2,2)= -dsin(-ypole)*diff
      ypdotmat(2,3)= dcos(-ypole)*diff
      ypdotmat(3,1)= 0.0d0
      ypdotmat(3,2)= -dcos(-ypole)*diff
      ypdotmat(3,3)= -dsin(-ypole)*diff

cd      write(*,60) ((ypdotmat(i,j),j=1,3),i=1,3)
cd 60    format(' dW/dYpdot partial matrix :',/,3(1x,3D22.14,/))

c
c compute ut1 and ut1dot (UT1-TAI) partial derivatives dS/dut1 dS/dut1_dot
c saved in matricies ut1utcmat, u1ucdotmat
c
      ut1taimat(1,1)= -dsin(sidtm)
      ut1taimat(1,2)= dcos(sidtm)
      ut1taimat(1,3)= 0.0d0
      ut1taimat(2,1)= -dcos(sidtm)
      ut1taimat(2,2)= -dsin(sidtm)
      ut1taimat(2,3)= 0.0d0
      ut1taimat(3,1)= 0.0d0
      ut1taimat(3,2)= 0.0d0
      ut1taimat(3,3)= 0.0d0

cd      write(*,70) ((ut1taimat(i,j),j=1,3),i=1,3)
cd 70    format(' dS/dut1 partial matrix :',/,3(1x,3D22.14,/))

      ut1taidotmat(1,1)= -dsin(sidtm)*diff
      ut1taidotmat(1,2)= dcos(sidtm)*diff
      ut1taidotmat(1,3)= 0.0d0
      ut1taidotmat(2,1)= -dcos(sidtm)*diff
      ut1taidotmat(2,2)= -dsin(sidtm)*diff
      ut1taidotmat(2,3)= 0.0d0
      ut1taidotmat(3,1)= 0.0d0
      ut1taidotmat(3,2)= 0.0d0
      ut1taidotmat(3,3)= 0.0d0

cd      write(*,80) ((ut1taidotmat(i,j),j=1,3),i=1,3)
cd 80    format(' dS/dut1dot partial matrix :',/,3(1x,3D22.14,/))
c
c From polar motion partials
c xpole:  dX/dxp =  P*N*S*dW/dXp*x
c
       call matmpy(xpolemat,sidten,dws_dxp,3,3,3)
cd      write(*,90) ((dws_dxp(i,j),j=1,3),i=1,3)
cd 90     format(' dWS/dXp :',/,3(1x,3D22.14,/))
       call pns(prec,rnut,dws_dxp,temp)
       call matmpy(temp,evec0(1,1),temp1,3,3,1)
cd       write(*,100)(temp1(j),j=1,3)
cd 100    format(' (dX/dXp)  :',/,3(1x,3D22.14,/))
       do 150 i=1,3
       polepart(i,1) = temp1(i)*casr
150    continue
c
c xpole rate dX/dxpdot =  P*N*S*dW/dxpdot*x
c
       call matmpy(xpdotmat,sidten,dws_dxpdot,3,3,3)
cd       write(*,160) ((dws_dxpdot(i,j),j=1,3),i=1,3)
cd 160    format(' dWS/dXpdot :',/,3(1x,3D22.14,/))
       call pns(prec,rnut,dws_dxpdot,temp)
       call matmpy(temp,evec0(1,1),temp2,3,3,1)
cd       write(*,170)(temp2(j),j=1,3)
cd 170    format(' (dX/dXpdot)  :',/,3(1x,3D22.14,/))
       do 180 i=1,3
       polepart(i,2) = temp2(i)*casr
180    continue
c
c ypole dX/dyp =  P*N*S*dW/dyp*x
c
       call matmpy(ypolemat,sidten,dws_dyp,3,3,3)
cd      write(*,190) ((dws_dyp(i,j),j=1,3),i=1,3)
cd 190    format(' dWS/dYp :',/,3(1x,3D22.14,/))
       call pns(prec,rnut,dws_dyp,temp)
       call matmpy(temp,evec0(1,1),temp3,3,3,1)
cd       write(*,200)(temp3(j),j=1,3)
cd200    format(' (dX/dYp)  :',/,3(1x,3D22.14,/))
       do 210 i=1,3
       polepart(i,3) = temp3(i)*casr
210    continue
c
c ypole rate dX/dypdot =  P*N*S*dW/dypdot*x
c
       call matmpy(ypdotmat,sidten,dws_dypdot,3,3,3)
cd      write(*,220) ((dws_dypdot(i,j),j=1,3),i=1,3)
cd 220    format(' dWS/dYpdot :',/,3(1x,3D22.14,/))
       call pns(prec,rnut,dws_dypdot,temp)
       call matmpy(temp,evec0(1,1),temp4,3,3,1)
cd       write(*,230)(temp4(j),j=1,3)
cd 230    format(' (dX/dYpdot)  :',/,3(1x,3D22.14,/))
       do 240 i=1,3
       polepart(i,4) = temp4(i)*casr
240    continue
c
c Form UT1-UTC partials
c
c ut1-utc  dX/dut1 = P*N*dS/dut1*W*x
c
       call matmpy(prot,ut1taimat,dsw_dut1,3,3,3)
cd       write(*,250) ((dsw_dut1(i,j),j=1,3),i=1,3)
cd250    format(' dSW/dut1 :',/,3(1x,3D22.14,/))
       call pns(prec,rnut,dsw_dut1,temp)
       call matmpy(temp,evec0(1,1),temp5,3,3,1)
cd       write(*,260)(temp5(j),j=1,3)
cd260    format(' (dX/dut1)  :',/,3(1x,3D22.14,/))
       do 270 i=1,3
       ut1part(i,1) = temp5(i)*radsec
270    continue
c
c ut1-utc rate dX/dut1dot = P*N*dS/dut1dot*W*x
c
       call matmpy( prot,ut1taidotmat,dsw_dut1dot,3,3,3)
cd      write(*,280) ((dsw_dut1dot(i,j),j=1,3),i=1,3)
cd 280   format(' dSW/dut1dot :',/,3(1x,3D22.14,/))
       call pns(prec,rnut,dsw_dut1dot,temp)
       call matmpy(temp,evec0(1,1),temp6,3,3,1)
cd       write(*,290)(temp6(j),j=1,3)
cd290    format(' (dX/dut1dot)  :',/,3(1x,3D22.14,/))
       do 300 i=1,3
       ut1part(i,2) = temp6(i)*radsec
300    continue
c
c      print*, 'polepart in eopart', polepart
      return
      end





