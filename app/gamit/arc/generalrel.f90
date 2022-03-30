    subroutine generalrel(J2, pos, vel, sunpos, angvel, gm, Ae, Me, Ie, c, grelacc)

! Emma-Kate Potter, 16 Sept 2010
! This program calculates the general relativistic contributions to acceleration

! The relativistic correction to the acceleration of an artificial Earth satellite taken from
! IERS conventions. 
! NOTE that the equation given in the IERS conventions 2003 documents has an error and the 
! updated draft documents (2010) are used instead.

     implicit none
     integer :: i
     real(kind=8) :: r, v, c, beta, gam, rs , Ie, Ae, Me
     real(kind=8) :: vdotv, rdotv, rdotJ 
     real(kind=8), dimension(3) ::  pos, vel, angvel, J, gm
     real(kind=8), dimension(3) ::  posS, velS
     real(kind=8), dimension(3) ::  rcrossv, vcrossJ, bigcross
     real(kind=8), dimension(3) ::  VScrossRS, RScrossv
     real(kind=8), dimension(3) ::  grelacc
     real(kind=8), dimension(6) ::  sunpos
! PT140205: general relativity accelerations due to J2
     real(kind=8), dimension(3) ::  a1,a3
     real(kind=8) :: J2,a2
     logical j2flag/.false./

!     print *,'Constants (c, gm, Ae, Me, Ie): ',c, gm, Ae, Me, Ie

!    CALCULATE distance from earth to satellite
     r = dsqrt(pos(1)*pos(1)+pos(2)*pos(2)+pos(3)*pos(3))
!     print *,'r: ',r
     v = dsqrt(vel(1)*vel(1)+vel(2)*vel(2)+vel(3)*vel(3))

!    ASSIGN VALUES to relativity constants, beta and gamma
     beta = 1.d0
     gam  = 1.d0

!    CALCULATE dot product of satellite velocity with itself (vdotv)
     vdotv = vel(1)*vel(1)+vel(2)*vel(2)+vel(3)*vel(3)
!    CALCULATE dot product of satellite position with velocity (rdotv)
     rdotv = pos(1)*vel(1)+pos(2)*vel(2)+pos(3)*vel(3)

!    CALCULATE Angular momentum of the earth per unit mass, J
     J = angvel*Ie/Me 
!     print*,'generalrel: J = ',J,dsqrt(j(1)*j(1)+j(2)*j(2)+j(3)*j(3))  ! |J| has a value of 981034718.66715....

!    CONVERT and assign Sun's position and velocity (m/day to m/s) vectors 
!     print *,'sunpos ',sunpos
     do i=1,3
       posS(i)=sunpos(i)
! PT140130: EK's code incremented sunpos by i+1 to get the velocities. Should it be i+3 ... ?
       velS(i)=sunpos(i+3)/(24.d0*60.d0*60.d0)
     enddo

!    CALCULATE the distance between the sun and the earth
     RS = dsqrt(posS(1)*posS(1)+posS(2)*posS(2)+posS(3)*posS(3))
!    CALCULATE the dot product of position and angular momentum J
     rdotJ = pos(1)*J(1)+pos(2)*J(2)+pos(3)*J(3)
!    CALCULATE the cross product of satellite position and velocity
     rcrossv = crossproduct(pos,vel)
!    CALCULATE the cross product of satellite velocity and the Earth's angular momentum
     vcrossJ = crossproduct(vel, J)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    CALCULATE the cross product of the Sun's velocity and the scaled Sun's position
     VScrossRS = crossproduct(velS, (-gm(3)*posS/(c*c*RS*RS*RS)))
!    CALCULATE the cross product of the above VScrossRS value and the satellite velocity
     bigcross = crossproduct(VScrossRS, vel)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE  the general relativity contribution to acceleration 
! DEBUG: according to Charlie Lineweaver, beta and gamma are related and, in special relativity, don't have the values of 1
!        beta = v/c   and gamma = 1/sqrt(1-beta^2)
!     beta = v/c
!     gam = 1.d0/dsqrt(1.d0-beta*beta)

     do i = 1, 3
       grelacc(i) = gm(1)/(c*c*r*r*r)*( (2.d0*(beta+gam)*gm(1)/r - gam*vdotv)*pos(i) + 2.d0*(1+gam)*(rdotv)*vel(i)) &    ! Schwartzchild term    (1e-8  m/s^2)
                     + (1.d0+gam) * gm(1)/(c*c*r*r*r) * ((3.d0/(r*r)*rcrossv(i)*rdotJ)+vcrossJ(i)) &                     ! Lens-Thirring effect  (1e-10 m/s^2)
                     + ((1.d0+2.d0*gam)*(bigcross(i)))                                                                   ! de Sitter effect      (1e-16 m/s^2)

! PT140205: add the acceleration perturbations caused by the Earth oblation (the J2 effect). Use
!           equations 4.2.23, 4.2.24 and 4.2.25 from Soffel (1989) "Relativity in Astrometry, Celestial
!           mechanics and Geodesy" 
       if( j2flag ) then
       if(i < 3)  then
         a1(i) = 2.d0*(beta+gam)*gm(1)/(c*c*r)*J2*(Ae/r)**2 * ( (pos(i)/r)*(2.d0-(9.d0*pos(3)**2/(r*r)) )*gm(1)/(r*r) )
       else
         a1(i) = 2.d0*(beta+gam)*gm(1)/(c*c*r)*J2*(Ae/r)**2 * ( (pos(i)/r)*(5.d0-(9.d0*pos(3)**2/(r*r)) )*gm(1)/(r*r) )
       endif
       a2    = 3.d0*(gam+1)   *gm(1)/(c*c*r)*J2*(Ae/r)**2 * ( (1.d0-5.d0*(pos(3)/r)**2) * (pos(1)*vel(1)+pos(2)*vel(2))/r  &
                                                             +(3.d0 - 5.d0*(pos(3)/r)**2) * pos(3)*vel(3)/r ) * v/r  
       if(i < 3)then
         a3(i) = -3.d0/2.d0*gam *gm(1)/(c*c)*J2*(Ae/r)**2 * ( (pos(i)/r)*(1.d0-(5.d0*pos(3)**2/(r*r)) )  *v*v/(r*r) )
       else
         a3(i) = -3.d0/2.d0*gam *gm(1)/(c*c)*J2*(Ae/r)**2 * ( (pos(i)/r)*(3.d0-(5.d0*pos(3)**2/(r*r)) )  *v*v/(r*r) )
       endif
!       grelacc(i) = grelacc(i) + a1(i) + a2 + a3(i)
! DEBUG
!       print*,'generalrel: i,a1(i)',i,a1(i),a2,a3(i) 
       endif
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! PT140130: try changing the sign of the general relativity acceleration. This actually increated the size of the
!           prefit residuals: magnitude of last peak in X increased from 0.44 to 2.26m. Using no genrel correction increased 
!           the prefit residuals from peak of 0.44 m to 1.35 m.
! PT140206: try scaling it
!     grelacc = -1.d0 * grelacc
 
     return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       contains
!      function to calculate the cross product of two vectors
       function crossproduct(a,b)
       implicit none
       real*8, dimension(3) :: a, b
       real*8, dimension(3) :: crossproduct

!      a × b = (a2b3 − a3b2) i + (a3b1 − a1b3) j + (a1b2 − a2b1) k
       crossproduct(1) = a(2)*b(3) - a(3)*b(2)
       crossproduct(2) = a(3)*b(1) - a(1)*b(3)
       crossproduct(3) = a(1)*b(2) - a(2)*b(1)

       end function crossproduct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     end

