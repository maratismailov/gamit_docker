      subroutine rep_S1S2Com( atides, slat, slon )

      implicit none

*     Routine to remove the Center of Mass S1/S2 signals
*     from the site S1/S2 load and replace with model from
*     fitting to 3-hr loading.u-strasbg.fr/ITRF/CM timeseries.
*     (differenced from CF version, atmib versions).

      real*8 pi
      parameter ( pi            = 3.1415926535897932D0 )

* PASSED
      real*4 atides(2,6) ! cos sin for U N E tide displacements
      real*4 slat, slon  ! Latitude and long of site (deg) needed
                         ! to rotate from XYZ to UNE

* LOCAL 
      real*4 dIERS_CoM(3,4)   ! Difference between IERS CoM model 
                         ! for XYZ dor S1 Cos/Sin and S2 Cos/Sin  
                         ! and Strasborg fit model (mm) (Based on
                         ! ~/TAH_docs/ATM_load/Chur.m)

      real*4 rot(3,3)    ! Rotation from XYZ to UNE (ie., dUNE = rot*dXYZ)
                         ! (NOTE: Switchs of last row to first from our
                         ! usual NEU rotation.
      real*4 dUNE(3)     ! Vector dUNE for dXYZ CoM component.

      integer*4 j,k,m,n  ! Loop counters for S1/S2 Cos/Sin and rot multiply
      integer*4 col      ! column number 


!     Note in entering data we need to enter by column 
!                        X       Y      Z
      data dIERS_CoM /  0.01,  0.17,   0.05,    ! Cos S1
     .                  0.22, -0.00,  -0.04,    ! Sin S1
     .                  0.13,  0.21,   0.11,    ! Cos S2
     .                  0.23,  0.05,   0.05 /   ! Sin S2 

****  Based on the slat and slon of the site, get the rotation 
*     XYZ to UNE
!      print *,'slat, slon ',slat, slon, pi
*     Up row
      rot(1,1) =  cos(slat*pi/180)*cos(slon*pi/180)
      rot(1,2) =  cos(slat*pi/180)*sin(slon*pi/180)
      rot(1,3) =  sin(slat*pi/180)
*     North row
      rot(2,1) = -sin(slat*pi/180)*cos(slon*pi/180)
      rot(2,2) = -sin(slat*pi/180)*sin(slon*pi/180)
      rot(2,3) = +cos(slat*pi/180)
*     East row
      rot(3,1) = -sin(slon*pi/180)
      rot(3,2) = +cos(slon*pi/180)
      rot(3,3) =  0.0d0

!     print *,'ROT ',rot(1,:),' R2 ',rot(2,:),' R3 ',rot(3,:)
!     print *,'dIERS_CoM ',dIERS_CoM
!     print *,'atides BS1 ',atides(1,:)
!     print *,'atides BS2 ',atides(2,:) 

****  Now multiply rot*dIERS_CoM to get NEU changes and apply
*     to UNE values.  Loop over S1/S2
      do j = 1, 2       ! S1/S2: First index is atides
         do k = 1, 2    ! Cos/Sin Compute (1 3 5 and 2 4 6)
             dUNE = 0.0
             col = (j-1)*2+k   ! Column in dIERS_CoM
             do m = 1,3     ! dUNE
                do n = 1,3  ! dXYZ of dIIERS 
                   dUNE(m) = dUNE(m) + rot(m,n)*dIERS_CoM(n,col)
                end do
             end do
!            print *,'Comp J, K, dUNE ',j,k,col, dUNE
*            Now assign to correct columns in altides
             if( k.eq.1 ) then
                atides(j,1) = atides(j,1) + dUNE(1)
                atides(j,3) = atides(j,3) + dUNE(2)
                atides(j,5) = atides(j,5) + dUNE(3)
             else
                atides(j,2) = atides(j,2) + dUNE(1)
                atides(j,4) = atides(j,4) + dUNE(2)
                atides(j,6) = atides(j,6) + dUNE(3)
             end if
         end do
      end do
!     print *,'atides AS1 ',atides(1,:)
!     print *,'atides AS2 ',atides(2,:) 

****  Thats all 
      end

                
      
