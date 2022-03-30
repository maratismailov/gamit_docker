CTITLE ETD_PARTIAL
 
      subroutine etd_partial

      implicit none 
 
 
c     routine to compute the extended earth tide partials
c     K1 frequency used for the diurnal tides and M2 is used
c     for the semi-diurnal tides
c MOD TAH 870616:  Converted units of partials to ps/m for
c     delays and (fs/s)/m for rates.
c     The ordering of the partials is
c     1 -- diurnal Radial in phase
c     2 -- diurnal Radial out of phase
c     3 -- diurnal Horizonal East in phase
c     4 -- diurnal Horizonal East out of phase
c     5 -- diurnal Horizonal north in phase
c     6 -- diurnal Horizonal north out of phase
c     7 -- semi-diurnal Radial in phase
c     8 -- semi-diurnal Radial out of phase
c     9 -- semi-diurnal Horizontal East in phase
c     10-- semi-diurnal Horizontal east out of phase
c     11-- semi-diurnal Horizontal north in phase
c     12-- semi-diurnal Horizontal north out of phase
c
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
      include '../includes/obs_apr.h'
      include '../includes/const_param.h'
 
 
 
*          i          - local site number
*          idex       - index for storing tidal partials
*          itide      - diurnal (1) and semi-diurnal (2) tides
*          j,s        - Counters for looping over XYZ and delay and
*                       and rate
 
      integer*4 i, idex, itide, j,s
 
*       colat         - co-latitude of site
*       cos_colat     - cosine of the co-latitude
*       phase         - phase of the normal mode
*       arg           - GST of observation (rad)

* MOD TAH 900217 Changed partials so that they will use the site position
*     partials, and replaced the plm with p20 calls explicitly.
 
      real*8 colat, phase, arg

*      x_to_l(3,3)  ! Rotations matrices to get from XYZ to
*               ! NEU (we need the transpose of thiis matrix)
*      du_in(3) ! Three compoenets of N,E and U motion in phase with 
*               ! potential
*      du_out(3)    ! Three components N,E, and U out of phase with
*               ! with potential
*      loc_coord(3)  ! Dummy local coordinates of site

      real*8 x_to_l(3,3), du_in(3), du_out(3), loc_coord(3)

*      cos_colat, sin_colat   ! Sine and cosine of co latitude
*      cos_phase, sin_phase   ! Sine and cosine of phase
*      dob_dl(3,2)   ! Parts of observables with respect u,e, and n
*                    ! for delay and rate

      real*8 cos_colat, sin_colat, cos_phase, sin_phase, dob_dl(3,2)
 
c.... Loop over two sites
      do i=1,2

*       First get the transformation from XYZ to NEU
        call XYZ_to_NEU(x_to_l, site_pos(1,site(i)), loc_coord)

*       Generate partials of observable with respect to u,e,and n
*       Parts of the deformation
*                     ! Loop delay and rate
        do s = 1,2 
*                     ! Loop over u, e and n
           do j = 1,3
              dob_dl(j,s) = 0.d0
           end do
           do j = 1,3
              dob_dl(1,s) = dob_dl(1,s) + x_to_l(3,j)*site_part(j,i,s)
              dob_dl(2,s) = dob_dl(2,s) + x_to_l(2,j)*site_part(j,i,s)
              dob_dl(3,s) = dob_dl(3,s) + x_to_l(1,j)*site_part(j,i,s)
           end do
        end do
 
c....   Loop over diurnal and semi diurnal tides
        do itide=1,2
 
c         calculate phase of mode (reference epoch JD 2451545.0)
          call gst_jd(epoch, arg)
          phase = (arg+loc_coord(2))*itide
          phase = mod(phase,2*pi)
 
c         calculate Legendre functions
          colat =   loc_coord(1)
          cos_colat = cos(colat)
          sin_colat = sin(colat)
          cos_phase = cos(phase)
          sin_phase = sin(phase)

*         Now compute the displacements.  We use different formulas
*         for diurnal and semidiurnal due cos and sine of time argument
*         and P21 versus P22
*                                    ! Diurnal tides
          if( itide.eq.1 ) then
*                                    ! North
             du_in(1) = -3.d0*cos(2.d0*colat)*sin_phase
*                                    ! East
             du_in(2) =  3.d0*sin_colat*cos_colat*cos_phase
*                                    ! Up
             du_in(3) =  3.d0*sin_colat*cos_colat*sin_phase

*            Now do out of phase for diurnal band (-cos phase)
*                                    ! North
             du_out(1) =  3.d0*cos(2.d0*colat)*cos_phase
*                                    ! East
             du_out(2) =  3.d0*sin_colat*cos_colat*sin_phase
*                                    ! Up
             du_out(3) = -3.d0*sin_colat*cos_colat*cos_phase

*         Now do the semidiurnal band
          else
             du_in(1)  = -6.d0*sin_colat*cos_colat*cos_phase
             du_in(2)  =  6.d0*sin(colat)**2*sin_phase
             du_in(3)  =  3.d0*sin(colat)**2*cos_phase
*            Out of phase
             du_out(1) = -6.d0*sin_colat*cos_colat*sin_phase
             du_out(2) =  6.d0*sin(colat)**2*cos_phase
             du_out(3) =  3.d0*sin(colat)**2*sin_phase
          end if

****      Now compute the partials with respect to delay and rate
*         This done by : 
*             d(delay)/d(amp) = d(delay)/d(X/Y/Z)*d(X/Y/Z)*d(amp)
*         Where the XYZ partials are stored in site_part and the
*         Rotation of neu to xyz is in x_to_l

          idex = (itide-1)*6

*                       ! Loop over delay and rate
          do s = 1,2

             etd_ext_part(1+idex,i,s) = dob_dl(1,s)*du_in(3)
             etd_ext_part(2+idex,i,s) = dob_dl(1,s)*du_out(3)
             etd_ext_part(3+idex,i,s) = dob_dl(2,s)*du_in(2)
             etd_ext_part(4+idex,i,s) = dob_dl(2,s)*du_out(2)
             etd_ext_part(5+idex,i,s) = dob_dl(3,s)*du_in(1)
             etd_ext_part(6+idex,i,s) = dob_dl(3,s)*du_out(1)
          end do
        end do
      end do
*
****  Thats all
      return
      end
 
 
 
