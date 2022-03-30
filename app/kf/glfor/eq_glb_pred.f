CTITLE EQ_GLB_PRED
 
      subroutine eq_glb_pred(cov_parm, sol_parm)

      implicit none 
 
*     Routine to update the parameter covariance matrix for
*     any earthquake related processes.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED variables
 
*   cov_parm(num_glb_parn, num_glb_parn)    - Parameter
*       - covarinace matrix
*   sol_parm(num_glb_parn)      - Solution vector
 
      real*8 cov_parm(num_glb_parn, num_glb_parn),
     .    sol_parm(num_glb_parn)
 
* LOCAL VARIABLES
 
*   i       - Loop counter
 
      integer*4 i
 
*   dt_eq   - Time difference (days) between earthquake
*           - and current experiment mid-point
 
      real*8 dt_eq
 
****  First scan the earthquakes to see if we need to anything
 
      do i = 1, num_eq
 
*         Get time difference from earthquake
          dt_eq = gepoch_expt - eq_epoch(i)
 
*         See if we need to do anything with coseismic displacement
          if( sign(1.d0,deltat)*dt_eq.gt.0 .and.
     .        .not.eq_co_applied(i)     ) then
 
*             We just passed the Earthquake epoch so apply the
*             covariance increment to all sites effected by the
*             earthquake
              call eq_app_cov(cov_parm,sol_parm,eq_apr_coseismic(1,i),
     .                deltat, i,'COSEISMIC')
              eq_co_applied(i) = .true.
 
          else if( dt_eq.lt.eq_dur(2,i) .and.
     .             dt_eq.gt.0.d0 ) then
 
*             Post-seimic deformation process
              call eq_app_cov( cov_parm, sol_parm, eq_mar_post(1,i),
     .                deltat, i,'POSTSEISMIC')
 
          else if(  -dt_eq. lt. eq_dur(1,i) .and.
     .               dt_eq.lt.0.d0  ) then
 
*             Pre-seismic deformation process.
              call eq_app_cov( cov_parm, sol_parm, eq_mar_pre(1,i),
     .                deltat, i,'PRESEISMIC')
          end if
      end do
 
****  Thats all
      return
      end
 
CTITLE EQ_APP_COV
 
      subroutine eq_app_cov( cov_parm, sol_parm, eq_noise,
     .                    dt, ne, type )
 
      implicit none 
 
*     Routine to apply the process noise to the solution
*     covariance matrix due to an earthquake.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED VARIABLES
 
*   ne          - Earthquake number
 
      integer*4 ne
 
*   cov_parm(num_glb_parn, num_glb_parn)    - Covariance matrix
*               - of estimated parameters
*   sol_parm(num_glb_parn)                  - Solution vector.
*   eq_noise(6)     - Process noise either as a sigma (if dt is
*                   - zero) or a random walk process.
*   dt              - Time difference for random walk (is zero
*                   - then coseimic displacement) (years)
 
 
      real*8 cov_parm(num_glb_parn, num_glb_parn),
     .    sol_parm(num_glb_parn), eq_noise(6), dt

*   type        - Indicates whether this is cosesimic, postseismic
*                 or preseimic (used when re-name feature hase been
*                 used

      character*(*) type
 
* LOCAL variables
 
*   i,j,k,l     - Loop counters
*   np          - Parametet number for start of site.
*   ps          - Site number for original site when covariances
*                 are copied.
 
      integer*4 i,k,l, np, ps
 
*   dist        - distance of site i from earthquake
*   rot_mat(3,3)    - Rotation maxtrix from NEU to XYZ (finally)
*   loc_coord(3)    - local coordinates of site (col-lat, long and
*                   - height)
*   cov_neu(3,3)    - Covariance matrix in NEU
*   cov_xyz(3,3)    - Covariance matrix in XYZ
*   temp_cov(3,3)   - Temporary matrix used in converting NEU
*                   - covariance matrix to xyz
*   total_eq_noise  - Total eq_noise.  If zero then we are
*                     not applying any process noise so we 
*                     don't need to do anything
 
      real*8 dist, rot_mat(3,3), loc_coord(3), cov_neu(3,3),
     .    cov_xyz(3,3), temp_cov(3,3), total_eq_noise

*   update_cov     - Set true if we need to update covariance.

      logical update_cov
 
***** See if we need do anything
      total_eq_noise = 0.d0
      do i = 1,6 
         total_eq_noise = total_eq_noise + eq_noise(i)
      end do

      if( total_eq_noise.eq.0 ) RETURN 
  
***** Start loop over the stations to see which ones are close
*     enough to the Earthquake to be effected.
 
      do i = 1, gnum_sites
          call eval_dist(eq_pos(1,ne), apr_val_site(1,1,i), dist)
 
          if( dist.le.eq_rad(ne) ) then

*             See if we should update the covariance matrix.  If
*             sites have not been renamed then we update, other
*             our actions depend on type
              if( .not.eq_rename(ne) ) then
                  update_cov = .true.
              else

*                 Need to look in more detail.  If this is post-
*                 seismic, only update renamed sites.
                  update_cov = .false.
                  if( type(1:2).eq.'PO' .and. 
     .                gsite_names(i)(7:8).eq.eq_codes(ne) ) then
                      update_cov = .true.
*                 If pre-seismic only update old site name
                  else if( type(1:2).eq.'PR' .and.
     .                gsite_names(i)(7:8).ne.eq_codes(ne) ) then
                      update_cov = .true.
                  else if( type(1:2).eq.'CO' ) then

*                     Coseismic.  This is more complicated. If
*                     we are going forward in time, then only
*                     do new station.  If backwards in time only
*                     the old station.
                      ps = 0
                      if( dt.gt.0.0 .and. 
     .                    gsite_names(i)(7:8).eq.eq_codes(ne) ) then
                          update_cov = .true.

*                         Find the original site (so we can copy 
*                         covariance elements
                          call find_preq('FORWARD',i,ne, ps)
                      end if
                      if( dt.lt.0.0 .and.
     .                    gsite_names(i)(7:8).ne.eq_codes(ne) ) then
                          update_cov = .true.
                          call find_preq('BACKWARD',i,ne, ps)
                      end if  
                      call copy_coveq(ps,i, cov_parm, sol_parm)
                  end if
              end if

*             Update the process noise if we have to.
              if( update_cov ) then 
*                 We are in spatial range of the Earthquake.
*                 Get the transformation from XYZ to NEU
                  call xyz_to_neu( rot_mat, apr_val_site(1,1,i),
     .                             loc_coord)
* MOD TAH 990802: Changed to xyz_to_geod call
* MOD TAH 030161: Changed back to NEU call
C                 call xyz_to_geod( rot_mat, apr_val_site(1,1,i),
C   .                             loc_coord)
 
*                 Transpose the matrix
                  do k = 1,2
                      call dvswp(rot_mat(k,k+1),3,
     .                           rot_mat(k+1,k),1, 3-k)
                  end do
 
*                 Set up the diagonal NEU matrix
                  do k = 1,3
                     do l = 1,3
                         cov_neu(k,l) = 0.d0
                     end do
                     if( type(1:2).ne.'CO')  then
*                        Markov process.  Add the constant and
*                        distance dependednt terms.
                         cov_neu(k,k) = (eq_noise(k) +
     .                          eq_noise(k+3)*(eq_depth(ne)/dist)**4)*
     .                          abs(dt)
                     else
*                        Co-seismic so apriori sigma passsed
                         cov_neu(k,k) = eq_noise(k)**2 +
     .                       (eq_noise(k+3)*(eq_depth(ne)/dist)**2)**2
                     end if
                  end do
 
                  call var_comp(rot_mat, cov_neu, cov_xyz,
     .                          temp_cov, 3,3,1 )

*                 Now apply to covarince matrix:
                  np = parn_site(1,1,i)
                  if( np.gt.0 ) then
                      write(*,120) eq_codes(ne),gsite_names(i),dt*365, 
     .                          (sqrt(cov_neu(k,k))*1000.0,k=1,3)
 120                  format('Process noise from ',a2,' for site ',
     .                        a,' being added: dt (days) ', F6.1,
     .                          ' SigNEU (mm)', 3F8.1)
                      do k = 0,2
                         do l = 0,2
                           cov_parm(np+k,np+l) =
     .                          cov_parm(np+k,np+l) + cov_xyz(k+1,l+1)
                         end do
                      end do
                  end if
              end if
          end if
      end do
 
****  Thats all
      return
      end
 
 
 
 
