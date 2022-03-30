
CTITLE write_glb_pos
 
      subroutine write_glb_pos( iout, options, use_sites,
     .                          cov_parm, sol_parm, descr )

      implicit none  
 
*     Routine to write out summary of the position adjs.  of the parameters from the
*     global solution.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
c     include '../includes/glorg_common.h'
*                                           ! To get nutation names
      include '../includes/globk_cmds.h'
      include '../includes/const_param.h'
 
*   icnt            - Counts number of components of site position
*                   - estimated
*   iel, jel, kel   - Lookup positions in the solution vector
*                   - and covariance matrix
*   i,j,k           - loop counters
*   iout            - Output lu
*   options         - Options for the output
*   use_sites       - Indicates which sites are used in the origin 
 
      integer*4 icnt, iel, i,j, iout, options, use_sites(*)
      integer*4 trimlen
 
*   covar(36)       - Six by six matrix needed for covariance
*                   - calculations
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix from
*                   - the solution
*   dt              - Time difference between site/source epoch
*                   - and current experiment
*   latitude, longitude - Lat and Long in degs
 
*   loc_coord(3)    - Local coordinates of the site (latitude,
*                   - longitude and radius)
*   NEU_covar(3,3)  - Covariance matrix for the NEU coordinates
*   pos_xyz_adj(3)  - adjustment to XYZ position
*   pos_NEU_adj(3)  - adjustment to NEU position
*   pos_xyz_fin(3)  - Final XYZ position
*   pos_NEU_fin(3)  - Final NEU position
*   pos_radec(2)    - Final positions for RA and Dec
*   rot_matrix(3,3) - Rotation matrix between XYZ and NEU
*   scr_real(10)    - Scratch area
*   sol_parm(num_glb_parn)  - Solution vector from
*                   - the solution
*   temp_covar(36)  - Temporay storage for COMP_VAR
 
      real*8 covar(36), cov_parm(num_glb_parn,num_glb_parn), dt,
     .    latitude, longitude, loc_coord(3), NEU_covar(3,3),
     .    pos_xyz_adj(3), pos_NEU_adj(3), pos_xyz_fin(3),
     .    pos_NEU_fin(3), rot_matrix(3,3), sol_parm(num_glb_parn),
     .    temp_covar(36)

*  sum_res(3)  - Sum of NE and U weighted residuals
*  sum_var(3)  - Sum of NE and U weighted residuals squares
*  sum_wgh(3)  - Sum of NE and U weights
*  wmean_neu(3) - Weighted Mean NE and U
*  wrms_neu(3) - Weighred RMS of NE and U
*  nrms_neu(3) - Normalized RMS of NE and U
*  num_var     - Number of values in sums
      real*8 sum_res(3), sum_var(3), sum_wgh(3), 
     .       wmean_neu(3), wrms_neu(3), nrms_neu(3)
      integer*4 num_var
 
*   descr -- Descriptor  string for run (list file name in glorg)

      character*(*) descr

*   orgchar     - * if site used in origin
 
      character*1 orgchar

      character*12 hsver
 
 
      logical kbit
 
 
 
      common covar, temp_covar
 
***** Output the header line
      write(iout,100) hsver(globk_version)
  100 format( /,' SUMMARY POSITION ESTIMATES FROM GLOBK Ver ',a,/,
     .          '    Long.       Lat.        dE adj.   dN adj.',
     .          '   dE +-     dN +-   RHO        dH adj.   dH +-  SITE',
     .        /,'    (deg)      (deg)         (mm)      (mm)  ',
     .          '    (mm)      (mm)               (mm)      (mm)' )

****  Clear statistics for computing RMS to reference sites
      do i = 1,3
        sum_res(i) = 0.d0
        sum_var(i) = 0.d0
        sum_wgh(i) = 0.d0
      end do
      num_var = 0
 
*     Site positions
      do i = 1, gnum_sites
*                             ! Indicates if any componets estimated
          icnt = 0
          dt   = (gepoch_out-site_epoch(i))/365.25d0
 
*         Do values first
*                             ! Loop over XYZ
          do j = 1,3
*                                              ! estimated
              if( parn_site(j,1,i).ne.0 ) then
                  iel  = parn_site(j,1,i)
                  icnt = icnt + 1
 
*                 Get final position of site
                  pos_xyz_adj(j) = sol_parm(iel)
                  pos_xyz_fin(j) = apr_val_site(j,1,i)    +
     .                             apr_val_site(j,2,i)*dt +
     .                             cont_nonsec(j,i)       +
     .                             pos_xyz_adj(j)
              else
                  pos_xyz_fin(j) = apr_val_site(j,1,i)    +
     .                             apr_val_site(j,2,i)*dt
*                             ! estimated
              end if
*                             ! Looping over XYZ
          end do
 
*****     Now see if we should do NEU components
*                                 ! Convert to NEU coordinates
          if( icnt.gt.0 ) then
              call rotate_geod(pos_xyz_adj, pos_neu_adj, 'XYZ','NEU',
     .                         pos_xyz_fin, loc_coord, rot_matrix)
 
*             Convert the local coordinates
              call loc_to_geod( loc_coord, pos_neu_fin )
 
*             Now compute the sigmas. Firstly get the elements of
*             covariance matrix we need
              call mov_cov(parn_site(1,1,i), parn_site(1,1,i), 3,
     .                     covar,3, cov_parm, num_glb_parn)
              call var_comp(rot_matrix, covar, NEU_covar, temp_covar,
     .                      3,3,1)
 
*                         ! Any local coordinates to be output
          end if
 
          latitude = 90.d0-loc_coord(1)*180/pi
          longitude = loc_coord(2)*180/pi
 
          if( kbit(use_sites,i) ) then
              orgchar = '*'
          else
              orgchar = ' '
          end if


*         Only write out value if position estimated.
          if( icnt.gt.0 ) then
              write(iout,200)   longitude, latitude,
     .                pos_neu_adj(2)*1.d3, pos_neu_adj(1)*1.d3,
     .                sqrt(abs(neu_covar(2,2)))*1.d3, 
     .                sqrt(abs(neu_covar(1,1)))*1.d3,
     .                neu_covar(1,2)/sqrt(abs(neu_covar(1,1)*
     .                                        neu_covar(2,2))),
     .                pos_neu_adj(3)*1.d3, 
     .                sqrt(abs(neu_covar(3,3)))*1.d3,
     .                gsite_names(i), orgchar
 200          format(2(1x,f10.5),1x,4(2x,f8.2),1x,f6.3,2x,
     .               2(2x,f8.2), 1x,a8,a1)

****          Accumulate statisttics on the velocity residuals for
*             the reference sites
              if( kbit(use_sites,i) ) then
                  num_var = num_var + 1
                  do j = 1,3
                     sum_res(j) = sum_res(j) + pos_neu_adj(j)/
     .                                         neu_covar(j,j)
                     sum_wgh(j) = sum_wgh(j) + 1.d0/neu_covar(j,j)
                     sum_var(j) = sum_var(j) + pos_neu_adj(j)**2/
     .                                         neu_covar(j,j)
                  end do
              end if

          end if
      end do

****  Compute and write out final statistics
      if( num_var.gt.0 ) then
          do j = 1, 3
             wmean_neu(j) = sum_res(j)/sum_wgh(j)
* MOD TAH 160510: Make WRMS about the MEAN 
             wrms_neu(j) = sqrt(sum_var(j)/sum_wgh(j)-wmean_neu(j)**2)
             nrms_neu(j) = sqrt((sum_var(j)-wmean_neu(j)**2*sum_wgh(j))
     .                                                  /num_var)
          end do
          write(iout,310) num_var, wrms_neu(2)*1000, 
     .                    wrms_neu(1)*1000, wrms_neu(3)*1000,
     .                    nrms_neu(2),nrms_neu(1),nrms_neu(3), 
     .                    descr(1:trimlen(descr))
          write(iout,320) num_var, wmean_neu(2)*1000,
     .                    wrms_neu(2)*1000/sqrt(float(num_var)),
     .                    wmean_neu(1)*1000, 
     .                    wrms_neu(1)*1000/sqrt(float(num_var)),
     .                    wmean_neu(3)*1000,
     .                    wrms_neu(3)*1000/sqrt(float(num_var)),
     .                    descr(1:trimlen(descr))


 310      format('POS STATISTICS: For ',i4,' RefSites WRMS ENU ',
     .            3(F6.2,1x), ' mm    NRMS ENU ',3(F6.2,1x),a)
 320      format('POS MEANS: For ',i4,' RefSites: ',
     .           'East ',F6.2,' +- ',F6.2,
     .           ' North ',F6.2,' +- ',F6.2,
     .           ' Up ',F6.2,' +- ',F6.2,' mm ',a)
      endif 
 
      return
      end
 
