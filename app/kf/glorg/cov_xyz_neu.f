CTITLE cov_xyz_neu
 
      subroutine cov_xyz_neu( direct, cov_parm, sol_parm )

      implicit none 
 
*     This routine will transform the complete variance-covariance
*     from with XYZ to NEU or from NEU back to XYZ
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED VARIABLES
 
* direct  - Direction of transformation.  If => 0 then from XYZ to
*           NEU; if <0 then NEU to XYZ
 
 
      integer*4 direct
 
*   cov_parm(num_glb_parn, num_glb_parn)    - Covariance matrix
*   sol_parm(num_glb_parn)          - Solution vector
 
 
      real*8 cov_parm(num_glb_parn, num_glb_parn),
     .    sol_parm(num_glb_parn)
 
*   i,j,k,l         - Loop counters
*   nq, np      - Parmeter number for rotating solution and
*               - covariance matrix
*   non_site    - First non-site position or rate parameter.  (Assummes
*                 site positions and rates are the first elements in the
*                 covariance matrix.
 
 
      integer*4 i,j,k,l, nq, np, non_site
 
*   rot_mat(3,3,max_glb_sites) - Rotation matrices (either from XYZ to
*                 NEU or NEU to XYZ depending on direct.
*   loc_coords(3) - Local coordinates of the site (not really used)
*   xyz_coords(3) - XYZ coordinates of site after accounting for
*                   velocity 
 
 
      real*8 rot_mat(3,3, max_glb_sites), loc_coords(3), dt,
     .       xyz_coords(3)

      non_site = 0
 
****  Now rotate the station coordinates from local frame to cartesian
*     XYZ.  Also find the first parameter after the site corrdinates.
*     These are not changed by the rotations
 
      do i = 1, gnum_sites
 
*         The transformation which rotates site coordinates from
*         XYZ to NEU
          dt = (gepoch_out-site_epoch(i))/365.25d0
          do j = 1,3
             xyz_coords(j) = apr_val_site(j,1,i) +
     .                       apr_val_site(j,2,i)*dt
          end do
          call xyz_to_neu( rot_mat(1,1,i), xyz_coords,
     .                     loc_coords )
 
*         If we are neu to xyz then transpose the matrix.
          if( direct.lt.0 ) then
              do j = 1,2
                 call dvswp(rot_mat(j,j+1,i),3,
     .                      rot_mat(j+1,j,i),1, 3-j)
              end do
          end if
          if( parn_site(3,1,i).ne.0 ) non_site = parn_site(3,1,i)+1
          if( parn_site(3,2,i).ne.0 ) non_site = parn_site(3,2,i)+1
      end do
 
***** Now rotate the solution vector and the covariance matrix for both
*     the value and the rate (This routine assumes that all coordintes
*     are estimated.
      do i = 1, gnum_sites
        do j = 1,2
 
*         Solution vector
          nq = parn_site(1,j,i)
          if( nq.ne.0 ) then
              call rotate_glsol( rot_mat(1,1,i), sol_parm(nq))
 
*****         Now rotate the covariance matrx for values and rates
              do k = i, gnum_sites
                  do l = 1, 2
                     np = parn_site(1,l,k)
                     if( np.gt.0 ) then
                         call rotate_glcov(rot_mat(1,1,i),
     .                       cov_parm(nq,np), num_glb_parn,
     .                       rot_mat(1,1,k) )
                     end if
                   end do
              end do
 
****          now do the remaining parameters in the solutions
              call rotate_runit( rot_mat(1,1,i), cov_parm(nq,non_site),
     .                           num_glb_parn, num_glb_parn-non_site+1)
           end if
         end do
      end do
 
 
****  Now fill in the symetric part of the matrix (This code actually
*     copies the diagonal as well but this should be no problem)
      do i = 1, num_glb_parn
          call dwmov(cov_parm(1,i),1, cov_parm(i,1), num_glb_parn,i)
      end do

c     write(*,100) direct, non_site
c100  format(' COV_XYZ_NEU: Direction ',i3,' non-site ', i4)
c     do i = 1, gnum_sites
c        do j = 1,3
c           np = parn_site(j,1,i)
c           if( np.gt.0 ) then
c               write(*,110) i,j, np, sol_parm(np), 
c    .                       sqrt(cov_parm(np,np))
c110            format(3i4, 2F20.3)
c           end if
c        end do
c     end do
 
 
****  Thats all.  We are now ready to generate the solution
      return
      end
 
 
CTITLE ROTATE_GLCOV
 
      subroutine rotate_glcov(R1, cov, dim, R2 )

      implicit none 
 
*     Routine to rotate the covariance matrix.  Cross element between
*     site 1 with rotation R1 and site site 2 with rotation R2
 
* PASSED VARIABLES
 
*   dim     - The row dimension of covaraince matrix
 
 
      integer*4 dim
 
*   R1(3,3), R2(3,3)    - The rotation matrices for each site
*   cov(dim,3)      - The cross covariance matrix (upper
*                   - right hand cornore pointed to)
 
 
 
      real*8 R1(3,3), R2(3,3), cov(dim,3)
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop couners
 
 
      integer*4 i,j,k
 
*   cov_R2T(3,3)    - Product of covariance matrix by R2 transposed
 
 
      real*8 cov_R2T(3,3)
 
****  Start.  First do:
*                      T
*     COV_R2T = Cov* R2
*
      do i = 1,3
          do j = 1,3
              cov_R2T(i,j) = 0.d0
              do k = 1,3
                  cov_R2T(i,j) = cov_R2T(i,j) + cov(i,k)*R2(j,k)
              end do
          end do
      end do
 
*     Now complete with R1 * cov_R2T
 
      do i = 1,3
          do j = 1,3
              cov(i,j) = 0.d0
              do k = 1,3
                  cov(i,j) = cov(i,j) + R1(i,k)*cov_R2T(k,j)
              end do
          end do
      end do
 
***** Thats all
      return
      end
 
CTITLE ROTATE_GLSOL
 
      subroutine rotate_glsol( R1, sol )

      implicit none 
 
*     Routine to apply the R1 solutiuon matrix to solution vector
*
* PASSED VARIABLES
 
*   R1(3,3)     - The rotation matrice for this site
*   sol(3)      - The solution vector (adjustment)
 
 
      real*8 R1(3,3), sol(3)
 
* LOCAL VARIABLES
 
*   i,j     - Loop couners
 
 
      integer*4 i,j
 
*   R1_sol(3)   - Product of R1*sol
 
 
 
      real*8 R1_sol(3)
 
****  Loop over each element
      do i = 1,3
          R1_sol(i) = 0.d0
          do j = 1,3
              R1_sol(i) = R1_sol(i) + R1(i,j)*sol(j)
          end do
      end do
 
****  Now return the new values
      do i = 1,3
          sol(i) = R1_sol(i)
      end do
 
***** Thats all
      return
      end
 
CTITLE ROTATE_RUNIT
 
      subroutine rotate_runit(R1, cov, dim, ic )

      implicit none 
 
*     Routine to compute the covariance between the rotated parameters
*     and those which remained unchanged in the transformation
 
* PASSED VARIABLES
 
*   dim     - The row dimension of covaraince matrix
*   ic      - The remaining number of parameters in the solution
 
 
      integer*4 dim, ic
 
*   R1(3,3)         - The rotation matrices for each site
*   cov(dim,ic)      - The cross covariance matrix (upper
*                   - right hand cornore pointed to)
 
 
 
      real*8 R1(3,3),cov(dim,ic)
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop couners
 
 
      integer*4 i,j,k
 
*   dcov(3)    - One 3 element column of the covariance matrix.  Must be
*                computed first since the values are needed during the
*                computation of each column
 
 
      real*8 dcov(3)
 
*     Multiple R1 by the columns of the covariance matrix
 
 
      do i = 1,ic
          do j = 1,3
             dcov(j) = 0.d0
             do k = 1, 3
                dcov(j) = dcov(j) + R1(j,k)*cov(k,i)
             end do
          end do
 
*         now copy the row since we are finished with it
          do j = 1,3
             cov(j,i) = dcov(j)
          end do
      end do
 
 
***** Thats all
      return
      end
 
