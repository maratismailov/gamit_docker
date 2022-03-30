CTITLE ROTATE_COV

      subroutine rotate_cov(R1, cov, dim, R2, dim2, ns )

*     Routine to rotate the covariance matrix.  Cross element between
*     site 1 with rotation R1 and site site 2 with rotation R2

* PASSED VARIABLES

*   dim     - The row dimension of covaraince matrix
*   dim2    - dimension of R2 (assumed square)
*   ns      - number of parmaters to rotate

      integer*4 dim, dim2, ns

*   R1(3,3), R2(dim2,dim2)    - The rotation matrices for each site
*   cov(dim,ns)      - The cross covariance matrix (upper
*                   - right hand cornore pointed to)


      real*8 R1(3,3), R2(dim2, dim2), cov(dim,ns)

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
          do j = 1,ns
              cov_R2T(i,j) = 0.d0
              do k = 1,ns
                  cov_R2T(i,j) = cov_R2T(i,j) + cov(i,k)*R2(j,k)
              end do
          end do
      end do

*     Now complete with R1 * cov_R2T

      do i = 1,3
          do j = 1,ns
              cov(i,j) = 0.d0
              do k = 1,3
                  cov(i,j) = cov(i,j) + R1(i,k)*cov_R2T(k,j)
              end do
          end do
      end do

***** Thats all
      return
      end


