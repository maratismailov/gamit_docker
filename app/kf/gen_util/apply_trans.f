CTITLE APPLY_TRANS
 
      subroutine apply_trans( parn, dim, sol_parm, gnum_sites,
     .                        rate_trans )

      implicit none 
 
 
*     Routine to apply a translation change
 
*   i,j,k       - Loop counter
*   iel         - position in sol_parm
*   dim         - Dim of parn
*   gnum_sites  - Number of sites
*   parn(3,dim,1)   - Parmeter numbers
 
      integer*4 i,j, iel, dim, gnum_sites, parn(3,dim,1)
 
*   corr        - Correction to rate
*   sol_parm(1) - Parameter estimates
*   rate_trans(3)   - Translation
 
      real*8  sol_parm(1), rate_trans(3)
 
*     Loop over sites
      do i = 1,gnum_sites
*                     ! Components
          do j = 1,3
              if( parn(j,1,i).ne.0 ) then
                  iel = parn(j,1,i)
                  sol_parm(iel) = sol_parm(iel) - rate_trans(j)
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
 
 
