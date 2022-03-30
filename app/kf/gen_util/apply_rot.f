CTITLE APPLY_ROT
 
      subroutine apply_rot( parn, dim, sol_parm, pmu_parts,
     .           gnum_sites, pmu_changes )

      implicit none 
 
 
 
*     Routine to apply a rotational change
 
*   i,j,k       - Loop counter
*   iel         - position in sol_parm
*   dim         - Dim of parn
*   gnum_sites  - Number of sites
*   parn(3,dim,1)   - Parmeter numbers
 
      integer*4 i,j, iel, dim, gnum_sites, parn(3,dim,1)
 
*   corr        - Correction to rate
*   sol_parm(1) - Parameter estimates
*   pmu_parts(3,3,1)  - PMU partials
*   pmu_changes(3)  - Changes to pmu values (mas)
 
      real*8 corr, sol_parm(1), pmu_parts(3,3,1), pmu_changes(3)
 
 
*     Loop over sites
      do i = 1,gnum_sites
*                     ! Components
          do j = 1,3
              if( parn(j,1,i).ne.0 ) then
                  iel = parn(j,1,i)
                  corr        = pmu_parts(1,j,i)*pmu_changes(1)
     .                        + pmu_parts(2,j,i)*pmu_changes(2)
     .                        + pmu_parts(3,j,i)*pmu_changes(3)
                  sol_parm(iel) = sol_parm(iel) - corr
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
 
