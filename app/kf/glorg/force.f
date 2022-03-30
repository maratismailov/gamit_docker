 
CTITLE force
 
      subroutine force( iout, cov_parm, sol_parm)

      implicit none 
 
*     This routine will loop over the force parameters and force
*     the listed parameters to be made equal adjustment.
 
* max_each_force - Max mnumber of forces at one time (set to one)
      integer*4 max_each_force
      parameter ( max_each_force = 1 ) 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glorg_common.h'
      include '../includes/orient.h'
 
* PASSED VARIABLES
 
*   cov_parm(num_glb_parn, num_glb_parn)    - Covariance matrix
*   sol_parm( num_glb_parn)     - Solution vector.
*   iout  - Output unit number
 
      real*8 cov_parm(num_glb_parn, num_glb_parn),
     .    sol_parm( num_glb_parn)

      integer*4 iout
 
* LOCAL
 
*   i       - Loop counter
*   ipivot(max_each_force)     - Pivot elements
 
      integer*4 i, ipivot(max_each_force)
 
*     The following arrays are used internally by force paramters.
*  dchi  - change in chi**2 for each condinition
*  tot_dchi - total change on chi**2
*  num_dchi - Number of condictions in chi**2 space 
*  val_prior - Value of parameter prior to forcing
*  sig_prior - Sigma prior to forcing
 
      real*8 cov_dcol(max_glb_parn, max_each_force),
     .    sol_dcol(max_each_force), dchi, tot_dchi, 
     .    avat(max_each_force,max_each_equate),
     .    equ_gn(max_glb_parn, max_each_force),
     .    scale(max_each_force), val_prior, sig_prior

      integer*4 num_dchi, np
 
***** Loop over all of the force paramss
      if( num_force.gt.0 ) then
          write(iout,100) num_force
 100      format(/' Forcing parameters : ',i4,' forces to be applied',/,
     .            '  #    dChi**2   Parameter        Forced to',10x,
     .            ' Force Sigma ',7x,
     .            ' From value ',10x,' +-')
      end if

      tot_dchi = 0.d0
      do i = 1, num_force

 
*         Now force the force.  Check to make sure that the sigma
*         if not already zero.
          if( cov_parm(param_force(i),param_force(i)).lt.1.d-16 ) then
              write(iout,110)  i, param_force(i)
 110          format(' ** WARNING ** Force number ',i4,
     .               ' (parameter ',i5,') has zero variance')
          else
              np = param_force(i)
              val_prior = sol_parm(np)
              sig_prior = sqrt(cov_parm(np,np))

              call force_parms( cov_parm, sol_parm, num_glb_parn,
     .            cov_dcol, sol_dcol, param_force(i), 1,  val_force(i),
     .            var_force(i), dchi, avat, equ_gn, ipivot, scale)

              tot_dchi = tot_dchi + dchi
              num_dchi = num_dchi + 1

              write(iout,120) i, dchi, param_names(param_force(i)), 
     .                        val_force(i), sqrt(var_force(i)),
     .                        val_prior, sig_prior
 120          format(I4,': ',F8.2,1x,A12,1x,2f20.6,1x,f20.6,1x,f15.6)
         end if
 
      end do

      write(iout, 140) tot_dchi/num_dchi, num_force, num_dchi
 140  format(' Total change in Chi**2/f is ',F8.2,' for ',i4,
     .       ' forces, and ',i4,' conditions')
      write(iout, 160) (sum_chi+tot_dchi)/(sum_chi_num+num_dchi),
     .                 nint(sum_chi_num+num_dchi)
 160  format(' Solution chi**2/f now ',f8.2,' with ',i7,
     .       ' degrees of freedom')

*     Now update the values for glout part
      sum_chi = sum_chi + tot_dchi
      sum_chi_num = sum_chi_num+num_dchi
 
****  Thats all
      return
      end
 
