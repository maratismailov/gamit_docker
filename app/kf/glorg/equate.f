 
CTITLE equate
 
      subroutine equate( iout, cov_parm, sol_parm)

      implicit none 
 
*     This routine will loop over the equate parameters and force
*     the listed parameters to be made equal adjustment.
 
 
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
*   ipivot(max_each_equate)     - Pivot elements
 
      integer*4 i,j, ipivot(max_each_equate)
 
*     The following arrays are used internally by equate paramters.
*  dchi  - change in chi**2 for each condinition
*  tot_dchi - total change on chi**2
*  num_dchi - Number of condictions in chi**2 space 
 
      real*8 cov_dcol(max_glb_parn, max_each_equate),
     .    sol_dcol(max_each_equate), dchi, tot_dchi, 
     .    avat(max_each_equate,max_each_equate),
     .    equ_gn(max_glb_parn, max_each_equate),
     .    scale(max_each_equate)

      integer*4 num_dchi
 
***** Loop over all of the equate paramss
      if( num_equates.gt.0 ) then
          write(iout,100) num_equates
 100      format(/' Equating parameters: ',i4,' Equates to be applied',
     .           /,'  #    Sigma    dChi**2   List of parameters')
      end if

      tot_dchi = 0.d0
      do i = 1, num_equates

*         Report status
          if( (i-1-int((i-1)/10)*10).eq.0 ) then
             write(*,110) i,num_equates
 110         format(i4,' of ',i4,' Equates completed')
          endif

*         Now force the equate
	  if( num_each_equate(i).gt.1 ) then
* MOD TAH 150223: Make sure none of the equate entries is zero:
              do j = 1, num_each_equate(i)
                if( param_equates(j,i).le.0 .or. 
     .              param_equates(j,i).gt. num_glb_parn) then
                    write(*,115) i, num_each_equate(i), 
     .                  param_equates(1:num_each_equate(i),i)   
 115                format('ERROR in equate parameter number: Equate ',
     .                      I5,' Number ',I5,' List ',
     .                      100I5,/,'Not completing')
                    exit
                 end if
              end do
              call equate_parms( cov_parm, sol_parm, num_glb_parn,
     .            cov_dcol, sol_dcol, param_equates(1,i), 
     .            num_each_equate(i), dchi, avat, equ_gn, ipivot, scale,
     .            eq_var(i))
              if( dchi.gt.0 ) then

                 tot_dchi = tot_dchi + dchi*(num_each_equate(i)-1)
                 num_dchi = num_dchi + num_each_equate(i) - 1
                 write(iout,120) i, sqrt(eq_var(i)), dchi, 
     .                         (param_names(param_equates(j,i)),
     .                         j=1,num_each_equate(i))
 120             format(I4,': ',F8.5,1x,F8.2,1x,8A12,:40(/,24x,8A12))
              else
                 write(iout,130) i, 
     .                         (param_names(param_equates(j,i)),
     .                         j=1,num_each_equate(i))
 130             format(I4,': ','Not completed    ',
     .                  1x,8A12,:40(/,24x,8A12))
              end if
	  end if
 
      end do

      write(iout, 140) tot_dchi/num_dchi, num_equates, num_dchi
 140  format(' Total change in Chi**2/f is ',F8.2,' for ',i4,
     .       ' equates, and ',i4,' conditions')
      write(iout, 160) (sum_chi+tot_dchi)/(sum_chi_num+num_dchi),
     .                  nint(sum_chi_num+num_dchi)
 160  format(' Solution chi**2/f now ',f8.2,' with ',i7,
     .       ' degrees of freedom')

*     Now update the values for glout part 
      sum_chi = sum_chi + tot_dchi
      sum_chi_num = sum_chi_num+num_dchi
 
****  Thats all
      return
      end
 
