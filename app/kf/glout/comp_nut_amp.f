      subroutine comp_nut_amp(parn, apr_val, apr_adj, 
     .         cov_parm, sol_parm,
     .         num_parn, nut_amp_est, nut_amp_adj, nut_amp_sig, output )

      implicit none 
 
*     Routine to compute the nutation ampltidues from the nutation
*     series coefficients.
 
*   num_parn    - Number of parameters estimated (total)
*   parn(4)     - Parmaeter numbers for long in, long out,
*               - obliquity in, obliquity out
 
      integer*4 num_parn, parn(4)
 
*   apr_val(4)      - Apriori values of nutaions in long and obl
*   apr_adj(4)      - Adjustments to apriori from val_nut_coeff/
*   nut_amp_est(4)  - Four values of the ampltiudes (a+ in,
*                   - a+ out, a- in, a- out)
*   nut_amp_adj(4)  - Four adjustment values
*   nut_amp_sig(4)  - Sigmas of the four values
 
*   cov_parm(num_parn, num_parn)    - Covarinace matrix
*   sol_parm(num_parn)              - Solution vector
 
      real*8 apr_val(4), apr_adj(4), nut_amp_est(4), nut_amp_adj(4),
     .    nut_amp_sig(4), cov_parm(num_parn, num_parn),
     .    sol_parm(num_parn)
 
*   output          - Indicates that there are values to be
*                   - output.
 
      logical output
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop counters
 
      integer*4 j
 
*   sine        - Sin of oblquity of ecliptic
*  nut_coe_adj(4)   - Adjustments to coefficients
*   nut_coe_sig(4)  - Sigmans of the coefficients
 
      real*8 sine, nut_coe_adj(4), nut_coe_sig(4)
 
      data sine / 0.39777716d0 /
 
 
****  START see if there is anything to output
      output = .false.
      do j = 1,4
          if( parn(j).gt.0 ) output = .true.
      end do
 
      if( .not.output ) RETURN
 
***** OK, we have values, not compute ampltitudes
 
      do j = 1,4
          if (parn(j).ne.0 ) then
              nut_coe_adj(j) = sol_parm(parn(j))
              nut_coe_sig(j) = sqrt(cov_parm(parn(j),parn(j)))
          else
              nut_coe_adj(j) = 0.d0
              nut_coe_sig(j) = 0.d0
          end if
      end do
 
*     Now compute the amplitudes and their sigmas
      nut_amp_adj(1) = -(nut_coe_adj(1)*sine+nut_coe_adj(3))/2
      nut_amp_adj(2) =  (nut_coe_adj(2)*sine-nut_coe_adj(4))/2
      nut_amp_adj(3) =  (nut_coe_adj(1)*sine-nut_coe_adj(3))/2
      nut_amp_adj(4) =  (nut_coe_adj(2)*sine+nut_coe_adj(4))/2
 
      nut_amp_est(1) = -((apr_val(1)+apr_adj(1))*sine
     .                 + (apr_val(3)+apr_adj(3)) )/2+nut_amp_adj(1)
      nut_amp_est(2) =  ((apr_val(2)+apr_adj(2))*sine
     .                 - (apr_val(4)+apr_adj(4)) )/2+nut_amp_adj(2)
      nut_amp_est(3) =  ((apr_val(1)+apr_adj(1))*sine
     .                 - (apr_val(3)+apr_adj(3)) )/2+nut_amp_adj(3)
      nut_amp_est(4) =  ((apr_val(2)+apr_adj(2))*sine
     .                 + (apr_val(4)+apr_adj(4)) )/2+nut_amp_adj(4)

      nut_amp_sig(1) = sqrt((nut_coe_sig(1)*sine)**2
     .                      +nut_coe_sig(3)**2)/2
      nut_amp_sig(2) = sqrt((nut_coe_sig(2)*sine)**2
     .                      +nut_coe_sig(4)**2)/2
      nut_amp_sig(3) = sqrt((nut_coe_sig(1)*sine)**2
     .                      +nut_coe_sig(3)**2)/2
      nut_amp_sig(4) = sqrt((nut_coe_sig(2)*sine)**2
     .                      +nut_coe_sig(4)**2)/2
 
****  Thats all
      RETURN
      end
