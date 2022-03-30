CTITLE COEFF_OC
 
      subroutine coeff_OC( type, indx, apr_val,  sol_obs )

      implicit none  
 
*     Routine to remove the apriori values from the coefficients
*     for either the extendended earth tides or the nutation series.
 
 
*   comp        - Series component
*   indx        - Indicates in or out of phase and the series term
*   phase       - Phase (1 for in phase and 0 for out of phase)
*   type        - Indicates type of parameter (33 for nutation)
 
      integer*4 comp, indx, phase, type
 
*   apr_val(2,1)  - Apriori values for the coefficients
*   sol_obs       - Solution vector read from solution
 
 
      real*8 apr_val(2,1), sol_obs
 
 
***** Get the phase and component
 
      comp  = mod( indx, 128)
      phase = (indx - comp)/128 + 1
 
      sol_obs = sol_obs - apr_val(phase,comp)
 
      return
      end
 
 
