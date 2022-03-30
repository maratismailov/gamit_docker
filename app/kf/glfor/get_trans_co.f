CTITLE GET_TRANS_COEFF
 
      subroutine get_trans_coeff ( row, col, coeff, deltat, coeff_in )

      implicit none  
 
*     Routine to get the transission coefficient to be used
*     for the next transission.  The main complication here concerns
*     whether the coefficient should be multipled by deltat,
*     abs(deltat) or not at all
 
*   col     - Column from the transission matrix
*   row     - Row of the transsion matrix
 
      integer*4 col, row
 
*   coeff_in    - Entry from trans_coeff to be used (may need to
*           - modified)
 
      real*4 coeff_in
 
*   coeff   - Transission coefficient to be used
*   deltat  - Epcoh change in years
 
      real*8 coeff, deltat
 
****  Now if the row does not equal the column then just multiple
*     by deltat
 
      if( row.ne.col ) then
          coeff = coeff_in*deltat
*                                 ! Two cases here: if coeff_in is
      else
*                                 ! -1 then just use as the value,
*                                 ! otherwise multiple by abs(deltat)
*                                 ! This later class is dissipation
*         These coefficients are meant to be be with SCALAR multiply
*         (DWSMY) rather than the pivot routines.  (Could be done
*         either way)
          if( coeff_in.ne.-1.0 ) then
              coeff = 1.d0 + coeff_in*abs(deltat)
          else
              coeff = 0.d0
          end if
      end if
 
***** Thats all
      return
      end
 
