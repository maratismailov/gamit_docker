CTITLE GLB_DOT_PART
 
      subroutine glb_dot_part( value, cov, step, part, pnt, indx)

      implicit none  
 
*     Routine to dot the partials with a covarince matrix when the
*     partials are in compressed form.
 
*   i       - Loop counter
*   indx    - Number of pairs of values in pnt
*   num     - Number of partials in current group
*   pel     - keeps track of position of the multiplication
*           - in the partials array
*   pnt(2,1)    - Pointers to parameter number and number of
*           - contiguous partials
 
*   step    - the step size to be used in going through cov
*           - =1 down a column
*           - = dim of matrix, across a row
 
      integer*4 i, indx, num, pel, pnt(2,1), step
 
*   iel     - Element number in cov to be used (used for convenience)
*           - (Note: we make this I*4 so that matrices > 180*180
*           - can be processed
      integer*4 iel
 
*   cov(1j)     - One column of covariance matrix
*   part(1)     - One compressed row of partials matrix
*   temp        - Temporary storage of result
*   value       - Result of the dot product
 
      real*8 cov(*), part(*), temp, value
 
 
*     Initialize value
      value = 0.d0
 
*     See if there are any partials for this observation
      if( indx.eq.0 ) RETURN
 
*     Loop over non-zero partials
      pel  = 1
      do i = 1,indx
 
*         Get the parameter number in cov
          iel = (pnt(1,i)-1)*step + 1
          num = pnt(2,i)
 
          call DWDOT( temp, cov(iel),step, part(pel),1, num)
 
*         ADD to total dot product
          value = value + temp
 
*         Update position in partial array
          pel = pel + num
      end do
 
****  Thats all
      return
      end
 
