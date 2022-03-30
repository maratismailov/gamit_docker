CTITLE CLEAR_PNT
 
      subroutine clear_pnt ( num, part_pnt )

      implicit none  
 
*     Routine to clear the pointers for the partial derivative matrices
 
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
 
*   i       - Loop counter
*   num     - Number of values to clear
*   part_pnt(2,max_glb_deriv,1) - Pointers to which partials should
*           - be used.
 
      integer*4 i, num, part_pnt(2,max_glb_deriv,1)
 
 
***** Loop over all local parameters
 
      do i = 1, num
          indx_pnt(i)     = 0
*                                 ! No parameter number
          part_pnt(1,1,i) = 0
*                                 ! zero number in group
          part_pnt(2,1,i) = 0
      end do
 
***** Thats all
      return
      end
 
 
