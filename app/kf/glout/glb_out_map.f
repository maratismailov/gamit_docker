CTITLE GLB_OUT_MAP
 
      subroutine glb_out_map

      implicit none  
 
*     Routine to map the use of the ema/vma area.  Currently space
*     is set aside for the covariance matrix and the solution vector
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'

* MOD TAH 190511: Introduced when number of parameters > 32767.
      integer*8 i8  !  i8 value of 1 to force I8 calculations
 
***** Add up the space we need
      i8 = 1
      icov_parm = istart_vma
      isol_parm = icov_parm + 2*(i8*num_glb_parn)*num_glb_parn
 
***** Thats all
      return
      end
 
