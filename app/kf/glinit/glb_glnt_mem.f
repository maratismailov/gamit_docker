CTITLE GLB_GLNT_MEM
 
      subroutine glb_glnt_mem( iexpt_names, isepoch_expts, iexpts_var,
     .            iapr_values, ema_data )

      implicit none 
 
*     This routine does the memory assignment for glinit.  We do
*     dynamically so that it can be de-assigned later.
 
      include '../includes/kalman_param.h'
 
* PASSED variables
 
*   iexpt_names     - Start of experiment names
*   isepoch_expts   - Start of epochs
*   iexpts_var      - Start of experiment variances
*   iapr_values     - Start of apriori values
 
*   ema_data(*)     - Place where addresses will be computed
*                   - from
 
      integer*8 iexpt_names, isepoch_expts, iexpts_var, iapr_values
      integer*4 ema_data(*)
 
* LOCAL variables
 
*   mallocg          - Allocates memory
*   memloc          - Location where memory is available
 
*   vma_i4wrd       - Number of i*4 words needed
*   vma_start       - Address of ema_data
 
 
      integer*4 vma_i4wrd
      integer*8 memassign
 
****  Start by computing the number of i*4 words needed
      vma_i4wrd = max_glb_sol*(sort_recl/4+4)+max_glbapr_sz*2+512
      isepoch_expts = memassign(vma_i4wrd,1, loc(ema_data))

      if( isepoch_expts.eq.0 ) then
          write(*,120) vma_i4wrd/(1024.0**2)*4
 120      format('*** DISASTER *** Not enough memory.  Try a larger',
     .                ' computer ',F8.2,' Mbytes needed')
          stop 'GLINIT: Not enough memory'
      end if
 
  
      iexpts_var    = isepoch_expts + 2*max_glb_sol
      iapr_values   = iexpts_var    + 2*max_glb_sol
      iexpt_names   = iapr_values   + 2*max_glbapr_sz + 128

****  Thats all
      return
      end
 
