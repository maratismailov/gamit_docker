CTITLE GLFOR_MEM
 
      subroutine glfor_mem( vma_data )
 
      implicit none 

*     Routine to assign the memory for the glfor run.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED variables
 
*   vma_data(1)     - Location in memory where we
*                   - start the addressing from
 
 
      integer*4 vma_data(1)
 
* Local Variables
 
*   vma_i8wrd       - Numner of bytes needed for this
*                   - run.
*   ferr            - Error returned from free (to free memory)
*   freeg            - Routine to free memory
*   memassign8      - Integer*8 routine to assign memory

* MOD TAH 150524: Introduced integer*8 version for large 
*     solutions >20000 parameters 
      integer*8 vma_i8wrd
      integer*8 memassign8

***** Start: set the starting address and get the amount
*         of memory needed.
c     if( istart_vma.ne.0 ) then
c         call freeg(loc(vma_data(istart_vma)))
c     end if

C     write(*,100) istart_vma, most_cparn_num
C100  format('MEMORY: Current vma start is ',i12,' Mapping for ',
C    .       i6,' parameters in cparn')

      if ( .not. glb_bak_soln ) then
          call glb_for_map( most_cparn_num, vma_i8wrd)
      else
          call glb_bak_map( most_cparn_num, vma_i8wrd)
      end if


*     Allocate the memory for the run 
      istart_vma = memassign8(vma_i8wrd,1,loc(vma_data))

      if( istart_vma.eq.0 ) then
          write(*,120) vma_i8wrd/(1024.0**2)*4
 120      format('*** DISASTER *** Not enough memory.  Try a larger',
     .                ' computer ',F12.2,' Mbytes needed')
          stop 'GLFOR: Not enough memory'
      end if
 
      write(*,140) vma_i8wrd/(1024.d0**2)*4
 140  format(' Allocating ',f12.2,' Mbytes of memory for run')
****  Thats all
      return
      end
 
 
 
