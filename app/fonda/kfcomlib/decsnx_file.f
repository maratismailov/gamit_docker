CTITLE DECSNX_FILE
 
      subroutine decsnx_file (unit, line )
 
*     Routine to decode the file  blocks from SINEX:
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit
 
*   line        - Line read from input
 
      character*(*) line
 
* LOCAL VARIABLES
 

****  Start Skip over the files blocks
      call snx_finbl(unit)
 
****  Thats all
      return
      end
 
 
