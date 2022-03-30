CTITLE DECODE_SNXB
 
      subroutine decode_snxb(unit, line, np, cov_parm, sol_parm)
 
      implicit none

*     This is the upper level routine that decodes the type of
*     SINEX block found.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
* unit   - Unit for reading the file
* np     - Number of parameters in this solution
 
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
* line   - Line read from file
 
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   indx        - Pointer in string
*   trimlen - Length of string
 
      integer*4 indx, trimlen
 
*   block_found - True if block found
  
      logical block_found
 
****  Start check which type of block this is:
 
      block_found = .false.
 
      indx = index(line, 'FILE')
      if( indx.eq.2 ) then
          block_found = .true.
          call decsnx_file(unit, line)
      end if
 
      indx = index(line, 'INPUT')
      if( indx.eq.2 ) then
          block_found = .true.
          call decsnx_input(unit, line)
      end if
 
      indx = index(line, 'SITE')
      if( indx.eq.2 ) then
          block_found = .true.
          call decsnx_site(unit, line)
      end if

      indx = index(line, 'SATELLITE')
      if( indx.eq.2 ) then
          block_found = .true.
          call decsnx_sat(unit, line)
      end if
 
      indx = index(line, 'SOLUTION')
      if( indx.eq.2 ) then
          block_found = .true.
          call decsnx_soln(unit, line, np, cov_parm, sol_parm) 
      end if
 
****  Check to see if block was found
      if( .not.block_found ) then
          write(*,500) line(1:trimlen(line))
 500      format('Sinex Block type not regonized',/,a)
 
          call snx_finbl(unit)
      end if
 
*     Thats all
      return
      end
 
