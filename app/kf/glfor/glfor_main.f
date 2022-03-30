      program glfor_main

      implicit none 

*     This is the programs that call glfor.  In globc the subroutine
*     glfor is called directly.  The argument passed here is to show
*     that glfor is being called as a MAIN program.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

      integer*4 pcontrol, cov_dcb(16), ema_data(1)

      common / progcon / pcontrol, cov_dcb
      common / globc_ema / ema_data

* ms_type  - Indicates that this is a program

      character*8 ms_type 

      ms_type = 'MAIN'
      call glfor(ms_type)
C     call glfor('MAIN')

      end 
