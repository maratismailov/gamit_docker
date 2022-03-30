CTITLE GLOUO
 
      subroutine glouo(ms_type, outfile)

      implicit none  
 
*     This segment decodes the runstring and opens the files
*     we will need
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* PASSED variables
* ms_type  - Sets the type of program call (MAIN for program)

      character*(*) ms_type, outfile
 
*   dcb(16)     - Dcb buffer for reading the global
*               - covariance matrix
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   ierr        - SEGLD error flag
*   iout        - Output LU
*   options     - options passed through runstring
 
      integer*4 dcb(16), ema_data(1), iout, options, ierr,
     .          num

*   kbit         - Check bit status
      logical kbit       
 
 
      common / progcon / iout, options
 
 
      common / globc_ema / ema_data
 
****  Decode the runstring
      iout = 100
      if( ms_type(1:4).eq.'MAIN' ) then
 
          call decode_glout_run( iout, options, glb_com_file,
     .                           glb_com_default)

      else
      
*         See if option give to erase file first.
          if( kbit(options,17) ) then
              call open_lu(iout,outfile,ierr,'unknown') 
          else
               call open_lu(iout,outfile,ierr,'append') 
          end if
          call report_error('IOSTAT',ierr,'open_lu',outfile,
     .            0,'GLOUT/GLOUO')
           if( ierr.ne.0 ) iout = 6
      end if

 
*     Now open the GLOBK common and read
      if( ms_type(1:4).eq.'MAIN' ) then 
          call rw_globk_common('R')
          istart_vma = 0
      end if

*     Now open and read the covariance matrix and solution vector
      if( ms_type(1:4).eq.'MAIN' ) then
          call glfor_mem( ema_data ) 
      end if

      call glb_out_map
*     This is needed on the HP
      inquire(file=glb_sol_file,iostat=ierr, number=num)
      if( num.gt.0 ) close(num)

      call rw_glb_covar('O', dcb, ema_data(icov_parm) )
      call rw_glb_covar('L', dcb, ema_data(icov_parm) )
 
***** Thats all
      return
      end
 
