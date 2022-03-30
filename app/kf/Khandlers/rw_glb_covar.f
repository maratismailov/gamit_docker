CTITLE RW_GLB_COVAR
 
      subroutine rw_glb_covar( option, dcb, cov )
 
      implicit none

*     Routine to open/creat/read/write the temporary copies
*     of the covariance matrices for the global Kalman filter.
* MOD TAH 190531: Updated to hande cov matrix with > 32767^2 values.
*     (Values with 128*nblks could be > 32767^2)
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
*                                         ! Needed for GLB_SOL_FILE
      include '../includes/globk_markov.h'
 
*   dcb(1)  - Dcb buffer
*   cov(1j) - Covariance matrix and solution vector to be
*           - operated on
*   ierr    - Error flag
*   nblks   - number of blocks in each covariance matrix
*   trimlen - HP function for length of string
 
      integer*4 dcb(16), cov(*), ierr, nblks, trimlen
 
*   jsize   - Size of the file to be created (words)
*   last_cov_rec    - Record number of the last covariance matrix
*           - in file
 
      integer*4 jsize, last_cov_rec
 
*   lenr  - Dummy for return of words read from file
*   lenr8 - I*8 version for readd8
 
      integer*4 lenr
      integer*8 lenr8
 
*   option  - Option takes on four values:
*           - C -- create file
*           - O -- open file
*           - R -- read the next covariance matrix
*           - P -- read the previous covariance matrix
*           - W -- write the covariance matrix
*           - L -- read last covariance matrix in the file
      character*(*) option
 
*   opt     - Upper case version if option
 
      character*1 opt
 
*   full_file_name  - Full name of the file with size and type
*           - added.
 
      character*132 full_file_name

      integer*8 i8  !  i8 value of 1 to force I8 calculations

      logical big   ! Set true 128*nbls will exceed I*4 variable

      data I8 / 1 /
 
 
***** Get the upper case version of the option 
      opt = option
      call casefold( opt )
 
*     Get number of blocks for each covariance matrix
 
C     nblks = 2*(num_glb_parn+1)*num_glb_parn/128 + 1
      nblks = 2*(num_glb_parn+i8)*num_glb_parn/128 + 1
* MOD TAH 190531: See if solution is "big"
      big = .false.
      if( num_glb_parn.gt.32767 ) big = .true. 
 
***** See what we want to do
 
      if( trimlen( glb_sol_file ).eq. 0 ) then
          glb_sol_file = glb_sol_default
      end if
 
*                             ! Create the file glb_sol_file
      if( opt.eq.'C' ) then
 
*         Purge file fisrt
          call clean_file( glb_sol_file )
 
*         Compute the file size
*                                     ! Save space for all experiments
          if( glb_bak_soln ) then
              jsize = num_glb_sol*nblks
          else
              jsize = nblks
          end if
 
*         Create file
          call FullFileName( glb_sol_file, 1, jsize, 128,
     .                       full_file_name)
          call FmpOpen( dcb, ierr, full_file_name, 'RWC', 1)
          call report_error('FmpOpen',ierr,'creat',full_file_name,
*                                             ! Kill if error
     .                      1,'RW_GLB_COVAR')
      end if
 
*                            ! open the file
      if( opt.eq.'O' ) then
 
          call FmpOpen( dcb, ierr, glb_sol_File, 'RWO', 1)
          call report_error('FmpOpen',ierr,'open',glb_sol_file,
*                                              ! Kill if error
     .                       1,'RW_GLB_COVAR')
      end if
 
*                             ! Read the next covariance matrix in the
      if( opt.eq.'R' ) then
*                             ! file
* MOD TAH 190531: See if "big"
          if ( big ) then
             call readd8(dcb,ierr, cov, (I8*128)*nblks, lenr8, 0)
          else   ! Old code 
             call readd(dcb,ierr, cov, 128*nblks, lenr, 0)
          endif 
          call report_error('VREAD',ierr,'read',glb_sol_file,
     .                      1,'RW_GLB_COVAR')
      end if
 
*                             ! Read the previous covariance matrix in
      if( opt.eq.'P' ) then
*                             ! the file.
          call eposn(dcb, ierr,-2*nblks,0 )
          call report_error('FMP',ierr,'back position',glb_sol_file,
     .                      1,'RW_GLB_COVAR')
          if ( big ) then
             call readd8(dcb,ierr, cov, (I8*128)*nblks, lenr8, 0)
          else
             call readd(dcb,ierr, cov, 128*nblks, lenr, 0)
          end if
          call report_error('VREAD',ierr,'read',glb_sol_file,
     .                      1,'RW_GLB_COVAR')
      end if
 
*                             ! Write the next covariance matrix in the
      if( opt.eq.'W' ) then
*                             ! file
* MOD TAH 190531: See if "big"
          if ( big ) then
             call writd8(dcb,ierr, cov, (I8*128)*nblks,0)
          else
             call writd(dcb,ierr, cov, 128*nblks,0)
          endif 

          call report_error('VWRIT',ierr,'writ',glb_sol_file,
     .                      1,'RW_GLB_COVAR')
      end if
 
*                             ! Read the last cov matrix
      if( opt.eq.'L' ) then
 
          if( glb_bak_soln ) then
              last_cov_rec = (num_glb_sol-1)*nblks + 1
          else
              last_cov_rec = 1
          end if
          call eposn(dcb, ierr,last_cov_rec,1)
          call report_error('FMP',ierr,'last position',glb_sol_file,
     .                      1,'RW_GLB_COVAR')
* MOD TAH 190531: See if "big"
          if ( big ) then
             call readd8(dcb,ierr, cov, (I8*128)*nblks, lenr8, 0)
          else
             call readd(dcb,ierr, cov, 128*nblks,lenr, 0)
          endif 
          call report_error('VREAD',ierr,'read',glb_sol_file,
     .                      1,'RW_GLB_COVAR')
      end if
 
***** Thats all
      return
      end
 
