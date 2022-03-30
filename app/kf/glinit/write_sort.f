CTITLE WRITE_SORT
      subroutine write_sort( expt_names, expts_var)

      implicit none 
 
*     This routine will write out the sorted list of experiment
*     names to a type 2 file with a 32 word (64 character) record
*     length.  This file may then be direct accessed for either
*     running forwards or backwards through the experiments.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
 
*   i                               - Loop counter
*   ierr                            - IOSTAT error
*   expt_names(max_glb_sol )        - Names of the global
*                                   - solutions
 
      integer*4  i, ierr
      character*(sort_recl)  expt_names(max_glb_sol )
 
*   expts_var(max_glb_sol)          - Variances to be given to
*                                     each experiment (to make chi**2
*                                     unity).
 
      real*8 expts_var(max_glb_sol)

      integer*4 trimlen
 
* MOD TAH 980519: Added reading and writing of forward and 
*   back chi**2/f to srt_file.
*   for_chi, bak_chi -- Forward and backwarsd chi**2/f
      real*4 for_chi, bak_chi
 
***** First purge and then create the sort file.
      if( trimlen(sort_file).eq.0 ) sort_file = glb_sort_default
      write(*,100) sort_file(1:trimlen(sort_file))
 100  format('Cleaning and creating sort_file ',a)
 
      call clean_file ( sort_file )
      open(200, file=sort_file, iostat=ierr, access='direct',
     .     status='unknown', recl=sort_recl+8+8)
      call report_error('IOSTAT',ierr,'open', sort_file, 1,
     .                  'WRITE_SORT')
 
***** Now write out the file names
      for_chi = -1
      bak_chi = -1
       
      do i = 1, num_glb_sol
          write(200,iostat=ierr,rec=i) expt_names(i), expts_var(i),
     .                                for_chi, bak_chi
          call report_error('IOSTAT',ierr,'writ',sort_file,1,
     .                      'WRITE_SORT')
      end do
 
***** Thats all
      close(200, iostat=ierr)
      call report_error('Close',ierr,'clos',sort_file, 0,
     .                  'WRITE_SORT')
 
      return
      end
 
