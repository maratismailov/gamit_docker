ctitle
 
      subroutine rw_covar(option, cov_parm)

      implicit none
 
c
c     Routine to read or write the covariance matrices from a solution
c     from a disk file.
c     If a back solution is to be run, then these matrices are wriiten
c     for each epoch of data in the solution.  Only one set is written
c     if there is no back solution.  The file will be used by 'OUTSL'
c     and 'BAKSL'.
c
C MOD TAH 870825 Replaced VREAD and VWRIT with I*4 length versions
C
c Include files
c -------------
*                                  ! The solvk parameter file
      include '../includes/kalman_param.h'
c
*                                  ! the solvk control common block
      include '../includes/solvk_cntl.h'
*                                  ! the Values block of KalObs file
      include '../includes/obs_values.h'
 
c Variables
c ---------
c option -- the option for the file manipulation.  These options are
c     'write' -- write the file, open or create if necessary
c     'read'  -- read the file at the next epoch of records
c     'read_last' -- read the last record of the file.
c cov_parm -- the covariance matrix of the solution
c
 
      character*(*) option
 
      real*8 cov_parm(par_num,par_num)
 
c
c
c Local variables
c ---------------
c file_not_open -- logical to indicate if file already open
c file_name -- the integer version of the file name
c file_cart -- the cartridge for the solution file
c jsize -- the size of the file (number of records and record length)
c     if the file needs to be created with ecreat.
c nbuf -- number of dcb buffer available.
c ifildc -- the sol_file dcb buffer.
c post_not_set -- indicates if the file has been positioned for reading
c last_cov_rec -- record number of the last covariance matrix in file
c
      logical file_not_open, post_not_set
 
c
*   ierr        - FMP error code.
*   lenr        - Number of words read (dummy)

      integer*4 nbuf, ifildc(144), ierr, lenr
c
      integer*4 jsize, last_cov_rec
c
c Scratch common
c --------------
c ident -- the idenification of the use of scratch common
c scr_data -- the scratch common area.  Set to the size for forsl
c     since this is the smallest scratch common available
c
 
      integer*4 ident, scr_data(scr_forsl)
 
*   full_file_name  - The full file name including the size if
*                   - file needs to be created.
      character*128 full_file_name
 
c
      common ident, scr_data
 
c
      save file_not_open, post_not_set
 
      data  file_not_open / .true. /,
     .      post_not_set  / .true. /
c
c.... Set the use of the scrarch common area
      ident = 6
c
c.... See if we need to open file
*                                 ! Try to open file
      if( file_not_open ) then
c
c....    Compute how large a dcb buffer we may use
         nbuf = (scr_forsl-17)/128
c
c....    Try to open file if we are reading it
*        Force FmpOpen to use next unit number.
         scr_data(1) = 0 
         call FmpOpen(scr_data,ierr, sol_file, 'rwo', nbuf)
c
c....    Check error
*                                                           !create file
         if( ierr.eq.-6 .and. option(1:5).eq.'write' ) then
c
c....       see how large the file will need to be
*                                  ! needs to be big enough for all epochs
            if( bck_soln ) then
               jsize = nblocks*total_epochs
*                                  ! just long enough for a single epoch
            else
               jsize = nblocks
            end if
c
c....       Try to create the file
 
            call FullFileName( sol_file, 1, jsize, 128,
     .                         full_file_name )
 
            write(icrt,100) full_file_name
 100        format(" Creating ",a)
 
            call FmpOpen(scr_data,ierr,full_file_name,'rwc', nbuf)
 
         end if
c
c....    Print any error message if ierr less than zero
         call report_error('FMGR',ierr,'open/creat',sol_file,1,
     .      'RW_COVAR')
c
c....    Save the DCB buffer and set file open
         call wmove(scr_data,ifildc,144)
         file_not_open = .false.
c
*                  ! file not open
      end if
c
c.... Copy the dcb buffer into the scratch area
      call wmove(ifildc, scr_data, 144)
c
c.... See if we are to write file or read file
*                                        ! write out covariance matrix
      if( option(1:5).eq.'write' ) then
c
c....    Use HP VWRIT routine to go from ema to disk
         call writd(scr_data, ierr, cov_parm, 128*nblocks,0)
         call report_error('FMGR',ierr,'writ',sol_file,1,'RW_COVAR')
c
*                                               ! read the file
      else
c
c....    See if we are reading backwards or if we wish to read last
c        covariance matrix.
c
         if( option(1:9).eq.'read_last' .or. post_not_set ) then
c                    ! position to last covariance matrix
c
c....       Get the record number of the last covariance matrix
            if( bck_soln ) then
               last_cov_rec = (total_epochs-1)*nblocks + 1
            else
               last_cov_rec = 1
            end if
c
c....       position the file
            call eposn(scr_data,ierr,last_cov_rec,1)
            post_not_set = .false.
c
*                                               ! position previous covar
         else
c
            call eposn(scr_data,ierr,-2*nblocks,0)
c
         end if
c
c....    Check the positioning error
         call report_error('FMGR',ierr,'position',sol_file,1,
     .      'RW_COVAR')
c
c....    Now read the covariance matrix
         call readd(scr_data,ierr,cov_parm,128*nblocks,lenr, 0)
c
c....    Check error
         call report_error('FMGR',ierr,'read',sol_file,1,'RW_COVAR')
c
*                       ! reading or writing the file
      end if
c
c.... Save the dcb before exiting
      call wmove(scr_data,ifildc,144)
c
      return
      end
 
 
