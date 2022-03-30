 
      program modear

      implicit none 
 
*     Program to model the phase and range measurements in an Earth
*     fixed frame taking as input a cfile and SP3 file.
*
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
 
      include 'modear.h'
 
* MAIN Program variables
 
*   len_run - Length of runstring
*   rcpar   - Function to read runstring
*   fmprename  - Function to rename file
*   ierr        - Rename error flag.
*   trimlen - Length of string 
*   null_terminate - Puts a null at the end of string
 
      integer*4 len_run, rcpar, fmprename, ierr, trimlen, null_terminate
 
****  Decode the runstring
      write(*,120)
 120  format(/' MODEAR: Earth fixed modeling of cfile')
      len_run = rcpar(1,in_cf)
      if( len_run.le.0 ) then
          call proper_runstring('modear.hlp','modear/cf name',1)
      end if
      write(*,130) in_cf(1:len_run)
 130  format(' Input cfile        : ',a)
 
      len_run = rcpar(2, sp3_file )
      if( len_run.le.0 ) then
          call proper_runstring('modear.hlp','modear/sp3 file',1)
      end if
      write(*,140) sp3_file(1:len_run)
 140  format(' SP3 Ephemeris file : ',a)
 
      len_run = rcpar(3, nav_file )
 
*     The nav file is optional if the SP3 file has clock in it
*     otherwize it is needed.
      if( len_run.le.0 ) then
          nav_file = 'NONE'
          len_run  = 4
      end if
      write(*,150) nav_file(1:len_run)
 150  format(' Navigation file    : ',a)
 
****  Now open and read the sp3 file and the Nav file
      call read_sp3
 
      if( nav_file(1:4).ne.'NONE' .and. svclk_OK ) then
         call read_nav
      end if
 
      call clean_sp3_clk

      if( .not.svclk_OK ) then
          write(*,200)
 200      format(' No clock information in sp3 file: Nav file ',
     .        'must be given')
          stop 'MODEAR: Nav file must be used with this SP3 file'
      end if
 
*     Now read the cfile to get the observables
      out_cf = in_cf
      in_cf  = out_cf(1:trimlen(out_cf)) // '.o'
      ierr = null_terminate(out_cf)
      ierr = null_terminate(in_cf)
     
      ierr = fmprename(out_cf, in_cf,'K')
      call report_error('RENAME',ierr,'renam',out_cf,1,'modear')
 
*     Read the cfile and compute the station clocks at each epoch.
*     These will be saved and directly used when the compute the
*     theoretical ranges and phases.
 
      call preread_cf
 
****  Now update the cfile.
      call update_cf
 
****  Thats all
      end
 
