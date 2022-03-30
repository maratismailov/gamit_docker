CTITLE DECODE_GLOBK_RUN
 
      subroutine decode_globk_run
 
      implicit none 

 
*     Routine to decode the globk runstring.  See main for desription
*     of contents of runstring
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   DecimalToInt    - Conversion from string to integer
*   ierr    - DecimalToInt error
*   len_run - Length of runstring
*   loglu   - HP function for user's Lu
*   rcpar   - HP function to read runstring
 
      integer*4 DecimalToInt, ierr, len_run, loglu, rcpar
 
*   runstring   - Part of runstring being decoded
 
 
      character*128 runstring
 
***** See if crt unit passed
 
      len_run = rcpar(1,runstring)
      if( len_run.gt.0 ) then
          crt_unit = DecimalToInt(runstring,ierr)
          call report_error('DecimalToInt',ierr,'decod',runstring,0,
*                                                   ! Kill below
     .                      'DECODE_GLOBK_RUN')
*                                                   ! Kill
          if( ierr.ne.0 ) call proper_runstring('globk.hlp','globk',-1)     
      else
          crt_unit = loglu(ierr)
      end if
 
*     See if printer unit passed
 
      len_run = rcpar(2,runstring)
      if( len_run.gt.0 ) then
          glb_prt_file = runstring
      else
          glb_prt_file = '6'
          prt_unit = 6
      end if
 
*     See if log_unit passed
 
      len_run = rcpar(3,runstring)
      if( len_run.gt.0 ) then
          glb_log_file = runstring
          if( len_run.eq.1 ) then
              log_unit = DecimalToInt(runstring,ierr)
          end if
      else
          glb_log_file = '6'
          log_unit = 6
      end if

*     Get the list file name (Contains the list of global files to be
*     included in this run)
 
      len_run = rcpar(4,runstring)
      if( len_run.gt.0 ) then
          list_file = runstring
      else
*                                                ! Kill
          call proper_runstring('globk.hlp','globk',-1)     
      end if
 
*     Get the markov command file name
 
      len_run = rcpar(5,runstring)
      if( len_run.gt.0 ) then
          glb_mar_file = runstring
      else
*                                                ! Kill
          call proper_runstring('globk.hlp','globk',-1)     
      end if
 
*     Now check the optional parameters

* MOD TAH 030607: Optional string for starts of lines
      len_run = rcpar(6,runstring)
      if( len_run.gt.0 ) then
          comopt = runstring
      else
          comopt = ' '
      end if 
      
*     Common file name
      len_run = rcpar(7,runstring)
      if( len_run.gt.0 ) then
          glb_com_file = runstring
      else
*                             ! Open routine will set default name
          glb_com_file = ' '
      end if
 
*     Sort file name
      len_run = rcpar(8,runstring)
      if( len_run.gt.0 ) then
          sort_file = runstring
      else
*                          ! Set default name
          sort_file = glb_sort_default
      end if
 
*     Sort direction.
      len_run = rcpar(9,runstring)
      if( len_run.gt.0 ) then
          sort_direction = DecimalToInt(runstring, ierr)
C         call report_error('DecimalToInt',ierr,'decod',runstring,
*                                                 ! Use default if error
C    .                      0,'DECODE_GLOBK_RUN')
C
*                                                ! error
          if( ierr.ne.0 ) then
              sort_direction = +1
          end if
      else
          sort_direction = +1
      end if
 
***** Thats all
      return
      end
 
