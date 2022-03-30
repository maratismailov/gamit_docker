CTITLE DECODE_GLINIT_RUN
 
      subroutine decode_glinit_run

      implicit none 
 
 
*     Routine to read the GLINIT runstring.  The form is
*
*     glinit, list_file, <glb_com_file>, <sort_file>, 
*             <sort_direction>, <eq_file> <glb_svs_file>
*
*     where all arguments are optional except for the name of the
*     list_file which gives the list of global solutions to be
*     used in this solution.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   DecimalToInt    - HP function for string to integer value
*   ierr            - Error flag for DecimalToInt
*   len_runstring   - Length of the runstring read
*   rcpar       - HP function for reading runstring
 
      integer*4 DecimalToInt, ierr, len_runstring, rcpar,
     .          trimlen
 
*   runstring   - Place to read runstring into
 
      character*64 runstring
 
****  Get the name of the list_file
 
      len_runstring = rcpar(1, list_file )
*                                         ! Error No name passed
      if( len_runstring.eq.0 ) then
          write(*,100)
  100     format(/' Incomplete runstring: You should use,',/,
     .            ' GLINIT, Input list, glb_com_file, sort_file,',
     .            ' sort_direction, eq_file',/,
     .            ' See GLOBK for usage (This is NOT a user program)')
          stop ' GLINIT Terminated: Incomplete runstring'
      end if
 
***** Now try out the optional arguments
 
      len_runstring = rcpar(2, glb_com_file )
      if( len_runstring.eq.0 ) glb_com_file = glb_com_default
 
      len_runstring = rcpar(3, sort_file )
      if( len_runstring.eq.0 ) sort_file = glb_sort_default
 
      len_runstring = rcpar(4, runstring)
*                                         ! Get sort direction
      if( len_runstring.gt.0 ) then
          sort_direction = DecimalToInt( runstring, ierr )
          call report_error('DecimalToInt',ierr,'decod',runstring,
*                                                  ! No Kill
     .                      0,'Decode_glinit_run')
      end if
 
      if( ierr.ne.0 .or. len_runstring.eq.0 ) sort_direction = +1
 
      len_runstring = rcpar(5, eq_inp_file(1) )
      if( len_runstring.le.0 ) eq_inp_file(1) = ' '
*     check to see if null terminated string.
      if( eq_inp_file(1)(1:1).eq.char(0) ) eq_inp_file(1) = ' '

      len_runstring = rcpar(6, glb_svs_file )
      if( len_runstring.le.0 ) glb_svs_file = ' '
*     check to see if null terminated string.
      if( glb_svs_file(1:1).eq.char(0) ) glb_svs_file = ' '

      if ( trimlen(glb_svs_file).gt.0 ) then
         make_svs_file = .true.
      else
         make_svs_file = .false.
      end if

***** Thats all
      return
      end
 
 
