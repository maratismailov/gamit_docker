CTITLE DECODE_GLOUT_RUN
 
      subroutine decode_glout_run( iout, options, glb_com_file,
     .                             glb_com_default )

      implicit none 
 
 
*     Routine to decode runstring
 
*   DecimalToInt    - Converts string to integer value
*   ierr            - Error during conversion
*   iout            - Ouput Lu
*   len_run         - Length of runstring
*   LogLu           - HP function for users LU
*   options         - Options for output
*   rcpar           - HP function to read runstring
*   lenfn           - Length of output file name
 
      integer*4 ierr, iout, len_run, LogLu, options,
     .    rcpar, indx, lenfn
     
*   kbit         - Check bit status
      logical kbit       
 
*   glb_com_file    - Name of the common file
 
      character*(*) glb_com_file, glb_com_default
 
*   runstring       - Runstring for integer values
*   outfile         - Output file name 
      character*128 runstring, outfile
 

***** Get the Lu for output      
      lenfn = rcpar(1, outfile)
 
*     Get options
      len_run = rcpar(2, runstring)
      if( len_run.gt.0 ) then
          indx = 0
          call decode_prt_opt(runstring, indx, options)
*                     ! Minimum output
      else
* MOD TAH 190603: Print help and quit if no runstring
          call proper_runstring('glout.hlp','glout',1)
          options = 0
      end if

!     write(*,120) len_run, options, trim(runstring)
!120  format('OPTIONS ',i4,1x,o12,1x,a)

      
* MOD TAH 970909: Moved file open to after runstring to see
*     if we need to erase first (ERAS option).      
      if( lenfn.gt.0 ) then
          iout = 200
          if( kbit(options,17) ) then
              call open_lu(iout, outfile, ierr, 'unknown')
          else
              call open_lu(iout, outfile, ierr, 'append')
          endif 
          call report_error('IOSTAT',ierr,'open_lu',outfile,
*                                                 ! Dont kill
     .                      0,'Decode_GLOUT_run')
          if( ierr.ne.0 ) then
              iout = LogLu( ierr )
          end if
*                     ! use user's terminal
      else
          iout = LogLu(ierr)
      end if
 
*     Get optional comfile name
      len_run = rcpar(3, glb_com_file)
      if( len_run.eq.0 ) then
            glb_com_file = glb_com_default
            write(*,100) glb_com_default  
 100        format(/,'**WARNING** No common file name passed using',
     .             ' default: ',a)
      end if
 
***** Thats all
      return
      end
 
