CTITLE DECODE_GLBAK_RUN
 
      subroutine decode_glbak_run( glb_com_file, glorg_command_file,  
     .                             com_orgopt )

      implicit none  
 
*     Routine to decode runstring for GLBAK
* MOD TAH 190525: Introduced when globk modified realize reference
* frame using GLORG apply_cond_full routines
 
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
*   glorg_command_file  - OPtional name of command file
*   outfile         - Output file name 
*   com_orgopt  - Optional string at start of lines 
 
      character*(*) glb_com_file, glorg_command_file,
     .              com_orgopt
 
*   runstring       - Runstring for integer values
      character*128 runstring
 
***** Get the com file name:  This must be passed.
      len_run = rcpar(1, glb_com_file )
      if( len_run.eq.0 ) then 
          write(*,100) 
 100      format('***DISASTER*** No common file passed')
          stop 'GLBAK: No common file name'
      end if


*     Get glorg_command file (if it passed) otherwize use the
*     file name in common file.
      len_run = rcpar(2, runstring)
      if( len_run.gt.0 ) then
          glorg_command_file = runstring
      else
          glorg_command_file = ' '
      end if
      

****  See if command read option passed (These are strings that
*     will be removed from the starts of lines).
      len_run = rcpar(3, com_orgopt)
 
***** Thats all
      return
      end
