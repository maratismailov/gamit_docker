 
      program glred

      implicit none 
 
*     Program which will schedule GLOBK for a series of single global
*     files which as given to it.  This program may be used to
*     generate global files with the source positions fixed from ones
*     which have all sources free.
*
*     The runstring is very similar to GLOBK's runstring:
*
*     CI> GLRED,crt,prt,log,input list,markov file,common,sort,list_mask
*
*     Where crt is user's LU
*           prt is print device (LU)
*           log is log unit number
*           input list is the name of the file containing the list of
*                  global files to be processed
*           markov file is the name of the markov control file to be
*                  passed to GLOBK.
*           common [file] is the name of the common file to be used
*                  (Can also be specified in the markov file)
*           sort [file] is the name of the sort file (can also be given
*                  in the markov file)
*           list mask is the mask to be used when the list file for GLOBK
*                  is generated.  Each list file contains the name of the
*                  global file to be processed by GLOBK.  The file name
*                  will be generated from the name of the current global
*                  file being processed.
*
 
      include '../includes/kalman_param.h'
 
*   crt         - User's LU number
*   DecimalToInt    - Converts string to integer
*   i,j         - Loop counters
*   ierr        - IOSTAT error flag
*   indx        - Pointer to position in string for Read_Line
*   iprm(5)     - Parameters returned from FmpRunProgram
*   len_run     - Length of runstring parameter
*   loglu       - Return users LU number
*   rcpar       - HP function to read runstring
*   trimlen     - String length function
*   num_in_gdl  - Number of files listed to the lsit file.
*   indx_save   - Saved value of index to see if + at end of line.
 
      integer*4 crt, DecimalToInt, i, ierr, indx, iprm(5), len_run,
     .    loglu, rcpar, trimlen, off_com, jerr, kerr, dumm, 
     .    num_in_gdl, indx_save
     
*   expts_var_read  - Variance to be given to the experiment.

      real*8 expts_var_read
      
*   Still_adding   - Logical to indicate that we are still 
*                    adding global files to the list. (Lines ending
*                    in + are added)
 
      logical still_adding
      
* MOD TAH 980519: Added explicit specification of diagonal 
*     scaling of matrices

* glb_diag -- Diagonal scaling factor in ppm
* glb_var  -- Complete matrix scaling

      real*8 glb_diag, glb_var
 
*   crt_string  - String containing CRT LU
*   log_string  - String containing LOG LU
*   prt_string  - String containing PRINT LU
 
      character*128 crt_string, log_string, prt_string

*   comopt      - Optional command line beginning string

      character*256 comopt
 
*   global_file - Name of the global file being processed
*   input_file  - Name of the file with the list of global files
*               - to be processed.
*   list_file   - Name of the list file for current global
*   list_mask   - Mask to be used to generate list file name
*   markov_file - Name of the markov file (Must be given)
*   common_file - Name of the global common file to used
*   sort_file   - Name of the sort file to be passed to GLOBK
 
      character*128 global_file, input_file, list_file, list_mask,
     .    markov_file, common_file, sort_file
 
*   line        - Line read from input file
 
      character*128 line
 
*   globk_run   - GLOBK run command line
 
      character*256 globk_run
      character*4 cdum
 
***** Start decoding the runstring
 
      crt_string = ' '
      prt_string = ' '
      log_string = ' '
      list_file  = ' '
      markov_file = ' '
      common_file = ' '
 
      sort_file   = ' '
 
*                                            ! Get CRT string for GLOBK
      len_run = rcpar(1, crt_string )
      if( len_run.gt.0 ) then
          crt = DecimalToInt( crt_string, ierr)
      end if
      if( len_run.eq.0 .or. ierr.ne.0 ) then
          crt = loglu(i)
      end if
 
*                                            ! Printer string
      len_run = rcpar(2, prt_string)
*                                            ! Log LU string
      len_run = rcpar(3, log_string)
*                                            ! Name of input file
      len_run = rcpar(4, input_file)
      if( len_run.eq.0 ) then
          call proper_runstring('glred.hlp','glred',1)
*                                            ! Report runstring and stop
      end if
 
*                                            ! Name of markov file
      len_run = rcpar(5, markov_file)
      if( len_run.eq.0 ) then
          call proper_runstring('glred.hlp','glred',1)
*                                            ! Report runstring and stop
      end if
 
*
      len_run = rcpar(6, comopt)             ! optional command line beginning
      if( len_run.eq.0 ) comopt = ' ' 
                                             ! Name of common file (optional)
      len_run = rcpar(7, common_file)
*                                            ! Name of sort file (optional)
      len_run = rcpar(8, sort_file)
 
*                                            ! List file mask (optional)
      len_run = rcpar(9, list_mask)
*                                            ! Use default
      if( len_run.eq.0 ) then
          list_mask = list_mask_default
      end if
 
***** Now loop over the input file, scheduling GLOBK to run on each of the
*     files
 
      open(100, file=input_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',input_file,0,'GLRED')
      if( ierr.ne.0 ) then
          call proper_runstring('glred.hlp','glred',1)
*                                            ! Report runstring and stop
      end if
 
***** Now loop over of the input file
 
      do while ( ierr.eq.0 )

          still_adding = .true.
          num_in_gdl = 0
          do while ( still_adding )
 
              read(100,'(a)',iostat=ierr) line
              jerr = ierr
              if( ierr.ne.0 ) still_adding = .false.
*                                             ! Get file name from
* MOD TAH 950106: Check file name to see if non-blank and does not
*             start with # or *
              if ( ierr.eq.0 .and. trimlen(line).gt.0 .and.
     .             line(1:1).ne.'#' .and. line(1:1).ne.'*' ) then
*                                                 ! line
                  indx = 1
                  call read_line( line, indx, 'CH', jerr, dumm,
     .                            global_file)

*                 Try to read the variance:
                  indx_save = indx
                  glb_var  = 1.d0
                  glb_diag = 0.d0
                  
                  call GetWord(line, cdum, indx)
                  if ( cdum(1:4).ne.'+   ' ) then
                      indx = indx_save
                      call read_line( line, indx, 'R8', kerr, 
     .                            glb_var, cdum )
                      if( kerr.ne.0 ) then
                          glb_var = 1.d0
                          if( index(line,'+').lt.indx_save ) 
     .                                  still_adding = .false.
                      else

* MOD TAH 980519:         see if diagonal passed
                          indx_save = indx
                          call GetWord(line, cdum, indx)
                          if ( cdum(1:4).ne.'+   ' ) then
                             indx = indx_save
                             call read_line( line, indx, 'R8', kerr, 
     .                                  glb_diag, cdum )
                             if( kerr.ne.0 ) then
                                glb_diag = 0.d0
                                if( index(line,'+').lt.indx_save ) 
     .                                    still_adding = .false.
                             else
* MOD TAH 980519:               In this case test against last thing
*                               read in line.     
                                if( index(line,'+').lt.indx ) 
     .                                    still_adding = .false.
                            endif
                         end if
                     end if
                  endif

*                 Compute the value of the variance scale to be written.
                  if( glb_diag.ne.0.d0 ) then
                      expts_var_read = -(glb_var + 
     .                              (1.d0+glb_diag/1.d6)/1000.d3)
                  else
                      expts_var_read = glb_var
                  end if

              else 
                  jerr = -1
              end if    
  
*                                     ! Schedule Globk
              if( jerr.eq.0 ) then
 
*                 Generate list file name
                  if( num_in_gdl.eq.0 ) then
                      list_file = list_mask
                      call wild_card( list_file, global_file )
 
*                     Create the list file
                      open(200,file=list_file, iostat=jerr, 
     .                         status='unknown')
                      call report_error('IOSTAT',jerr,'open',list_file,
     .                                  0,'GLRED')
                  end if
 
              end if
 
*                                 ! Only continue is no errors
              if( jerr.eq.0 ) then
 
*                 Write the global file name into the list file
                  num_in_gdl = num_in_gdl + 1
                  write(200,'(a,1x,f25.16)', iostat=ierr) 
     .                 global_file(1:trimlen(global_file)), 
     .                 expts_var_read
              end if
          end do
              
          close(200)
 
*             Now schedule GLOBK, build up runstring
          if( num_in_gdl.gt. 0 ) then
              globk_run = 'globk ' //
     .            crt_string (1:max(1,trimlen(crt_string  ))) // ' ' //
     .            prt_string (1:max(1,trimlen(prt_string  ))) // ' ' //
     .            log_string (1:max(1,trimlen(log_string  ))) // ' ' //
     .            list_file  (1:max(1,trimlen(list_file   ))) // ' ' //
     .            markov_file(1:max(1,trimlen(markov_file ))) // ' ' //
     .            comopt     (1:max(1,trimlen(comopt      ))) // ' ' //
     .            common_file(1:max(1,trimlen(common_file ))) // ' ' //
     .            sort_file  (1:max(1,trimlen(sort_file   )))

              write(*,'(a)') globk_run(1:trimlen(globk_run))
 
              call execute( globk_run, iprm, 1, 100, off_com)
 
*             Now purge the list file since we do not need it
              open(200,file=list_file, iostat=ierr, status='old')
              close(200, status='delete', iostat= ierr)
              call report_error('IOSTAT',ierr,'clos',list_file,0,
     .                          'GLRED')
          end if
 
      end do
 
***** Thats all
      close(100)
      end
