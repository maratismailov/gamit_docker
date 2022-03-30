CTITLE READ_GLB_MARKOV
      subroutine read_glb_markov

      implicit none
 
 
*     Routine to read the markov control file for GLOBK.  See
*     GLOBK::HELP for list of the commands and there arguments.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/globk_cmds.h'
 
*   iel         - Command number from get_cmd
*   ierr        - IOSTAT error
*   indx        - Pointer to current position in the command
*               - buffer
 
      integer*4 iel, ierr, indx, i, jerr

*   values -- Dummy entry for readline

      real*8 values

* MOD TAH 980517: Added source command to allow a second command
*     file to be read

*   push_unit -- Pushed unit number
*   curr_unit -- Current unit number
*   len_comopt -- length of comopt string

      integer*4 push_unit, curr_unit, len_comopt, trimlen

*   new_cmd_file -- Name of new command file
      character*256 new_cmd_file
 
*   glinit_run  - Indicates that GLINIT has been run
*   process_line - Set true if current line from file is to be
*     processed (either starts with blank or comopt string).
 
      logical glinit_run, process_line, updated 
 
*   buffer      - Line read from markov file
 
      character*256 buffer
 
***** Say that glinit needs to be run
      glinit_run = .false.
 
***** Open the command file
      push_unit = 0
      curr_unit = 99
      num_eqfiles = 0
 
      open(curr_unit, file=glb_mar_file, iostat=ierr, status='old')
*                                                            ! Kill
      call report_error('IOSTAT',ierr,'open',glb_mar_file,1,
     .                  'READ_GLB_MARKOV')
 
***** Now read the file, and decode the commands until we reach EOF
 
      ierr = 0
      make_svs_file = .false.
* MOD TAH 180402: Set defualt to use PRN names for GPS rather than
*     new GNSS names (will be changed later)
      use_prnn = .true.

      len_comopt = trimlen(comopt)
      decnum = 1
      decoff = 1

      do while ( ierr.eq.0 )
 
          read(curr_unit,'(a)', iostat=ierr) buffer

* MOD TAH 980517: If there is an error see if we can pop the unit
*         stack
          if( ierr.ne.0 .and. push_unit.ne. 0  ) then
              close(curr_unit)
              curr_unit = push_unit
              push_unit = 0
              read(curr_unit,'(a)', iostat=ierr) buffer
          endif
              
*                                 ! Check error
          if( ierr.ne.-1) then
              call report_error('IOSTAT',ierr,'read',glb_mar_file,
     .                          0,'READ_GLB_MARKOV')
          end if
 
*         Process command if the error was zero
          if( ierr.eq.0 ) then
 
              indx = 1
 
*             Check for comment (nonblank in firt charcater)
*                                             ! Process
              process_line = .false.
              if( len_comopt.gt.0 ) then
                  call decode_comopt(buffer, comopt, updated) 
              end if 
              if( buffer(1:1).eq.' ' ) process_line = .true.
               
              if( process_line ) then
 
                  call get_cmd( buffer, glb_commands,
     .                max_glb_commands, iel, indx )
 
*                                             ! Command found
* MOD TAH 980517: See if SOURCE command used.  If so handle locally
* MOD TAH 180616: Increased command number of SOURCE command after the additon
*                 of the USE_PRNN command.  Change from 78 to 79.
                  if( iel.eq. 79 ) then
* MOD TAH 070105: Added fatal is source used in source file
                      if( push_unit.ne.0 ) then
                         call report_stat('FATAL','globk',
     .                         'read_glb_mar',buffer,
     .                         'SOURCE COMMAND in SOURCE file',0)
                      end if

                      call read_line(buffer,indx,'CH', jerr, values,
     .                    new_cmd_file)
                      call wild_card( new_cmd_file, list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
                      call wild_date( new_cmd_file, wild_mjd, 'N' )

                      push_unit = curr_unit
                      curr_unit = 100
                      open(curr_unit, file=new_cmd_file, iostat=ierr,
     .                     status='old' )
                      call report_error('IOSTAT',ierr,'open',
     .                     new_cmd_file, 0, 'Source command file')
                      if( ierr.ne.0 ) then
                          ierr = 0
                          curr_unit = push_unit
                          push_unit = 0
                      end if
                  else if( iel.gt.0 ) then
                      call process_glb_command( buffer, indx, iel,
     .                                          glinit_run )
                  end if
*                                             ! Not a comment
              end if
*                                             ! No file reading error
          end if
*                                             ! Until EOF or error found
      end do
 
***** Make sure that GLINIT was run, If not run it now
 
      if( .not. glinit_run ) then
          call run_glinit

      end if

* MOD TAH 200416: Add the GGVersion information gdescription variable
      call add_GGV( gdescription ) 

*     Do the final check on whether we should a site (based on
*     number of times it is used)
      do i = 1, gnum_sites
         if( times_used(i).le.0 ) call sbit(guse_site,i,0)
* MOD TAH 000901: Check for sites named _XCL, delete these sites
         if( gsite_names(i)(5:8).eq.'_XCL' ) call sbit(guse_site,i,0)

      end do
 
****  Thats all, close command file and start run
      close(curr_unit)
      return
      end
 
