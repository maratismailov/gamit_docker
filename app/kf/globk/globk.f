      program globk
 
      implicit none 
 
*     This program is the main part of the global Kalman filter
*     software.
*     The operation of this program is to read the command file,
*     compute the parameter list and the state transission pointers,
*     and to schedule the other programs in this suite:
*     GLINIT: Sorts the experiments into time order, and get the list
*             of site and source names and positions.
*     GLFOR : Read the sorted global files, form the partials, and does
*             the forward solution.
*     GLBAK : Similar to GLFOR, but does the back solution as well.
*     GLOUT : Outputs the solution to the print device.
*
*     The runstring for GLOBK is:
*     CI> GLOBK,<crt>,<prt>,<log>,list,markov,<commom>,<sort>,<direct>
*
*     where crt  is users terminal (optional, default 1)
*           prt  is print device   (optional, default crt)
*           log  is log device (one line per input experiment) (optional,
*                   default no log output)
*           list is name of file containing the list of experiments to
*                   be processed.  One experiment name per line.
*           markov is the name of the command file containing the GLOBK
*                   commands
*           common is the name of the common file to be used (default
*                   GLBCOM::53)
*           sort is the name of the file with the time sorted list of
*                   experiments (Default GLBSOR::53)
*           direct is the direction of the sorting, +1 for forward, -1
*                   for backward. (Default +1).  The back solution will
*                   come out in the opposite order.
*
*
*                                     09:42 AM MON., 3 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/globk_cmds.h'
c      external globk_cmd_bd
 
*   trimlen     - Gets length of string
*   jerr, fmppurge - Error return and function to delete files
*     at the end of runs.

      integer*4 trimlen, jerr, fmppurge

*   crt_file    - String for name of unit 6

      character*4 crt_file
      
*   kbit        - Logical to test bit status
      logical kbit 
 
***** Start, Decode the runstring. If incomplete then the program will
*     stop

      call globk_header( 6, globk_version)

      glb_com_file = ' '
      call decode_globk_run
 
 
*     Now read the command file.  Note: GLINIT will be run when a
*     command after those which allow input of the runstring options
*     is executed.
      call report_stat('status','globk','main', list_file,
     .                 'Start: Processing list file',0)
      call read_glb_markov

*     Update the apriori's as given the in the apriori files
      call get_apr_positions

*     Now get the parameter list.  This assigns the parameters to be
*     estimated to specific parameter numbers.  It also forms the
*     state transission pointers ie. pointers to which parameters have
*     non unit transission elements
 
      call glb_param_list
 
*     Now run the forward global filter

      call report_stat('status','globk','main',' ',
     .                 'Running forward Kalman filter',0)
      call run_glfor
 
*     If we are doing a back solution then run now
 
      if( trimlen(glb_out_file).gt.0 ) then
          call report_stat('status','globk','main',glb_out_file,
     .                     'Saving combined global file',0)
          call run_glsave
      end if
 
      if ( glb_bak_soln ) then
          call report_stat('status','globk','main',glb_bak_file,
     .                     'Running back Kalman filter',0)
          call run_glbak
      end if
 
*     Now output the solution to each output device
      crt_file = '6'  

*     See if user has not suppressed the output to crt  
      if( .not.kbit(crt_opts,18) ) then
          call report_stat('status','globk','main',' ',
     .                 'Outputting filter run to screen',0)
          call run_glout( crt_file, crt_opts )
      end if

*     Check also if output suppressed.     
      if( glb_prt_file(1:1).ne.'6' .and.
     .    .not.kbit(prt_opts,18)) then
          call report_stat('status','globk','main',glb_prt_file,
     .                 'Outputting filter run to file',0)
          call run_glout( glb_prt_file, prt_opts )
      end if

*     See if glorg is to be run.
      if( trimlen(glr_cmd_file).gt.0 ) then

*         See if we have output file name.  If not then generate it
*         from the print file name
          if( trimlen( glb_org_file ).eq. 0 ) then
              call gen_org_name( glb_org_file, glb_prt_file, 'org' )
          end if
          call report_stat('status','globk','main',glb_org_file,
     .                 'Running glorg to output file',0)
          call run_glorg( glr_cmd_file, glb_org_file, org_opts )
          write(*,*) 'GLORG_FINISHED'
      end if

* MOD TAH 050524: See if scracth files are to be deleted
      if( del_scratch ) then
         print *,'DELETING COM, SOL, SRT and MAKE_SVS FILES'
         jerr = 0
         if( trimlen(glb_com_file).gt.0 ) jerr = FmpPurge(glb_com_file)
         call report_error('IOSTAT',jerr,'delet',glb_com_file,0,'globk') 
         if( trimlen(glb_sol_file).gt.0 )jerr = FmpPurge(glb_sol_file)
         call report_error('IOSTAT',jerr,'delet',glb_sol_file,0,'globk') 
         if( trimlen(sort_file).gt.0 )jerr = FmpPurge(sort_file)
         call report_error('IOSTAT',jerr,'delet',sort_file,0,'globk') 
         if( make_svs_file ) then
            jerr = FmpPurge(glb_svs_file)
            call report_error('IOSTAT',jerr,'delet',glb_svs_file,
     .                      0,'globk') 
         endif
      endif

 
***** Thats all.  Now interpret the results.
      end
 
*     Include the block data in the main program since some compilers/linkers (specifically OSX) will not correctly 
*     load the initialized block data routines from library archives.
      include 'globk_cmd_bd.f' 
 
