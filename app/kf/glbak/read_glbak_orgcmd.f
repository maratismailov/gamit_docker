CTITLE READ_GLBAK_ORCMD
 
      subroutine read_glbak_orgcmd(iout)  
 
      implicit none
 
*     Routine to read the commands from the glorg command file
*     for use in GLOBAK.  Only the frame stabilization commands
*     are read.  com_orgopt is used for start_line options.
*
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
 
*   i,j,k       - Loop counters
*   ierr        - IOSTAT errors
*   iel         - command number
*   indx        - index in READLINE
*   jerr        - Error return from decode_param
 
 
      integer*4 i, j, ierr, iel, indx, jerr,
     .          trimlen

      real*8 val_save, var_save   ! Force Value and variance read from line

      logical kbit, process_line

* MOD TAH 980517: Added reading multiple command files
* push_unit  -- Unit number of pushed file
* curr_unit  -- Current number

       integer*4 push_unit, curr_unit
       integer*4 def_conts(max_glb_site_wrds)

* new_cmd_file -- Name of new command file

       character*256 new_cmd_file
 
*   buffer      - Line read from various files
*   name        - Name for plate or site
 
      character*150 buffer
      character*8 name
      character*4 cdum
 
*   values(6)   - Position and rate read from apriori file
*   sol_parm(1) - Solution vector
 
      real*8 values(6), sol_parm(num_glb_parn)

* Variables for checking aprioris on equates
*  ref  - Reference site number
*  apr  - Apriori value

      integer*4 ref
      real*8 apr

      logical updated  ! For com_orgopt line

      integer*4 iout   ! Output unit
 
***** Start open the command file.  Check to see if command file.
      if( trimlen(glorg_command_file).eq.0 ) RETURN
      push_unit = 0
      curr_unit = 99
      open(curr_unit, file=glorg_command_file, iostat=ierr, 
     .       status='old' )
      call report_error('IOSTAT',ierr,'open',glorg_command_file,
     .    1, 'read_glbak_orgcmd')

      glr_cmd_file = glorg_command_file

*     Firstly clear the adjustments to the parameters
      do i = 1, num_glb_parn
          parm_change(i) = 0
      end do
      do i = 1, gnum_sites
          call sbit( cov_sites, i,1 )
      end do
      do i = 1, max_glb_site_wrds
         def_conts(i) = -1
      end do 

      cnd_pos_bits = 0
      cnd_rat_bits = 0
      cnd_hgt_var(1) = 10.d0
      cnd_hgt_var(2) = 10.d0

      use_ratio(1) = 3.d0
      use_ratio(2) = 3.d0

*     Set the default number of stabalization iterations to
*     be 1 (current number used to glorg).
      num_stab_iter = 4
      stab_nsig  = 4.d0
*     Put the default half way between constant and site 
*     dependant weighting.
      stab_rel   = 0.5d0

*     Set the minium RMS and DH sif diffeenrces
      stab_min_rms(1)  = 0.003
      stab_min_rms(2)  = 0.003
      stab_min_dh(1)   = 0.005
      stab_min_dh(2)   = 0.005
      stab_min_dne(1)  = 0.0005
      stab_min_dne(2)  = 0.0001 

*     Set default to estimation translation when plate euler poles
*     are estimated

      do i = 1,7
         cond_var(i,1) = 0.d0
         cond_var(i,2) = 0.d0
      end do

****  Scan for sites that have zero epoch and non-zero apriori velocity.
*     (Site not in the apriori files will often have zero epoch and this
*     is not a problem if the velocity is zero).
      do i = 1, gnum_sites
         if( site_epoch(i).eq.0 .and. kbit(guse_site,i) ) then
*            If apriori velocity is zero, then changing the epoch is OK  
             if ( apr_val_site(1,2,i).eq.0 .and. 
     .            apr_val_site(2,2,i).eq.0 .and. 
     .            apr_val_site(3,2,i).eq.0      ) then
                 site_epoch(i) = gepoch_out
             else
                 write(*,210) gsite_names(i), (apr_val_site(j,2,i),
     .                    j=1,3)
                 if( iout.ne.6 )
     .           write(iout,210) gsite_names(i), (apr_val_site(j,2,i),
     .                        j=1,3)
 210             format('** WARNING ** Site ',a8,' has zero reference ',
     .               'epoch, Apriori Velocity ',3F6.3,' m/yr, Setting',
     .               ' epoch to solution epoch')
                 site_epoch(i) = gepoch_out
             end if
          end if
      end do

 
*     Now loop over the command file
 
      ierr = 0
      jerr = 0

      do while ( ierr.eq.0 )
 
          read(curr_unit,'(a)', iostat=ierr ) buffer

* MOD TAH 980517: See if we can pop back to a previous unit
*         read
          if( ierr.ne.0 .and. push_unit.ne. 0 ) then
              close(curr_unit)
              curr_unit = push_unit
              push_unit = 0
              read(curr_unit,'(a)', iostat=ierr ) buffer
          end if

****      See if we should process or not
          process_line = .false.
          indx = 1

          call decode_comopt( buffer, com_orgopt, updated)

          if( buffer(1:1).eq.' ' .and. trimlen(buffer).gt.0 ) 
     .                                    process_line = .true.
          if( ierr.ne.0 ) process_line = .false.

          if( process_line ) then
              call get_cmd( buffer, glorg_commands,
     .            num_glorg_commands, iel, indx ) 
 
*****         Skip most of the glorg commands and process only those
*             associated with frame realization.
 
*             ! USE_SITES or STAB_SIT command
* MOD TAH 190526: Save the stab_sites in cov_sites array so that it does
*             mot over the use_site array (needed for backwards solution)
              if( iel.eq.5 .or. iel.eq.23 ) then
                  call decode_stab_opt( buffer(indx+1:), gsite_names, 
     .                gnum_sites, cov_sites, def_conts,  list_file, 
     .                gepoch_start, gepoch_end, guse_site)
 
              end if
 

*             ! POS_ORG: Sets translation/rotation/scale to used/
              if( iel.eq.10 ) then
                  call casefold(buffer)
                  cnd_pos_bits = 0
                  call decode_option(buffer, org_types, 7, 
     .                               cnd_pos_bits, -1)
              end if  

*             ! RATE_ORG: Unlikely to be used in back solution.
              if( iel.eq.11 ) then
                  call casefold(buffer)
                  cnd_rat_bits = 0
                  call decode_option(buffer, org_types, 7,
     .                               cnd_rat_bits, -1)
              end if  

*                                     ! CND_HGTV
              if( iel.eq.12 ) then
                  call read_line(buffer, indx, 'R8',jerr,
     .                           cnd_hgt_var(1), cdum)
                  call report_error('IOSTAT',jerr,'decod',
     .                           buffer,0,'READ_GLORG_COMMAND') 
                  if( jerr.ne.0 ) cnd_hgt_var(1) = 1.d0
                  call read_line(buffer, indx, 'R8',jerr,
     .                           cnd_hgt_var(2), cdum)
                  if( jerr.ne.0 ) cnd_hgt_var(2) = 
     .                                         cnd_hgt_var(1)

*                 See if the ratio values have been given as well.
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values(1), cdum)
                  if( jerr.eq.0 ) use_ratio(1) = values(1)
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values(1), cdum)
                  if( jerr.eq.0 ) use_ratio(2) = values(1)
              end if


*             ! COND_SIG: Sigmas on refeence frame realization,
              if( iel.eq.17 ) then
*                 Get Rotation sigma
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values(1), cdum)
                  do i = 1,7
                     cond_var(i,1) = values(1)**2
                     cond_var(i,2) = values(1)**2
                  end do
*                 See if translation sigma passed
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values(1), cdum)
                  if( jerr.eq.0 ) then
                      do i = 4,6
                         cond_var(i,1) = values(1)**2
                         cond_var(i,2) = values(1)**2
                      end do
                  end if
*                 See if Scale sigma given
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values(1), cdum)
                  if( jerr.eq.0 ) then
                      cond_var(7,1) = values(1)**2
                      cond_var(7,2) = values(1)**2
                  end if
         
*                 Get rotation rate sigma
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values(1), cdum)
                  do i = 1,7
                     cond_var(i,2) = values(1)**2
                  end do
*                 See if translation rate sigma passed
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values(1), cdum)
                  if( jerr.eq.0 ) then
                      do i = 4,6
                         cond_var(i,2) = values(1)**2
                      end do
                  end if
*                 See if Scale rate sigma given
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values(1), cdum)
                  if( jerr.eq.0 ) then
                      cond_var(7,2) = values(1)**2
                  end if
              end if

*                                    ! STAB_ITE
*             Get the number of iterations to use in the stabalization
*             and the 'nsig' editing condition for removing sites.
              if( iel.eq.18 ) then
                  call read_line(buffer, indx, 'I4',jerr,
     .                           num_stab_iter, cdum)
                  call read_line(buffer, indx, 'R8',jerr,
     .                           stab_rel , cdum)
                  if( stab_rel.lt.0 ) stab_rel = 0.d0
                  if( stab_rel.gt.1 ) stab_rel = 1.d0
                  call read_line(buffer, indx, 'R8',jerr,
     .                           stab_nsig, cdum)
              end if

*                                   ! STAB_MIN
*             Sets the minium values for dh sig difference and
*             rms fit (Read in order dh_sig [pos], RMS [pos],
*             dh_sig [Rate], RMS [rate]
              if( iel.eq.19 ) then
                  call read_line(buffer, indx, 'R8',jerr,
     .                           stab_min_dh(1) , cdum)
                  call read_line(buffer, indx, 'R8',jerr,
     .                           stab_min_rms(1) , cdum)
                  call read_line(buffer, indx, 'R8',jerr,
     .                           stab_min_dh(2) , cdum)
                  if( jerr.ne.0 ) stab_min_dh(2) = stab_min_dh(1)
                  call read_line(buffer, indx, 'R8',jerr,
     .                           stab_min_rms(2) , cdum)
                  if( jerr.ne.0 ) stab_min_rms(2) = stab_min_rms(1)
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values , cdum)
                  if( jerr.eq.0 ) stab_min_dne(1) = values(1)
                  call read_line(buffer, indx, 'R8',jerr,
     .                           values , cdum)
                  if( jerr.eq.0 ) stab_min_dne(2) = values(1)
                  
              end if
 
*                                   ! SOURCE
*             Get the name of new command file to source to (cannnot
*             contain another source command in current version)
              if( iel.eq.20 ) then
                  call read_line(buffer,indx,'CH', jerr, values,
     .                new_cmd_file)
                  call wild_card(new_cmd_file, list_file)
                  push_unit = curr_unit
                  curr_unit = 100
                  open(curr_unit, file=new_cmd_file, iostat=ierr, 
     .                 status='old' )
                  call report_error('IOSTAT',ierr,'open',new_cmd_file,
     .                 0, 'Source command file')

*                 Go back to original unit if error opening file
                  if( ierr.ne.0 ) then
                      ierr = 0
                      curr_unit = push_unit
                      push_unit = 0
                  end if
              end if

 
          end if
 
      end do

***** Thats all
      close(curr_unit)
      return
      end
