CTITLE READ_GLORG_COMMANDS
 
      subroutine read_glorg_commands ( sol_parm, options )
 
      implicit none
 
*     Routine to read the commands for glorg.  This routine also
*     updates the arprioi values if these are to be changed.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
 
*   i,j,k       - Loop counters
*   ierr        - IOSTAT errors
*   iel         - command number
*   indx        - index in READLINE
* MOD TAH 150130: Made param_num an array to allow multiple sites to be 
*   returned.  Add num_pn for number in array
*   param_num   - Parameter number(s) decoded with decode_param
*   jerr        - Error return from decode_param
*   iwc         - Position of wild card character
 
 
      integer*4 i, j, ierr, iel, indx, jerr,
     .          jel, kel, jndx, trimlen, options, iwc

      integer*4 param_num(max_eq), num_pn

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

      integer*4 iout
      common / progcon / iout
 
***** Start open the command file
      push_unit = 0
      curr_unit = 99
      open(curr_unit, file=glorg_command_file, iostat=ierr, 
     .       status='old' )
      call report_error('IOSTAT',ierr,'open',glorg_command_file,
     .    1, 'Read_glorg_commands')

      glr_cmd_file = glorg_command_file

*     Firstly clear the adjustments to the parameters
      do i = 1, num_glb_parn
          parm_change(i) = 0
      end do
      do i = 1, gnum_sites
          call sbit( use_sites, i,1 )
          call sbit( cov_sites, i,1 )
          plate_number(i) = 0
          assign_number(i) = 0
      end do
      do i = 1, max_glb_site_wrds
         def_conts(i) = -1
      end do 
      num_equates = 0
      equate_loc  = .true.  
      num_force   = 0
      first_eqf   = .false.
      num_plates = 0

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

      glr_sol_file = ' '

*     Set default to estimation translation when plate euler poles
*     are estimated
      PlateTrans = .true.

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
 
*****         Now process the commands
 
*                                     ! Aprioris
              if( iel.eq.2 ) then
                  call read_line(buffer,indx,'CH', jerr, values,
     .                glorg_apr_file)
                  call wild_card( glorg_apr_file, list_file)
                  call update_glorg_aprs( jerr ) 

* MOD TAH 970905: Add the file name to list of names of apr files
                  num_apr_files = num_apr_files + 1
                  if( num_apr_files.le.max_apr_files .and.
     .                jerr.eq.0 ) then
                     glb_apr_file(num_apr_files) = 
     .                  glorg_apr_file(1:trimlen(glorg_apr_file)) //
     .                               ' (glorg)'
                  end if 
              end if
*                                     ! Plate definition
              if( iel.eq.3 ) then
*                 see if we have this plate
                  call casefold(buffer)
                  call getword(buffer, name, indx)
                  jndx = 1
                  call get_cmd( name  , plate_names,
     .                 num_plates, jel, jndx )
                  if( jel.lt.1 .and. name(1:4).ne.'NONE' ) then
                      num_plates = num_plates + 1
                      if( num_plates.gt.max_plates-1 ) then
                          write(*,310) max_plates,
     .                          buffer(1:trimlen(buffer))
 310                      format('**WARNING** Too many plates',
     .                           ' Max allowed is ',i3,/,
     .                           ' IGNORING ENTRY:',a)
                          name = ' '
                      else
                          jel = num_plates
                          plate_names(jel) = name
                      end if
                  elseif ( name(1:4).eq.'NONE' ) then
                      jel = 0
                  end if

****              now pull the names of the sites
                  do while (trimlen(name).gt.0 )
                     call getword(buffer, name, indx)
                     jndx = 1
                     call get_cmd(name, gsite_names,
     .                    gnum_sites, kel, jndx)
* MOD TAH 040327: Added ALL option and allowed wild cards
                     if( kel.gt.0 .and.kel.ne. 999999) then
                         if( parn_site(1,2,kel).gt.0 ) then
                             plate_number(kel) = jel
                         end if
                     elseif ( kel.eq.999999 ) then   ! ALL option
                         do j = 1, gnum_sites
                            if( parn_site(1,2,j).gt.0 ) then
                               plate_number(j) = jel
                            end if
                         end do
                     else    ! See if we can match on wild card
                         iwc = 0
                         iwc = index( name,'*')
                         if( iwc.eq.0 ) iwc = index(name,'@')
                         if( iwc.gt.0 ) then
*                           OK: Found wild card, match all names up
*                           to this point
                            call casefold(name)
                            do j = 1, gnum_sites
                               if( name(1:iwc-1).eq.
     .                            gsite_names(j)(1:iwc-1) .and.
     .                            parn_site(1,2,j).gt.0 ) then
                                  plate_number(j) = jel
                               endif
                            end do
                         end if
                     end if
* END WILD CARD MOD
                  end do
              END IF

* MOD TAH 980612: Introdced Assign_p command to assign site to
*             plate.
*                                     ! Plate assignment (ASSIGN_P)
              if( iel.eq.21 ) then
*                 see if we have this plate
                  call casefold(buffer)
                  call getword(buffer, name, indx)
                  jndx = 1
                  call get_cmd( name  , plate_names,
     .                 num_plates, jel, jndx )
                  if( jel.lt.1 .and. name(1:4).ne.'NONE' ) then
                      num_plates = num_plates + 1
                      if( num_plates.gt.max_plates ) then
                          write(*,310) max_plates,
     .                          buffer(1:trimlen(buffer))
                          name = ' '
                      else
                          jel = num_plates
                          plate_names(jel) = name
                      end if
                  elseif ( name(1:4).eq.'NONE') then
                     jel = 0
                  end if

****              now pull the names of the sites
                  do while (trimlen(name).gt.0 )
                     call getword(buffer, name, indx)
                     jndx = 1
                     call get_cmd(name, gsite_names,
     .                    gnum_sites, kel, jndx)
* MOD TAH 040327: Added ALL option and allowed wild cards
                     if( kel.gt.0 .and.kel.ne. 999999) then
                         if( parn_site(1,2,kel).gt.0 ) then
                             assign_number(kel) = jel
                         end if
                     elseif ( kel.eq.999999 ) then   ! ALL option
                         do j = 1, gnum_sites
                            if( parn_site(1,2,j).gt.0 ) then
                               assign_number(j) = jel
                            end if
                         end do
                     else    ! See if we can match on wild card
                         iwc = 0
                         iwc = index( name,'*')
                         if( iwc.eq.0 ) iwc = index(name,'@')
                         if( iwc.gt.0 ) then
*                           OK: Found wild card, match all names up
*                           to this point
                            call casefold(name)
                            do j = 1, gnum_sites
                               if( name(1:iwc-1).eq.
     .                            gsite_names(j)(1:iwc-1) .and.
     .                            parn_site(1,2,j).gt.0 ) then
                                  assign_number(j) = jel
                               endif
                            end do
                         end if
                     end if
* END WILD CARD MOD
                   end do
              END IF
*                                     ! USE_SITES or STAB_SIT command
              if( iel.eq.5 .or. iel.eq.23 ) then
                  call decode_stab_opt( buffer(indx+1:), gsite_names, 
     .                gnum_sites, use_sites, def_conts,  list_file, 
     .                gepoch_start, gepoch_end, guse_site)
 
              end if
 
*                                     ! COV_SITES
              if( iel.eq.6) then
                  call casefold(buffer)
                  call decode_option( buffer, gsite_names, gnum_sites,
     .                cov_sites, -1)
 
              end if
 
*                                     ! EQUATE
              if( iel.eq.7 ) then
                  call decode_equates( buffer, indx, num_equates,
     .                 num_each_equate, param_equates, eq_var, 
     .                 max_equates, max_each_equate, 'EQUA', options )
              end if

*                                     ! LOCAL_EQ
              if( iel.eq.8 ) then
                 call GetWord( buffer, cdum, indx)
                 call casefold(cdum)
                 if( cdum(1:1).eq.'N' ) then
                     equate_loc = .false.
                 else
                     equate_loc = .true.
                 end if
              end if 

*                                     ! FORCE 
              if( iel.eq.9 ) then
                  ref = 0
                  call decode_param( buffer, indx, param_num, num_pn, 
     .                               jerr, ref, apr, options) 
                  if( jerr.eq.0 .and. num_pn.gt.0 ) then
                      do i = 1, num_pn 
                         num_force = num_force + 1
                         param_force(num_force) = param_num(i) 
* MOD TAH 150130: When multiple sites; only read the force value once 
                         if( i.eq.1 )  call read_line(buffer,indx,'R8', 
     .                                      jerr, val_save, cdum)
                         val_force(num_force) = val_save

* MOD TAH 000302:        See if the sigma of the force has been passed
* MOD TAH 150130:        Only read for first entry
                         if( i.eq.1 ) call read_line(buffer,indx,'R8', 
     .                                     jerr, var_save, cdum)
                         if( jerr.eq.0 ) then
                            var_force(num_force) = var_save**2 
                         else
                             var_force(num_force) = 0.d0
                         end if 
                      enddo                    
                  end if
              end if

*                                     ! POS_ORG
              if( iel.eq.10 ) then
                  call casefold(buffer)
                  cnd_pos_bits = 0
                  call decode_option(buffer, org_types, 7, 
     .                               cnd_pos_bits, -1)
              end if  

*                                     ! RATE_ORG
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

*                                     ! EQ_DIST
              if( iel.eq.13 ) then
                  call decode_eq_dist( buffer, indx, options )
              end if

*                                     ! UNEQUATE
              if( iel.eq.14 ) then
                  call decode_unequates( buffer, indx, num_equates,
     .                 num_each_equate, param_equates,
     .                 max_equates, max_each_equate )
              end if

*                                     ! FIRST_EQF
              if( iel.eq.15 ) then
                 first_eqf  = .true.
              end if 

*                                     ! CONSTRAI
              if( iel.eq.16) then
                  call decode_equates( buffer, indx, num_equates,
     .                 num_each_equate, param_equates, eq_var,
     .                 max_equates, max_each_equate, 'CONS', options )
              end if

*                                     ! COND_SIG
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

*                                     ! OUT_SOL
*             Command to allow glorg to write out a new solution 
*             file so that it can be saved with glsave. 
              if( iel.eq.22 ) then
                  call read_line(buffer,indx,'CH', jerr, values,
     .                 glr_sol_file)
                  cdum = glr_sol_file
*                 If none passed as name; remove now
                  call casefold(cdum)
                  if( cdum(1:4).eq. 'NONE' ) glr_sol_file = ' '
                  call wild_card(glr_sol_file, list_file)
                  
              end if

*             NOTES: Command 23 (stab_sit is the same as use_site
*             and therefore is no entry here:
*             if( iel.eq.23 ) -- Treated with iel=5.

*             See if we are turing off plate translation
              if( iel.eq.24 ) then
                  PlateTrans = .false.
              end if

* MOD TAH 150528: Check for EQ_4CHAR command
              if( iel.eq.25 ) then 
                  call decode_eq_4char( buffer, indx, options )
              endif 
  
          end if
 
      end do

*     Now that we have finished reading all the commands, update
      call glorg_upd_apr ( sol_parm ) 

***** Thats all
      close(curr_unit)
      return
      end

CTITLE glorg_upd_apr

      subroutine glorg_upd_apr ( sol_parm )

      implicit none

*     Routine to apply the corrections for the changes to the
*     apriori values.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'

      real*8 sol_parm(num_glb_parn)

      integer*4 i
 
***** Update the solution vector for any changes which have been made
      do i = 1, num_glb_parn
          sol_parm(i) = sol_parm(i) - parm_change(i)
      end do

      return
      end

CTITLE decode_equates
 
      subroutine decode_equates( buffer, indx, num, each, params,
     .                           var, max_eq, max_each, type, options )

      implicit none

*     Routine to decode the parameter numbers to equated in this
*     solution.

*   indx    - pointer to position to string
*   num     - Number of equates upto this point
*   max_eq  - Maxiumum number of equates allowed
*   max_each    - Maximum number allowed for each each equate
*   each(max_eq) - Number of equates for of all equates
*   params(max_each, max_eq)    - Parameter numbers to be equated.
*   options - Option that set if we should update the aprioris
 
      integer*4 indx, num, max_eq, max_each, each(max_eq),
     .    params(max_each, max_eq), options

*   var     - Variance to be assigned to equate

      real*8 var(max_eq)
 
*   buffer  - Line read from command file.
*   type    - Type of constraint:  EQUA - Equate variance = 0
*                                  CONS - Constraint (sigma given)
 
      character*(*) buffer, type
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error (look for -1 to tell end of list)
*   nt, np  - Temporary pointers to number of parameters read
*   ref     - Reference site number
 
      integer*4 ierr, nt, ref
* MOD TAH 150130: Made np and array with nnp values
      integer*4 np(max_eq), nnp, i

      real*8 sigma, apr
 
*   cdum    - Dummy string
 
      character*4 cdum

*     Start looping over the input buffer
 
*     Increment number of equates found
      num = num + 1
 
*     See if we have to many
      if( num.gt.max_eq ) then
          write(*,*) ' Too many EQUATE lines given.  Max is ',max_eq
          RETURN
      end if

*     See if constraint (read sigma first)
      if( type(1:2).eq.'CO' ) then
          call read_line(buffer,indx,'R8', ierr, sigma, cdum)
          var(num) = sigma**2
      else
          var(num) = 0.d0
      end if

 
*     initialzie the number of equates for this equate.
      nt = 0
      ierr = 0
      ref = 0
      do while ( nt.lt.max_each .and. ierr.eq.0 )
c         call read_line( buffer, indx, 'I4', ierr, np, cdum)
          call decode_param(buffer,indx, np, nnp, ierr, ref, apr, 
     .                      options)
 
*         If no error reading the value say it
          if( ierr.eq.0 .and. nnp.gt.0 ) then
              do i = 1, nnp
                 nt = nt + 1
                 params(nt,num) = np(i)
              end do 
          end if
 
      end do

*     Now sort the list of equates
      call esort( nt, params(1,num) )
 
*     See if we ended up with any equates from this line.  If we didnot
*     decrement the total count
      if( nt.gt.1 ) then
          each(num) = nt
      else
          num = num - 1
      end if
 
****  Thats all
      return
      end
 
CTITLE UPDATE_GLORG_APRS
 
      subroutine update_glorg_aprs( ierr )
 
      implicit none
 
*     Routine to read an aprori file and update the site coordinates and
*     velocities.  The changes here are added to the aprioris and
*     later will be subtracted from the adjustments.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
 
 
*   i,j,k   - Loop counters
*   ierr,jerr   - IOSTAT errors
*   iel,jel     - Site number
*   indx    - Line index
 
 
      integer*4 i,j, ierr,jerr, iel,jel, indx
 
*   buffer  - Line read from file
 
 
      character*150 buffer, cpbuff
 
*   dt          - Time change for old site epoch
*   dtn         - Time change for new site epoch
*   values(3,2) - Values read from apriori file
*   new_epoch   - Epoch read from aproori file
*   new_jd      - JD corresponding to new_epoch
 
 
      real*8 dt, dtn, values(3,2), new_epoch, new_jd

* Copies needed so that we can update the non-secular parameters
* (Note: there are some problems with this for merged solutions)
* save_num_nonsec -- Original number of non-secular terms
* save_param_nonsec(2,max_nonsec) -- Parameter pointers
      integer*4 save_num_nonsec, save_param_nonsec(2,max_nonsec)
* save_val_nonsec(8,max_nonsec) -- Actual parameters of non-secular 
*     terms
      real*8 save_val_nonsec(8,max_nonsec)

* cval -- String read from buffer
      character*8 cval

****  OK, Before we start, save any non-secular terms already used 
      save_num_nonsec = num_nonsec
      do i = 1, num_nonsec 
         do j = 1,2 
            save_param_nonsec(j,i) = param_nonsec(j,i)
         end do
         do j = 1,8
            save_val_nonsec(j,i) = apr_val_nonsec(j,i)
         end do
      end do

****  Clear the number of non-secular terms so that new values
*     are read
      num_nonsec = 0
 
****  OPen the apriori file
      open(101, file= glorg_apr_file, status='old', iostat=ierr)
 
      call report_error('IOSTAT',ierr,'open', glorg_apr_file,
     .    0,'update_glorg_apr')

      if( ierr.ne.0 ) glorg_apr_file = ' '
 
*     Now loop over the data file
      ierr = 0
      do while ( ierr.eq.0 )
          read(101,'(a)', iostat=ierr ) buffer

* MOD TAH 0807224: See if reference frame name is passed
          if( (buffer(1:16).eq.'+REFERENCE_FRAME' .or.
     .         buffer(1:16).eq.'+REFERENCE FRAME')  .and.
     .        ierr.eq.0                            ) then
              cpbuff = buffer(17:)
              call trimlead(cpbuff)
              reference_frame = cpbuff
          endif


          if( ierr.eq.0 .and. buffer(1:1).eq.' ' ) then
 
*             See if site name
              indx = 1
              call get_cmd(buffer, gsite_names, -gnum_sites, iel,
     .                         indx)
 
*             If match get the rest
              if( iel.gt.0 ) then

*                 Show that we have updated these coordinates.
                  call sbit(gapr_updated,iel,1)

*                 Get the values for position and velocity.
                  call multiread(buffer, indx, 'R8', jerr, values,
     .                            cval,6)
                  call read_line(buffer, indx, 'R8', jerr, new_epoch,
     .                           cval)
                  call decyrs_to_jd( new_epoch, new_jd)
 
                  dt = (gepoch_out - site_epoch(iel))/365.25d0
                  dtn= (gepoch_out - new_jd)/365.25d0
 
*                                 ! Update position
                  do i = 1,3
                      if( parn_site(i,1,iel).ne.0 ) then
                          jel = parn_site(i,1,iel)
                          parm_change(jel) = parm_change(jel) +
     .                        (values(i,1)+values(i,2)*dtn)-
     .                        (apr_val_site(i,1,iel) +
     .                         apr_val_site(i,2,iel)*dt)
                      end if
                  end do
 
                  do i = 1,3
                      if( parn_site(i,2,iel).ne.0 ) then
                          jel = parn_site(i,2,iel)
                          parm_change(jel) = parm_change(jel)+
     .                           (values(i,2) -
     .                            apr_val_site(i,2,iel))
                      end if
                  end do
 
                  site_epoch(iel) = new_jd
 
                  do i = 1,3
                      do j = 1,2
                          apr_val_site(i,j,iel) = values(i,j)
                      end do
                  end do
              else
* MOD TAH 991110: See if the extended apriori position models have
*                 been given.
                  indx = 1
                  call GetWord(buffer, cval, indx)
                  call casefold(cval)     
                  if( cval(1:8).eq.'EXTENDED' ) then
*                     OK, we have extended non-secular model entry
                      call decode_nonsec(buffer, indx)
                  end if          
              end if
          end if
      end do

****  OK, we have finished reading the file.  Now update the effects
*     of the nonsecular terms.  This needs to be carefully done since
*     the only non-secular terms that can be updated are log terms when
*     logs are estimated in the solution.  (There is an issue here with
*     constraints as well since only loosely constrained parameters can
*     have their apriori values changed.)
      call merge_nonsec(save_num_nonsec, save_param_nonsec,
     .                  save_val_nonsec)

*     Update the apr_val_log values 
      do i = 1, gnum_sites
         call get_nonlog(i,apr_val_log(1,i))
      end do
 
***** Thats all
      if( ierr.eq.-1 ) ierr = 0
      close(101)
      return
      end

CTITLE decode_unequates
 
      subroutine decode_unequates( buffer, indx, num, each, params,
     .                           max_eq, max_each )

      implicit none
 
*     Routine to remove equates from list of equates.  Used mainly
*     to remove parameters added with the eq_dist command.

*   indx    - pointer to position to string
*   num     - Number of equates upto this point
*   max_eq  - Maxiumum number of equates allowed
*   max_each    - Maximum number allowed for each each equate
*   each(max_eq) - Number of equates for of all equates
*   params(max_each, max_eq)    - Parameter numbers to be equated.
 
      integer*4 indx, num, max_eq, max_each, each(max_eq),
     .    params(max_each, max_eq)
 
*   buffer  - Line read from command file.
 
      character*(*) buffer
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error (look for -1 to tell end of list)
*   nt, np  - Temporary pointers to number of parameters read
*   ref     - Reference site number
*   apr     - Apriorui value (not used)
 
      integer*4 ierr, nt, i, j, k, n, ref
* MOD TAH 150130: Made np and array with nnp valiue
      integer*4 np(max_eq), nnp 
      real*8 apr
 

*     initialzie the number of equates for this equate.
      nt = 0
      ierr = 0
      do while ( nt.lt.max_each .and. ierr.eq.0 )
*         Set ref = 0 for each call, so that apr will not
*         be checked since this is an un-equate.
          ref = 0
          call decode_param(buffer,indx, np, nnp, ierr, ref, apr, 0)
 
*         If no error reading the value say it
          if( ierr.eq.0 .and. nnp.gt.0 ) then
              do i = 1, num
                 nt = each(i)
                 do j = 1, nt
* MOD TAH 150130: Implemented array of parameters for multiple sites
                    do n = 1, nnp ! Loop over list 
                       if( params(j,i).eq.np(n) ) then

*                          remove this parmater from list
                           do k = j+1, nt
                              params(k-1,i) = params(k,i)
                           end do
                           each(i) = each(i) - 1
                       end if
                    enddo
                 end do
              end do
          end if
      end do

****  Thats all
      return
      end

CTITLE DECODE_EQ_DIST

      subroutine decode_eq_dist( line, indx, options )

      implicit none

*     This routine will make up a list of equates for all sites within
*     a specified distance of each other

* MOD TAH 060321: Added eq_dist for logs 

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
 
* PASSED VARIABLES
 
*   indx    - Current position in string (updated)
*   options - option that sets if we should update aprioris
 
      integer*4 indx, options 
 
*   line    - Line read from input
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   trimlen - Length of string
*   jndx    - Local positioning in string
*   jerr    - Local IOSTAT error on reading value
*   iel, jel, kel   - General element numbers from get_cmd
*   ne, nt  - Short verions of num_equates and num_each_equate
*   equated(max_glb_site_wrds) - Bit set when site is used in
*             equate.
 
      integer*4  jndx, jerr,  jel, kel, i,j, ne,
     .          nt, equated(max_glb_site_wrds)

*   bl  - Baseline lengths between sites
*   dapr - Apriori value difference

      real*8 bl, dapr

      real*8 eq_sigma, var  ! Sigma of equate (m) and variance (m^2)

*   kbit  - Test bit setting
      logical kbit
 
*   next_word   - Next word in string.
*   types(12)   - Types of parameters
 
      character*8 next_word, types(15), cdum

*   iout        - Output unit number

      integer*4 iout

* xhi_OK -- Set true if this site is not effected by being named to 
*     end if _xhi
      logical xhi_OK

      common / progcon / iout
 
      data types / 'XPOS    ', 'YPOS    ', 'ZPOS    ',
     .            'XDOT    ', 'YDOT    ', 'ZDOT    ',
     .            'NPOS    ', 'EPOS    ', 'UPOS    ',
     .            'NDOT    ', 'EDOT    ', 'UDOT    ',
     .            'NLOG    ', 'ELOG    ', 'ULOG    '  /

 
****  Get the next word from the string

      call read_line(line, indx, 'R8', jerr, eq_dist, cdum)
      call report_error('IOSTAT',jerr,'decod', line, 0, 
     .                  'DECODE_EQ_DIST')
      if( jerr.eq.0 ) then

          call GetWord( line, next_word, indx)
          call casefold( next_word )
          jndx = 1
          call get_cmd(next_word, types, 15, jel, jndx )
*                                 ! OK, compute parameter number
*                                 ! for POS and DOT options (logs are
*                                 ! done below).
* MOD TAH 150528: See if constraint sigma included
          call read_line(line, indx, 'R8', jerr, eq_sigma, cdum)
          if( jerr.eq.0 ) then
              var = eq_sigma**2
          else
              var = 0.0d0
          endif

          if( Jel.gt.0 .and. Jel.le.12 ) then
*                                                 ! Set back to XYZ
              if( jel.gt.6 ) jel = jel - 6

*             Now see if position or rate and set type.
              if( jel.gt.3 ) then
                  jel = jel - 3
                  kel = 2
              else
                  kel = 1
              end if

*             Now loop over all site combinations getting the equates
              do i = 1, max_glb_site_wrds
                 equated(i) = 0
              end do
 
              ne = num_equates

              do i = 1, gnum_sites-1
                 xhi_OK = .true.
                 if ( jel.eq.3 .and. 
     .                gsite_names(i)(5:8).eq.'_XHI' ) xhi_OK = .false.

                 if( .not. kbit(equated,i) .and. 
     .                parn_site(jel, kel,i).ne. 0 .and.
     .                xhi_OK ) then
                    nt = 0
                    do j = i+1, gnum_sites
                       bl = sqrt(
     .                   (apr_val_site(1,1,i)-apr_val_site(1,1,j))**2+
     .                   (apr_val_site(2,1,i)-apr_val_site(2,1,j))**2+
     .                   (apr_val_site(3,1,i)-apr_val_site(3,1,j))**2)
                      if ( jel.eq.3 .and. 
     .                     gsite_names(j)(5:8).eq.'_XHI' ) 
     .                                                xhi_OK = .false.

                       if( bl.le.eq_dist .and. 
     .                     parn_site(jel,kel,j).ne.0 .and.
     .                     xhi_OK ) then

*                          Add this condition.  If nt is zero then
*                          this is first
                           if( nt.eq.0 ) then
                               nt = nt + 1
                               if( ne+1.gt.max_equates ) then
                                   write(*,100) max_equates
 100                               format('***ERROR*** Too many',
     .                                    ' constraints is EQ_DIST',/,
     .                                 16x,'Max Allowed ',i6)
                               else
                                   ne = ne + 1
                                   param_equates(nt,ne) =
     .                                   parn_site(jel, kel,i)
                                   eq_var(ne) = var
                               end if
                           end if
                           if( parn_site(jel, kel,j).gt.0 ) then
                               if( nt+1.gt. max_each_equate ) then
                                   write(*,150) max_each_equate,
     .                                   gsite_names(j)
  150                              format('***ERROR*** To many',
     .                                 ' equates needed to add',
     .                                 ' site ',a8,/,
     .                               16x,'Problem in EQ_DIST')
                               else
                                   nt = nt + 1
                                   param_equates(nt,ne) =
     .                                  parn_site(jel, kel,j)
                                   call sbit(equated,j,1)
                                   num_each_equate(ne) = nt
                               end if
                           end if

* MOD TAH 980929: Check to see if apriori values agree.
                           dapr = apr_val_site(jel,kel,i) -
     .                            apr_val_site(jel,kel,j)
                           if( abs(dapr).gt.1.d-4 ) then
* MOD TAH 030109: Call routine to update aprioris
                               call eqfixa(apr_val_site, site_epoch, 
     .                             parn_site, parm_change, 
     .                             gsite_names, jel, kel, i,j,
     .                             iout, options, gepoch_out)
                           end if
                       end if
                    end do
                 end if
              end do
*         Check for log
          elseif( Jel.gt.12 .and. Jel.le.15 ) then

*             Convert Jel to 1,2,3 for XYZ
              Jel = Jel - 12
*             Now loop over all site combinations getting the equates
              do i = 1, max_glb_site_wrds
                 equated(i) = 0
              end do
 
              ne = num_equates

              do i = 1, gnum_sites-1

                 if( .not. kbit(equated,i) .and. 
     .                parn_log(jel,i).ne. 0 ) then
                    nt = 0
                    do j = i+1, gnum_sites
                       bl = sqrt(
     .                   (apr_val_site(1,1,i)-apr_val_site(1,1,j))**2+
     .                   (apr_val_site(2,1,i)-apr_val_site(2,1,j))**2+
     .                   (apr_val_site(3,1,i)-apr_val_site(3,1,j))**2)

                       if( bl.le.eq_dist .and. 
     .                     parn_log(jel,j).ne.0 ) then

*                          Add this condition.  If nt is zero then
*                          this is first
                           if( nt.eq.0 ) then
                               nt = nt + 1
                               if( ne+1.gt.max_equates ) then
                                   write(*,200) max_equates
 200                               format('***ERROR*** Too many',
     .                                 ' constraints is EQ_DIST LOG',/,
     .                                 16x,'Max Allowed ',i6)
                               else
                                   ne = ne + 1
                                   param_equates(nt,ne) =
     .                                   parn_log(jel,i)
                                   eq_var(ne) = var
                               end if
                           end if
                           if( parn_log(jel,j).gt.0 ) then
                               if( nt+1.gt. max_each_equate ) then
                                   write(*,250) max_each_equate,
     .                                   gsite_names(j)
  250                              format('***ERROR*** To many',
     .                                 ' equates needed to add',
     .                                 ' site ',a8,/,
     .                               16x,'Problem in LOG EQ_DIST')
                               else
                                   nt = nt + 1
                                   param_equates(nt,ne) =
     .                                  parn_log(jel,j)
                                   call sbit(equated,j,1)
                                   num_each_equate(ne) = nt
                               end if
                           end if

* MOD TAH 980929: Check to see if apriori values agree.
                           dapr = apr_val_log(jel,i) -
     .                            apr_val_log(jel,j)
                           if( abs(dapr).gt.1.d-4 ) then
* MOD TAH 030109: Call routine to update aprioris
                               call eqfixlog(jel, i, j, iout, options)
                           end if
                       end if
                    end do
                 end if
              end do

          end if
      end if

*     Save the total number of equates
      num_equates = ne
C     write(*,*) 'EQ_DIST List of equates'
C     do ne = 1, num_equates
C         write(*,*) 'Equates', ne, nt, ' Params ',
C    .         (param_equates(i,ne), i=1,num_each_equate(ne))  
C     end do 


****  Thats all.
      return
      end
 
CTITLE SET_PARAM_NAMES
 
      subroutine set_param_names

      implicit none
 
*     Routine to set the names of parameters (mainly for output
*     with the equate commands in glorg)
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
 
* Local variables
 
*   i,j,k   - Loop counters
*   iel     - Parameter number
 
      integer*4 i,j,k, iel
 
*   site_xyzn(3,2)  - Site parameter codes for XYZ
*   site_neun(3,2)  - Site parameter codes for NEU
 
      character*(4) site_xyzn(3,2), site_neun(3,2), site_logn(3)
 
      data site_xyzn / ' XP ',' YP ',' ZP ',' XD ',' YD ',' ZD ' /
      data site_neun / ' NP ',' EP ',' UP ',' ND ',' ED ',' UD ' /
      data site_logn / ' NL ',' EL ',' UL ' /
 
****  Start by just putting in the numeric values of the parameters.
*     Then we will add the names of those we are interested in (can
*     be updated later).
 
      do i = 1, num_glb_parn
          write(param_names(i), 120) i
 120        format(i12)
      end do
 
****  Now add the parameter names
      do i = 1, gnum_sites
 
*         Loop over XYZ or NEU
          do j = 1, 3
*             Loop over position and velocity
              do k = 1,2
                  iel = parn_site(j,k,i)
                  if( iel.gt.0 ) then
                      if( equate_loc ) then
                          param_names(iel) = gsite_names(i) //
     .                        site_neun(j,k)
                      else
                          param_names(iel) = gsite_names(i) //
     .                        site_xyzn(j,k)
                      end if
                  end if
              end do
          end do
*         Check the earthquake log parameters
          do j = 1, 3
             iel = parn_log(j,i)
             if( iel.gt.0 ) then
                param_names(iel) = gsite_names(i) //
     .                        site_logn(j)
             end if
          end do
 
      end do
 
***** Thats all for the moment
      return
      end
 
CTITLE Check_use_site

      subroutine check_use_site( i, OK ) 

      implicit none

*     Routine to check that position has been estimated for the 
*     use_sites command.
* MOD TAH 030923: Added chech to make sure that _XPS sites are not
*     included.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'

      integer*4 i, j
      logical OK, kbit

***** Scan over the sites in use_Sites and check the parn_site
*     array
*     See if site really used
      OK = .true.
      do j = 1, 3
         if( parn_site(j,1,i).eq.0 ) OK = .false.
      end do

*     Check the name for _XPS sites
      if( gsite_names(i)(5:6).eq.'_X' ) OK = .false.


***** Thats all
      return
      end

CTITLE EQFIXA 

      subroutine eqfixa(apr_val_site, site_epoch, parn_site, 
     .                  parm_change, gsite_names,jel, kel, i,j,
     .                  iout, options, gepoch_out)

      implicit none

*     Routine to update or report differences in aprioris during
*     equates.

* PASSED VARIABLES
*  jel - component number 1,2 or 3
*  kel - Type number 1 position, 2 velocity
*  i   - Reference site number
*  j   - Second site number, its aprioris are updated.
*  iout - Output unit numbers
*  options - Output options for glorg, if bit 21 set then update is
*         done.
*  parn_site - Parameter numbers associates with site pos/vel estimates

      integer*4 jel, kel, i, j, iout, options
      integer*4 parn_site(3,2,*)

*  apr_val_site - Apriori coordinates of sites
      real*8 apr_val_site(3,2,*)
*  site_epoch - Epoch to which apriori coordinates referr
      real*8 site_epoch(*), gepoch_out

*  parm_change  - Changes in the parameter values when the aprioris are
*     updated
      real*8 parm_change(*)
      character*(*) gsite_names(*)   ! Names of sites

* LOCAL VARIABLES
      integer*4 k, np
      real*8 dapr, td  ! change in apriori value and total change
      real*8 dt, dtn   ! Differences in epoch between old and new
      logical kbit ! Checks status of bits
      logical OK   ! Remains true while OK to continue

      character*9 ctype(2)
      character*2 ptype(3)
      character*4 units(2)

      data ctype / 'Position', 'Velocity' /
      data units / 'm   ','m/yr' /

      data ptype / 'X ','Y ','Z' /

****  See if we should fix or just report the difference in apriori

      dapr = apr_val_site(jel,kel,i) - apr_val_site(jel,kel,j)
      OK = .true.
      if ( kbit(options,21) ) then
*         Fix the difference.  For kel=1 position only do this if less
*         than 1 meter separation
          td = sqrt((apr_val_site(1,1,i)-apr_val_site(1,1,j))**2+
     .              (apr_val_site(2,1,i)-apr_val_site(2,1,j))**2+
     .              (apr_val_site(3,1,i)-apr_val_site(3,1,j))**2)
          if( td.gt.1.d0 ) then
              if ( jel.eq.1 ) 
     .        write(iout,120) ctype(kel), gsite_names(i), gsite_names(j)
 120          format('**WARNING** Positions for ',a,' equated sites ',
     .            a8,' and ',a8,' differ by more than 1 m',/,
     .           '            NO UPDATE OF POSITION')
              OK = .false.
          endif
          if( OK ) then
*            Update both postion and velocity
             dt = (gepoch_out - site_epoch(j))/365.25d0
             dtn = (gepoch_out - site_epoch(i))/365.25d0
             site_epoch(j) = site_epoch(i)
*            Change the position
             do k = 1, 3
                dapr = (apr_val_site(k,1,i)+apr_val_site(k,2,i)*dtn)-
     .                 (apr_val_site(k,1,j)+apr_val_site(k,2,j)*dt)
                apr_val_site(k,1,j) = apr_val_site(k,1,i)
                np = parn_site(k, 1,j)
                if( np.gt.0 )
     .              parm_change(np) = parm_change(np)+dapr
             end do
             do k = 1, 3
                dapr = apr_val_site(k,2,i)-apr_val_site(k,2,j)
                apr_val_site(k,2,j) = apr_val_site(k,2,i)
                np = parn_site(k, 2,j)
                if( np.gt.0 )
     .             parm_change(np) = parm_change(np)+dapr
 
             end do
             write(iout,140) 'Position and Velocity', gsite_names(j),
     .                        gsite_names(i)
 140         format('FIXA option set: Updating ',a,' of ',
     .               a,' from ',a)
           else if( kel.eq.2 ) then
              do k = 1, 3
                 dapr = apr_val_site(k,kel,i) - apr_val_site(k,kel,j)
                 apr_val_site(k,kel,j) = apr_val_site(k,kel,i)
                 np = parn_site(k, kel,j)
                 parm_change(np) = parm_change(np)+dapr
              end do
              write(iout,140) ctype(kel), gsite_names(j),gsite_names(i)
           end if
      else 
           write(iout,220) ptype(jel),ctype(kel),
     .          gsite_names(i), gsite_names(j),
     .          dapr, units(kel)
 220      format('**WARNING** FIXA option not set: ',a,a,' of ',
     .            a8,' and ',a8,' differ by ',F10.4,1x,a)
      endif

****  Thats all
      return
      end

CTITLE DECODE_STAB_OPT
 
      subroutine decode_stab_opt( buffer, cont_types, num_types, conts,
     .                          default_cont, list_file, gepoch_start,
     .                          gepoch_end, guse_site )
 

      implicit none
 
*     Routine to read a buffer line from the markov file and decode
*     the contribution types and set the corresponding bits in the
*     'conts' array.
*     If the 'RESET' command is given then conts is set equal to
*     the default_cont.
*     The type may be proceeded with a minus sign to force a
*     contribution to be turned off. (An optional plus sign can
*     also be used, but no sign defaults to plus).

* MOD TAH 030314: Added kalman_param include to get dimensioning
      include '../includes/kalman_param.h'
* MOD TAH 030921: Added wild card feature in decode.  Strings of the form
*     XXXX_* or XXXX_@ will set (or reset) all bits assocaited with strings
*     that match XXXX_ (position of * or @ set the length compared).
* MOD TAH 110816: Added form SITEA/SITEB/SITEC to allow hiarchical 
*     choice of reference sites depending on which sites are available.

* MOD TAH 961206:  Added feature to ignore anything after ! or # in line.
*
*                                         08:47 PM WED., 18 Feb., 1987
*
 
*   bit_state   - state of the bit to be set (set to zero if
*               - minus sign preceeds the type)
*   conts(1)    - The contributions bit map.  A bit is set for
*               - each contribution to be set
*   default_cont(1) - The default contribution pattern.  Set when
*               - the RESET command is given.
*   dummy       - dummy value used as place holder in READ_LINE
*   i           - Loop counter
*   iel         - index in cont_types which matches the next
*               - command in buffer.
*   ierr        - Error flag from READ_LINE.  Error generated
*               - by either by IOSTAT error on read, or -1
*               - generated by reaching the end of the string.
*   indx        - Index used to keep track of where we are
*               - the buffer (used by READ_LINE)
*   next_cont_len   - length of the next_cont string found
*   num_types   - number of contributions allowed. Gives the
*               - length cont_types array.  This value is also
*               - used to compute number of words in the bit
*               - mapped conts when the default value needs to
*               - be assigned.
*   num_words   - number of words in Default_cont (computed from
*               - num_types)
*   trimlen     - HP function to return length of a string
 
      integer*4 bit_state, conts(*), default_cont(*), dummy, i, iel,
     .    ierr, indx, next_cont_len, num_types, num_words, trimlen,
     .    guse_site(max_glb_site_wrds)

* gepoch_start, gepoch_end - Start and end dates of solution
      real*8 gepoch_start, gepoch_end
 
*   buffer      - line read from the Markov file.  Should contain
*               - the list of contributions to be applied or
*               - reset (if minus sign used)
*   cont_types(1)   - the array of types which are set.  This
*               - array must have the entries in the same order
*               - which the bits will be set.
*  list_file    - Name of the gdl file
 
      character*(*) buffer, cont_types(*), list_file
 
*   next_cont   - The next contribution string from buffer.  One
*               - extra character is added to allow for plus or
*               - minus sign.
* MOD TAH 110816: Added form SITEA/SITEB/SITEC to allow hiarchical 
*                 choice.
 
      character*128 next_cont
      character*16  split_cont

* date(5) -- Date read for converting range
* curr_set -- Conts set in this command
* cuur_unset -- Conts unset by this command

      integer*4 date(5), curr_set(max_glb_site_wrds), kndx, 
     .          curr_unset(max_glb_site_wrds), jndx, jerr, iwc

      integer*4 k, split_cont_len
      
      real*8 sectag, jd
      logical apply, kbit, done, site_OK 

      character*4 cdum
  
 
***** Loop through the buffer getting all the contrubution names
*     until we run out of entries or an error (ierr.ne.0)
 
      ierr = 0
*                 ! Start at the beginning of the line
      indx = 1
      do i = 1, max_glb_site_wrds
         curr_set(i) = 0
         curr_unset(i) = 0
      end do
      apply = .true.
 
      do while ( ierr.eq.0 )
 
*         Get the next string in buffer
          call read_line( buffer, indx, 'CH', ierr, dummy, next_cont)
          call casefold(next_cont)
          if( next_cont(1:1).eq.'#' .or. 
     .        next_cont(1:1).eq.'!' ) ierr = -1

*         See if start of restrictions on list
          if( next_cont(1:4).eq.'R   ') ierr = +1
 
*         Only process if ierr is zero (no error)
          if( ierr.eq.0 ) then
 
*             strip + or - sign from start of next_cont and set
*             bit_state acordingly
              call get_bit_state( next_cont, bit_state)
 
*             Now see if we can find this contribution, if we cant
*             ignore entry unless it is RESET.
 
*                             ! Set not found value
              iel = -1

* MOD TAH 110816: Start spliting the site names based on separation
*             by / (allows hiarchial listing of frame sites) 
              jndx = 1     ! Start of string
              kndx = jndx
              done = .false.
              
              do while ( .not. done )
                 k = index(next_cont(jndx+1:),'/')
                 if( k.gt.0 ) then
                     jndx = jndx + k
                 else
                    jndx = trimlen(next_cont)+1
                    done = .true.
                 end if

                 split_cont = next_cont(kndx:jndx-1)
                 kndx = jndx + 1
                 split_cont_len = min( len(cont_types(1)),
*                                                            ! We use this
     .                               trimlen( split_cont ) )
*                                    ! quantity so that shortened version
*                                    ! of the contributions can be used.
*                                    ! The first match will be used.
* MOD TAH 030921: Check for wild cards in the split_cont string
                 iwc = 0
                 iwc = index(split_cont,'*')
                 if( iwc.eq.0 ) iwc = index(split_cont,'@')
                 if( iwc.gt.0 ) split_cont_len = 
     .                                       min(split_cont_len,iwc-1)
* MOD TAH 030921: Removed check on finding element so that all are searched
*                when there is a wild card.  (Short strings will now find all
*                as well)
                 i = 0
                 do while ( i.lt.num_types )
 
                     i = i + 1
                     if( (cont_types(i)(1:split_cont_len) .eq.
     .                   split_cont(1:split_cont_len)) ) then  ! Make sure site is used
*                        See if site OK to use
                         call check_use_site(i, site_OK)
                         if( site_OK ) then
                            iel = i
                            done = .true.   ! Find site so don't use others
* MOD TAH 030921:           Set the bit as we find each element
                            if( bit_state.eq.1 ) then
                                call sbit(curr_set,iel,1)
                            else
                                call sbit(curr_unset,iel,1)
                            endif
                         end if
                     endif
                 end do
              end do
 
*             If we found next_cont then set bit in cont according to
*             bit state
 
*                                     ! Found
              if( iel.gt.0 ) then
* MOD TAH 030921: Moved the code below upto the above loops over
*                 num_types so that wild cards can be used.
C                 if( bit_state.eq.1 ) then
C                     call sbit(curr_set,iel,1)
C                 else
C                     call sbit(curr_unset,iel,1)
C                 endif
*                                     ! Not found, check for reset
              else
*                                     ! or clear
*                                                       ! Reset to
                  num_words = (num_types+31)/32
                  if( next_cont(1:6).eq.'RESET ' ) then
*                                     ! default value.
                      do i = 1,num_words
                          conts(i) = default_cont(i)
                      end do
*                                     ! RESET cont
                  end if
                  if( next_cont(1:4).eq.'ALL ' ) then
                      do i = 1, num_words
                         conts(i) = -1
                      end do 
                  end if
*                                                      ! Clear values
                  if( next_cont(1:6).eq.'CLEAR ' .or.
     .                next_cont(1:5).eq.'NONE '      ) then
                      do i = 1,num_words
                          conts(i) = 0
                      end do
*                                     ! CLEAR conts
                  end if
*                                     ! contribution found
              end if
*                                     ! no error reading buffer
          end if

*         If ierr = 1, then start reading the restriction range
          if( ierr.eq.1 ) then
*             See if a list file name restriction
              jndx = indx
              call read_line( buffer, indx, 'CH', ierr, dummy, 
     .                        next_cont)
*             See if this is number or string
              call check_num(next_cont, jerr)
              if( jerr.ne.0 ) then
                  call get_bit_state(next_cont, bit_state)
*                 See if we make decision based on this
                  jndx = index(list_file,
     .                         next_cont(1:trimlen(next_cont)))
                  if( bit_state.eq.1 ) then
*                     list file should contain this string
                      if( jndx.eq.0 ) apply = .false.
                  else
                      if( jndx.eq.1 ) apply = .false.
                  endif
                  jndx = indx
              endif

*             Now read the start time needed
              indx = jndx
              do i = 1, 5
                 call read_line(buffer,indx,'I4',jerr,date(i),cdum)
              end do
              sectag = 0.0d0
              call ymdhms_to_jd(date,sectag,jd)
              if( gepoch_start.lt.jd .and. apply ) apply = .false.
*             If we still have entries in the string check the end date
              if( jerr.eq.0 ) then
                 do i = 1, 5
                    call read_line(buffer,indx,'I4',jerr,date(i),cdum)
*                   See if first value is there
                    if( i.eq.1 .and. jerr.ne.0 ) then
                        date(1) = 2100
                    endif
                 end do
                 sectag = 0.0d0
                 call ymdhms_to_jd(date,sectag,jd)
                 if( gepoch_end.gt.jd .and.apply ) apply = .false.
             endif
          endif 
*                                     
      end do               ! looping until end of string

****  Now see if we should apply these lines
      if( apply ) then
         do i = 1, max_glb_site_wrds*32
            if( kbit(curr_set,i) ) call sbit(conts,i,1)
            if( kbit(curr_unset,i) ) call sbit(conts,i,0)
         end do
      else
         write(*,320) list_file(1:trimlen(list_file)),
     .                buffer(1:trimlen(buffer))
 320     format('For ',a,' Stabilization not using: ',a)
      end if

***** Thats all
      return
      end
 
CTITLE eqfixlog 

      subroutine eqfixlog(jel, ref,upd, iout, options)

      implicit none

*     Routine to update or report differences in aprioris logarithm
*     terms during equates

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'

* PASSED VARIABLES
*  jel - component number 1,2 or 3
*  ref - Reference site number
*  upd - Second site number, its aprioris are updated.
*  iout - Output unit numbers
*  options - Output options for glorg, if bit 21 set then update is
*         done.
*  parn_site - Parameter numbers associates with site pos/vel estimates

      integer*4 jel, ref, upd, iout, options

* LOCAL VARIABLES
      integer*4 i,j, k, np
      logical kbit ! Checks status of bits
      real*8 dapr    ! Difference in NEU apriori
     .,      old_log(3), new_log(3)  ! old and new NEU of log term changed
     .,      old_xyz(3), new_xyz(3)  ! old and new Log term comnverted to XYZ

      character*2 ptype(3)

      data ptype / 'N ','E ','U' /

****  See if we should fix or just report the difference in apriori

      dapr = apr_val_log(jel,ref) - apr_val_log(jel,upd)
      if( kbit(options,21) ) then
          np = parn_log(jel,upd)
          if( np.gt.0 ) 
     .       parm_change(np) = parm_change(np)+dapr


****      Save the in and out values because we need to remove the 
*         old either NEU contribution from XYZ values before adding
*         in the new value
          do i = 1,3
            old_log(i) = 0.d0
            new_log(i) = 0.d0
          end do
          old_log(jel) = apr_val_log(jel,upd)
          new_log(jel) = apr_val_log(jel,ref)

          apr_val_log(jel,j) = apr_val_log(jel,i)
          write(iout,140) ptype(jel), gsite_names(upd),gsite_names(ref),
     .            dapr
 140      format('FIXA option set: Updating ',a,' Log of ',
     .               a,' from ',a,' by ',F7.4,' m')

****      Now save the new apriori back in the non-secular parameters
          do k = 1, num_nonsec
             if( param_nonsec(1,k).eq.upd .and. 
     .           param_nonsec(2,k).eq. 4       ) then   ! type 4 is log

****             Convert the old values to XYZ and remove 
                 call nonsec_convert('TOXYZ',1,old_log,
     .                         old_xyz, apr_val_site(1,1,upd))
*                Remove these from the apr_val_nonsec entries
                 do j = 1,3
                    apr_val_nonsec(2+j,k) = apr_val_nonsec(2+j,k) -
     .                      old_xyz(j)
                 end do
*                Compute new term and add
                 call nonsec_convert('TOXYZ',1,new_log,
     .                         new_xyz, apr_val_site(1,1,upd))
*                Remove these from the apr_val_nonsec entries
                 do j = 1,3
                    apr_val_nonsec(2+j,k) = apr_val_nonsec(2+j,k) +
     .                      new_xyz(j)
                 end do
             end if
          end do

      else 
           write(iout,220) ptype(jel),
     .         gsite_names(ref), gsite_names(upd),
     .         dapr
 220      format('**WARNING** FIXA option not set: ',a,' Log of ',
     .        a8,' and ',a8,' differ by ',F7.4,' m')
      endif

****  Thats all
      return
      end

CTITLE DECODE_EQ_4CHAR

      subroutine decode_eq_4char( line, indx, options )

      implicit none

*     This routine will make up a list of equates for all sites with the
*     same 4-character IDs

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
 
* PASSED VARIABLES
 
*   indx    - Current position in string (updated)
*   options - option that sets if we should update aprioris
 
      integer*4 indx, options 
 
*   line    - Line read from input
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   trimlen - Length of string
*   jndx    - Local positioning in string
*   jerr    - Local IOSTAT error on reading value
*   iel, jel, kel   - General element numbers from get_cmd
*   ne, nt  - Short verions of num_equates and num_each_equate
*   equated(max_glb_site_wrds) - Bit set when site is used in
*             equate.
 
      integer*4  jndx, jerr,  jel, kel, i,j, ne,
     .          nt, equated(max_glb_site_wrds)

*   dapr - Apriori value difference

      real*8 dapr

      real*8 eq_sigma, var  ! Sigma of equate (m) and variance (m^2)

*   kbit  - Test bit setting
      logical kbit
 
*   next_word   - Next word in string.
*   types(12)   - Types of parameters
 
      character*8 next_word, types(15), cdum

*   iout        - Output unit number

      integer*4 iout

* xhi_OK -- Set true if this site is not effected by being named to 
*     end if _xhi
      logical xhi_OK

      common / progcon / iout
 
      data types / 'XPOS    ', 'YPOS    ', 'ZPOS    ',
     .            'XDOT    ', 'YDOT    ', 'ZDOT    ',
     .            'NPOS    ', 'EPOS    ', 'UPOS    ',
     .            'NDOT    ', 'EDOT    ', 'UDOT    ',
     .            'NLOG    ', 'ELOG    ', 'ULOG    '  /

 
****  Get the next word from the string

      call GetWord( line, next_word, indx)
      call casefold( next_word )
      jndx = 1
      call get_cmd(next_word, types, 15, jel, jndx )
*                             ! OK, compute parameter number
*                             ! for POS and DOT options (logs are
*                             ! done below).
* MOD TAH 150528: See if constraint sigma included
      call read_line(line, indx, 'R8', jerr, eq_sigma, cdum)
      if( jerr.eq.0 ) then
          var = eq_sigma**2
      else
          var = 0.0d0
      endif

      if( Jel.gt.0 .and. Jel.le.12 ) then
*                                             ! Set back to XYZ
          if( jel.gt.6 ) jel = jel - 6

*         Now see if position or rate and set type.
          if( jel.gt.3 ) then
              jel = jel - 3
              kel = 2
          else
              kel = 1
          end if

*         Now loop over all site combinations getting the equates
          do i = 1, max_glb_site_wrds
             equated(i) = 0
          end do
 
          ne = num_equates

          do i = 1, gnum_sites-1
             xhi_OK = .true.
             if ( jel.eq.3 .and. 
     .            gsite_names(i)(5:8).eq.'_XHI' ) xhi_OK = .false.

             if( .not. kbit(equated,i) .and. 
     .            parn_site(jel, kel,i).ne. 0 .and.
     .            xhi_OK ) then
                nt = 0
                do j = i+1, gnum_sites
                  if ( jel.eq.3 .and. 
     .                 gsite_names(j)(5:8).eq.'_XHI' ) 
     .                                            xhi_OK = .false.

                   if( gsite_names(i)(1:4).eq. gsite_names(j)(1:4) 
     .                 .and. parn_site(jel,kel,j).ne.0 .and.
     .                 xhi_OK ) then

*                      Add this condition.  If nt is zero then
*                      this is first
                       if( nt.eq.0 ) then
                           nt = nt + 1
                           if( ne+1.gt.max_equates ) then
                               write(*,100) max_equates
 100                           format('***ERROR*** Too many',
     .                                ' constraints is EQ_4CHAR',/,
     .                             16x,'Max Allowed ',i6)
                           else
                               ne = ne + 1
                               param_equates(nt,ne) =
     .                               parn_site(jel, kel,i)
                               eq_var(ne) = var
                           end if
                       end if
                       if( parn_site(jel, kel,j).gt.0 ) then
                           if( nt+1.gt. max_each_equate ) then
                               write(*,150) max_each_equate,
     .                               gsite_names(j)
  150                          format('***ERROR*** To many',
     .                             ' equates needed to add',
     .                             ' site ',a8,/,
     .                           16x,'Problem in EQ_4CHAR')
                           else
                               nt = nt + 1
                               param_equates(nt,ne) =
     .                              parn_site(jel, kel,j)
                               call sbit(equated,j,1)
                               num_each_equate(ne) = nt
                           end if
                       end if

* MOD TAH 980929: Check to see if apriori values agree.
                       dapr = apr_val_site(jel,kel,i) -
     .                        apr_val_site(jel,kel,j)
                       if( abs(dapr).gt.1.d-4 ) then
* MOD TAH 030109: Call routine to update aprioris
                           call eqfixa(apr_val_site, site_epoch, 
     .                         parn_site, parm_change, 
     .                         gsite_names, jel, kel, i,j,
     .                         iout, options, gepoch_out)
                       end if
                   end if
                end do
             end if
          end do
*     Check for log
      elseif( Jel.gt.12 .and. Jel.le.15 ) then

*         Convert Jel to 1,2,3 for XYZ
          Jel = Jel - 12
*         Now loop over all site combinations getting the equates
          do i = 1, max_glb_site_wrds
             equated(i) = 0
          end do
 
          ne = num_equates

          do i = 1, gnum_sites-1

             if( .not. kbit(equated,i) .and. 
     .            parn_log(jel,i).ne. 0 ) then
                nt = 0
                do j = i+1, gnum_sites

                   if( gsite_names(i)(1:4).eq. gsite_names(j)(1:4) 
     .                 .and. parn_log(jel,j).ne.0 ) then

*                      Add this condition.  If nt is zero then
*                      this is first
                       if( nt.eq.0 ) then
                           nt = nt + 1
                           if( ne+1.gt.max_equates ) then
                               write(*,200) max_equates
 200                           format('***ERROR*** Too many',
     .                             ' constraints is EQ_4CHAR LOG',/,
     .                             16x,'Max Allowed ',i6)
                           else
                               ne = ne + 1
                               param_equates(nt,ne) =
     .                               parn_log(jel,i)
                               eq_var(ne) = var
                           end if
                       end if
                       if( parn_log(jel,j).gt.0 ) then
                           if( nt+1.gt. max_each_equate ) then
                               write(*,250) max_each_equate,
     .                               gsite_names(j)
  250                          format('***ERROR*** To many',
     .                             ' equates needed to add',
     .                             ' site ',a8,/,
     .                           16x,'Problem in LOG EQ_4CHAR')
                           else
                               nt = nt + 1
                               param_equates(nt,ne) =
     .                              parn_log(jel,j)
                               call sbit(equated,j,1)
                               num_each_equate(ne) = nt
                           end if
                       end if

* MOD TAH 980929: Check to see if apriori values agree.
                       dapr = apr_val_log(jel,i) -
     .                        apr_val_log(jel,j)
                       if( abs(dapr).gt.1.d-4 ) then
* MOD TAH 030109: Call routine to update aprioris
                           call eqfixlog(jel, i, j, iout, options)
                       end if
                   end if
                end do
             end if
          end do

      end if

*     Save the total number of equates
      num_equates = ne
C     write(*,*) 'EQ_4CHAR List of equates'
C     do ne = 1, num_equates
C         write(*,*) 'Equates', ne, nt, ' Params ',
C    .         (param_equates(i,ne), i=1,num_each_equate(ne))  
C     end do 


****  Thats all.
      return
      end
