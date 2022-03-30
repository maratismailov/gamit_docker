CTITLE PROCESS_GLB_COMMAND
 
      subroutine process_glb_command( buffer, indx, iel, glinit_run)

      implicit none
 
 
*     Routine to process the globk command found in the command file.
*     The remaining arguments may need to be read from the command
*     buffer
*     The commands are processed in two parts.  The commands which
*     can be passed through the runstring are in one group, the
*     remaining commands are in the other.  Once any command from the
*     second group is executed, the GLINIT program is run (only once)
*     so that we can get the site and source names.
*
*     NOTE: the GLOBK common file is actually created in GLINIT.
*     Therefore some information is saved temporarily here so that
*     it can be added to the common.
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/globk_cmds.h'
 
*   iel         - Command number
*   ierr        - IOSTAT/READ_LINE error
*   indx        - pointer to current position in command buffer
*   is          - Site or source number from command
*   num_period  - Period number for coefficients
*   num_name    - Name of the coefficient to be estimated
*   num_etd_site - Number of etd site to be estimated.
*   num_term    - Specific coefficient to be estimated
*   gnum_temp   - Temporary version of the the number of sites
*                 used for global tides when there are more
*                 global sites than the number allowed in
*                 solutiuon (max_etd_sites)
*   num_miss    - Number of missing arguments.  Reported only for
*                 some commands.
*   num_use     - Sets number of times a site must a used.
*   iwc         - Position of wild card in name
 
      integer*4 iel, ierr, indx, is, i,j, num_period, num_name, 
     .          num_etd_site, num_term, gnum_temp, trimlen, 
     .          num_use, date(5), iwc
      real*8    sectag

* k, jndx  -- Counter and test position in string
* eol      -- Character number for ! or # (buffer cleared after this
*             position in string

      integer*4 k, jndx, eol, jerr
      character*256 word
      character*8 name

*   OK          - Indicates that etd coefficient site is OK
*   glinit_run  - Indicates that glinit has been run
 
      logical glinit_run, OK, kbit, done
 
*   vals(100)    - Upto 100 values read from command line (before
*               - transferr to ema)
*   decon      - Deconstrain when free_log used
*   dtl, dpdt  - Log time and log value
*   logsig     - log sigma to be used to make position sigma
*                equal to decon constraint.
 
      real*8 vals(100), decon, dtl, dpdt, logsig 
      real*4 val4    ! Dummy real*4 for read_line
 
*   cval        - Dummy character string for read_line
 
      character*16 cval

*   glb_opt_types(3) -- Options for outputing global files
      character*8 glb_opt_types(3)
 
*   buffer      - Buffer read from command file
 
      character*(*) buffer

      data glb_opt_types / 'SITES   ','EOP     ','ORBITS  ' /
 
*     Use computed goto to get to command.  Split this into two parts
*     so that glinit can be run if need be.  Because of the way computed
*     are handled in FTN77 (compared to FTN4/66) we need to explicitily
*     check if iel is > than our range or if it is past the point where
*     glinit should be run. (In FTN77 if iel is > then last label then
*     first label is executed where as in earlier versions of fortran
*     the last would have been executed
 
***** Any command past iel = 6 will cause GLINIT to run, if it has
*     not been run.
 
*     Check to see if glinit run yet
*                                                    ! Run now
      if( (iel.gt.7 .and. iel.ne.92) .and. .not. glinit_run ) then
          call run_glinit
          glinit_run = .true.
          wild_mjd = (gepoch_start + gepoch_end)/2 - 2 400 000.5d0
          write(*,40) wild_mjd
40        format('Using ',F10.3,' MJD for wild date')
      end if
      if( (iel.le.7 .or. iel.eq.92) .and. glinit_run ) then
          write(*, 50) buffer(1:max(1,trimlen(buffer)))
 50       format('**ERROR** The command ',a,/,
     .           '          can not be used at this time.  Move to',
     .           ' the top of the globk command file'/,
     .           '  Command being ignored')
          call report_stat('FATAL','GLOBK','read_glb_mar',glb_mar_file,
     .           'Commands out of order',0)
          RETURN
      end if
  
*     Now make sure command in range.  This should trap ALL command if
*     used
* MOD TAH 980517: Added SOURCE command to source to a new file (can only
*     be used in to level).  This is command number 78 and is processed
*     in read_glb_markov rather than in here.
* MOD TAH 030611: Increased from 90 to 92 for the decimate and glb_opt commands.
* MOD TAH 080802: Increased 95 to 96 when FREE_LOG command added

      if( iel.gt.97 ) iel = 98 

* MOD TAH 981020: Clear any part of the buffer after a ! or # symbol
      eol = index( buffer,'!' )
      if( eol.eq.0 ) eol = index(buffer,'#')
      if( eol.gt.0 ) then
          buffer(eol:) = ' '
      end if

***** Now continue processing the commands
      goto (  100,  200,  300,  350,  375,  400,  450,
     .                                500,  600,  700,  800,  900, 1000,
     .       1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
     .       2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
     .       3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
     .       4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000,
     .       5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000,
     .       6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000,
     .       7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900, 8000,
     .       8100, 8200, 8300, 8400, 8500, 8600, 8700, 8800, 8900, 9000,
     .       9100, 9200, 9300, 9400, 9500, 9600 ) iel
 
****  COM_FILE: Get common file name
  100 continue
          call read_line(buffer, indx,'CH',ierr, vals, glb_com_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( glb_com_file, list_file)
*                                             ! If error use default
          if( ierr.ne.0 ) glb_com_file = ' '
          return
 
****  SRT_FILE: Get sort file name
  200 continue
          call read_line(buffer, indx,'CH',ierr, vals, sort_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( sort_file , list_file)
*                                          ! If error use default
          if( ierr.ne.0 ) sort_file = ' '
          return
 
****  LST_FILE: Get list file name
  300 continue
          call read_line(buffer, indx,'CH',ierr, vals, list_file)
          call report_error('READ_LINE',ierr,'read',buffer,1,
*                                                  ! Kill, no data
     .                      'PROCESS_GLB_COMMAND')
          return

****  EQ_FILE: Get the earthquake file name 
  350 continue
          num_eqfiles = num_eqfiles + 1
          if( num_eqfiles.gt.max_eqfiles) then
              call report_stat('fatal','globk','process_glb',
     .           'eq_files','Too many eq-files',max_eqfiles)
          end if
          call read_line(buffer, indx,'CH',ierr, vals,
     .          eq_inp_file(num_eqfiles))
          call report_error('READ_LINE',ierr,'read',buffer,1,
*                                                  ! Kill, no data
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( eq_inp_file(num_eqfiles), list_file)
          cval = eq_inp_file(num_eqfiles)
          call casefold(cval)
* MOD TAH 160324: Removed RESET so that name can be used as RESET command
          if( cval(1:4).eq.'NONE' ) then
              num_eqfiles = 0
              eq_inp_file(1) = ' '
          endif
          return

****  MAKE_SVS: Tells glinit to make the svs_file.
  375 continue
          make_svs_file = .true.
          call read_line(buffer, indx,'CH',ierr, vals, glb_svs_file)
* MOD TAH 190617: See if NONE is the name which will reset the command.
*         (Useful if options turn in svs_file e.g.
          call report_error('READ_LINE',ierr,'read',buffer,1,
     .                      'make_svs/process_glb_command')
          if( glb_svs_file(1:5).eq.'NONE ' ) then
*            Reset the command
             make_svs_file = .false.
             glb_svs_file = ' '
             write(*,*) 'Resetting make_svs_file command'
          else
*            Normal processing
             call wild_card( glb_svs_file, list_file)

* MOD TAH 970317: See if option argument passed to zero the
*            radiation parametrers.
             call read_line(buffer,indx,'CH',ierr, vals, cval)
             call casefold(cval)
             if( ierr.eq.0 .and. cval(1:1).eq.'Z' ) then
                 glb_svs_file(trimlen(glb_svs_file)+1:) = '_Z'
             else  
                 glb_svs_file(trimlen(glb_svs_file)+1:) = '_A'
             end if 
          endif 
          return

****  SORT_DIR: Get sort direction
  400 continue
          call read_line(buffer, indx,'I4',ierr, sort_direction, cval)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          if( ierr.ne.0 ) sort_direction = +1
          return

****  USE_PRNN <Y/N>: Sets whether to used PRN names rather than
*     GNSS versions of satellite names
 450  continue
          call GetWord( buffer, word, indx )
          if( word(1:1).eq.'Y' .or. word(1:1).eq.'y' ) then
             use_prnn = .true.
          else
             use_prnn = .false.
          endif
          return

****  APR_FILE:  Get name of file with apriori site/source positions. 
  500 continue

*         Increment the number and make sure we don't have too many
          num_apr_files = num_apr_files + 1
*         Leave one file name space for glorg
          if( num_apr_files.ge. max_apr_files ) then
              write(*,520)  max_apr_files, buffer
 520          format('**ERROR** Maximum number of apriori files ',
     .               'allowed exceeded (',i2,' allowed)',/,
     .               '          Ignoring this command: ',a)
          else
              i = num_apr_files
              call read_line(buffer, indx,'CH',ierr,vals, 
     .                       glb_apr_file(i))
              call report_error('READ_LINE',ierr,'read',buffer,0,
     .                          'PROCESS_GLB_COMMAND')
              call wild_card( glb_apr_file(i), list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
              call wild_date( glb_apr_file(i), wild_mjd, 'N' )
* MOD TAH 091006: See of NONE or RESET passed
              if( glb_apr_file(i)(1:4).eq.'NONE' .or.
     .            glb_apr_file(i)(1:5).eq.'RESET' ) then
                  num_apr_files = 0
              endif

              return
          end if
 
****  SOL_FILE:  Get name of solution scratch file
  600 continue
          call read_line(buffer, indx,'CH',ierr, vals, glb_sol_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( glb_sol_file, list_file )
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
          call wild_date( glb_sol_file, wild_mjd, 'N' )
*                                             ! If error use default
          if( ierr.ne.0 ) glb_sol_file = ' '
          return
 
****  BAK_FILE: Get name of back file
  700 continue
          call read_line(buffer, indx,'CH',ierr,vals, glb_bak_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          if( ierr.eq.0 ) then
*             Check for wild card
              call wild_card( glb_bak_file, list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
              call wild_date( glb_bak_file, wild_mjd, 'N' )
              call set_use( glb_bak_file, glb_bak_soln )
          end if
          return
 
****  APR_AXO: Get apriori sigmas for axis offest
  800 continue
 
*         Get name of site
          call get_cmd( buffer, gsite_names, gnum_sites, is,indx)
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 2)
          call assign_ema_val4( vals, 2, apr_axo, is, gnum_sites)
          return
 
****  BAK_PRTS: Sets the sites to printed in the back file.  If the
*               CLEAR option is used, these will be only sites 
*               printed.  Otherwise they will be added to the 
*               default lists.
  900 continue
 
*         Get names of the sites
          if( index(buffer,'CLEAR').gt.0 ) clear_bo_site = .true.
          call decode_option( buffer, gsite_names, gnum_sites,
     .                        bak_out_site, -1)
          return
 
****  APR_GAM: Get apriori sigmas for gamma
 1000 continue
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 1)
          call assign_ema_val4( vals, 1, apr_gamma, 1, 1)
          return
 
****  APR_NANG: Get apriori sigmas for nutation angles, and seasonal model
 1100 continue
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 8)
          call assign_ema_val4( vals, 8, apr_nut_ang, 1, 1)
          return
 
****  APR_NCOE: Get apriori sigmas for nutation series coefficients
 1200 continue
 
*         Get name of nutation components
          call get_cmd( buffer, nut_names, max_nut_coeff, is,indx)
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 2)
          call assign_ema_val4( vals, 2, apr_nut_coeff, is,
     .                          max_nut_coeff)
          return
 
****  APR_RAO: Get apriori sigmas for right ascension orgin.
 1300 continue
 
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 1)
          call assign_ema_val4( vals, 1, apr_rao , 1, 1)
          return
 
****  APR_SITE: Get apriori sigmas for site positions and rates
 1400 continue
 
*         Get name of site
          call get_cmd( buffer, gsite_names, gnum_sites, is,indx)
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 6)
          call assign_ema_val4( vals, 6, apr_site, is, gnum_sites)
          return
 
****  APR_SOUR: Get apriori sigmas for source positions and rates
 1500 continue
 
*         Get name of source
          call get_cmd( buffer, gsource_names, gnum_sources,
     .                      is,indx)
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 4)
          call assign_ema_val4( vals, 4, apr_source, is, gnum_sources)
          return
 
****  APR_TID : Get apriori sigmas for h,l, and lag
 1600 continue
 
*         Get name of site
          call get_cmd( buffer, gsite_names, gnum_sites, is,indx)
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 3)
          call assign_ema_val4( vals, 3, apr_tid , is, gnum_sites)
          return
 
****  AMR_TRAN: Get apriori sigmas for site origin translation
 1700 continue
 
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 3)
          call assign_ema_val4( vals, 3, apr_tran, 1, 1)
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 3)
          call assign_ema_val4( vals, 3, mar_tran, 1, 1)
          do i = 1,3
             apr_tran(i,2) = 0.0
             mar_tran(i,2) = 0.0
          end do
          return
 
****  APR_WOB : Get apriori sigmas for wobble values, rates and seasonal
*               terms
 1800 continue
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 4)
          call assign_ema_val4( vals, 4, apr_wob , 1, 1)
          return
 
****  APR_UT1 : Get apriori sigmas for UT1 values, rates and seasonal
*               terms
 1900 continue
 
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 2)
          call assign_ema_val4( vals, 2, apr_ut1 , 1, 1)
          return
 
****  MAR_AXO: Get variances for markov axis offest
 2000 continue
 
*         Get name of site
          call get_cmd( buffer, gsite_names, gnum_sites, is,indx)
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 2)
          call assign_ema_val4( vals, 2, mar_axo, is, gnum_sites)
          return
 
****  MAR_NANG: Get variances for markov nutation angles and seasonal model
 2100 continue
 
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 8)
          call assign_ema_val4( vals, 8, mar_nut_ang, 1, 1)
          return
 
****  MAR_SITE: Get variances for markov site positions and rates
 2200 continue
 
*         Get name of site
          call get_cmd( buffer, gsite_names, gnum_sites, is,indx)
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 6)
          call assign_ema_val4( vals, 6, mar_site, is, gnum_sites)
          return
 
****  MAR_SOUR: Get variances for markov source positions and rates
 2300 continue
 
*         Get name of source
          call get_cmd( buffer, gsource_names, gnum_sources,
     .                      is,indx)
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 4)
          call assign_ema_val4( vals, 4, mar_source, is, gnum_sources)
          return
 
 
****  MAR_WOB : Get variances for markov wobble values, rates and seasonal
*               terms
 2400 continue

          call multiread( buffer, indx, 'R8', ierr, vals, cval, 4)
          call assign_ema_val4( vals, 4, mar_wob , 1, 1) 
          return
 
****  MAR_UT1 : Get variances for markov UT1 values, rates and seasonal
*               terms
 2500 continue

          call multiread( buffer, indx, 'R8', ierr, vals, cval, 2)
          call assign_ema_val4( vals, 2, mar_ut1 , 1, 1)
          return
 
****  VAL_AXO: Get apriori values for axis offest, rate and epoch
 2600 continue
 
*         Get name of site
          call get_cmd( buffer, gsite_names, gnum_sites, is,indx)
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 3)
          call assign_ema_val8( vals, 2, apr_val_axo, is, gnum_sites)
 
*         Get epoch for rate
          call decyrs_to_jd( vals(3), vals(4) )
          call assign_ema_val8( vals(4), 1, axo_epoch, is, gnum_sites)
 
          return
 
****  VAL_ETD: Get apriori values for extended earth coefficients
 2700 continue
 
*         Get name of earthtide components
          call get_cmd( buffer, etd_names, max_etd_coeff, is,indx)
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 2)
          call assign_ema_val8( vals, 2, apr_val_etd_coeff,
     .                          is, max_etd_coeff)
          return
 
****  VAL_NCOE: Get apriori values for nutation series coefficients
 2800 continue
 
*         Get name of nutation components
          call get_cmd( buffer, nut_names, max_nut_coeff, is,indx)
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 2)
          call assign_ema_val8( vals, 2, apr_val_nut_coeff, is,
     .                          max_nut_coeff)
          return
 
****  VAL_WOB : Get apriori values for wobble values, rates and seasonal
*               terms
 2900 continue
 
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 8)
          call assign_ema_val8( vals, 8, apr_val_wob , 1, 1)
          return
 
****  VAL_UT1 : Get apriori values for UT1 values, rates and seasonal
*               terms
 3000 continue
 
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 6)
          call assign_ema_val8( vals, 6, apr_val_ut1 , 1, 1)
          return
 
***** COMP_RES : Compute postfit residials option.  If this command given
*                compute_glb_res is set true
 3100 continue
          compute_glb_res = .true.
          return
 
****  OUT_PMU : Get the name of the output file for polar motion and UT1
 3200 continue
          call read_line( buffer, indx, 'CH',ierr,vals, pmu_table_file)
          call wild_card( pmu_table_file, list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
          call wild_date( pmu_table_file, wild_mjd, 'N' )
          return
 
****  OUT_NUT : Get the name of the output file for nutation values
 3300 continue
          call read_line( buffer, indx, 'CH',ierr,vals, nut_table_file)
          call wild_card( nut_table_file, list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
          call wild_date( nut_table_file, wild_mjd, 'N' )
          return
 
****  OUT_GLB : Get the name of the output global solution file (BINARY)
 3400 continue
          call read_line( buffer, indx, 'CH',ierr,vals, glb_out_file)
          call wild_card( glb_out_file, list_file )
****      See if none or NONE passed
          cval = glb_out_file
          call casefold(cval)

* MOD TAH 190624: Added wild_date option (I will increment in back solution
*         if user does not use wild_date features).
*         Only call wild_date if replace name is not used.
          if( cval(1:4).ne.'NONE' ) then
              call wild_date( glb_out_file, wild_mjd, 'I' )
          endif

          if ( cval(1:4).eq.'NONE' ) glb_out_file = ' '

*         See if glorg output option passed
          call read_line( buffer, indx, 'CH',ierr,vals, cval)
          call casefold(cval)
          if( index(cval,'GLORG').gt.0 ) then
*             This name can be overwritten latter when the glorg command
*             if read
              glr_sol_file = 'GLORG.SOL'
          end if
 
          return
 
***** OUT_APR : Get the name for output of final site and source positions
 3500 continue
          call read_line( buffer, indx, 'CH',ierr,vals, apr_table_file)
          call wild_card( apr_table_file, list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
          call wild_date( apr_table_file, wild_mjd, 'N' )
          return
 
***** BAK_OPTS: Get back solution options
 3600 continue
          call decode_prt_opt(buffer, indx, bak_opts)
C         call read_line( buffer, indx, 'I4',ierr,bak_opts,cval)
          return
 
***** FCN_PER : Get FCN period
 3700 continue
          call read_line (buffer, indx, 'R4',ierr,nut_period, cval)
          return
 
***** DESCRIPT : Description for solution
 3800 continue
          read(buffer(indx:),'(a)', iostat=ierr) gdescription
          call report_error('IOSTAT',ierr,'read',buffer,0,
     .                      'PROCESS_GLOBK_COM')
          return
 
***** CRT_OPT  : OPtions for output to crt
 3900 continue
          call decode_prt_opt(buffer, indx, crt_opts)
          return
 
 
***** PRT_OPT  : OPtions for output to crt
 4000 continue
          call decode_prt_opt(buffer, indx, prt_opts)
          return
 
***** USE_SITE : Select sites to be used in the solution
 4100 continue
          call casefold(buffer)
          call decode_option( buffer, gsite_names, gnum_sites,
     .                        guse_site, -1)
          return
 
***** USE_SOUR : Select sources to be used in the solution
 4200 continue
          call casefold(buffer)
          call decode_option( buffer, gsource_names, gnum_sources,
     .                        guse_source, -1)
          return

***** APR_UANG : Apriori sigma for diur/semi ut1 estimates.
 4300 continue
          call sub_char( buffer(indx:),'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread(buffer,indx, 'R4',ierr, apr_eor_ut1, cval,4)
          return

***** APR_XYAN : Apriori sigma for diur/semi xy estimates.
 4400 continue
          call sub_char( buffer(indx:),'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread(buffer,indx, 'R4',ierr, apr_eor_xy, cval,6)
          return

***** APR_ETD  : Apriori sigma for diur/semi etd estimates.
 4500 continue
          call sub_char( buffer(indx:),'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')

*         Get the site name
          call get_cmd(buffer, gsite_names, gnum_sites, is, indx )
          call multiread(buffer,indx, 'R8',ierr, vals, cval,12 )
          if( is.eq.999999 ) then
              do i = 1, gnum_sites
                  do j = 1,12
                      apr_eor_etd(j,i) = vals(j)
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = 1,12
                 apr_eor_etd(j,is) = vals(j)
             end do
          end if
          return

***** APR_SVS : apriori sigmas for the SV orbits
 4600 continue
          call sub_char( buffer(indx:),'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')

*         Get the svs name
          call get_cmd(buffer, gsvs_names, gnum_svs, is, indx )
*         Read each type separately and check if we actually used
*         parameter (type_svs_elem set in glinit)
* MOD TAH 981020: Made input more flexible:  Rad arguments may have
*         R to set them all.

*         Work along a line see what we have 
          do j = 1, max_svs_elem
             vals(j) = 0.d0
          end do
          j = 0
          done = .false.
          do while ( j.lt. max_svs_elem .and. .not.done )    
             call GetWord( buffer, word, indx )
             call casefold(word) 

*            See if we have reached end of line
             if( trimlen(word).eq.0 ) then
                 done = .true.
             else

*               See what we have:
                call check_num( word, jerr )
                jndx = 1 

                if( jerr.eq.0 ) then   ! Value is numeric
                    j = j + 1
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                else  if ( index(word,'R').gt.0 ) then
*                   Remaining radiation parameters given, get the value
                    call sub_char(word,'R',' ')
                    j = j + 1 
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                    do k = j, max_svs_elem-3
                         vals(k) = vals(j)
                    end do
                    j = max_svs_elem-3
                else if ( index(word,'A').gt.0 ) then
*                   Remaining radiation parameters given, get the value
                    call sub_char(word,'A',' ')
                    j = j + 1 
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                    do k = j,max_svs_elem 
                      vals(k) = vals(j)
                    end do
                    done = .true.
                else 
                    write(*,4610) word(1:max(1,trimlen(word))),
     .                            'APR_SVS'
 4610               format('**WARNING** Error decoding ',a,
     .                     ' in ',a,' command')
                end if
             end if
          end do

*         Now check all values are needed                  
          do j = 1, max_svs_elem
              if ( .not.kbit(type_svs_elem,j) ) then
                   vals(j) = 0.d0
              end if
          end do

*         Now copy the values to the correct satellite entries
          if( is.eq.999999 ) then
              do i = 1, gnum_svs
                  do j = 1,max_svs_elem
                      apr_svs(j,i) = vals(j)
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = 1,max_svs_elem
                 apr_svs(j,is) = vals(j) 
             end do
          end if
          return

***** MAR_UANG : mariori sigma for diur/semi ut1 estimates.
 4700 continue
          call sub_char( buffer(indx:),'F', '0')
          call sub_char( buffer(indx:), 'f', '0')
          call multiread(buffer,indx, 'R4',ierr, mar_eor_ut1, cval,4)
          return

***** MAR_XYAN : mariori sigma for diur/semi xy estimates.
 4800 continue
          call sub_char( buffer(indx:),'F', '0')
          call sub_char( buffer(indx:), 'f', '0')
          call multiread(buffer,indx, 'R4',ierr, mar_eor_xy, cval,6)
          return

***** MAR_ETD  : mariori sigma for diur/semi etd estimates.
 4900 continue
          call sub_char( buffer(indx:),'F', '0')
          call sub_char( buffer(indx:), 'f', '0')

*         Get the site name
          call get_cmd(buffer, gsite_names, gnum_sites, is, indx )
          call multiread(buffer,indx, 'R8',ierr, vals, cval,12 )
          if( is.eq.999999 ) then
              do i = 1, gnum_sites
                  do j = 1,12
                      mar_eor_etd(j,i) = vals(j)
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = 1,12
                 mar_eor_etd(j,is) = vals(j)
             end do
          end if
          return

***** MAR_SVS : mariori sigmas for the SV orbits
 5000 continue
          call sub_char( buffer(indx:),'F', '0')
          call sub_char( buffer(indx:), 'f', '0')

*         Get the svs name
          call get_cmd(buffer, gsvs_names, gnum_svs, is, indx )

* MOD TAH 981020: Made input more flexible:  Rad arguments may have
*         R to set them all.

*         Work along a line see what we have 
          do j = 1, max_svs_elem
             vals(j) = 0.d0
          end do
          j = 0
          done = .false.
          do while ( j.lt. max_svs_elem .and. .not.done )    
             call GetWord( buffer, word, indx )
             call casefold(word)

*            See if we are end of line
             if ( trimlen(word).eq.0 ) then
                done = .true.
             else
*               See what we have:
                call check_num( word, jerr )
                jndx = 1 

                if( jerr.eq.0 ) then   ! Value is numeric
                    j = j + 1
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                else if ( word(1:1).eq.'!' .or. word(1:1).eq.'#' ) then
                    done = .true.
                else if ( index(word,'R').gt.0 ) then
*                   Remaining radiation parameters given, get the value
                    call sub_char(word,'R',' ')
                    j = j + 1 
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                    do k = j, max_svs_elem-3
                         vals(k) = vals(j)
                    end do
                    j = max_svs_elem-3
                   else if ( index(word,'A').gt.0 ) then
*                   Remaining radiation parameters given, get the value
                    call sub_char(word,'A',' ')
                    j = j + 1 
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                    do k = j,max_svs_elem 
                      vals(k) = vals(j)
                    end do
                    done = .true.
                else 
                    write(*,4610) word(1:max(1,trimlen(word))),
     .                            'MAR_SVS'
                end if
             endif
          end do

*         Now check all values are needed:
          do j = 1, max_svs_elem
              if( .not.kbit(type_svs_elem,j) ) then
                  vals(j) = 0.d0
              end if
          end do

****      Now copy values to correct entries
          if( is.eq.999999 ) then
              do i = 1, gnum_svs
                  do j = 1,max_svs_elem
                      mar_svs(j,i) = vals(j)
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = 1,max_svs_elem
                 mar_svs(j,is) = vals(j)
             end do
          end if
          return

****  APR_UCOE : Aprioro sigmas for the UT1 diurn/semi coefficients
 5100 continue

*         Get the period first.
          call get_cmd(buffer, coeff_periods, max_coeff_periods,
     .                     num_period, indx )

*         Get the coeffiencent type
          call get_cmd(buffer, ut1_names, max_ut1_names,
     .                     num_name, indx )

*         Compute the coeffiecient number and read the values
          num_term = (num_period-1)*max_ut1_names + num_name
          call multiread(buffer,indx, 'R8',ierr, vals, cval,2 )
          apr_ut1_coeff(1,num_term) = vals(1)
          apr_ut1_coeff(2,num_term) = vals(2)

          return

***** APR_XYCO : Aprioro sigmas for the xy  diurn/semi coefficients
 5200 continue
*         Get the period first.
          call get_cmd(buffer, coeff_periods, max_coeff_periods,
     .                     num_period, indx )

*         Get the coeffiencent type
          call get_cmd(buffer, xy_names, max_xy_names,
     .                     num_name, indx )

*         Compute the coeffiecient number and read the values
          num_term = (num_period-1)*max_xy_names + num_name
          call multiread(buffer,indx, 'R8',ierr, vals, cval,2 )
          apr_xy_coeff(1,num_term) = vals(1)
          apr_xy_coeff(2,num_term) = vals(2)

          return

***** APR_ECOE : Aprioro sigmas for the etd diurn/semi coefficients
 5300 continue

*         Get the period first.
          call get_cmd(buffer, coeff_periods, max_coeff_periods,
     .                     num_period, indx )

*         Get the coeffiencent type
          call get_cmd(buffer, etd_names, max_etd_names,
     .                     num_name, indx )

*         Now get the site number
          call get_cmd(buffer, gsite_names, gnum_sites,
     .                     num_etd_site, indx )

*         Now check if value is in range.  If it is get values, otherwize
*         warn user and ignore
          OK = .false.
          if( num_etd_site.lt. max_etd_sites ) OK = .true.
          if( num_etd_site.eq.999999 .and. 
     .        gnum_sites.le.max_etd_sites ) OK = .true.

*         if we have too many sites then just set the maximum number
*         allowed and print warning.
          gnum_temp = gnum_sites
          if( num_etd_site.eq.999999 .and.
     .        gnum_sites.gt.max_etd_sites ) then
              gnum_temp = max_etd_sites
              OK = .true.
          end if

          if( OK ) then

*             Compute the coeffiecient number and read the values
              num_term = (num_period-1)*max_etd_names + num_name
              call multiread(buffer,indx, 'R8',ierr, vals, cval,2 )

              if( num_etd_site.eq.999999 ) then
                  do i = 1, gnum_temp 
                      apr_etd_coeff(1,num_term,i) = vals(1)
                      apr_etd_coeff(2,num_term,i) = vals(2)
                  end do
              else
                  apr_etd_coeff(1,num_term,num_etd_site) = vals(1)
                  apr_etd_coeff(2,num_term,num_etd_site) = vals(2)
              end if
          else

*             Report that there is a problem  
              write(*,'(2a,/,a,i3,a,/,a,/,a,a8,a)')
     .     ' *** WARNING *** Too many global sites for '
     .                , 'site dependent coefficients'
     .     ,' Limit is ',max_etd_sites,' sites. '
     .     ,' If you want global estimate, do not use ALL for site name'
     .     ,' Use name of site 1 instead (',gsite_names(1),')'
c           above replaced the following to avoid DEC OSF4 compiler warning
c               write(*,'('' *** WARNING *** Too many global sites for'',
c     .              '' site dependent coefficients'',/,
c     .              '' Limit is '',i3,'' sites.  If you want global'',
c     .              '' estimate, do not use ALL for site name,'',/,
c     .              '' use name of site 1 instead ('',a8,'')'')')
c     .              max_etd_sites, gsite_names(1)
          end if
          return

***** SVS_FILE : Name of the satellite emphemeris file.  Only read if
*     make_svs command have not been used.
 5400 continue 
          if( .not.make_svs_file ) then
              call read_line(buffer, indx,'CH',ierr,vals, glb_svs_file)
              call wild_card( glb_svs_file, list_file)  
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
              call wild_date( glb_svs_file, wild_mjd, 'N' )
            
              call report_error('READ_LINE',ierr,'read',buffer,0,
     .                          'PROCESS_GLB_COMMAND')
          end if
          return

***** GLB_TIDES : Indicates the use of glolbal tides
 5500 continue
          call GetWord( buffer, cval,indx )
          call casefold(cval)
          if( cval(1:1).eq.'Y' ) then
              glb_glb_tides = .true.
          else
              glb_glb_tides = .false.
          end if
          return

***** IN_NUT   : Name of the nutation series coefficient file
 5600 continue
          call read_line(buffer, indx,'CH',ierr,vals, nut_inp_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( nut_inp_file, list_file)
          return

***** IN_PLAN  : Name of the planetary series coefficient file
 5700 continue
          call read_line(buffer, indx,'CH',ierr,vals, plan_inp_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( plan_inp_file, list_file)
          return

***** IN_PMU   : Name of the polar motion UT1 tabular file
 5800 continue
          call read_line(buffer, indx,'CH',ierr,vals, pmu_inp_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( pmu_inp_file, list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
          call wild_date( pmu_inp_file, wild_mjd, 'N' )
          return

***** IN_SD    : Name of the Short period rotation file    
 5900 continue
          call read_line(buffer, indx,'CH',ierr,vals, sd_inp_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( sd_inp_file, list_file)
          return

****  APR_NEU : Get apriori sigmas for site NEU and rates
 6000 continue

*         Get name of site
          call get_cmd( buffer, gsite_names, gnum_sites, is,indx)

*         Do not need a special fix type here since all XYZ will be
*         turned on by any of NEU.  Therefore just put in 0 for F or f
* MOD TAH 950116: Change to -1 and we set to zero later as with apr_site
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 6)
          call assign_ema_val4( vals, 6, apr_neu , is, gnum_sites)
          return

****  MAR_NEU : Get variances for markov site NEU and rates
 6100 continue

*         Get name of site
          jndx = indx
          call get_cmd( buffer, gsite_names, gnum_sites, is,indx)
          call multiread( buffer, indx, 'R8', ierr, vals, cval, 6)
 
* MOD TAH 030607: See if a wild card character has been passed as part
*         of the name.
          if( is.gt.0 ) then     ! Regular name
              call assign_ema_val4( vals, 6, mar_neu , is, gnum_sites)
          else                   ! See if wild card
              call GetWord(buffer,name,jndx)
              iwc = 0
              iwc = index( name,'*')
              if( iwc.eq.0 ) iwc = index(name,'@')
              if( iwc.eq.0 ) iwc = trimlen(name)+1 ! MOD TAH 131016: Short name as wild
                                    ! as well
              if( iwc.gt.0 ) then
*                 OK: Found wild card, match all names up to this point
                  call casefold(name)
                  do i = 1, gnum_sites
                     if( name(1:iwc-1).eq.gsite_names(i)(1:iwc-1)) then
                         call assign_ema_val4( vals, 6, mar_neu , 
     .                                         i, gnum_sites)
                         if( 1.eq.2 )
     .                   write(*,6110) name, i, gsite_names(i)
 6110                    format('+ WILD CARD: Matching ',a8,' to ',
     .                          'site ',i4,1x,a8)
                     end if
                  end do
              end if
          end if

          return

***** SVS_MARF : Name of the file containing the markov statistics
*         for satellite orbits.
 6200 continue

*         Get file name
          call read_line(buffer, indx,'CH',ierr, vals, svs_mar_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
*                                             ! If error use default
          call wild_card( svs_mar_file, list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
          call wild_date( svs_mar_file, wild_mjd, 'N' )

          if( ierr.ne.0 ) svs_mar_file = ' '
          return

***** APR_RAD : apriori sigmas for the Satellite Radiation 
*         parameters only.
 6300 continue
          call sub_char( buffer(indx:),'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')

*         Get the svs name
          call get_cmd(buffer, gsvs_names, gnum_svs, is, indx )
*         Read each type separately and check if we actually used
*         parameter (type_svs_elem set in glinit)
          done = .false.
          j = 6
          do while ( j.lt. max_svs_elem-3 .and. .not.done )    
             call GetWord( buffer, word, indx )
             call casefold(word)

*            See if we are end of line
             if ( trimlen(word).eq.0 ) then
                done = .true.
             else
*               See what we have:
                call check_num( word, jerr )
                jndx = 1 

                if( jerr.eq.0 ) then   ! Value is numeric
                    j = j + 1
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                else  if ( index(word,'R').gt.0 ) then
*                   Remaining radiation parameters given, get the value
                    call sub_char(word,'R',' ')
                    j = j + 1 
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                    do k = j, max_svs_elem-3
                         vals(k) = vals(j)
                    end do
                    j = max_svs_elem-3
                else 
                    write(*,4610) word(1:max(1,trimlen(word))),
     .                            'APR_RAD'
                end if
             endif
          end do
* MOD TAH 190610: Copied test from apr_svs to see values are needed.
*         Now check all values are needed                  
          do j = 1, max_svs_elem
              if ( .not.kbit(type_svs_elem,j) ) then
                   vals(j) = 0.d0
              end if
          end do

****      Now save the values
          if( is.eq.999999 ) then
              do i = 1, gnum_svs
                  do j = 7,max_svs_elem - 3
                      apr_svs(j,i) = vals(j)
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = 7,max_svs_elem - 3
                 apr_svs(j,is) = vals(j) 
             end do
          end if
          return

***** MAR_RAD : markov Process noise for the satellite
*         radition parameters only
 6400 continue
          call sub_char( buffer(indx:),'F', '0')
          call sub_char( buffer(indx:), 'f', '0')

*         Get the svs name
          call get_cmd(buffer, gsvs_names, gnum_svs, is, indx )
*         Read each type separately and check if we actually used
*         parameter (type_svs_elem set in glinit)
          done = .false.
          j = 6
          do while ( j.lt. max_svs_elem-3 .and. .not.done )    
             call GetWord( buffer, word, indx )
             call casefold(word)

*            See if we are end of line
             if ( trimlen(word).eq.0 ) then
                done = .true.
             else
*               See what we have:
                call check_num( word, jerr )
                jndx = 1 

                if( jerr.eq.0 ) then   ! Value is numeric
                    j = j + 1
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                else  if ( index(word,'R').gt.0 ) then
*                   Remaining radiation parameters given, get the value
                    call sub_char(word,'R',' ')
                    j = j + 1 
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                    do k = j, max_svs_elem-3
                         vals(k) = vals(j)
                    end do
                    j = max_svs_elem-3
                else 
                    write(*,4610) word(1:max(1,trimlen(word))),
     .                            'MAR_RAD'
                end if
             endif
          end do

* MOD TAH 190610: Copied test from apr_svs to see values are needed.
*         Now check all values are needed                  
          do j = 1, max_svs_elem
              if ( .not.kbit(type_svs_elem,j) ) then
                   vals(j) = 0.d0
              end if
          end do

*         Now save the values 
          if( is.eq.999999 ) then
              do i = 1, gnum_svs
                  do j = 7,max_svs_elem - 3
                      mar_svs(j,i) = vals(j)
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = 7,max_svs_elem - 3
                 mar_svs(j,is) = vals(j)
             end do
          end if
*
*     MAX_CHII -- Sets the maxiumum chi**2 increment allowed
*     in the solutions, and optionally the max coordinate difference
*     and max rotation allowed allowed before removed apriori.
*     (Upto 3 arguments can be passed, m and mas for last two).
 6500 continue

          call read_line(buffer, indx,'R4',ierr,max_chi_inc, cval)
          call read_line(buffer, indx,'R8',ierr,vals, cval)
          if( ierr.eq.0 ) max_prefit_diff = vals(1)
          call read_line(buffer, indx,'R8',ierr,vals, cval)
          if( ierr.eq.0 ) max_eop_rot = vals(1)
          
          return
          
*     USE_POS -- Sets the use of sites by latitude and longitude
*     range.  Form of command is:
*     use_pos <+/-> <LL Lat> <LL Long> <UR Lat> <UR Long>
*     where LL is lower left, and UR is upper right cornors
*     of box.  Maybe followed by use_site command to put sites 
*     back.
 6600 continue

          call edit_use_pos(buffer, indx)
          return

*     USE_NUM -- Sets the use site array based on the number of
*     times a site has been used.
 6700 continue
          call read_line(buffer, indx,'I4',ierr,num_use, cval)
          do i = 1, gnum_sites 
             if( times_used(i).lt.num_use ) then
                 call sbit(guse_site, i, 0 )
             end if
          end do
          return

*     APR_TRAN -- Sets the translation and rate of change apriori
*     sigmas (m in XYZ and m/yr XYZ)
 6800 continue
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread(buffer, indx,'R4',ierr,apr_tran, cval,6)
          return

*     MAR_TRAN -- Sets the translation and rate of change markov
*     variances (m**2/yr in XYZ and (m/yr)**2/yr XYZ)
 6900 continue
          call sub_char( buffer(indx:), 'F', '0')
          call sub_char( buffer(indx:), 'f', '0')
          call multiread(buffer, indx,'R4',ierr,mar_tran, cval,6)

*         See if the optional arguments have been passed (i.e. White
*         noise process
          do i = 1, 3
             call read_line(buffer, indx,'R4',ierr,mar_tran(i,3),cval)
             if( ierr.ne.0 ) mar_tran(i,3) = 0.d0
          end do
          return

*     APR_SCAL -- Sets the scale and scale rate of change apriori
*     sigmas (ppb and ppb/yr)
 7000 continue
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread(buffer, indx,'R4',ierr,apr_scale, cval,2)
          return
         
*     MAR_SCAL -- Sets the scale and scale rate of change markov
*     variances (ppb**2/yr and (ppb/yr)**2/yr)
 7100 continue
          mar_scale(3) = 0 
          call sub_char( buffer(indx:), 'F', '0')
          call sub_char( buffer(indx:), 'f', '0')
          call multiread(buffer, indx,'R4',ierr,mar_scale, cval,2)
* MOD TAH 130716: See if white noise variarance has been passed
          call read_line(buffer,indx,'R4',ierr,val4, cval)
          if( ierr.eq.0 ) mar_scale(3) = val4
          return
         
***** ORG_OPT  : OPtions for glorg run     
 7200 continue
          call decode_prt_opt(buffer, indx, org_opts)
          return
 
****  ORG_CMD : Get command files name
 7300 continue
          call read_line(buffer, indx,'CH',ierr, vals, glr_cmd_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
*                                             ! If error use default
          call wild_card( glr_cmd_file, list_file)
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
          call wild_date( glr_cmd_file, wild_mjd, 'N' )
          if( ierr.ne.0 ) glr_cmd_file = ' '
          return

****  ORG_OUT : Output file name for glorg run
 7400 continue
          call read_line(buffer, indx,'CH',ierr, vals, glb_org_file)
          call report_error('READ_LINE',ierr,'read',buffer,0,
     .                      'PROCESS_GLB_COMMAND')
          call wild_card( glb_org_file, list_file)
*                                             ! If error use default
* MOD TAH 190624: Added wild_date option (N will stop name being incremented).
          call wild_date( glb_org_file, wild_mjd, 'N' )
          if( ierr.ne.0 ) glb_org_file = ' '
          return


***** NO_DIRCP: Stops the direct copy of the covariance matrices 
 7500 continue

          no_direct_copy = .true.

          return

***** SOURCE: This command is executed in the main reading routine and
*     we should never end up here.  Print warning if we do.
 7600 continue
          write(*,7605) buffer(1:40)
 7605     format('***WARNING*** SOURCE command interpretation error',/,
     .           'First 40 characters of Command line ',a40)
          return

***** RAD_RESE: Tells program if radiation parameters should be reset
*     with each new arc.  The default is YES.
 7700 continue
          call GetWord( buffer, cval,indx )
          call casefold(cval)
          if( cval(1:1).ne.'N' ) then
              rad_reset = .true.
          else
              rad_reset = .false.
          end if
          return

***** APR_SVAN: Gets the apriori sigmas for the satellite antenna
*     offsets (may also be given in the apr_svs command)
 7800 continue

          call sub_char( buffer(indx:),'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')

*         Get the svs name
          call get_cmd(buffer, gsvs_names, gnum_svs, is, indx ) 

*         Read each type separately and check if we actually used
*         parameter (type_svs_elem set in glinit)
          done = .false.
          j = max_svs_elem-3
          do while ( j.lt. max_svs_elem .and. .not.done )    
             call GetWord( buffer, word, indx )
             call casefold(word)

*            See if we are end of line
             if ( trimlen(word).eq.0 ) then
                done = .true.
             else
*               See what we have:
                call check_num( word, jerr )
                jndx = 1 

                if( jerr.eq.0 ) then   ! Value is numeric
                    j = j + 1
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                else  if ( index(word,'A').gt.0 ) then
*                   Remaining radiation parameters given, get the value
                    call sub_char(word,'A',' ')
                    j = j + 1 
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                    do k = j, max_svs_elem
                         vals(k) = vals(j)
                    end do
                else 
                    write(*,4610) word(1:max(1,trimlen(word))),
     .                            'APR_SVAN'
                end if
             endif
          end do

*         Now save the values 
          if( is.eq.999999 ) then
              do i = 1, gnum_svs
                  do j = max_svs_elem-2, max_svs_elem
                      apr_svs(j,i) = vals(j)
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = max_svs_elem-2, max_svs_elem
                 apr_svs(j,is) = vals(j)
             end do
          end if

          return

***** MAR_SVAN: Sets the Markov process for the offsets.  The process
*     noise is random walk
 7900 continue

          call sub_char( buffer(indx:),'F', '0')
          call sub_char( buffer(indx:), 'f', '0')

*         Get the svs name
          call get_cmd(buffer, gsvs_names, gnum_svs, is, indx ) 

*         Read each type separately and check if we actually used
*         parameter (type_svs_elem set in glinit)
          done = .false.
          j = max_svs_elem-3
          do while ( j.lt. max_svs_elem .and. .not.done )    
             call GetWord( buffer, word, indx )
             call casefold(word)

*            See if we are end of line
             if ( trimlen(word).eq.0 ) then
                done = .true.
             else
*               See what we have:
                call check_num( word, jerr )
                jndx = 1 

                if( jerr.eq.0 ) then   ! Value is numeric
                    j = j + 1
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                else  if ( index(word,'A').gt.0 ) then
*                   Remaining radiation parameters given, get the value
                    call sub_char(word,'A',' ')
                    j = j + 1 
                    call read_line(word,jndx,'R8',jerr,vals(j), cval)
                    do k = j, max_svs_elem
                         vals(k) = vals(j)
                    end do
                else 
                    write(*,4610) word(1:max(1,trimlen(word))),
     .                            'MAR_SVAN'
                end if
             endif
          end do

*         Now save the values 
          if( is.eq.999999 ) then
              do i = 1, gnum_svs
                  do j = max_svs_elem-2, max_svs_elem
                      mar_svs(j,i) = vals(j)
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = max_svs_elem-2, max_svs_elem
                 mar_svs(j,is) = vals(j)
             end do
          end if


          return 

***** VAL_SVAN: Sets the values of the antenna offsets by satellite
*     PRN number.  These values also appear in the ephemeris file.  The
*     values given over-write those in the ephemeris file.
 8000 continue

*         Get the svs name
          call get_cmd(buffer, gsvs_names, gnum_svs, is, indx )

*         Now read the values from the line
          do j = 1, 3
              call GetWord(buffer, word, indx)
              call check_num(word, jerr)
              jndx = 1
              if ( jerr.eq.0 ) then
                 call read_line(word,jndx,'R8',jerr,vals(j+3), cval)
              end if
              if( jerr.eq.0 ) then
                 vals(j) = vals(j+3)
              else
                 vals(j) = -999.d0
              endif
          end do

*         Now copy to final locations
          if( is.eq.999999 ) then
              do i = 1, gnum_svs
                  do j = 1,3
                      if( vals(j).ne. -999.d0 ) then
                          apr_val_ant(j,i) = vals(j)
                      end if
                  end do
             end do
          else if ( is.gt.0 ) then
             do j = 1,3
                 if( vals(j).ne. -999.d0 ) then
                     apr_val_ant(j,is) = vals(j)
                 end if
             end do
          end if

          return

***** MUL_PMU : Sets parameters for the multiple day estimates of 
*     multiday pmu parmeters
 8100 continue 

*         Decode the options; Number of values per block
          mul_pmu_opt = 0
          call read_line(buffer,indx,'I4',jerr,num_mul_pmu, cval)
*         Spacing of values (days)
          call read_line(buffer,indx,'R8',jerr,spacing_mul_pmu, cval)
*         Start Epoch: See if type start-epoch passed
          call GetWord(buffer, word, indx)
          call casefold(word)
          call check_num( word, jerr )
          if ( jerr.eq.0 ) then   ! Date passed
               jndx = 1
               call read_line(word, jndx,'I4',jerr,date(1), cval)
               call multiread(buffer, indx,'I4',jerr,date(2), cval,4)
               sectag = 0.d0
               call ymdhms_to_jd(date, sectag, start_mul_pmu)
               call sbit(mul_pmu_opt,2,1)
          else if ( index(word,'GPSW').gt.0 ) then
*              Start at gpsw of first experiment +- the following
*              number of days
               call sbit(mul_pmu_opt,3,1)
               call read_line(buffer,indx,'R8',jerr,start_mul_pmu,cval)
          else if( index(word,'START').gt.0 ) then
*              MOD TAH 120411: Uses the epoch of the first experiment in
*              list to be processed (gepoch_start)
               call read_line(buffer,indx,'R8',jerr,vals,cval)
               start_mul_pmu = gepoch_start + vals(1)  ! Offset from start 
                                    ! offset is normally 0.5
               call sbit(mul_pmu_opt,2,1)
          end if

*         See if anything else on line
          jerr = 0
          do while ( jerr.eq.0 )
              call GetWord(buffer, word, indx)
              call casefold(word)
              if( index(word,'IND').gt.0 ) then
*                 Set the estimates to treated independently
                  call sbit(mul_pmu_opt,1,1)
              else if( index(word,'NOUT').gt.0 ) then 
                  call sbit(mul_pmu_opt,4,1)
              else if( index(word,'WARN').gt.0 ) then
                  call sbit(mul_pmu_opt,5,1)
              else if( index(word,'PROP').gt.0 ) then
                  call sbit(mul_pmu_opt,6,1)
              end if
              if( trimlen(word).eq.0 ) jerr = -1
          end do
          return

***** APP_PTID: Sets parameters for applying the pole tide to hfile
*     results
 8200 continue

*         Read the list of h-file names from the command buffer
          jerr = 0
          num_ptide_hfs = 0
          ptide_opt = 1
          mean_pole_def = 'IERS10'   ! Default mean pole model.
         
          do while ( jerr.eq.0 )
             eol  = 0
             call GetWord(buffer, word, indx)
* MOD TAH 140405: See if +-OPT option passed as first argument
* MOD TAH 200220: See if Mean pole specified.  Default is IERS10
*         with IERS20 new option for ITRF2020.
             if( trimlen(word).eq.0 ) jerr = -1
             if( num_ptide_hfs.eq.0 ) then  ! Check OPT
                 cval = word
                 call casefold(cval)
                 if( index(cval,'OPT').eq.1 .or. 
     .               index(cval,'OPT').eq.2 .and.
     .               trimlen(cval).le.4 ) then
                     if( cval(1:1).eq.'+' .or. cval(1:1).eq.'O') then
                         call sbit(ptide_opt,2,1) 
                     elseif( cval(1:1).eq.'-' ) then 
* MOD TAH 200220: Changed remove bit to 4 from 3.  3 now used for SE poleitde
                         call sbit(ptide_opt,4,1)
                     endif
                     eol = 1  ! Set so that option will not
                              ! interpretted as hfile code.
                 endif
* MOD TAH 202020: See if Solid-Earth pole should be removed or added if
*                not already applied
                 if( index(cval,'SEPT').eq.1 .or. 
     .               index(cval,'SEPT').eq.2 .and.
     .               trimlen(cval).le.5 ) then
                     if( cval(1:1).eq.'+' .or. cval(1:1).eq.'S') then
                         call sbit(ptide_opt,1,1) 
                     elseif( cval(1:1).eq.'-' ) then 
* MOD TAH 200220: Set remove bit to bit 3 for SOLID-Earth pole tide.
                         call sbit(ptide_opt,3,1)
                     endif
                     eol = 1  ! Set so that option will not
                              ! interpretted as hfile code.
                 endif
* MOD TAH 200220: Check for IERSxx models. Name must be correct length
                 if( index(cval,'IERS').eq.1 .and. 
     .              trimlen(cval).eq. 6 ) then 
                    if( cval(5:6).ne.'10' .and. cval(5:6).ne.'20' .and.
     .                  cval(5:6).ne.'96'  ) then
                        write(*,8220) cval
 8220                   format('Mean Pole ',a,' is not valid: ',
     .                         ' Ingoring')
                     else
                         mean_pole_def = cval
                     endif
                     eol = 1
                 end if  
* MOD TAH 200717: Added REP option to report changes being made.
                 if( index(cval,'REP').eq.1 .and.
     .               trimlen(cval).le.3 ) then
                     call sbit(ptide_opt,16,1) 
                     eol = 1
                 endif
                                    
             endif

             if( trimlen(word).gt.0 .and. eol.eq.0 ) then
                 if( num_ptide_hfs.ge. max_ptide_hfs ) then
                     write(*,8250) word
 8250                format('**WARNING*** To many hfile names ',
     .                      'given in APP_PTID command. ',
     .                      'Ignoring ',a)
                 else
                     num_ptide_hfs = num_ptide_hfs + 1
                     ptide_hfs(num_ptide_hfs) = word
                 endif
             elseif( eol.eq.0 ) then
                 if( num_ptide_hfs.eq.0 ) then
                     num_ptide_hfs = 1
                     ptide_hfs(num_ptide_hfs) = 'ALL'
                 end if
             endif
          end do
          if( num_ptide_hfs.eq.0 ) then
              num_ptide_hfs = 1
              ptide_hfs(num_ptide_hfs) = 'ALL'
          endif

          return

***** APR_ROT :  Rotation parameter apriori sigmas.  Values are
*     given as XYZ rotation sigma and then XYZ rotation rates.   
 8300 continue
          call sub_char( buffer(indx:), 'F', '-1')
          call sub_char( buffer(indx:), 'f', '-1')
          call multiread(buffer, indx,'R4',ierr,apr_rot, cval,6)
          return

***** MAR_ROT : Rotation markov process noise.  Values are given
*     XYZ rotation rw process, IRW process, and White noise levels.
 8400 continue
          call sub_char( buffer(indx:), 'F', '0')
          call sub_char( buffer(indx:), 'f', '0')
          call multiread(buffer, indx,'R4',ierr,mar_rot, cval,6)

*         See if the optional arguments have been passed (i.e. White
*         noise process
          do i = 1, 3
             call read_line(buffer, indx,'R4',ierr,mar_rot(i,3),cval)
             if( ierr.ne.0 ) mar_rot(i,3) = 0.d0
          end do
          return 

***** EP_TOL  : Allows user to specify tolerances for match
*     between epochs for EOP parameters and Orbital ephemerides
 8500 continue 

          call read_line(buffer, indx,'R8',ierr, vals, cval)
          if( vals(1).ne.-1.d0 ) tol_mul_pmu = vals(1)
          call read_line(buffer, indx,'R8',ierr, vals, cval)
          if( vals(1).ne.-1.d0 ) tol_svs_eph = vals(1)
          RETURN
 
***** IRW_MOD : Let's user set which Itergrated random walk 
*     model to use:  The 'OLD' option is for backward compatablity.
 8600 continue
          call getword(buffer,cval, indx)
          call casefold(cval)
          if( cval(1:1).eq.'O' ) old_irw = .true.
          RETURN

***** SIG_NEU : Command to allow noise to be added to the
*     components of individual sites.  Format
*     sig_neu  <site_name> [hf_code] <sig N> <sig E> <sig U> <Start> <End>
 8700 continue

*         Check we don't have too many renames
          if( num_ss+1.gt.max_ss ) then
              write(*,8710) max_ss,buffer(1:trimlen(buffer))
 8710         format('**WARNING** Max site_neu entries exceeded.',
     .               i5,' Max allowed',/,
     .               'Ignoring ',a)
              RETURN
          end if

*****     Start pulling entries off the lines
          num_ss = num_ss + 1
          call getword(buffer, ss_codes(num_ss),indx)
          call casefold(ss_codes(num_ss))
* MOD TAH 041118: If @ appears at the end of the name, remove
*         so that it will match names
* MOD TAH 050429: Added check that @ is at end of name (need because
*         @_GXi is treated differently.  Replaced ALL will @.
          if( ss_codes(num_ss)(1:4).eq.'ALL ' ) 
     .        ss_codes(num_ss) = '@'
          jndx = index(ss_codes(num_ss),'@')
* MOD TAH 050509: Fix: Only replace the @ with a blank if it is not
*         the first character.  This way eg., VILL@ will be replaced
*         by VILL which will match all VILL_xxx names. @ by itself
*         is explicitly handled in glfor/remove_param.f. @_XXX is
*         also explicitly handled.
          if( jndx.gt.1 .and. jndx.eq.trimlen(ss_codes(num_ss)) ) then 
              ss_codes(num_ss)(jndx:) = ' '
          endif

*         See if hfile name restriction given
          jndx = indx
          call getword(buffer, cval, indx )
*         See if number or string
          call check_num( cval, jerr )
          if( jerr.ne.0 ) then
              ss_hfiles(num_ss) = cval
          else
              ss_hfiles(num_ss) = ' '
              indx = jndx
          end if
 
****      Read the sigmas associated with this line
          do j = 1,3
              call read_line(buffer,indx,'R4',jerr,
     .            ss_sig(j,num_ss),cval)
              if( jerr.ne.0 ) ss_sig(j,num_ss) = 0.d0
          end do

*****     Now start reading the oprional parts of the command.
*         Start time for the site_sig (ymdhms)
          do j = 1,5
              call read_line(buffer,indx,'I4',ierr,date(j),cval)
              if( ierr.ne.0 ) date(j) = 1
          end do
*         If there is an error, default the start to 1900
          if( ierr.ne.0 ) date(1) = 1900
 
*         Process the starting data
          sectag = 0.0d0
          call ymdhms_to_jd(date, sectag, ss_times(1,num_ss))
 
****      Get the end time:
          do j = 1,5
              call read_line(buffer,indx,'I4',ierr,date(j),cval)
              if( ierr.ne.0 ) date(j) = 1
          end do
*         If there is an error, default the end to 2100
          if( ierr.ne.0 ) date(1) = 2100
 
*         Process the starting data
          sectag = 0.0d0
          call ymdhms_to_jd(date, sectag, ss_times(2,num_ss))

*         Thats all
          RETURN

***** UNI_WGHT -- Option to unify weight for stations depending on 
*     the number of times they are used in different networks on 
*     the same day.  Only implemented for single-day combinations
 8800 continue
          call getword(buffer,cval, indx)
          call casefold(cval)
          if( cval(1:1).eq.'Y' ) then
              uni_wght = .true. 
          else
              uni_wght = .false.
          endif
          RETURN

***** DECIMATE -- Select step size and offset
 8900 continue
          call read_line(buffer,indx,'I4',ierr,decnum, cval)
          if( decnum.le.0 ) decnum = 1
          call read_line(buffer,indx,'I4',ierr,decoff, cval)
          if( decoff.le.0 .or. ierr.ne.0 ) decoff = 1

          return

***** GLB_OPT -- Decodes the options for saving the combined global
*         files.
 9000 continue
          call decode_option( buffer, glb_opt_types, 3, glb_out_opt, 0)
          return

****  APR_ATM -- Sets apriori sigmas for atmospheric delay estimates
 9100 continue

*         Get name of site.  This command allows the name COMMON which
*         will only set the apriori for sites that are used more than
*         once.  Get this option first
          jndx = indx
          call GetWord(buffer,cval, jndx)
          call casefold(cval)
          if( cval(1:7).eq.'COMMON ' ) then
*            OK COMMON form used.  Get the apriori sigma and assign
             call read_line( buffer, jndx, 'R8', ierr, vals, cval)
*            Now find the sites used more than once and assign values
             do i = 1, gnum_sites
           	if( times_used(i).gt.1 ) then
           	    call assign_ema_val4( vals, 1, apr_atm , i, 
     .     					    gnum_sites)
           	endif
             end do
          else      ! Looks like a regular use of the command
             call get_cmd( buffer, gsite_names, gnum_sites, is,indx)

             call read_line( buffer, indx, 'R8', ierr, vals, cval)
             call assign_ema_val4( vals, 1, apr_atm , is, gnum_sites)
          endif

          return


***** MAR_ATM -- Set Markov process noise for atmospheric delay 
*     estimates
 9200 continue

*         Get name of site.  This command allows the name COMMON which
*         will only set the apriori for sites that are used more than
*         once.  Get this option first
          jndx = indx
          call GetWord(buffer,cval, jndx)
          call casefold(cval)
          if( cval(1:7).eq.'COMMON ' ) then
*            OK COMMON form used.  Get the apriori sigma and assign
             call read_line( buffer, jndx, 'R8', ierr, vals, cval)
*            Now find the sites used more than once and assign values
             do i = 1, gnum_sites
           	if( times_used(i).gt.1 ) then
           	    call assign_ema_val4( vals, 1, mar_atm , i, 
     .     					    gnum_sites)
           	endif
             end do
          else      ! Looks like a regular use of the command
             call get_cmd( buffer, gsite_names, gnum_sites, is,indx)

             call read_line( buffer, indx, 'R8', ierr, vals, cval)
             call assign_ema_val4( vals, 1, mar_atm , is, gnum_sites)
          endif


          return

***** DEL_SCRA -- Sets option to delete scratch files, (com, sol and srt)
*         when globk completes.  NOTES: glorg can not be run afterwards
*         if this option is set.  Command must be used as DEL_SCRA Y
 9300 continue
          call getword(buffer,cval, indx)
          call casefold(cval)
          if( cval(1:1).eq.'Y' ) then
              del_scratch = .true. 
          else
              del_scratch = .false.
          endif
          RETURN

***** FREE_LOG -- Loosens the constraints on earthquake log tersm that
*     have been previously estimates
 9400 continue
*         Since earthquakes have read, simply replace the constant 
*         log constraints
*         See if constraint has been passed. Default is 10cm
          decon = 0.1d0
          call read_line(buffer, indx, 'R8', ierr, vals, cval)
          if( ierr.eq.0 ) decon = vals(1)

          do i = 1, num_eq
             do j = 1,3
                if( eq_log_sig(j,i).gt.0 .or. 
     .              eq_log_sig(j+3,i).gt.0) then
*                   Set the value of the log so that site position sigma is
*                   no more that specified amount
                    if( decon.gt.0 ) then 
                       logsig = decon
                       dtl = (gepoch_start-eq_epoch(i))/eq_log_tau(i)
                       if( dtl.gt.0 ) then
                           dpdt = log(1+dtl)
                           if( dpdt.gt.1.d0 ) then
                              logsig = decon/dpdt
                           endif
                       endif
                    else
                       logsig = 0.d0
                    endif
                    eq_log_sig(j,i) = logsig   ! Set sigma to +-0.1 m
                    eq_log_sig(j+3,i) = 0.d0   ! Set distance term to zero
                end if
             end do
          end do

          RETURN

***** APP_MODL: Command to apply and remove loading models in solution
 9500 continue

*         See what models are requested: + to apply, - to remove if already
*         applied.  Options are ATML and HYDR
          call sbit(appload_mod,1,1)   ! Set bit to show command used
          call decode_option(buffer, load_types, 3, appload_mod, 1)
          RETURN

         
***** Place for ALL option (iel = 999999) Not valid in this case
*     THIS CODE SHOULD BE ALWAYS LAST IN THE SUBROUTINE

 9600 continue
 
          write(crt_unit,9610)
 9610     format(' ALL option not valid for as a command or',
     .           ' use of dummy command')
          return
 
****  Thats all
      end
 
CTITLE EDIT_USE_POS
 
      subroutine edit_use_pos( buffer, indx)

      implicit none
 
*     Routine to set guse_site based on the geophraphic
*     location of the site.
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* Passed variables
 
*   indx    - Current position in string
 
      integer*4 indx
 
 
      character*(*) buffer
 
* Local variables
 
*   ierr    - IOSTAT error
*   opt     - Set to 1 for sbit if + used, 0 is - given.
*   i       - Loop counter
 
      integer*4 ierr, opt, i
 
*   geod_pos(3) - Geodetic position, co-lat, long and height
*           - (rads and m)
*   rot_mat(3,3)    - Not used (returned by XYZ_to_GEOD)
*   llrange(4)  - LL Lat and long, and UR Lat and Long
*               - (read in degrees, converted to co-lat and
*               - rads)
*   tmp         - USed to switch latitudes if out of order
 
      real*8 geod_pos(3), rot_mat(3,3), llrange(4), tmp
 
*   ctype       - + or - character depending on it to use or
*           - not to use.
*   cdum    - Dummy character for multiread
 
      character*4 ctype, cdum
 
*   set     - Set true if point in box and should be set.
 
 
      logical set
 
***** Start see if + or - option used.
      call Getword( buffer, ctype, indx)
      if( ctype(1:1).eq.'+' ) then
          opt = 1
      else if( ctype(1:1).eq.'-' ) then
          opt = 0
      else
          call report_stat('Warning','globk','edit_use_pos',
     .        buffer,'+ or - argument expected',0)
          return
      end if
 
***** Now read the position arguments
      call multiread(buffer, indx,'R8',ierr, llrange, cdum,4)
 
****  Check and convert llrnage
      if( llrange(1).gt.llrange(3) ) then
          tmp = llrange(3)
          llrange(3) = llrange(1)
          llrange(1) = tmp
      end if
 
*     Now convert to co-latitude and radians
      llrange(1) = pi/2 - llrange(1)/rad_to_deg
      llrange(3) = pi/2 - llrange(3)/rad_to_deg
 
      llrange(2) = llrange(2)/rad_to_deg
      if( llrange(2).lt.0 ) llrange(2) = llrange(2) + 2*pi
      llrange(4) = llrange(4)/rad_to_deg
      if( llrange(4).lt.0 ) llrange(4) = llrange(4) + 2*pi

****  Now scan over sites and see if they are in this range
 
      do i = 1, gnum_sites
 
          call XYZ_to_GEOD( rot_mat, apr_val_site(1,1,i), geod_pos)
 
*         Now see if in box
          set = .false.
          if( geod_pos(1).le.llrange(1) .and.
     .        geod_pos(1).ge.llrange(3)       ) then
 
*             Inside latitude range, see if in longitude range.  Here
*             We need to be careful of Greenwich crossing
              if( llrange(2).gt. llrange(4) ) then
 
*                 We are crossing long=0 so do test in two parts
                  if( geod_pos(2).ge. llrange(2) .or.
     .                geod_pos(2).le.llrange(4)         ) then
                      set = .true.
                  end if
              else
****              Normal check
                  if( geod_pos(2).ge. llrange(2) .and.
     .                geod_pos(2).le.llrange(4)         ) then
                      set = .true.
                  end if
              end if
          end if
 
*****     See we should set the guse_site entry
          if( set ) call sbit(guse_site, i, opt)
      end do
 
***** Thats all
      return
      end
 
