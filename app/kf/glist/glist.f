      program glist

      implicit none 
 
*     This program will sort the list of input global solution
*     files into ascending or descending order; and build up the
*     list of aprori site,source positions and other information
*     about the solution input parameters.
*
*     A Table is then produced which shows the occurence of the sites
*     and sources.
*
*     This program is derived from GLINIT.
*                                  8:56 AM  TUE., 18  AUG., 1987
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
 
*   dumm       - Dummy place holder
*              - from read_line
*   i,j        - Loop counter
*   ierr       - IOSTAT error
*   indx       - Positioner for Read_line.
*              - Keeps track of where we
*              - are in solution
*   irm(5)     - 5 parameters to be returned
*   jerr       - Read_line error
*   kerr       - Decodeing error
*   most_sites, most_sources, most_svss - Most site/sources/svss 
*                used in any one experiment.
 
      integer*4 dumm, ierr, indx,  jerr, kerr, trimlen,
     .          most_sites, most_sources, most_svss
 
*   cr                              - ASCII carriage return
 
      character*1 cr
 
*   buffer                          - Line read from list file.
 
      character*80 buffer
 
*   sepoch_expts( max_glb_sol )     - Epochs of each solution.
*                                   - List will be sorted based
*                                   - on these values
*   apr_values(max_glbapr_sz)           - Storage for covariance
*                                   - matrix and solution vector

      character*(sort_recl) expt_names(max_glb_sol), out_gdl
 
      real*8 sepoch_expts( max_glb_sol ), apr_values(max_glbapr_sz),
     .       expts_var(max_glb_sol),
     .       expts_var_read, expts_diag_read, val8, glb_var

      real*8 crunjd, grunjd, sectag  ! Run time of current hfile being
              ! read (crunjd), run time of latest solution in this
              ! epoch of data (grunjd) and sectag for conversion.
              ! Values are used to see if SV antenna offsets should
              ! by updated. 

      integer*4 max_len, i

      character*8 cdum

* MOD TAH 031112: Added + character for experiments within 6-hours of
*     each other
      character*8 plus_char

* MOD TAH 110529: Added option to outout all renames (not just those 
*     used in current hfiles).  Set by adding :A to end of eq-file names
      logical report_allrn

* MOD TAH 140327: Added report_mods option
      logical report_mods 

C     common / globc_ema / expt_names, sepoch_expts, apr_values

* MOD TAH 980513: added external for globk_cmd_bd
c      external globk_cmd_bd
 
***** Start decode the runstring for this program
 
      cr = char( 13 )
      out_gdl = ' '
      report_allrn = .false.
      eq_reset = .false.
      report_mods = .true.
      guse_site(:) = -1

      call decode_glist_run(out_gdl, report_allrn) 
*     Announce program
      write(*, 50)
   50 format(/'GLIST: Contents of GLOBK Binary H-files',/,
     .        '---------------------------------------',/,
     .        ' Starting to read input data list')
      if( trimlen(out_gdl).gt.0 ) then
          write(*,60) out_gdl(1:trimlen(out_gdl))
  60      format('Sorted List written to ',a)
      end if

****  Read any earthquake file given
      call read_eq_file
* MOD TAH 100128: Added sort of earthquakes and resolving ?PS names
      call sort_eq
      call sort_unf_rn
 
***** Open the input list and loop over it
 
      open( 100, file=list_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open', list_file,1,'glist')
*                                                   ! Kill if error
 
***** Now loop over the input list, and save epochs, and apriori
*     information
 
      ierr = 0
      num_glb_sol = 0

      most_sites = 0
      most_sources = 0
      most_svss = 0

      do while ( ierr.eq.0 )
 
*                                              ! Line from file
          read(100,'(a)', iostat=ierr) buffer
*         Don't worry about reading experiment variance. Not used here.
          expts_var_read = 1.0
          expts_diag_read = 0.0
*                                             ! Get file name from
          if ( ierr.eq.0 .and. buffer(1:1).ne.'#' .and. 
     .         buffer(1:1).ne.'*' ) then
*                                             ! line
              indx = 1
              call read_line( buffer, indx, 'CH', jerr, dumm,
     .                        glb_inp_file)

* MOD TAH 991121: See if scale and diagonal passed
              call read_line( buffer, indx, 'R8', kerr, val8,
     .                        cdum)
              if( kerr.eq.0 ) expts_var_read = val8
              call read_line( buffer, indx, 'R8', kerr, val8,
     .                        cdum)
              if( kerr.eq.0 ) then
                  expts_diag_read = val8
                  expts_var_read = -(expts_var_read + 
     .                          (1.d0+expts_diag_read/1.d6)/1000.d3)
              endif
                  
 
*****         Now open this global and get the epoch information
              if( jerr.eq.0 ) then
 
                  call fmpopen(cglb_dcb,kerr,glb_inp_file,'rwo',1)
                  call report_error('FmpOpen',kerr,'open',
     .                 glb_inp_file,0,'GLSAVE')
 
                  if( kerr.eq.0 ) then
                      call rw_glb_header('R', kerr)

* MOD TAH 970610: Check to see if hfile valid.
                      if( kerr.ne.0 ) then
                         if( kerr.eq.2001 ) then
                             write(*,*) '**WRONG BYTE ORDER** Run swaph' 
                         end if
                         call report_error('RW_GLB_HEADER',kerr,'read',
     .                       glb_inp_file,0,'GLINIT') 
                      end if 
                  end if
 
* MOD TAH 131111: Compute the run-time JD for this experiment
                  sectag = crun_time(6)
                  call ymdhms_to_mjd(crun_time,sectag, crunjd )
                  
*                                         ! we read header OK, so continue
                  if( kerr.eq.0 ) then
 
                      num_glb_sol = num_glb_sol + 1
                      if( num_glb_sol.gt.max_glb_sol ) then
                          call report_stat('fatal','glist','main',
     .                         ' ','Too many experiments',max_glb_sol)
                      endif
                      expts_var(num_glb_sol) = expts_var_read

                      most_sites = max(most_sites, cnum_sites)
                      most_sources = max(most_sources, cnum_sources)
                      most_svss = max(most_svss, cnum_svs )
 
                      write(*,100) num_glb_sol, 
     .                      glb_inp_file(1:trimlen(glb_inp_file)),
     .                      cepoch_expt
  100                 format(' Global ',i4,2x,a,' Epoch ',F13.3)
 
                      call save_epoch(num_glb_sol, cepoch_expt,
     .                      glb_inp_file, expts_var_read,
     .                      sepoch_expts, expt_names, expts_var) 
 
*****                 Get the apriori name and position information and
*                     save this as well
 
                      call save_apr_inf( apr_values, crunjd, 
     .                                   grunjd )
*                     call check_apr_inf( apr_values )
                      call solution_inf

*                                     ! Global open OK
                  end if
*                                     ! Name read OK from buffer
              end if
*                                     ! List file read OK
          end if
*                                     ! Looping over the input file
      end do
 
***** Now sort experiment epochs, and create/write the type 2 file
*     with the sort listed
 
      call sort_epochs( expt_names, sepoch_expts, expts_var)
 
* MOD TAH 1001228: Now sort the earthquakes and site renames

***** Now sort the names of the site and sources

* MOD TAH 021005: Make sure we have aprioris
      call check_apr_crd
 
      call sort_apriori
 
      close ( 100  )
      
****  See if output sorted list
      if( trimlen(out_gdl).gt.0 ) then
          open(200,file=out_gdl,status='new',iostat=kerr)
          call report_error('IOSTAT',kerr,'open',out_gdl,0,'glist')
          if( kerr.eq.0 ) then
              write(200,210) list_file(1:trimlen(list_file))
 210          format('# Sorted list from ',a)
              max_len = 0
              do i = 1, num_glb_sol
                 if( trimlen(expt_names(i)).gt.max_len) then
                     max_len = trimlen(expt_names(i))
                 end if
              end do
              do i = 1, num_glb_sol
                 plus_char = ' '
                 if( i.lt. num_glb_sol ) then
                     if( sepoch_expts(i+1)-sepoch_expts(i).lt.0.25 )
     .                   plus_char = '+'
                 endif

                 if( expts_var(i).ne.1 ) then
                     if( expts_var(i).lt.0 ) then
                         glb_var = expts_var(i)
                         expts_diag_read = (abs(glb_var*1.d3) - 
     .                          int(abs(glb_var*1.d3)))*1000.d0
                         expts_var_read = abs(glb_var) - 
     .                          expts_diag_read/1000.d3
                         expts_diag_read = (expts_diag_read-1)*1000.d3
                         write(200,220) expt_names(i)(1:max_len),
     .                                  expts_var_read, expts_diag_read,
     .                                  plus_char(1:1)
 220                     format(a,1x,F10.3,1x,F6.2,1x,a1)
                     else 
                         write(200,225) expt_names(i)(1:max_len),
     .                                  expts_var(i), plus_char(1:1)
 225                     format(a,1x,F10.3,1x,a1)
                     end if
                 else
                     write(200,230) expt_names(i)(1:max_len),
     .                                  plus_char(1:1)
 230                 format(a,1x,a1)
                 end if
              end do
          end if 
          close(200)
      end if   
 
***** Produce summary 
      call sum_site_source( glb_com_file, expt_names, most_sites,
     .     most_sources, most_svss, apr_values, report_allrn, 
     .     report_mods )
 
***** Thats all
      end
 
CTITLE DECODE_GLIST_RUN
 
      subroutine decode_glist_run( out_gdl, report_allrn )
 
      implicit none 
 
*     Routine to read the GLIST runstring.  The form is
*
*     GLIST, list_file, out_put_file , <sort_file>, <sort_direction>.
*
*     where all arguments are optional except for the name of the
*     list_file which gives the list of global solutions to be
*     used in this solution.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

*   out_gdl -- Optional output gdl

      character*(*) out_gdl

      logical report_allrn  ! Set true when :A added eq-file name 
                ! and causes all renames to be reported).
 
*   DecimalToInt    - HP function for string to integer value
*   ierr            - Error flag for DecimalToInt
*   len_runstring   - Length of the runstring read
*   rcpar       - HP function for reading runstring
 
      integer*4 DecimalToInt, ierr, len_runstring, rcpar, i
      integer*4 is, js   ! Start and stop of name separators
      integer*4 trimlen
 
*   runstring   - Place to read runstring into
 
      character*64 runstring
      character*2048 eqent  ! Name of earthuake files
      character*2048 apren  ! Name of apr  files
 
****  Get the name of the list_file
 
      len_runstring = rcpar(1, list_file )
*                                         ! Error No name passed
      if( len_runstring.eq.0 ) then
          call proper_runstring('glist.hlp','glist',1)
      end if
 
***** Now try out the optional arguments
 
      len_runstring = rcpar(2, glb_com_file )
      if( len_runstring.eq.0 ) glb_com_file = '6'
 
      len_runstring = rcpar(3, runstring)
*                                         ! Get sort direction
      if( len_runstring.gt.0 ) then
          sort_direction = DecimalToInt( runstring, ierr )
          call report_error('DecimalToInt',ierr,'decod',runstring,
*                                                 ! No Kill
     .                      0,'Decode_glist_run')
      end if
 
      if( ierr.ne.0 .or. len_runstring.eq.0 ) sort_direction = +1

****  See if earthquake file passed:
      len_runstring = rcpar(4, eqent)
c MOD SCM 1/26/2010: Need to set num_eqfiles = 1, else the eq_file 
c     is not read by read_eq_file. 
* MOD TAH 1/26/2010: Fixed to allowed multiple names
      num_eqfiles = 0
      if( len_runstring.eq.0 .or. eqent(1:4).eq.'none') then
         eq_inp_file(1) = ' '
      else     ! See if should split names about
         is = 1
         js = 1
         do while ( js.lt. len_runstring+1 )
            js = js + 1
            if( eqent(js:js).eq.':' .or.
     .          eqent(js:js).eq.'+' .or.
     .          eqent(js:js).eq.'=' .or.
     .          eqent(js:js).eq.' ' ) then
               num_eqfiles = num_eqfiles + 1
               if( num_eqfiles.gt.max_eqfiles ) then
                   write(*,120) max_eqfiles, eqent(1:len_runstring)
 120               format('**FATAL** Too many eathquake files. ',
     .                    'Max allowed 'i3,/,'Names ',a)
                   stop 'GLIST: Too many earthquake files'
               endif
* MOD TAH 130221: Allow for RESET option to set names to GPS
               if( eqent(is:js-1).eq.'RESET' ) then
                  eq_reset = .true.
                  num_eqfiles = num_eqfiles - 1
               else
                  eq_inp_file(num_eqfiles) = eqent(is:js-1)
               end if
               is = js + 1
            end if
         end do
* MOD TAH 110529: Check last name
         if( eq_inp_file(num_eqfiles).eq.'A ' ) then
             report_allrn = .true.
             num_eqfiles = num_eqfiles - 1
         end if

         write(*,140) num_eqfiles,(i,eq_inp_file(i)
     .                    (1:trimlen(eq_inp_file(i))),i=1,num_eqfiles) 
 140     format('There are ',i2,' earthquake files ',/,
     .        30('# ',i2,' File : ',a,/,:))
         if( report_allrn ) write(*,'(a)') 
     .                     'Renames will be reported :A option'

      endif

****  Get name of output gdl file.
      len_runstring = rcpar(5, out_gdl)
      if( len_runstring.eq.0 ) out_gdl = ' '
      
****  See if apriori file is passed
      len_runstring = rcpar(6, apren)
      num_apr_files = 0
      if( len_runstring.gt.0 ) then
* MOD TAH 140723: See if multiple apriori files passed
         is = 1
         js = 1
         do while ( js.lt. len_runstring+1 )
            js = js + 1
            if( apren(js:js).eq.':' .or.
     .          apren(js:js).eq.'+' .or.
     .          apren(js:js).eq.'=' .or.
     .          apren(js:js).eq.' ' ) then
               num_apr_files = num_apr_files + 1
               if( num_apr_files.gt.max_apr_files ) then
                   write(*,220) max_apr_files, apren(1:len_runstring)
 220               format('**FATAL** Too many aproiri files. ',
     .                    'Max allowed 'i3,/,'Names ',a)
                   stop 'GLIST: Too many apriori coordinate files'
               endif
               glb_apr_file(num_apr_files) = apren(is:js-1)
               is = js + 1
            end if
         end do
      end if      

* MOD TAH 140723: Add globk_command file name as option to get
*     the use_site, use_pos, use_num options.
      glb_mar_file = ' '
      len_runstring = rcpar(7, glb_mar_file)
      js = index(glb_mar_file,':')
      if( js.gt.0 ) then
          comopt = glb_mar_file(js+1:)
          glb_mar_file(js:) = ' '
      else
          comopt = ' '
      endif 

***** Thats all
      return
      end
  
CTITLE SUM_SITE_SOURCE
 
      subroutine sum_site_source( out_file, expt_names, most_sites,
     .           most_sources, most_svss, apr_values, report_allrn, 
     .           report_mods )
 
      implicit none 
 
*     Routine to produce a summary of the sites and source usage by
*     experiment from the global files.
 
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
*   glb_name(32)    - Name of input glb_file
*   i,j         - Loop counter
*   ierr        - IOSTAT error
*   expt_names(32,1)   - Names of experiments in this global
*   trimlen     - HP function for length of string
*   used(max_glb_sources)   - Counter for number of times used
*   times(2,max_glb_sources) - Start and stop times for this site/source
*   most_sites, most_sources - Most used at any one time
 
      integer*4 i, ierr, trimlen,
     .    used(max_glb_sources+max_glb_sites), unit, 
     .    most_sites, most_sources, most_svss

      real*8 times(2,max_glb_sources+max_glb_sites)
      real*8 apr_values(*)

      character*(*) expt_names(max_glb_sol)
 
*   out_file    - name of the output file or LU
 
      character*(*) out_file

      logical report_allrn 
      logical report_mods    ! Generated report with models used
 
***** Open the output file
      unit = 100 
      call open_lu(unit, out_file, ierr, 'unknown' )
      call report_error('IOSTAT',ierr,'open',out_file,1,
     .                  'SUM_SITE_SOURCE/OUTFILE')
 
*     Write out the header line
      write(unit,100, iostat=ierr) 'SITE',
     .    list_file(1:max(1,trimlen(list_file)))
 100  format(/,'GLIST: Contents of GLOBK Binary H-files',/,
     .         '---------------------------------------',//,
     .         'Summary of ',a,' occurences in ',a)

      if( num_eqfiles.gt.0 ) then 
         write(unit,120) num_eqfiles,(i,eq_inp_file(i)
     .                 (1:trimlen(eq_inp_file(i))),i=1,num_eqfiles) 
 120     format('There are ',i2,' earthquake files ',/,
     .        30('# ',i2,' File : ',a,/,:))
         if( report_allrn ) write(unit,'(a)') 
     .                  'Renames will be reported :A option'
      endif

      if( num_apr_files.gt.0 ) then
         write(unit,140) num_apr_files,
     .                   (i,trim(glb_apr_file(i)),i=1,num_apr_files) 
 140     format('There are ',i2,' apriori coordinate files ',/,
     .        30('# ',i2,' File : ',a,/,:))
      endif

 
*     Now loop over sort file
 
*     Now loop over the experiments
 
      call clear_used( used, times, gnum_sites)

* MOD TAH 140723: See if globk command file passed
      if( trimlen(glb_mar_file).gt.0 ) then
          write(unit,150) trim(glb_mar_file), trim(comopt)
          if( unit.ne.6 ) write(*,150) trim(glb_mar_file), 
     .         trim(comopt)
 150      format('READING USE Commands from ',a,' with option ',a)
          call read_use_cmds
      endif

*     Clear guse_site if no data
      do i = 1, gnum_sites
         if( times_used(i).eq. 0 ) call sbit(guse_site,i,0)
      end do

      
* MOD TAH 000902: See if apr file passed
      if( num_apr_files.gt.0 ) then
          write(*,*) 'Reading apriori positions '
          call get_apr_positions

* MOD TAH 140722: Test output of the apriori coordinates
          write(unit, 210) (trim(glb_apr_file(i)),i=1,num_apr_files)
 210      format('APRIORI COORDINATES FROM ',10(a,1x))
          write(unit, 215)
 215      format('A      # Site Name    JD Ref Epoch          X (m)',
     .           '            Y (m)            Z (m)        Vx (m/yr)',
     .           '    Vy (m/yr)    Zy (m/yr)     USED') 
         do i = 1,gnum_sites
             write( unit, 220) i, gsite_names(i), site_epoch(i),
     .            apr_val_site(:,1,i), apr_val_site(:,2,i),
     .            times_used(i)
 220         format('A ',I6,1x,a8,1x,F16.2,2x,3(F16.5,1x),3(F12.5,1x),
     .              ' # ',i6)
          enddo
          if( num_nonsec.gt.0 ) then 
               write(unit, 250) (trim(glb_apr_file(i)),
     .                           i=1,num_apr_files)
 250           format('EXTENDED APRIORI COORDINATES FROM ',10(a,1x))
               write(unit,255)
 255           format('E      # Site Name  Type         JD Ref Epoch',
     .                '     Tau (d)     XI (m)     YI (m)     ZI (m)',
     .                '     XQ (m)     YQ (m)     ZQ (m)')
               do i = 1, num_nonsec
                  write(unit,260) i, gsite_names(param_nonsec(1,i)),
     .                nonsec_types(param_nonsec(2,i)), 
     .                apr_val_nonsec(:,i)
 260              format('E ',I6,1x,a8,3x,a,1x,F16.2,1x,F12.4,
     .               6(F10.5,1x))
               end do
           end if

      end if      
 
      call write_header( unit,'sites', gsite_names, gnum_sites)
      
      do i = 1, num_glb_sol
 
 
*                         ! Get the name of the global file
          glb_inp_file = expt_names(i)
 
          call fmpopen(cglb_dcb,ierr,glb_inp_file,'ro',1)
          call report_error('FmpOpen',ierr,'open',glb_inp_file,0,
     .                      'SUM_SITE_SOURCE/GLB_INP_FILE/SITES')
 
*                                 ! Continue
          if( ierr.ge.0 ) then
              call rw_glb_header('R',ierr)

*             MOD TAH 150825: Get the PRN to SVN information
              call read_svinf_rec

*             Read the names block
              call rw_names_block('R')

              call eq_name_change('NO')
              call check_apr_inf( apr_values )
              
*             Now accumulate and write the information
              call accum_write( unit, used, times, ltog_sites,
     .                          cnum_sites, cepoch_expt, 
     .                          glb_inp_file, most_sites,'SITE') 
          end if
 
          call fmpclose( cglb_dcb,ierr)
      end do
 
*     Write summary
! MOD TAH 150829: Removed since now reduntant with write_sum_pos
!     call write_sum(unit, used, times, gsite_names, gnum_sites)

*     Write the updated summary entries
      call write_sum_pos( unit, used, times)

*     Write the renames
      if( report_allrn ) then 
          call report_eq( unit , 'ALL')
      else
          call report_eq( unit , 'UNIQ')
      end if 

****  Now do the sources

*     SKIP IF NO SOURCES

      IF( GNUM_SOURCES.GT.0 ) THEN
 
         write(unit,100, iostat=ierr) 'SOURCE',
     .     list_file(1:max(1,trimlen(list_file)))
 
         call clear_used( used, times, gnum_sources)
 
         call write_header( unit, 'sources',  gsource_names, 
     .                      gnum_sources)
 
         do i = 1, num_glb_sol
 
*                                ! Get global name from EMA
             glb_inp_file = expt_names(i)
 
             call fmpopen(cglb_dcb,ierr,glb_inp_file,'ro',1)
             call report_error('FmpOpen',ierr,'open',glb_inp_file,0,
     .                      'SUM_SITE_SOURCE/GLB_INP_FILE/SOURCES')
 
*                                    ! Continue
             if( ierr.ge.0 ) then
                 call rw_glb_header('R',ierr)
 
*                MOD TAH 150825: Get the PRN to SVN information
                 call read_svinf_rec

*                Read the names block
                 call rw_names_block('R')
 
*                Now accumulate and write the information
                 call accum_write( unit, used, times, 
     .                          ltog_sources, cnum_sources,
     .                          cepoch_expt, glb_inp_file,most_sources,
     .                          'SOURCE')
             end if
             call Fmpclose(cglb_dcb, ierr)
 
         end do


*        Write summary
         call write_sum(unit, used, times, gsource_names, gnum_sources)

      END IF


*     Summarize satellites

      IF( gnum_svs. GT.0 ) THEN
 
         write(unit,100, iostat=ierr) 'SATELLITES',
     .     list_file(1:max(1,trimlen(list_file)))
 
         call clear_used( used, times, gnum_svs)
 
         call write_header( unit, 'satellites',  gsvs_names, 
     .                      gnum_svs)
 
         do i = 1, num_glb_sol
 
*                                ! Get global name from EMA
             glb_inp_file = expt_names(i)
 
             call fmpopen(cglb_dcb,ierr,glb_inp_file,'ro',1)
             call report_error('FmpOpen',ierr,'open',glb_inp_file,0,
     .                      'SUM_SITE_SOURCE/GLB_INP_FILE/SVS')
 
*                                    ! Continue
             if( ierr.ge.0 ) then
                 call rw_glb_header('R',ierr)
 
*                MOD TAH 150825: Get the PRN to SVN information
                 call read_svinf_rec

*                Read the names block
                 call rw_names_block('R')
 
*                Now accumulate and write the information
                 call accum_write( unit, used, times, 
     .                          ltog_svs, cnum_svs,
     .                          cepoch_expt, glb_inp_file,most_svss,
     .                          'SVS')
             end if
             call Fmpclose(cglb_dcb, ierr)
 
         end do


*        Write summary
         call write_sum(unit, used, times, gsvs_names, gnum_svs)

      END IF

* MOD TAH 140327: Add block of output to write models used in each hfile.
      if( report_mods ) then

         write(unit,420) 
 420     format(/,'MODEL REPORT and USAGE',/,
     .            'YYYY MM DD    Dur Vers Sites  SVS    Parms',
     .            ' DRYZ WETZ DMAP WMAP  IonSrc   MagField Frame ',
     .            '    Prec      SRP    Time  OctalGAMIT MOD  ',
     .            'OctalLoad  SDEOP   ETIDE   OCEANTIDE  ATMLOAD S1/S2',
     .            '     HYDROLD   NUT     GRAV     ERAD    ANTTh ',
     .            '    H-file')
!2010 04 10    1.0    59   32      800 J2000    IAU76    BERNE    GPST        23600013     201403 IERS10   IERS03   FES2004E  NCEP CM ECMWFCMT          IAU00    EGM08    NCLE1    ANTB    

         do i = 1, num_glb_sol
 
*                                ! Get global name from EMA
             glb_inp_file = expt_names(i)
 
             call fmpopen(cglb_dcb,ierr,glb_inp_file,'ro',1)
             call report_error('FmpOpen',ierr,'open',glb_inp_file,0,
     .                         'SUM_SITE_SOURCE/MODELS')
   
*                                    ! Continue
             if( ierr.ge.0 ) then
                 call rw_glb_header('R',ierr)
                 call sub_null(cspeopmod)
                 call sub_null(cetidemod)
                 call sub_null(cotidemod) 
                 call sub_null(coatmlmod)
                 call sub_null(catmtdmod)
                 call sub_null(chydromod)
                 call sub_null(cgnut)
                 call sub_null(cggrav)
                 call sub_null(ceradmod)
                 call sub_null(cantradmod)
* MOD TAH 140403: Added ATM model lines
                 call sub_null(cdryzen)
                 call sub_null(cwetzen)
                 call sub_null(cdrymap)
                 call sub_null(cwetmap)
                 call sub_null(cionsrc)
                 call sub_null(cmagfield)

*                Write out line: 
                 call write_models(unit)
 
             end if
             call Fmpclose(cglb_dcb, ierr)
 
         end do
      end if
      write(unit,'(1x)')

 
***** Thats all
      close(100)
 
      return
      end
 
CTITLE CLEAR_USED
 
      subroutine clear_used( used, times, gnum)
 
      implicit none 
 
*     Clear the counters for the number of times a site or source
*     is used
 
*   gnum    - number of values to be cleared
*   i       - Loop counters
*   used(1) - Counters for number of times used
*   times(2,gnum) = Start and stop years
 
      integer*4 gnum, i, used(gnum)
      real*8 times(2,gnum)
 
      do i = 1, gnum
          used(i) = 0
          times(1,i) = 2100.0
          times(2,i) = 1900.0
      end do
 
****  Thats all
      return
      end
 
CTITLE WRITE_HEADER
 
      subroutine write_header( unit, type,  names, gnum)
 
      implicit none 
 
*     routine to write out the header
 
*   gnum        - Number of items in names
*   i,j         - Loop counter
*   ierr        - IOSTAT error
*   pos         - Current position in line
*   unit        - Unit for output
 
      integer*4 gnum, unit
 
*   names(1)    - List of names for header
*   type        - Type of quantity (sites or sources)
 
      character*(*) names(1), type
 
 
*     Now loop over the characters in the names and write them down
*     the page

      write(unit,100) type, gnum, type
 100  format('Use of ',a,' for all ',i4,1x,a) 
 
***** Thats all
      return
      end
 
CTITLE ACCUM_WRITE
 
      subroutine accum_write( unit, used, times, ltog, gnum, epoch, 
     .                        glb_file, maxn, type )

      implicit none 
 
*     Routine to add next line to file, with stars in the site
*     or source name position

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

      character*(*) type   ! type of list (SITE,SVS,SOURCE).
 
*   date(5)     - Date of current experiment
*   gnum        - Number of elements in ltog
*   i,j         - Loop counters
*   ierr        - IOSTAT error
*   ltog(1)     - Converstion for local site/source number to
*               - Global site source number
*   pos         - Position for star in buffer
*   trimlen     - Length of string routine
*   unit        - Unit for output
*   used(1)     - Counter for number of times used
*   maxn         - Maximum number of entries per line
*   numline - number of lines to be written excluding first and last
*   startlast, lenlast - Start of last line and its length.
 
      integer*4 date(5), gnum, i,j, ierr, ltog(1), pos, trimlen, unit,
     .    used(1), maxn, numline, startlast, lenlast

*   sort_used(max_glb_sites) - Sorted list of sites used (this is same size
*         for sources as well as sites
*   lenl   - Length of line to be written.

      integer*4 sort_used(max_glb_sites+max_glb_sources), lenl
 
*   epoch       - Epoch of this experiment
*   nepoch      - Nearest integer day to epoch
*   sec_tag     - seconds tag of epoch
*   times(2,*)  - Start and stop years for this site
*   dec_yrs     - Deciminal years
 
      real*8 epoch, nepoch, sec_tag, times(2,*), dec_yrs
 
*   line        - Line to be written out
*   glb_file    - Name of global file
 
*     character*16192 line
      character*32760 line
      character*(*) glb_file 

      logical kbit
 
****  get date of experiment
 
      nepoch = anint( epoch )
      call JD_to_YMDHMS( nepoch, date, sec_tag )
      call jd_to_decyrs( epoch, dec_yrs)
 
*     Start building line
* MOD TAH 150829: Keep 4-digit year
*     date(1) = date(1) - 1900
*     if ( date(1).ge. 100 ) date(1) = date(1) - 100
      write(line,100, iostat=ierr) (date(j),j=1,3)
 100  format(1x,I4,2i3.2)
 
*     Now put in the stars for those which have used
 
c     do i = 1, gnum
c         used(ltog(i)) = used(ltog(i)) + 1
c         pos = 10 + ltog(i)
c         line(pos:pos) = '*'
c     end do

*     Make list of used values and sort the list
* MOD TAH 140723: Added test to see if site should be used or not
*     based on globk command file.
      j = 0
      do i = 1, gnum
         if( kbit(guse_site,ltog(i)) ) then
             j = j + 1
             sort_used(j) = ltog(i)
         end if
      end do
*     Reset gnum based on used sites
!     gnum = j
   
C     call esort( gnum,  sort_used)

      do i = 1,gnum
          used(ltog(i)) = used(ltog(i)) + 1
          times(1,ltog(i)) = min(times(1,ltog(i)), dec_yrs)
          times(2,ltog(i)) = max(times(2,ltog(i)), dec_yrs)
*         Only output to the number of entries when guse_site 
*         is considered. 
          if( i.le.j ) then
             pos = 15 + 5*(i-1)
             write(line(pos:), 120) sort_used(i)
 120         format(i4)
          endif
      end do

      do i = j+1, maxn
          pos = 15 + 5*(i-1)
          write(line(pos:), 120) -1      
      end do

      pos = trimlen(line)
C     write(line(pos+2:),'(a)') glb_file(1:trimlen(glb_file))
      line(pos+2:) = trim(glb_file) // ' ' // type
c         
C     write(unit,'(a)', iostat=ierr )
C    .        line(1:trimlen(line))
      lenl = trimlen(line)

*     get the number of lines excluding the first and last
C     numline = (lenl-267)/255
*     Get length of last line
C     lenlast = lenl - 267 - numline*255
C     startlast = numline*255 + 267 + 1
C     if( lenlast.le.0 .and. numline.eq.1 ) then
C         write(unit,'(a)',iostat=ierr) line(1:lenl)
C     else
C         write(unit,'(a)',iostat=ierr) line(1:267)
C         do i = 1, numline
C             write(unit,'(12x,a)',iostat=ierr) 
C    .            line((i-1)*255+268:i*255+267)
C         end do
C         if( lenlast.gt.0 ) 
C    .    write(unit,'(12x,a)',iostat=ierr) 
C    .        line(startlast:startlast+lenlast)
C     end if
C     write(unit,140, iostat=ierr) line(1:128),
C    .            (line((i-1)*116+129:i*116+128),
C    .                   i = 1,(lenl-128)/116+1)
C140  format(a128,:/,100(12x,a116,:/))
 
      write(unit,'(a)') line(1:lenl)

****  Thats all
      return
      end
 
CTITLE WRITE_SUM
 
      subroutine write_sum( unit, used, times, names, gnum)
 
      implicit none 
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

*     Write out summary of the occurences of sites or sources
 
*   gnum        - number of values
*   i,j         - Loop counters
*   iel         - Look up value
*   ierr        - IOSTAT error
*   end_pass    - Number for any pass
*   last_pass   - number in last pass
*   passes      - Number of passes to print all
*   pos         - Current position on line
*   unit        - UNit for output
*   used(1)     - Counter for number of times used
 
      integer*4 gnum, i,j, iel, ierr, end_pass, last_pass, passes, pos,
     .    unit, used(1)

*   times(2,*)  - Start and stop year for each site/source

      real*8 times(2,*)
 
*   names(1)    - Names to be listed
 
      character*(*) names(1)
 
*   line        - Each line of output
 
      character*256 line
      logical line_out  ! Set true when line to be output
      logical kbit

      pos = 0

***** Output the number of occurrances
* MOD TAH 150829: Removed double column format and expanded year resoultions
      write(unit,100, iostat=ierr)
  100 format(/,' SUMMARY of site/svs/source occurences ',/,
     .         '   G# Name      Used Start      End       Dur (yr)')
 
      do j = 1, gnum
          if( used(j).gt.0 .and. kbit(guse_site,j))  then
*            Site used so output
             write(line,120) j, names(j), used(j),
     .             times(1,j), times(2,j),
     .             times(2,j)-times(1,j)
  120        format(1x,i4,' ', a,1x,i5,1x,f8.3,' - ',f8.3,1x,f7.3)
             write(unit,'(a)',iostat=ierr) trim(line)
         endif
      end do
 
C     write(unit,100, iostat=ierr)
C 100 format(/,' SUMMARY of site occurences ',/,
C    .         '   G#  Site     Used  Start   End    Dur  ',
C    .         '   G#  Site     Used  Start   End    Dur') 
C
*     We will put 2 values per line
C     passes = gnum/2 + 1
C     last_pass = mod( gnum,2)

* MOD TAH 140723: Removed output sites that were not used
C     iel = 0
C     j = 0
C     line = ' '
C     line_out = .false.
C     do while ( iel.lt. gnum ) 
C         iel = iel + 1
C         if( used(iel).gt.0 .and. kbit(guse_site,iel))  then
*            Site used so add to list
C            j = j + 1
C            if( j.eq.2 )  line_out = .true.
C            if( j.gt.2 ) then
C               j = 1
C               line = ' '
C            end if
C            pos = (j-1)*42 + 2
C            write(line(pos:),120) iel, names(iel), used(iel),
C    .             times(1,iel), times(2,iel),
C    .             times(2,iel)-times(1,iel)
C 120        format(i4,'. ', a,1x,i4,1x,f6.1,'-',f6.1,1x,f5.2)
C        else
C            if( times(1,iel).eq.2100.0 ) times(1,iel) = 0
C            if( times(2,iel).eq.1900.0 ) times(2,iel) = 0
C        endif
C        if( line_out ) then
C            write(unit,'(a)',iostat=ierr) line(:pos+42)
C            line_out = .false.
C        endif
C     end do
*     Make sure last line output when ony 1 entry
C     if( j.eq.1 ) write(unit,'(a)',iostat=ierr) line(:pos+42)


C     do i = 1, passes
C
C         if( i.eq.passes ) then
C             end_pass = last_pass
C         else
C             end_pass = 2
C         end if
C
C         line = ' '
*                                 ! Build line
C         do j = 1, end_pass
C             pos = (j-1)*42 + 2
C             iel = (i-1)*2  + j
C             write(line(pos:),120) iel, names(iel), used(iel),
C    .              times(1,iel), times(2,iel),
C    .              times(2,iel)-times(1,iel)
C 120         format(i4,'. ', a,1x,i4,1x,f6.1,'-',f6.1,1x,f5.2)
C         end do
C
C         write(unit,'(a)',iostat=ierr)
C    .        line(:pos+42)
C     end do
 
****  Thats all
      return
      end
 
 
CTITLE WRITE_SUM_POS
 
      subroutine write_sum_pos( unit, used, times)
 
      implicit none 
 
*     Write out summry of the occurences and positions of sites
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'

*   i,j         - Loop counters
*   used(1)     - Counter for number of times used
 
      integer*4 i, iel, ierr, unit, used(*)
      logical kbit 

*   times(2,*)  - Start and stop year for each site/source

      real*8 times(2,*)

* Additional variables for lat and long
      real*8 rot_mat(3,3), geod_pos(3), long, lat
 
***** Output the number of occurrances
 
      write(unit,100, iostat=ierr)
  100 format(/'. SUMMARY of site position and occurences P',/,
     .        '.  Long      Lat          Ht     #  First   Last ',
     .        '      Dur. Name        Seq  GNum P',/,
     .        '.   deg      deg          km                     ',
     .        '      yrs                         P')

 
      iel = 0 
      do i = 1, gnum_sites

*        Get the geodetic longitude and latitude of site
         call XYZ_to_GEOD( rot_mat, apr_val_site(1,1,i), geod_pos)

         long = geod_pos(2)*180.d0/pi
         lat  = (pi/2-geod_pos(1))*180.d0/pi

         if( times(1,i).gt.0 .and. times(2,i).gt. 0 .and.
     .       kbit(guse_site,i) ) then
             iel = iel + 1
             write(unit,120) long, lat, geod_pos(3)/1000.d0, used(i),
     .                 times(1,i), times(2,i), times(2,i)-times(1,i),
     .                 gsite_names(i), iel, i
 120         format(1x,f9.4,1x,f8.4,1x,f9.4, i5,1x,f8.3,1x,f8.3,f7.3,1x,
     .             a8,1x,i5,1x,i5,' P')
         end if
      end do

****  Thats all
      return
      end

CTITLE CHECK_APR_INF

      subroutine check_apr_inf( apr_values )

      implicit none 

*     Routine to check the aprioris and in this file versus ones in 
*     previous files.  Should be called save_apr_inf since this
*     routine will infact force the values to the same.

*     Due to the structure of the program, the apriori codes are
*     re-read here.

*     Routine only checks the station positions

* INCLUDES
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'

*   apr_codes(128)  - One record worth of apriori codes
*   i,j,k           - Loop counters
*   ierr            - File manipluation error routines
*   indx            - Index for the parameter being decoded eg. could
*                   - site number, source number etc.
*   len_read        - Length from readd.
*   type            - Type of parameter being decoded. See
*                   - GLB_HEADER_DEF.FTNI for descriptions.

 
      integer*4 apr_codes(128), i,j, ierr, indx, len_read, type

* k -- Loop counter
* doy -- Day of year
* yr  -- Year
* sec -- Integer seconds into start of day
* date(5) -- Date
* trimlen -- Length of string routine.
* pos     -- Element number from apriori values.

      integer*4 k, doy, yr, sec, date(5), trimlen, pos

 
*   apr_values(*)    - apriori values read from file (stored
*                    - temporaryily in ema). Use large size to
*                    - convince compiler to used I*4)
*   dXYZ(3) -- Position differences (m)
*   sectag  -- Seconds part of date

      real*8 apr_values(*), dxyz(3), sectag

      logical first_call
    
      data  first_call / .true. /
 
***** Loop over the apriori code records and get the code values.
*     These are then decoded and the corresponding values are
*     read from apriori values block
      do i = 1, 3
         dxyz(i) = 0.d0
      end do

***** Read all of the apriori values into the ema area.
      if( first_call ) then 
         write(*,110)
 110     format('Difference in the h-file aproiri coordinates of ',
     .          '>0.5 m either other h-files or apriori coordinate',
     .          ' files',/,
     .          'DAPR Site Name     dX (m)      dY (m)      dZ (m)',
     .          ' YYYY MM DD HR MN  DOY     H-file apriori XYZ ',
     .          'coordinates    H-file')
          first_call = .false.
      endif

      call readd(cglb_dcb , ierr, apr_values, 128*cnum_apr_vals,
     .           len_read, crec_apr_vals)
      call report_error('VREAD',ierr,'read','apriori values',
     .                  0,'GET_APRIORIS')
      if( ierr.ne.0 ) RETURN

      do i = 1, cnum_apr_types
 
          call readd(cglb_dcb, ierr, apr_codes, 128, len_read,
     .               crec_apr_types+i-1)
 
*****     Now start decoding and saving the corresponding values.
          do j = 1, 128
 
*             Only decode if a valid value
              if( j+(i-1)*128.le. cnum_apr_codes ) then
                  call decode_code ( apr_codes(j), type, indx )
 
*                 Now if these values match the previous ones
*                 (zero apriori means that it has not been set
*                 yet.
                  pos = (i-1)*128+j
                  if( type.eq.7 .and. 
     .                apr_val_site(1,1,ltog_sites(indx)).ne.0 
     .                                                    ) then
                      dxyz(1) = apr_val_site(1,1,ltog_sites(indx)) -
     .                         apr_values(pos)
                  else if ( type.eq.8  .and. 
     .                apr_val_site(2,1,ltog_sites(indx)).ne.0 
     .                                                   ) then
                      dxyz(2) = apr_val_site(2,1,ltog_sites(indx)) -
     .                          apr_values(pos)
                  else if ( type.eq.9  .and. 
     .                apr_val_site(3,1,ltog_sites(indx)).ne.0 
     .                                                   ) then
                      dxyz(3) = apr_val_site(3,1,ltog_sites(indx)) -
     .                          apr_values(pos)

*                    Now that we have the Z-value see if all the
*                    components are OK.
                     if( abs(dxyz(1)).gt.0.5 .or. 
     .                   abs(dxyz(2)).gt.0.5 .or.
     .                   abs(dxyz(3)).gt.0.5       ) then 

*                        The apriori's differ so warn the user.  
*                        Convert the time
                         call jd_to_ymdhms(cepoch_expt, 
     .                                     date, sectag)
                         call jd_to_yds(cepoch_expt, yr, doy, sec)
*                        Now write out the warning
C                        write(*,220) 
C    .                      gsite_names(ltog_sites(indx)),
C    .                      dxyz, date, doy,
C    .                      (apr_val_site(k,1,ltog_sites(indx)),
C    .                                         k=1,3),
C    .                       glb_inp_file(1:trimlen(glb_inp_file))
                         write(*,220) 
     .                      gsite_names(ltog_sites(indx)),
     .                      dxyz, date, doy, apr_values(pos-2:pos),
     .                       glb_inp_file(1:trimlen(glb_inp_file))
 220                     format('DAPR ',a8,3F12.2,
     .                          I5,4i3,2x,i3.3,2x,3F12.2,1x,a)
*                        Now zero dXYZ for next site
                         do k = 1,3
                            dXYZ(k) = 0.d0
                         end do
                     
                     end if
                  end if
              end if
          end do
      end do

****  Thats all
      return
      end

CTITLE WRITE_MODELS 
 
      subroutine write_models( unit )
 
      implicit none 


* INCLUDES
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'

* PASSED
      integer*4 unit  ! Output unit number

* LOCAL
      integer*4 date(5), ierr

      real*8 sectag
      real*8 dur   ! Duration of data in h-file

      call jd_to_ymdhms(cepoch_expt, date,sectag)
      dur = cepoch_end-cepoch_start
* MOD TAH 200429: Increased number of digits in version number/
      write(unit,120,iostat=ierr) date(1:3), dur, cglb_vers/100., 
     .    cnum_sites, cnum_svs, cnum_parn,
     .    cdryzen, cwetzen, cdrymap,  cwetmap, cionsrc,  cmagfield,
     .    cgframe, cgprec, cgsrpmod, cgtime, cgamit_mod,
     .    cload_mod, cspeopmod, cetidemod,
     .    cotidemod, coatmlmod, catmtdmod, chydromod, cgnut, cggrav, 
     .    ceradmod, cantradmod, trim(glb_inp_file)
 120  format(I4,1x,I2.2,1x,I2.2,1x,F6.1,1x,F4.2,1x,I5,1x,I4,1x,I8,1x,
     .    4(a4,1x),1x,2(a8,1x), 
     .    4(a8,1x),1x,o10,1x,o10,1x, 10(a8,1x),1x,a)
 



      end 


CTITLE check_apr_crd

      subroutine check_apr_crd

      implicit none 

*     Routine to try to replace zeros in apriori coordinates

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'

      integer*4 i,j,k   ! Loop counters
      logical matched   ! Set true when site names matched


*     Loop over sites
      do i = 1, gnum_sites
         if( sqrt(apr_val_site(1,1,i)**2+
     .            apr_val_site(2,1,i)**2+
     .            apr_val_site(3,1,i)**2).lt.6.d6 ) then
*            See if we can find a site with the same first
*            4 character
             matched = .false.
             do j = 1, gnum_sites
                if( i.ne.j .and. .not.matched .and.
     .              gsite_names(i)(1:4).eq.gsite_names(j)(1:4) .and.
     .              sqrt(apr_val_site(1,1,j)**2+
     .                   apr_val_site(2,1,j)**2+
     .                   apr_val_site(3,1,j)**2).gt.6.d6 ) then

****               Copy site coords
                   do k = 1,3
                      apr_val_site(k,1,i) = apr_val_site(k,1,j)
                   end do
                   if( gsite_names(i)(6:8).ne.'GPS' )
     .             write(*,120) gsite_names(j), gsite_names(i)
 120               format('* Copying coordinate of ',a,' to ',a)
                   matched = .true.
                end if
             end do
             if( .not.matched ) then
                write(*,140) gsite_names(i)
 140            format('* No coordinates found for ',a)
             end if
          end if
       end do

*****  Thats all
       return
       end          

CTITLE READ_USE_CMDS

      subroutine read_use_cmds

      implicit none

*     Routine to read use_site, use_num and use_pos commands from a globk
*     command file and sites the bits in guse_site accordingly

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'
      include '../includes/globk_cmds.h'

* LOCAL VARIABLES
*   iel         - Command number from get_cmd
*   ierr        - IOSTAT error
*   indx        - Pointer to current position in the command
*               - buffer
 
      integer*4 iel, ierr, indx, i, jerr

*   values -- Dummy entry for readline

      real*8 values

*   push_unit -- Pushed unit number
*   curr_unit -- Current unit number
*   len_comopt -- length of comopt string
*   num_use     - Sets number of times a site must a used.

      integer*4 push_unit, curr_unit, len_comopt, trimlen, num_use

* eol      -- Character number for ! or # (buffer cleared after this
*             position in string

      integer*4 eol

*   new_cmd_file -- Name of new command file
      character*256 new_cmd_file
 
*   process_line - Set true if current line from file is to be
*     processed (either starts with blank or comopt string).
 
      logical process_line, updated 

*   buffer      - Line read from markov file
      character*256 buffer

*   cval        - Dummy character string for read_line
      character*16 cval

*     Open the commnand file
      push_unit = 0
      curr_unit = 99
      num_use = 0    ! Use all
      len_comopt = trimlen(comopt)
      open(curr_unit,file=glb_mar_file,iostat=ierr,status='old')
      call report_error('IOSTAT',ierr,'open',glb_mar_file,1,
     .                  'glist/read_use_cmds')

****  Now loop over over file
      do while ( ierr.eq.0 ) 
         read(curr_unit,'(a)',iostat=ierr) buffer

*        See if in source mode
          if( ierr.ne.0 .and. push_unit.ne. 0  ) then
              close(curr_unit)
              curr_unit = push_unit
              push_unit = 0
              read(curr_unit,'(a)', iostat=ierr) buffer
          endif

         if( ierr.eq.0  ) then
!             Start to decode command
              indx = 1

              if( len_comopt.gt.0 ) then
                  call decode_comopt(buffer, comopt, updated) 
              end if 
              if( buffer(1:1).eq.' ' ) process_line = .true.
              if( process_line ) then
 
                  call get_cmd( buffer, glb_commands,
     .                max_glb_commands, iel, indx )
 
*                                             ! Command found
*                 See if SOURCE command used.  If so handle locally
* MOD TAH 210402: Fixed indexing bug due to added pre-glinit command: 78>79
                  if( iel.eq. 79 ) then
*                     Fatal if source used in source file
                      if( push_unit.ne.0 ) then
                         call report_stat('FATAL','globk',
     .                         'glist/read_use_cmds',buffer,
     .                         'SOURCE COMMAND in SOURCE file',0)
                      end if

                      call read_line(buffer,indx,'CH', jerr, values,
     .                    new_cmd_file)
                      call wild_card( new_cmd_file, list_file)

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
*                     Process actual commands.
*                     call process_glb_command( buffer, indx, iel,
*    .                                          glinit_run )
*                     Only process the use_xxx commands
                      eol = index( buffer,'!' )
                      if( eol.eq.0 ) eol = index(buffer,'#')
                      if( eol.gt.0 ) then
                           buffer(eol:) = ' '
                      end if
*
***** USE_SITE : Select sites to be used in the solution
* MOD TAH 210402: Fixed indexing bug due to added pre-glinit command: 43>44
                      if ( iel .eq. 44 ) then
                          call casefold(buffer)
                          call decode_option( buffer, gsite_names, 
     .                         gnum_sites,guse_site, -1)
***** USE_POS -- Sets the use of sites by latitude and longitude
* MOD TAH 210402: Fixed indexing bug due to added pre-glinit command: 68>69
                      else if( iel .eq. 69 ) then    
                          call edit_use_pos(buffer, indx)
* MOD TAH 210402: Fixed indexing bug due to added pre-glinit command: 60>70
                      else if( iel.eq. 70 ) then 
***** USE_NUM -- Sets the use site array based on the number of
*     times a site has been used.
                          call read_line(buffer, indx,'I4',jerr,
     .                                            num_use, cval)
                      endif    ! One if the USE_ commands
                  end if       ! If command #  greater than zero.
              end if     ! Process line
          endif          ! Read error is zero
      end do             ! Looping over command file
      close(curr_unit)

****  Finally apply any number limit 
      if( num_use.gt.1 ) then
         do i = 1, gnum_sites 
             if( times_used(i).lt.num_use ) then
                call sbit(guse_site, i, 0 )
             end if
         end do
      endif

      end

*     Include the block data in the main program since some compilers/linkers (specifically OSX) will not correctly 
*     load the initialized block data routines from library archives.
      include '../globk/globk_cmd_bd.f' 
      
