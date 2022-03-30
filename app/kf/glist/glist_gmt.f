      program glist
 
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
 
 
*   dumm                            - Dummy place holder
*                                   - from read_line
*   i,j                             - Loop counter
*   ierr                            - IOSTAT error
*   indx                            - Positioner for Read_line.
*                                   - Keeps track of where we
*                                   - are in solution
*   irm(5)                          - 5 parameters to be returned
*   jerr                            - Read_line error
*   kerr                            - Decodeing error
*   most_sites, most_sources - Most site/sources used in any one
*          experiment.
 
      integer*4 dumm, ierr, indx,  jerr, kerr, trimlen,
     .          most_sites, most_sources
 
*   cr                              - ASCII carriage return
 
      character*1 cr
 
*   buffer                          - Line read from list file.
 
      character*80 buffer
 
*   sepoch_expts( max_glb_sol )     - Epochs of each solution.
*                                   - List will be sorted based
*                                   - on these values
*   apr_values(max_glbapr_sz)           - Storage for covariance
*                                   - matrix and solution vector

      character*128 expt_names ( max_glb_sol )
 
      real*8 sepoch_expts( max_glb_sol ), apr_values(max_glbapr_sz),
     .       expts_var(max_glb_sol), expts_var_read
 
      real*8 crunjd, grunjd, sectag  ! Run time of current hfile being
              ! read (crunjd), run time of latest solution in this
              ! epoch of data (grunjd) and sectag for conversion.
              ! Values are used to see if SV antenna offsets should
              ! by updated. 
C     common / globc_ema / expt_names, sepoch_expts, apr_values

* MOD TAH 980513: added external for globk_cmd_bd
c      external globk_cmd_bd
 
***** Start decode the runstring for this program
 
      cr = char( 13 )
 
      call decode_glist_run

****  Read any earthquake file given
      call read_eq_file
* MOD TAH 100128: Added sort of earthquakes and resolving ?PS names
      call sort_eq
      call sort_unf_rn
 
***** Print out header
 
      write(*, 50)
   50 format(/' GLIST: Summarize global solution contents',//,
     .        ' Starting to read input data list')
 
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

      do while ( ierr.eq.0 )
 
*                                              ! Line from file
          read(100,'(a)', iostat=ierr) buffer
*         Don't worry about reading experiment variance. Not used here.
          expts_var_read = 1.0
 
*                                             ! Get file name from
          if ( ierr.eq.0 .and. buffer(1:1).ne.'#' .and. 
     .         buffer(1:1).ne.'*' ) then
*                                             ! line
              indx = 1
              call read_line( buffer, indx, 'CH', jerr, dumm,
     .                        glb_inp_file)
 
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
                      most_sites = max(most_sites, cnum_sites)
                      most_sources = max(most_sources, cnum_sources)
 
                      write(*,100) num_glb_sol, 
     .                             glb_inp_file(1:trimlen(glb_inp_file))
  100                 format(' Global ',i4,2x,a)
 
                      call save_epoch(num_glb_sol, cepoch_expt,
     .                      glb_inp_file, expts_var_read,
     .                      sepoch_expts, expt_names, expts_var) 
 
*****                 Get the apriori name and position information and
*                     save this as well
 
                      call save_apr_inf( apr_values , crunjd, 
     .                                   grunjd )
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
 
***** Now sort the names of the site and sources
 
      call sort_apriori
 
      close ( 100  )
 
***** Produce summary
 
      call sum_site_source( glb_com_file, expt_names, most_sites,
     .                      most_sources )
 
***** Thats all
      end
 
CTITLE DECODE_GLIST_RUN
 
      subroutine decode_glist_run
 
 
*     Routine to read the GLIST runstring.  The form is
*
*     GLIST, list_file, out_put_file , <sort_file>, <sort_direction>.
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
 
      integer*4 DecimalToInt, ierr, len_runstring, rcpar, i
      integer*4 is, js   ! Start and stop of name separators
      integer*4 trimlen

*   runstring   - Place to read runstring into
 
      character*64 runstring
      character*2048 eqent  ! Name of earthuake files
 
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
c MOD SCM 1/26/2010: Need to set num_eqfiles = 1, else the eq_file is not read by read_eq_file
      num_eqfiles = 0
      if( len_runstring.eq.0 .or. eqent(1:4).eq.'none') then
         eq_inp_file(1) = ' '
      else     ! See if should split names about
         is = 1
         js = 1
         do while ( js.lt. len_runstring+1 )
            js = js + 1
            if( eqent(js:js).eq.':' .or.
     .          eqent(js:js).eq.'=' .or.
     .          eqent(js:js).eq.' ' ) then
               num_eqfiles = num_eqfiles + 1
               if( num_eqfiles.gt.max_eqfiles ) then
                   write(*,120) max_eqfiles, eqent(1:len_runstring)
 120               format('**FATAL** Too many eathquake files. ',
     .                    'Max allowed 'i3,/,'Names ',a)
                   stop 'GLIST: Too many earthquake files'
               endif
               eq_inp_file(num_eqfiles) = eqent(is:js-1)
               is = js + 1
            end if
         end do
         write(*,140) num_eqfiles,(i,eq_inp_file(i)
     .                    (1:trimlen(eq_inp_file(i))),i=1,num_eqfiles) 
 140     format('There are ',i2,' earthquake files ',/,
     .        30('# ',i2,' File : ',a,/,:))

      endif
 
***** Thats all
      return
      end
 
 
CTITLE SUM_SITE_SOURCE
 
      subroutine sum_site_source( out_file, expt_names, most_sites,
     .           most_sources )
 
 
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
     .    most_sites, most_sources

      real*8 times(2,max_glb_sources+max_glb_sites)

      character*(*) expt_names(max_glb_sol)
 
*   out_file    - name of the output file or LU
 
      character*(*) out_file
 
***** Open the output file
      unit = 100 
      call open_lu(unit, out_file, ierr, 'unknown' )
      call report_error('IOSTAT',ierr,'open',out_file,1,
     .                  'SUM_SITE_SOURCE')
 
*     Write out the header line
 
      write(unit,100, iostat=ierr) 'SITE',
     .    list_file(1:max(1,trimlen(list_file)))
  100 format(/' Summary of ',a,' occurences in ',a)
 
*     Now loop over sort file
 
*     Now loop over the experiments
 
      call clear_used( used, times, gnum_sites)
      do i = 1, gnum_sites
          times_used(i) = 0
      end do
 
      call write_header( unit,'sites', gsite_names, gnum_sites)
 
      do i = 1, num_glb_sol
 
 
*                         ! Get the name of the global file
          glb_inp_file = expt_names(i)
 
          call fmpopen(cglb_dcb,ierr,glb_inp_file,'ro',1)
          call report_error('FmpOpen',ierr,'open',glb_inp_file,0,
     .                      'SUM_SITE_SOURCE')
 
*                                 ! Continue
          if( ierr.ge.0 ) then
              call rw_glb_header('R',ierr)

*             MOD TAH 150825: Get the PRN to SVN information
              call read_svinf_rec

*             Read the names block
              call rw_names_block('R')

              call eq_name_change('NO')
 
*             Now accumulate and write the information
              call accum_write( unit, used, times, ltog_sites,
     .                          cnum_sites, cepoch_expt, 
     .                          glb_inp_file, most_sites,
     .                          gsite_names ) 
          end if
 
          call fmpclose( cglb_dcb,ierr)
      end do
 
*     Write summary
      call write_sum(unit, used, times, gsite_names, gnum_sites)

*     Write the updated summary entries
      call write_sum_pos( unit, used, times)

*     Write the renames
      call report_eq( unit, 'UNIQ' )

****  Now do the sources
*     No longer do radio sources (not needed)

***** Thats all
      close(100)
 
      return
      end
 
CTITLE CLEAR_USED
 
      subroutine clear_used( used, times, gnum)
 
 
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
 100  format(' Use of ',a,' for ',i4,1x,a) 
 
***** Thats all
      return
      end
 
CTITLE ACCUM_WRITE
 
      subroutine accum_write( unit, used, times, ltog, gnum, epoch, 
     .                        glb_file, maxn, names )
 
 
*     Routine to add next line to file, with stars in the site
*     or source name position

      include '../includes/kalman_param.h'
 
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
 
      integer*4 date(5), gnum, i, ierr, ltog(1), trimlen, unit,
     .    used(1), maxn

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
*   names       - Names of sites 

      character*8 names(*) 
      character*2096 line
      character*(*) glb_file 
      integer*4 icount
      data icount /0/
      save icount

      icount=icount+1
 
****  get date of experiment
 
      nepoch = anint( epoch )
      call JD_to_YMDHMS( nepoch, date, sec_tag )
      call jd_to_decyrs( epoch, dec_yrs)
 
*     Start building line
      date(1) = date(1) - 1900
 
*     Make list of used values and sort the list
      do i = 1, gnum
         sort_used(i) = ltog(i)
      end do

      call esort( gnum,  sort_used)

      do i = 1,gnum
          used(ltog(i)) = used(ltog(i)) + 1
          times(1,ltog(i)) = min(times(1,ltog(i)), dec_yrs)
          times(2,ltog(i)) = max(times(2,ltog(i)), dec_yrs)
          write(unit,510, iostat=ierr ) icount,sort_used(i),
     .              names(sort_used(i))
 510      format (1x,i4,1x,i5,2x,a8)
      end do

      write(unit,520,iostat=ierr) icount,glb_file
 520  format (1x,'GLBFILE:',1x,i4,1x,a)
      lenl = trimlen(line)

****  Thats all
      return
      end
 
CTITLE WRITE_SUM
 
      subroutine write_sum( unit, used, times, names, gnum)
 
 
*     Write out summry of the occurences of sites or sources
 
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

      pos = 0

***** Output the number of occurrances
 
      write(unit,100, iostat=ierr)
  100 format(/' SUMMARY of occurences ')
 
*     We will put 2 values per line
      passes = gnum/2 + 1
      last_pass = mod( gnum,2)
 
      do i = 1, passes
 
          if( i.eq.passes ) then
              end_pass = last_pass
          else
              end_pass = 2
          end if
 
          line = ' '
*                                 ! Build line
          do j = 1, end_pass
              pos = (j-1)*42 + 2
              iel = (i-1)*2  + j
              if( times(1,iel).eq.2100.0 ) times(1,iel) = 0
              if( times(2,iel).eq.1900.0 ) times(2,iel) = 0
              write(line(pos:),120) iel, names(iel), used(iel),
     .              times(1,iel), times(2,iel),
     .              times(2,iel)-times(1,iel)
  120         format(i4,'. ', a,1x,i4,1x,f6.1,'-',f6.1,1x,f5.2)
          end do
 
          write(unit,'(a)',iostat=ierr)
     .        line(:pos+42)
      end do
 
****  Thats all
      return
      end
 
 
CTITLE WRITE_SUM_POS
 
      subroutine write_sum_pos( unit, used, times)
 
 
*     Write out summry of the occurences and positions of sites
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'

*   i,j         - Loop counters
*   used(1)     - Counter for number of times used
 
      integer*4 i, iel, ierr, unit, used(*)

*   times(2,*)  - Start and stop year for each site/source

      real*8 times(2,*)

* Additional variables for lat and long
      real*8 rot_mat(3,3), geod_pos(3), long, lat
 
***** Output the number of occurrances
 
      write(unit,100, iostat=ierr)
  100 format(/'. SUMMARY of site position and occurences P',/,
     .        '.  Long      Lat         Ht     #  First   Last ',
     .        '   Dur.  Name     Seq P',/,
     .        '.   deg      deg         km                     ',
     .        '   yrs                P')

 
      iel = 0 
      do i = 1, gnum_sites

*        Get the geodetic longitude and latitude of site
         call XYZ_to_GEOD( rot_mat, apr_val_site(1,1,i), geod_pos)

         long = geod_pos(2)*180.d0/pi
         lat  = (pi/2-geod_pos(1))*180.d0/pi

         if( times(1,i).gt.0 .and. times(2,i).gt. 0 ) then
             iel = iel + 1
             write(unit,120) long, lat, geod_pos(3)/1000.d0, used(i),
     .                 times(1,i), times(2,i), times(2,i)-times(1,i),
     .                 gsite_names(i), iel
 120         format(1x,f9.4,1x,f8.4,1x,f8.4, i5,1x,f7.2,1x,f7.2,f6.2,1x,
     .             a8,1x,i4,' P')
         end if
      end do

****  Thats all
      return
      end
 
*     Include the block data in the main program since some compilers/linkers (specifically OSX) will not correctly 
*     load the initialized block data routines from library archives.
      include '../globk/globk_cmd_bd.f' 
      
