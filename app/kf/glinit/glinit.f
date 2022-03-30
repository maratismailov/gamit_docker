      subroutine GLINIT(ms_type) 

      implicit none 
 
 
*     This program will sort the list of input global solution
*     files into ascending or descending order; and build up the
*     list of aprori site,source positions and other information
*     about the solution input parameters.
*
*     This program is normally run from GLOBK before any processing
*     of the Markov file is carried out.
*
*                                         11:22 PM SAT., 25 Jul., 1987
* MOD TAH 890228 Forced EmaSize to be relocated to root so that size
*     could be returned at true value always.  EmaSize sometimes fails
*     possibly due to its check value randomly appearing in memory.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/globk_cmds.h' 

c      external globk_cmd_bd

* ms_type   - Character string which indicates that glfor is a program
*          (MAIN) or a subroutine (SUBR)

      character*(4) ms_type

* MAIN PROGRAM VARIABLES 

*   dumm                            - Dummy place holder
*                                   - from read_line
*   i,j                             - Loop counter
*   ierr                            - IOSTAT error
*   expt_names(max_glb_sol )        - Names of the global
*                                   - solutions
*   indx                            - Positioner for Read_line.
*                                   - Keeps track of where we
*                                   - are in solution
*   jerr                            - Read_line error
*   kerr                            - Decodeing error
 
      integer*4 dumm, i, j,ierr, indx, jerr, kerr, lerr, fmpclose
 
*   cr                              - ASCII carriage return
 
      character*1 cr
 
*   buffer                          - Line read from list file.
 
      character*80 buffer
 
*   sepoch_expts( max_glb_sol )     - Epochs of each solution.
*                                   - List will be sorted based
*                                   - on these values
*   expts_var(max_glb_sol)          - Variances to be given to
*                                     each experiment (to make chi**2
*                                     unity).
*   apr_values( max_glbapr_sz )   - Storage for covariance
*                                   - matrix and solution vector
*   iexpt_names  - Start of names areas
*   isepoch_expts - Start of epochs area.
*   iexpts_var   - Start of expts_var area
*   iapr_values  - Start of apriori area 
*   
 
      integer*4 ema_data(1), trimlen
      integer*8 iexpt_names, isepoch_expts, iexpts_var, iapr_values

*   num_ephem  -  Number of satellites with current ephemeris elements
*   ephem_ep(max_ephem) - Epochs of current ephemeris
*   prev_ephem  - Epoch of previous ephermis elements

* The following values are to allow B1950 and J2000 values to be mixed.
*   ephem_frame(max_ephem) - Frames associated with each of the ephemeris
*                 epochs to be output
*   ephem_prec(max_ephem)  - Precession modes assocated with each ephemeris
*   ephem_nut(max_ephem)  - nutation modes assocated with each ephemeris
*   svs_frame(max_glb_svs) - Frames for each IC value saved
*   svs_prec(max_glb_svs)  - Precession modes associated with each IC value saved
*   svs_nut(max_glb_svs)  - Precession modes associated with each IC value saved
*   lc -- number of line from experiment list
*   inc -- Set true if experiment to be included
*   firstday -- JD of the first day of data in the hfiles
*   oneday -- Used to see if _XPS should be converted to _XCL (done if more 
*     than one day of data).

      integer*4 num_ephem, lc
      real*8 ephem_ep(max_ephem), prev_ephem, firstday
      logical inc, oneday, kbit
      
* MOD TAH 980519: Added explicit specification of diagonal 
*     scaling of matrices

* glb_diag -- Diagonal scaling factor in ppm
* glb_var  -- Complete matrix scaling

      real*8 glb_diag, glb_var

      character*8 ephem_frame(max_ephem), ephem_prec(max_ephem),
     .            ephem_nut(max_ephem),
     .            svs_frame(max_glb_svs), svs_prec(max_glb_svs),
     .            svs_nut(max_glb_svs)
     
 
      common / globc_ema / ema_data

*  expts_var_read  - Value for experiment varaiance read from file
*                    (defaults to 1 if no value present)

      real*8 expts_var_read

      real*8 crunjd, grunjd, sectag  ! Run time of current hfile being
              ! read (crunjd), run time of latest solution in this
              ! epoch of data (grunjd) and sectag for conversion.
              ! Values are used to see if SV antenna offsets should
              ! by updated. 

      character*8 cdum
 
***** Start decode the runstring for this program
      cr = char( 13 )

      if( ms_type.eq.'MAIN' ) then 
          call decode_glinit_run
      end if

***** Initialize some global parameters
 
      call global_init

****  Assign the memory for this run
      call glb_glnt_mem( iexpt_names, isepoch_expts, iexpts_var, 
     .                   iapr_values, ema_data )

      istart_vma = isepoch_expts

***** Read the earthquake list file

      call read_eq_file
      call sort_eq
      call sort_unf_rn
 
***** Open the input list and loop over it
 
      open( 102, file=list_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open', list_file,1,'GLINIT')
*                                                   ! Kill if error

***** See if we are to make the SVS_FILE
      if( make_svs_file ) then
          open(201, file=glb_svs_file, iostat=ierr, status='unknown')
      end if
 
***** Now loop over the input list, and save epochs, and apriori
*     information
 
      ierr = 0
      num_glb_sol = 0
      gnum_comb = 0
      most_cparn_num = 0
      num_ephem = 0
      prev_ephem = 0.d0
* MOD TAH 030223: Do better bookeeping on start and stop times
      gepoch_start = 1.d30
      gepoch_end = 0.d0 
      grunjd = 0.d0

* MOD TAH 980514: Set the used status on renames and eq's to zero.
*     Bit will be set when entry used.
      do i = 1, max_rn_wrd 
         rn_used(i) = 0
      end do
      do i = 1, max_eq_wrd
         eq_used(i) = 0
      end do

* MOD TAH 030607: New features
      glb_out_opt = 0   ! Makes all estimated parametres output
      lc = 0

* MOD TAH 030924: Implementation of _XPS rename feature.  If there
*     is more than one day between the hfiles (solution times) then
*     the XPS will be excluded by conversion to _XCL renames.  Other
*     wise the XPS will be treated at any other name.
      oneday = .true.
      firstday = 0.d0
      do while ( oneday .and. ierr .eq.0 ) 
          read(102,'(a)', iostat=ierr) buffer
 
*                                             ! Get file name from
* MOD TAH 950106: Check file name to see if non-blank and does not
*         start with # or *
          inc = .false.
          if ( ierr.eq.0 .and. trimlen(buffer).gt.0 .and.
     .         buffer(1:1).ne.'#' .and. buffer(1:1).ne.'*' ) then
*                                             ! line
              indx = 1
              call read_line( buffer, indx, 'CH', jerr, dumm,
     .                        glb_inp_file)

*****         Now open this global and get the epoch information
              if( jerr.eq.0 ) then
 
                  call fmpopen(cglb_dcb,kerr,glb_inp_file,'rwo',0)
                  call report_error('FmpOpen',kerr,'open',
     .                 glb_inp_file,0,'GLSAVE')
 
                  if( kerr.eq.0 ) then
                      call rw_glb_header('R', kerr)

* MOD TAH 970610: Check if header OK.
                      if( kerr.ne.0 ) then
                         if( kerr.eq.2001 ) then
                             write(*,*) '**WRONG BYTE ORDER** Run swaph' 
                         end if
                         call report_error('RW_GLB_HEADER',kerr,'read',
     .                       glb_inp_file,0,'GLINIT') 
                      end if 
                  end if

*                 Check the time of the solution
                  if( firstday.eq.0.d0 ) firstday = cepoch_expt
                  if( abs(cepoch_expt-firstday).gt.1.d0 ) 
     .                                          oneday = .false.
              endif
          endif
      enddo

****  Now see if we should change the _XPS to _XCL
      if( .not.oneday ) then
          write(*,*) 'More than 1-day of data, converting _XPS to _XCL'
          do i = 1, num_renames
             if( rn_codes(2,i)(6:8).eq.'XPS' ) then
                 rn_codes(2,i)(6:8) = 'XCL'
             endif
          end do
      end if

****  Now read the hfiles proper
      close(102)
      open( 102, file=list_file, iostat=ierr, status='old')
      ierr = 0

      do while ( ierr.eq.0 )
 
*                                              ! Line from file
          read(102,'(a)', iostat=ierr) buffer
 
*                                             ! Get file name from
* MOD TAH 950106: Check file name to see if non-blank and does not
*         start with # or *
          inc = .false.
          if ( ierr.eq.0 .and. trimlen(buffer).gt.0 .and.
     .         buffer(1:1).ne.'#' .and. buffer(1:1).ne.'*' ) then
*                                             ! line
              lc = lc + 1
*             Check to see if should be included.  (Must be after the
*             offset if it is to be included).

              if( mod((lc-decoff),decnum).eq.0 .and.
     .            lc.ge.decoff ) then
                  inc = .true.
              end if
              if( .not.inc ) then 
cd               write(*,120) buffer(1:trimlen(buffer))
cd 120             format('Decimate skipping: ',a)
              end if
          end if

          if( inc ) then
              indx = 1
              call read_line( buffer, indx, 'CH', jerr, dumm,
     .                        glb_inp_file)

*             Try to read the variance:
              glb_var = 1.d0
              glb_diag = 0.d0
              call read_line( buffer, indx, 'R8', kerr, glb_var,
     .                        cdum )
              if( kerr.ne.0 ) then
                  expts_var_read = 1.d0
              else

* MOD TAH 980519: see if diagonal passed
                  call read_line( buffer, indx, 'R8', kerr, glb_diag,
     .                        cdum )
                  if( kerr.ne.0 ) then
                       glb_diag = 0.d0
                       if( glb_var.eq.0.d0 ) glb_var = 1.d0
                       expts_var_read = glb_var
                  else
                       if( glb_diag.gt.0.d0 ) then
                           expts_var_read = -(glb_var + 
     .                          (1.d0+glb_diag/1.d6)/1000.d3)
                       else
                           expts_var_read = glb_var
                       endif
                  endif
              endif
               
*****         Now open this global and get the epoch information
              if( jerr.eq.0 ) then
 
                  call fmpopen(cglb_dcb,kerr,glb_inp_file,'rwo',0)
                  call report_error('FmpOpen',kerr,'open',
     .                 glb_inp_file,0,'GLSAVE')
 
                  if( kerr.eq.0 ) then
                      call rw_glb_header('R', kerr)

* MOD TAH 970610: Check if header OK.
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
                  
                  
****              See if the ephemeris time has changed.  If it has
*                 can we are making the svs_file, update now
                  if( make_svs_file .and. kerr.eq.0 ) then
                     if( csvs_epoch.ne.prev_ephem ) then
                        call upd_svs_file(201, num_ephem, ephem_ep, 
     .                        ephem_frame, ephem_prec, ephem_nut,
     .                        svs_frame, svs_prec, svs_nut)

*     
                         num_ephem = 1
                         prev_ephem = csvs_epoch
                         ephem_ep(num_ephem) =  cepoch_expt

*                        Clear the orbital elements so that the new
*                        values for the next arc will be set.
                         do i = 1,max_glb_svs
                            do j = 1,max_svs_elem
                               apr_val_svs(j,i) = 0.d0
                            end do
                         end do
                         
*                        Save the frame and preccesion mode
                         if( ichar(cgframe(1:1)).eq.0 ) 
     .                                                cgframe = 'J2000'
                         if( ichar(cgprec(1:1)).eq.0 ) cgprec = 'IAU76'
                         if( ichar(cgnut(1:1)).eq.0 ) cgnut = 'IAU80'
                         call check_ascii(cgframe)
                         call check_ascii(cgprec)
                         call check_ascii(cgnut)

                         ephem_frame(num_ephem) = cgframe
                         ephem_prec(num_ephem) = cgprec
                         ephem_nut(num_ephem) = cgnut
                         
                      else
* MOD TAH 131111: Only increment number of entries if the experiment
*                 day has changed.  This is to handle multi-day orbit
*                 arcs.
                         if( cepoch_expt.ne.ephem_ep(num_ephem) ) then
                             num_ephem = num_ephem + 1
                             ephem_ep(num_ephem) =  cepoch_expt

*                            Save the frame and preccesion mode
                             if( ichar(cgframe(1:1)).eq.0 ) 
     .                                                cgframe = 'J2000'
                             if( ichar(cgprec(1:1)).eq.0 ) 
     .                                                 cgprec = 'IAU76'
                             if( ichar(cgnut(1:1)).eq.0 ) 
     .                                                 cgnut = 'IAU80'
                             call check_ascii(cgframe)
                             call check_ascii(cgprec)
                             call check_ascii(cgnut)

                             ephem_frame(num_ephem) = cgframe
                             ephem_prec(num_ephem) = cgprec
                             ephem_nut(num_ephem) = cgnut
                          end if
                       end if
                  end if

*                                         ! we read header OK, so continue
                  if( kerr.eq.0 ) then
 
                      num_glb_sol = num_glb_sol + 1
* MOD TAH 160613: Make sure we don't exceed limit (max_glb_sol)
                      if( num_glb_sol.gt.max_glb_sol ) then
                          write(*,200) max_glb_sol
 200                      format('ERROR: Maximum number of global',
     .                           ' input files (',I6,') exceeded',/
     .                           'Increase max_glb_sol in',
     .                           ' kf/includes/kalman_param.f')
                          call report_stat('FATAL','GLOBK','glinit',
     .                         list_file,'Too many global files. Max ',
     .                         max_glb_sol)
                      end if

*                     Accumulate number of solution records that are for
*                     combined solutions.
                      gnum_comb = gnum_comb + cnum_comb
 
                      write(*,210) num_glb_sol
  210                 format(' Global ',i4,$)
                      if( mod(num_glb_sol,10).eq.0 ) then 
                          write(*,'(1x)')
                      end if
 
*                     Save the frame and prec mode for these values. (Only
*                     save if the apriori value is zero (and therefore not
*                     set yet).
                      do i = 1, max_glb_svs
                         if( apr_val_svs(1,i).eq.0.d0 ) then
                             call check_ascii(cgframe)
                             call check_ascii(cgprec)
                             call check_ascii(cgnut)
                             svs_frame(i) = cgframe
                             svs_prec(i) = cgprec
                             svs_nut(i)  = cgnut

                         end if
                      end do 
                      
****                  Save the experiment epoch, name, and variance
* MOD TAH 030223: Keep track of the start and end epochs
                      if( cepoch_start.ne.0 ) then
                          if ( cepoch_start.lt.gepoch_start ) then
                               gepoch_start = cepoch_start
                          end if
                          if( cepoch_end.gt. gepoch_end ) then
                               gepoch_end = cepoch_end
                          end if
                      end if


                      call save_epoch(num_glb_sol, cepoch_expt,
     .                      glb_inp_file, expts_var_read,
     .                      ema_data(isepoch_expts),
     .                      ema_data(iexpt_names), ema_data(iexpts_var))
       
*****                 Get the apriori name and position information and
*                     save this as well

*                     Now get and save the apriori values from the
*                     binary hfile. 
                      call save_apr_inf( ema_data(iapr_values), crunjd, 
     .                                   grunjd)
                      call solution_inf

*****                 Keep track of the maximum number of local
*                     parameters needed.
                      if( cnum_parn.gt.most_cparn_num ) 
     .                    most_cparn_num = cnum_parn
 
****                  Close the file
                      lerr = fmpclose(cglb_dcb)

*                                     ! Global open OK
                  end if
*                                     ! Name read OK from buffer
              end if
*                                     ! List file read OK
          end if
*                                     ! Looping over the input file
      end do


****  See if we have more epehemeris element to write out
      if( make_svs_file ) then
          call upd_svs_file(201, num_ephem, ephem_ep, 
     .                      ephem_frame, ephem_prec, ephem_nut, 
     .                      svs_frame, svs_prec, svs_nut)
          close(201)
      end if
 
***** Now sort experiment epochs, and create/write the type 2 file
*     with the sort listed
      write(*,'(1x)')
 
      call sort_epochs(ema_data(iexpt_names), ema_data(isepoch_expts),
     .                 ema_data(iexpts_var) )
 
***** Now sort the names of the site and sources
* MOD TAH 051202: Before sorting, make sure we have all the apriori
      call check_all_crd
 
      call sort_apriori
 
***** Set all sites to be used by default
      do i = 1, gnum_sites

*         if there is no data on a site due to renames then
*         set not to use
          if( times_used(i).gt.0 ) then
              call sbit( guse_site,i,1)
          else
              call sbit( guse_site,i,0)
          end if

* MOD TAH 130212: See if any sites are _X in name, if so remove
          if( .not. oneday .and. gsite_names(i)(5:6).eq.'_X' .and.
     .        kbit(guse_site,i) ) then
              write(*,'(a,1x,a)') 'More than 1-day, removing site ',
     .              gsite_names(i)
              call sbit( guse_site,i,0)
              times_used(i) = 0
          end if


      end do
 
      do i = 1, gnum_sources
          call sbit( guse_source,i,1)
      end do
 
 
***** Now write out the sorted list of experiments
 
      call write_sort(  ema_data(iexpt_names), ema_data(iexpts_var))
 
***** Finally create the GLOBK common with the information in it
 
*                               ! Will create file if needed.
      if( ms_type(1:4).eq.'MAIN' ) then
          call rw_globk_common('W')
      end if

***** Thats all
      close ( 102 )
      
      return
      end
 
      subroutine debug_out(type,num_glb_sol,
     .                      sepoch_expts, 
     .                      expt_names,
     .                      expts_var )

      implicit none 
     
      include '../includes/kalman_param.h'
     
      integer*4 num_glb_sol
      real*8 sepoch_expts(*), expts_var(*)
      character*(sort_recl) expt_names(*)
      character*(*) type
      
      integer*4 i
      
      write(*,100) type, num_glb_sol
 100  format('DEBUG: Location ',a,' Number glb sol ',i5)
      do i = 1, num_glb_sol
          write(*,150) i, sepoch_expts(i), expts_var(i),
     .                 expt_names(i)
 150      format(i4,F20.3,f10.3,1x,a64)
      end do
      
      return
      end     
