      subroutine glsave(ms_type,option)

      implicit none  
 
*     This program will save a global file from a GLOBK run in the
*     form needed for GLOBK.  See GLB_HEADER_DEF.FTNI for descripion
*     of the format of global files.  This program reads the sol_file
*     from GLOBK, collects the needed header information, compresses
*     the non-used portions of the parameters from the matrices and saved
*     the information to disk.
*
*     The parameter values are saved with the apriori values back into
*     the adjustments i.e., total station, source positions, nutation
*     angles and pole position/UT1-AT values are saved.
*
* MOD TAH 190619: Introduced new ms_type which allows frame realized
*     GLX files to be saved during back solutions.
*
*                                    August 25, 1987.
*
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
 
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* PASSED variables

* ms_type - sets if called from main progrom (MAIN) or
*           subroutine (SUBR)
*         - BACK option set to save GLX file during back solution.
* option  - Option which if set to NO will not remap the memory
*           space.
      character*(4) ms_type, option
 
*   dcb(16)         - DCB buffer used to read solution file
*   ema_data(max_vma_space) - Space for the covariance matrix and
*               - solution vector in ema.
*   ierr            - File reading error flag
*   len_runstring   - Length of the runstring paramter
*   rcpar           - Hp function to read runstring
*   FmpClose       - HP1000 compartability function to close file
*   nr             - Current count in runsting (in case -M option oassed)
 
      integer*4 dcb(16), ema_data(1), len_runstring,
     .    rcpar, FmpClose, ierr, trimlen, num, nr

*   midpt_option - Set true if mid-point time option selected (-M option)
*   nomidpt -- Set true to turn off midpt option that could be (-N option)
*             embedded in a com file.

      logical midpt_option, nomidpt, kbit

* MOD TAH 190511: Introduced when number of parameters > 32767.
      integer*8 i8  !  i8 value of 1 to force I8 calculations

      character*128 runstring

      integer*4 i,j  ! Loop counters
 
      common / globc_ema / ema_data
 
***** START, decode the runstring
      i8 = 1
      midpt_option = .false.
      nomidpt = .false.
      back_type = .false.   ! Indicates saving during back solution.

      if( ms_type(1:4).eq.'MAIN' ) then 
          len_runstring = rcpar(1, glb_com_file )
          if( len_runstring.le.0 ) then
              call proper_runstring('glsave.hlp','glsave',1)
          end if

*         See of -M option passed which will make the h-file refer its
*         epoch to center time MOD TAH 070824
          if( glb_com_file(1:2).eq.'-M' .or. 
     .        glb_com_file(1:2).eq.'-m') then
*             Mid point option passed
              midpt_option = .true.
              nr = 1
*             Now get the actual com file name
              len_runstring = rcpar(nr+1, glb_com_file )
          elseif( glb_com_file(1:2).eq.'-N' .or. 
     .            glb_com_file(1:2).eq.'-n') then
*             Mid point option passed
              nomidpt = .true.
              nr = 1
*             Now get the actual com file name
              len_runstring = rcpar(nr+1, glb_com_file )
          else
              nr = 0
          endif
              
 
*         Now read the SOLVK common
          call rw_globk_common('R')
          istart_vma = 0

*****     See if the name of the global file passed in the
*         runstring
          len_runstring = rcpar(nr+2, runstring)
          if( len_runstring.gt.0 ) glb_out_file = runstring
          len_runstring = rcpar(nr+3, runstring)
          if( len_runstring.gt.0 ) gdescription = runstring
          len_runstring = rcpar(nr+4, runstring)
          if( len_runstring.gt.0 ) then
              glb_sol_file = runstring
          end if        

          log_unit = 6

      else
          if( trimlen(glr_sol_file).gt.0 ) then
              glb_sol_file = glr_sol_file
          end if
          if( ms_type.eq.'BACK' ) then
*             Generate the name of glb_out_file based on wild
*             card input.  The original name is saved in the
*             subroutine so that it can be re-mapped multiple
*             times.
              call wild_date(glb_out_file, gepoch_expt, 'Inc')
              back_type = .true.
          endif
          write(*,110) trim(glb_out_file),
     .                 trim(glb_sol_file), ms_type 
 110      format('Creating ',a,' using ',a,' MS_TYPE ',a)
      end if
 
*     Get the solution covariance matrix
C     write(*,*) 'START_VMA (Initially)       ', istart_vma
      if( option(1:2).ne.'NO' ) then 
          call glfor_mem( ema_data )
      end if
      
C     write(*,*) 'START_VMA (Memory allocated)', istart_vma
      icov_parm = istart_vma 
C     isol_parm = icov_parm +  2*num_glb_parn*num_glb_parn
      isol_parm = icov_parm +  2*(i8*num_glb_parn)*num_glb_parn

* MOD TAH 1906020: Only read files if this is SUBR/MAIN call.
*     FOR BACK type, we will use copies that are in memory.
      if( ms_type.ne.'BACK' ) then     
         if( num.gt.0 ) close(num)
 
         call rw_glb_covar('O', dcb, ema_data(icov_parm))
         call rw_glb_covar('L', dcb, ema_data(icov_parm))
      endif

C*DEBUG
C     print *,'Calling list_soln ',num_glb_parn
C     call list_soln(num_glb_parn,ema_data(icov_parm), 500)

* MOD TAH 070824: Get the status of mid-point epoch output
      if( kbit(prt_opts,29) .or. kbit(org_opts,29) ) then
          midpt_option = .true.
      endif

      if( midpt_option .and. .not.nomidpt ) then
          gepoch_out = (gepoch_start+gepoch_end)/2.d0
          write(*,120) gepoch_out
 120      format('GLSAVE: Reference global file to mid-point epoch ',
     .           'JD: ',F15.5)
      else
          gepoch_out = gepoch_expt
          write(*,140) gepoch_out
 140      format('GLSAVE: Reference global file to end epoch ',
     .           'JD: ',F15.5)
      endif

* MOD TAH 070824: Check on the PMU status and change of epochs
*     do not match (associated with midppoint epoch option)
      if( abs(gepoch_expt-gepoch_out).gt.1.d-3 ) then
         call pmu_midpt
      endif
    
***** Make sure we have an output file name
      if( trimlen(glb_out_file).eq.0 ) then
          write(*,220) 
 220      format('GLSAVE: No name for output binary hfile. ',
     .           'Add name to runstring')
          stop 'GLSAVE: No  name for output binary hfile'
      endif
***** Now generate the parameter codes for both the apriori values
*     and the solved for parameters.  We also determine the number of
*     first parameter to be saved.
      call gw_glb_codes
 
***** Now create and write the header records
 
      call gw_glb_header
 
*     Now create and write the names information
 
      call gw_glb_names
 
*     Now create and write the solution description.  Need to pass
*     covariance matrix in to get ratio of sigma to apriori sigma.
* MOD TAH 190620: Pass in the type of solution.  type BACK will use
*     current global file
      if( ms_type.eq.'BACK' ) then
         call gw_description('BACK')
      else
         call gw_description('NORM')
      endif

*     The full names are determined in gw_description, so write
*     them out now. 
      call gw_glb_full

*     Now create and write the codes for the apriori and the estimated
*     parameters
      call gw_codes
 
*     create and save the apriori site and source positions
      call gw_aprioris

* MOD TAH 981104: Save the multi-day parameter epochs
      call gw_mul_ep
 
*     Save the solution and covariance matrix
* MOD TAH 190620: If this is BACK type, then use icov_sav_parm
      if( ms_type.eq.'BACK' ) then
         call gw_soln( ema_data(icov_sav_parm) )
      else
         call gw_soln( ema_data(icov_parm) )
      endif

      
*     Save the apriori covariances used in the analysis
      call gw_cons

*     Close file
      ierr = FmpClose( cglb_dcb )
      call report_error('FmpClose',ierr,'clos', glb_out_file,
     .                   0,'GLSAVE')

* MOD TAH 961028: Close the sol_file (opened again later)
      ierr = FmpClose( dcb )


***** Thats all
      end
