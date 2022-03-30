CTITLE glb_upd_apr

      subroutine glb_upd_apr( epoch, sepoch, do_pmu, sol_obs, do_rn )

      implicit none

*     This routine will update the apriori's for satellite orbits
*     and later polar motion while glfor and glbak are running.
*
* MOD TAH 0709226: Added SV epoch for update associated with 
*     referring solutions to mid-point epoch.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED VARIABLES

*   epoch       - Time at which the elements are required
*   sepoch      - Time at which satellite elements are required
*                 (during filter same as epoch, in output can be 
*                 different)
*   sol_obs(cnum_parn)  - Solution vector.  For renaming sites we
*                 may update this vector

      real*8 epoch, sepoch, sol_obs(cnum_parn)

*   do_pmu      - If true then pole position/UT1 will be updated
*                 This is always done expect in the glout/glorg call.
* MOD TAH 190531: Use do_pmu to test if orbit frame should be checked.
*                 For the orbit frame test to work we need a h-file to
*                 have been read and this does not happen in glorg.
*                 Only causes problem whrn back solution is run with
*                 glorg (190520 update).
*   do_rn       - If true update the site names.  Not done in 
*                 glout
*   do_ptide    - If true then pole tide will be applied.

      logical*4 do_pmu, do_rn, do_ptide
      logical chg_load(2)   ! Set true id we need to change load
                 ! application for atm or hydrology
 
*   ierr        - Fmpread error
*   trimlen     - Length of line
*   date(5)     - Date read from apriori file
*   indx,jndx   - Pointer for position in string
*   is          - sv number in list of svs
 
      integer*4  ierr, date(5), indx, jndx, is, trimlen, j, jerr, i
      
*   frame_pos   - Sum of positions of possible frame choices in svs
*                 line.  Used to see if the values are actually
*                 on the line
*   frame_match - Position of frame matched to current file (0 if
*                 there is no match)
*   prec_match  - Same as frame_match but for the precession model
*   lpt         - Length of the ptide_hfs substr name
*   nut_pos,nut_match  - Position and match of nutation on line to cgnut

      integer*4  frame_pos, frame_match, prec_match, lpt, 
     .           nut_pos, nut_match    

*   found       - Logical which stays set until we have found some
*                 entries and then it is used for exit
*   frame_OK    - Set true if the frame and precession matches or
*                 if their is no frame information in the svs_file
*   kbit        - Check bit

      logical found, frame_OK, kbit

* NOD TAH 200220: Added for more complete treatment of pole tide
      logical cstatus(3) !  Three entries 
                         ! (1) status known or not (SINEX vs h-files)
                         ! (2) Solid Earth pole-fide applied
                         ! (3) Ocean pole-tide applied
      character*6 cmpole ! Mean pole model (string to match names 
                         ! used in mean_pole.f used in curent file

*   sectag      - second tag of date
*   tepoch      - Test julian date to see if with in 0.5 days of the
*                 center of our experiment.
*   eepoch      - Experimenemt epoch truncated to an hour in the same
*                 way as the ephemeris epoch is truncated.

      real*8 sectag, tepoch, eepoch 

*   cprn       - PRN name read from file.
*   cval       - dummy for multiread
*   cptide     - Test string for pole-tide applicatiom

      character*8 cprn, cval, cptide

*   line        - Line read from file

      character*512 line

      character*128 message  ! Message for pole-tide warning

***** Try to open the apriori satellite position file

      if( trimlen(glb_svs_file).gt.0 .and. gnum_svs.gt.0 ) then

          open(105, file=glb_svs_file, status='old', iostat=ierr )
          call report_error('IOSTAT',ierr,'open', glb_svs_file, 0,
     .                      'glb_upd_apr')

*         if there was no open error, then try to read the values from
*         the file.
* MOD TAH 980217: Get experiment epoch truncated to and hour (same as
*         ephemeris entry).
* NOM TAH 070926: Get ephemeris at sepoch
          call jd_to_ymdhms(sepoch, date, sectag)
          date(5) = 0
          sectag  = 0.d0
          call ymdhms_to_jd(date, sectag, eepoch)
          
          found = .false.

          if( ichar(cgnut(1:1)).eq. 0 ) cgnut = 'IAU80'

          do while ( ierr.eq.0 )

              read(105,'(a)', iostat=ierr ) line
              if( line(1:1).eq.' ' .and. ierr.eq.0 .and. 
     .            trimlen(line).gt.0 )   then

*                 See if the date on this line matchs experimentr
                  indx = 1
                  call multiread(line, indx, 'I4', ierr, date, cval, 4)
                  date(5) = 0
                  sectag = 0.d0
                  call ymdhms_to_jd( date, sectag, tepoch)
* MOD TAH 980217: Short term fix to close orbits; set offset to 0.01 days
*                 (under an hour; effectively must match hour exactly).
* MOD TAH 980217: Use truncated eepoch for comparison.
* MOD TAH 981215: Now user settable tolerance.  Default is 0.1 days                  
                  if( abs(tepoch-eepoch).le. tol_svs_eph ) then
                      found = .true.

*                     we have found a time with 12 hours of the center
*                     of our experiment.  Find the SV Name
                      call getword(line, cprn, indx )
                      jndx = 1
                      call get_cmd(cprn, gsvs_names, gnum_svs, is, 
     .                                 jndx)
     
*****                 See if the Frame and Precession information is
*                     at the ends of the lines.  If it is then only
*                     use these values if the the same as the current
*                     binary h-file.
                      frame_pos = index(line,'B1950') + 
     .                            index(line,'J2000')
* MOD TAH 190531: Added do_pmu test (true for glfor and glbak but
*                 not for glout and glorg).
                      if( frame_pos.gt.0 .and. do_pmu) then

*                         OK, frame information on line, see if it
*                         matches the current file
                          frame_match = index(line,cgframe(1:5))
                          prec_match  = index(line,cgprec(1:5))
                          if( frame_match.gt.0 .and. 
     .                        prec_match.gt.0         ) then
                              frame_OK = .true.
                          else
                              frame_OK = .false.
                          end if
* MOD TAH 080128: Added check on nutation frame
*                         Now check to see if nutation matchs
                          nut_pos = index(line,'IAU80') +
     .                              index(line,'IAU00')
                          if( frame_OK .and. nut_pos.gt.0 ) then
*                             Check nutation consisentency
                              nut_match = index(line,cgnut(1:5))
                              if( nut_match.eq.0 ) then
                                 frame_OK = .false.
                              endif
                          endif
                          
                      else
*                         This must be old svs_file so assume OK and
*                         see what happens.
                          frame_OK = .true.
                      end if

*                     See if we have found the epoch and frame we
*                     nned.  
                      if( is.gt.0 .and. frame_OK ) then
                           call multiread(line, indx, 'R8', ierr,
     .                          apr_val_svs(1,is), cval, 9 )
                           do j = 10, max_svs_elem
                               call read_line(line, indx, 'R8', jerr,
     .                                apr_val_svs(j,is), cval)
                               if( jerr.ne.0 ) apr_val_svs(j,is) = 0
                           end do

****                       See if we should replace the satellite antenna
*                          offsets
                           do j = 1, 3
                              if( apr_val_ant(j,is).ne.-999.d0 ) then
                                  apr_val_svs(max_svs_elem-3+j,is) =
     .                               apr_val_ant(j,is)
                              endif
                           end do 
                      end if
                  else
*                     If we have already found svs and now time is
*                     wrong set ierr so that we get out
* MOD TAH 940620: Commented line below so that we will keep reading.  
*                 This ensures that ephemeris elements will be found
*                 if they are in the file
C                     if( found ) ierr = -1
                  end if
              end if
         end do

         close(105)
      end if
***** See if we need to read the pmu table file
      if( trimlen(pmu_inp_file).ne.0 .and. do_pmu ) then
          call read_pmu_inp  
      end if

***** See if we have any site movements associated with renaming
*     sites
      if( do_rn ) then
          call rn_app_dpos( sol_obs )
      end if

****  See if we should apply the pole tide correction.  We check
*     do_rn because this routine is also called by the output
*     programs.  For these we do not want to apply the corrections
*     again.
* MOD TAH 030517: If the bytes need swapping, then value should 
*     be > 2**30 since the EOP bits will be in the top byte.  Therefore 
*     if greater than 2**30, swap the bytes, otherwise leave untouched. 
* MOD TAH 110512: Removed the test below because no longer valid (only
*     applied to old swaph'd files.
!     if( cgamit_mod.gt. 2**30 )  then
!         call swap_bytes(4, cgamit_mod, 1)
!         write(*,505) cgamit_mod
!505      format('Swapping CGAMIT_MOD: Octal Value ',o11.11)
!     end if 
      do_ptide = .false.  ! Set not to do anything
      if( kbit(ptide_opt,1) .and. do_rn ) then

* MOD TAH 200220: Check status to see what we can do
          call ptd_status( cstatus, cmpole, cgamit_mod )
          if( .not.cstatus(1) ) then
*            We don't know status
             write(message,510) cgamit_mod
 510         format('Pole-tide changed requested but status ',
     .              'unknown: GAMIT_MOD o',o11.11)

             call report_stat('WARNING','GLOBK','glb_upd_apr',
     .              '',message,0)

          else   ! See if other critera are met.

*            See if this hfile is one we need to apply to
             do i = 1, num_ptide_hfs
                lpt = trimlen(ptide_hfs(i))
                cptide = ptide_hfs(i)
                call casefold(cptide)
                if( index(glb_inp_file,ptide_hfs(i)(1:lpt)).gt.0 .or.
     .              cptide(1:3).eq.'ALL' ) do_ptide = .true.
             end do
          endif

* MOD TAH 200220: Re-worked logic with fewer calls.
*         Start with SE-pole-tide (solid earth).
          if( do_ptide ) then
*             ptide_opt in common control what is done with
*             model application.
              call app_ptide ( sol_obs, cstatus(2), cmpole )
*             See if ocean pole tide to be applied
* MOD TAH 200717: Added test to see if removing as well OPT as well
              if( kbit(ptide_opt,2) .or. kbit(ptide_opt,4)) then   ! Apply ocean ptide 
                 call app_optd ( sol_obs, cstatus(3), cmpole )
              endif
              ptd_updated = .true.
          end if
       endif

* MOD TAH 130418: See if there is a change in the loads to be applied
*     Only change load values if APP_MODL command has been used
* MOD TAH 160303: Load corrections should only be applied when we
*     doing the rename (do_rn).  In glorg, a summary sol_obs array is
*     passed and this code should be executed.

      if( kbit(appload_mod,1) .and. do_rn ) then   ! Command has been used, 
*                                                   see what the status is
*         Test atmospheric load status (see if change in model)
          chg_load(1) = kbit(appload_mod,2) .neqv. kbit(cload_mod, 9)
          chg_load(2) = kbit(appload_mod,3) .neqv. kbit(cload_mod,25)
          if( chg_load(1) .or. chg_load(2) ) then
*            Update the load values in the solution vector
             call change_load( sol_obs )
          endif
      end if
 
***** Thats all
      return
      end 

CTITLE LOAD_UPD

      logical function load_upd()

      implicit none

*     Function to test if we need to change the status of loading applied
*     to the input globk solution
     
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* LOCAL VARIABLES

      logical kbit

***** Test each part
      load_upd = .false.
      if( kbit(appload_mod,1) ) then  ! Command to change load used, so 
*              see if changes are needed (ie., move to apply or remove 
*              loads if currently in other state).
*         Atmospheric load (when different then state is being 
*         changed) 
          if( kbit(appload_mod,2) .neqv. kbit(cload_mod, 9) ) 
     .                                           load_upd = .true.
*         Hydrologic load 
          if( kbit(appload_mod,3) .neqv. kbit(cload_mod,25) ) 
     .                                           load_upd = .true.
      endif

***** If function returns true, then we need to read load values from 
*     the sinf_rec records (Site information records) and save in the
*     gatmload/ghydload records; Note these are by global site number
*     so need to use ltog_sites array for mapping

      RETURN
      end



