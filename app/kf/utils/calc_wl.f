CTITLE calc_WL
      program calc_WL
 
*     This program is based on Calc_ion and will compute the 
*     VLBI equivalnent of the WL obsersable.
 
*     CI> calc_WL <KalObs1> <KalObs2> 
*
*     where <KalObs1> and <KalObs2> are the KalObs files corresponding
*     to the different frequency bands.
*
*     CI> calc_WL KalObs_X.kal KalObs_S.kal
*
*     tells the program to use the two files KalObs_X.kal and KalObs_S.kal
*     to calculate the the WL observable.
*
*     Throughout this program, one of the files is referred to as the
*     X-Band file, and the other as the S-Band.  The assignment of these
*     labels is for ease of variable-handling only;  the frequencies that
*     the files correspond to is irrelevant since the frequency information
*     is stored in the KalObs files.  Thus the two files could be from
*     X- and K-Band data bases, for example.
 
 
      include '../includes/kalman_param.h'
      include 'calc_wl.h'
 
*       len_KalObs_s    - Length of the S-Band KalObs file runstring
*   ,   len_KalObs_x    - Length of the X-Band KalObs file runstring
*   ,   len_obs_code    - Length of the <obs_code> runstring
*   ,   rcpar           - HP runstring function
 
      integer*4 len_KalObs_s, len_KalObs_x, len_obs_code, rcpar
 
*       help_file       - The name of the help file for this
*                       -   program
 
      character*64 help_file
 
      data help_file / 'calc_wl.hlp' /
 
***** Initialize the runstring variables
      sx_KalObs(1) = ' '
      sx_KalObs(2) = ' '
 
***** Recover the runstrings
      len_KalObs_s = rcpar(1,sx_KalObs(1))
      len_KalObs_x = rcpar(2,sx_KalObs(2))
 
***** Check that all three runstrings were given.  If not, Print help file
*         and abort
      if (len_KalObs_s * len_KalObs_x .eq. 0)
     .    call proper_runstring(help_file,1,1)
 
***** Calculate the ionospheric delays
      call calc_WL_delays
 
      end
 
CTITLE calc_WL_DELAYS
 
      subroutine calc_WL_delays
 
*     J.L. Davis                   4:32 PM  THU., 16  APR., 1987
 
*     General driver routine for calc_WL
 
 
*       i                   - Loop counter
 
      integer*4 i
 
***** Loop over files (1=S, 2=X)
      do i = 1, 2
 
*****     Read the files into EMA
          call read_into_ema(i)
 
      end do
 
***** Calculate ionospheric contributions to observables
      call calc_WL_obs
 
      end
 
CTITLE READ_INTO_EMA
 
      subroutine read_into_ema(f)
 
*     J.L. Davis                   5:34 PM  THU., 16  APR., 1987
*
*     Reads though the Fth KalObs file and stores the pertinent
*     iformation into EMA.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
      include 'calc_wl.h'
 
*       error                   - Error flag
*   ,   f                       - Frequency (1 or 2)
*   ,   i                       - Loop counter
 
      integer*4 error, f, i
 
 
      character*13 prog
 
      data prog / 'READ_INTO_EMA' /
 
***** Open the KalObs file
      call open_KalObs(sx_KalObs(f),error)
 
***** Check for error
      call report_error('FmpOpen',error,'open',sx_KalObs(f),1,prog)
 
***** Read the KalObs header
      call rw_KalObs_header('R',error)
 
***** Check for error
      call report_error('FMP',error,'read',sx_KalObs(f),1,prog)
 
***** Check that there are not too many observations
      call check_num_obs(max_obs,sx_KalObs(f))
 
***** Output the header
      call out_header(6,sx_Kalobs(f))
 
***** Store the header information
      call store_header(f)
 
***** Loop through KalObs file
      do i = 1, num_obs
 
*****     Read a data record
          call rw_KalObs_block('R','data',site,error,i)
 
*****     Check for error
          call report_error('FMP',error,'read',sx_KalObs(f),1,prog)
 
*****     Store observation in EMA
          call store_obs(f,i)
 
      end do
 
***** Close this KalObs
      call close_KalObs(error)
 
***** Check for error
      call report_error('FMP',error,'clos',sx_KalObs(f),1,prog)
 
      END
 
CTITLE STORE_HEADER
 
      subroutine store_header(f)
 
*     J.L. Davis                   5:50 PM  THU., 16  APR., 1987
*
*     Stores the header information for the Fth KalObs file
*
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include 'calc_wl.h'
 
*       f               - The file index
 
      integer*4 f
 
***** Store the number of observations
C     call MoveWords(num_obs,sx_num_obs(f),1)
      sx_num_obs(f) = num_obs 
 
***** Move the site names
      call MoveChars(site_names,sx_site_names(1,f),num_sites)
 
***** Move source names
      call MoveChars(source_names,sx_source_names(1,f),num_sources)
 
      end
 
CTITLE STORE_OBS
 
      subroutine store_obs(f,ob)
 
*     J.L. Davis                   6:07 PM  THU., 16  APR., 1987
*
*     Routine to store the pertinent information from a KalObs
*     file record into EMA
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
      include 'calc_wl.h'
 
*       bit_off         - Bit offset in ion_flag
*   ,   f               - KalObs index
*   ,   ob              - KalObs record number
*   ,   type            - Loop index
 
      integer*4 bit_off, f, ob, type
 
***** Data flag
      sx_data_flag(ob,f) = data_flag
 
***** FRNGE quality code
      sx_FRNGE_code(ob,f) = ichar(FRNGE_qual(2:2))
      if(FRNGE_qual(2:2).lt.'5' .or.
     .   FRNGE_qual(2:2).gt.'A' ) then
         sx_data_flag(ob,f) = sx_data_flag(ob,f) + 1
      end if
 
***** KalObs record number
      sx_KalObs_rec(ob,f) = ob
 
***** Ion flag
      sx_ion_flag(ob,f) = ion_flag
 
***** Site numbers
      sx_site(1,ob,f) = site(1)
      sx_site(2,ob,f) = site(2)
 
***** Source number
      sx_source(ob,f) = source
 
***** Sigmas
*                                            ! Group delay
      sx_sigma(1,ob,f) = db_sigma(1)
*                                            ! Phase delay
      sx_sigma(2,ob,f) = db_sigma(2)
*                                            ! Phase delay rate
      sx_sigma(3,ob,f) = db_sigma(4)
 
***** Epoch
      sx_epoch(ob,f) = epoch
 
***** O-C
      sx_omc(1,ob,f) = observable(1) + dble(num_grp_amb) * grp_amb
*                                                         ! Group delay
*    .               - db_theoretical(1)
      sx_omc(2,ob,f) = observable(2) + dble(num_phs_amb) * phs_amb
*                                                         ! Phase delay
*    .               - db_theoretical(2)
*                                                         ! Phase delay rate
      sx_omc(3,ob,f) = observable(4)
*     sx_omc(3,ob,f) = observable(4) - db_theoretical(4)
 
***** Effective frequencies
      call effective_freqs(ob,f)
 
***** Initialize to no ion and no mathcing ion
      do type = 1, 2
 
*****     Caluclate bit offset in ion_flag for this type
          bit_off = 6 * (type - 1)
 
*****     Set the bits if we are doing this type
          if (calc_WL_cont(type)) then
 
*****         Turn off "ionosphere available" bit
              call sion(sx_ion_flag(ob,f),bit_off+4,0)
 
*****         Turn on "no matching data" bit
              call sion(sx_ion_flag(ob,f),bit_off+2,1)
 
          end if
 
*****     If this is the second file, initialize the corresponding record
*             Number
          if (f .eq. 2) rec_file_1(ob) = 0
 
      end do
 
      end
 
CTITLE calc_WL_OBS
 
      subroutine calc_WL_obs
 
*     J.L. Davis                   2:05 PM  FRI., 17  APR., 1987
*
*     Calculates ionospheric contribution to observables
*
 
      include '../includes/kalman_param.h'
      include 'calc_wl.h'
 
*       krec1, krec2            - Record numbers, files 1 and 2
*   ,   last_krec2              - Remember the previous krec2
*   ,   match_obs               - Routine for finding matching records
 
      integer*4 krec1, krec2, last_krec2, match_obs
 
***** Initialize record pointer for searching through records
      last_krec2 = 0
 
***** Loop through first file's observations
*
*     Write out the header:
      write(*,120)
 120  format('  Rec   Site 1  2 Source   WL     +-    Dphase X   +-',/,
     .       '                           (cycles)       (cycles)')
      do krec1 = 1, sx_num_obs(1)
 
*****     Find a match for this record in the second file's records
          krec2 = match_obs(krec1,last_krec2)
 
*****     Did we find a match?
          if (krec2 .eq. 0) then
 
*****         We did not find a match.  Report.
              call report_no_match(6,krec1)
 
          else
 
*****         Store the file 1 record in the pointer
              rec_file_1(krec2) = krec1
 
*****         We found a match.  Calculate ionospheric delay for this
*                 observation
              call calc_observation(krec1,krec2)
 
*****         Set the KalObs flags
              call set_KalObs_flags(krec1,krec2)
 
*****         Remember last record number
              last_krec2 = krec2
 
          end if
 
      end do
 
      end
 
CTITLE CALC_OBSERVATION
 
      subroutine calc_observation(krec1,krec2)
 
*     J.L. Davis                   2:13 PM  FRI., 17  APR., 1987
*
*     Calculates the ionospheric contribution for the observation with
*     record number KREC1 from file #1 and KREC2 from file #2
*
*     The basic ionospheric formula used is
*
* I(type,1) = Kappa(type) * [tau(type,1) - tau(type,2)] / (1 - beta(type))
*
*     where <type> is the type of observable (1=group delay, 2=phase delay,
*     3 = phase delay rate), <tau(type,i)> is the O-C for observable
*     of type <type> at frequency <i>, and <beta(type)> is given by
*
*     beta(type) = [f_eff(type,1) / f_eff(type,2)]^2
*
*     where <f_eff(type,i)> is the effective frequency for the observable
*     of type <type> at frequency <i>.  The definition of <I> is
*     slightly different for the different types.  The definitions are:
*
*     type 1 (group delays):  I(1,1) = The ion contribution to the group
*                                      delay at frequency 1 as cal-
*                                      culated from the group delays.
*     type 2 (phase delays):  I(2,1) = The ion contribution to the group
*                                      delay at frequency 1 as calculated
*                                      from the phase delays.
*     type 3 (phase delay rates):  I(3,1) = The ion contribution to the
*                                           phase delay rates at frequency
*                                           1 as calculated from the
*                                           phase delay rates.
*
*     Because of these inconsistent definitions, we have introduced the
*     variable <kappa>, seen above in the formula for <I>.  For the
*     group delays and phase delay rates (types 1 and 3), <kappa>=1,
*     because the observables are being used to determine the ion
*     contribution for themselves.  For <type>=2, we have that the
*     phase delays are being used to determine the ion contribution to
*     the GROUP delays, and so we must account for the fact that the
*     effective frequencies for the group and phase delays for a given
*     frequency band are different, and that the contribution for the
*     group delays has a different sign than for the phases.  Thus we
*     have
*
*     kappa(2) = - [f_eff(2,1) / f_eff(1,1)]^2
*
*     N.B.  The term "inconsistent" used above refers only to the
*           need for special treatment in this routine.  This
*           inconsistency is necessary for consistency in other
*           parts of the Kalman filter software.
*
 
 
      include '../includes/kalman_param.h'
      include 'calc_wl.h'
 
*       krec1                   - Record number from first file
*   ,   krec2                   - Record number from second file
*   ,   type                    - Loop variable (gr-ph-pdr)
 
      integer*4 krec1, krec2, type, cor
 
*       ak_factor               - Temp variable
*   ,   beta                    - Temp variable
*   ,   delta_omc               - Temp variable
*   ,   kappa                   - Temp variable
*   ,   sum_sq_sig              - Temp variable
 
 
      real*8 ak_factor, beta, delta_omc, kappa, sum_sq_sig
      real*8 A,B,C,D, AX, AS, wlgn, fxp, fsp, fxg, fsg, fx,fs,
     .       swl, txg, txp, tsg, tsp, E, F, G, pxo, pxg, dpx,
     .       spx, bx, bs, tcg, tcp, betp
 
***** Loop over all three observable types: 1 = group delay; 2 = phase
*         delay; 3 = phase delay rate.
      do type = 1, 3
 
*****     Beta is relates the ionospheric delay at the two frequencies
*             for this observable type
          beta = (sx_eff_frq(type,krec1,1)/sx_eff_frq(type,krec2,2))**2
 
*****     Kappa only matters for the phase delays.  It relates the ion-
*             ospheric delay at the phase delay effective frequency to
*             the phase delay at the group delay effective frequency.
          if (type .eq. 2) then
 
              kappa = -(sx_eff_frq(2,krec1,1)/sx_eff_frq(1,krec1,1))**2
 
          else
 
              kappa = 1
 
          end if
 
*****     Calculate the constant product to be used for the delays and
*             the sigmas
          ak_factor = kappa / (1.0D+00 - beta)
 
*****     Calculate the difference in (O-C)'s for this observation
          delta_omc = sx_omc(type,krec1,1) - sx_omc(type,krec2,2)
 
*****     Compute the ionospheric delay for this epoch at the frequency
*             of the first file
          sx_ion_cont(type,krec1) = ak_factor * delta_omc
 
*****     Calculate the sum of the sigmas squared
          sum_sq_sig = sx_sigma(type,krec1,1) ** 2
     .               + sx_sigma(type,krec2,2) ** 2
 
*****     Compute the sigma for this delay
          sx_ion_sigma(type,krec1) = abs(ak_factor) * sqrt(sum_sq_sig)
 
      end do

***** Now compute the WL data type/
      fxg = sx_eff_frq(1,krec1,1)
      fsg = sx_eff_frq(1,krec2,2)
      fxp = sx_eff_frq(2,krec1,1)
      fsp = sx_eff_frq(2,krec2,2)
      fx = sx_phs_ref(1)
      fs = sx_phs_ref(2)
      beta = (fxg/fsg)**2
      betp = (fxp/fsp)**2
      txg = sx_omc(1,krec1,1) 
      tsg = sx_omc(1,krec2,2) - txg
      txp = sx_omc(2,krec1,1) - txg
      tsp = sx_omc(2,krec2,2) - txg
      txg = 0.d0

      tcg = (tsg-beta*txg)/(1.d0-beta)
      tcp = (tsp-betp*txp)/(1.d0-betp)

      A = -beta*(fx-fs)/(1.d0-beta)
      B = -(fx-fs)/(1.d0-beta)
      C = (fx*(fxg/fxp)**2 - fsp*(fxg/fsp)**2)/(1.d0-beta)

      AX =   A - C
      AS = -(B - C)

c     tc = (sx_omc(1,krec2,2) - sx_omc(1,krec1,1)*beta)/(1-beta)
c     tg = sx_ion_cont(1,krec1)
c     wlg = tc*(sx_phs_ref(1)-sx_phs_ref(2)) - tg*
c    .      (  (sx_eff_frq(1,krec1,1)/sx_eff_frq(2,krec1,1))**2*
c    .                    sx_phs_ref(1)   -
c    .         (sx_eff_frq(1,krec1,1)/sx_eff_frq(2,krec2,2))**2*
c    .                    sx_phs_ref(1)  )
      wlgn = AX*txg + AS*tsg

c     wlp =  sx_omc(2,krec1,1) * sx_phs_ref(1) -
c    .       sx_omc(2,krec2,2) * sx_phs_ref(2)   
      wlp =  txp*fx - tsp*fs
      WL  =  (wlp - wlgn) * 1.d-6
      swl = sqrt( (ax*sx_sigma(1,krec1,1))**2 + 
     .            (as*sx_sigma(1,krec2,2))**2  )*1.d-6

*     Now compare X-band phase delay
      E = -fx*beta/(1.d0-beta)
      F = fx/(1.d0-beta)
      G = fx*(fxg/fxp)**2/(1.d0-beta)

      bx = e - g
      bs = f + g

      pxg = bx*txg + bs*tsg
      pxo = txp*fx
      dpx = (pxo - pxg)*1.d-6
      spx = sqrt( (bx*sx_sigma(1,krec1,1))**2 + 
     .            (bs*sx_sigma(1,krec2,2))**2  )*1.d-6

      write(*,200) 'WL ', krec1, sx_site(1,krec1,1), 
     .           sx_site(2,krec1,1), sx_source(krec1,1), WL, swl,
     .           dpx, spx, tcp-tcg,
     .          cor(sx_data_flag(krec1,1),sx_data_flag(krec2,2)),'E'
 200  format(a3,i5,' SSS ',3I3, 2F8.2,1x, F12.2, F8.2,F9.3, o12,1x,a)
 
      end
 
CTITLE MATCH_OBS
 
 
      integer*4 function match_obs(krec1,last_krec2)
 
*     J.L. Davis                   4:13 PM  FRI., 17  APR., 1987
*
*     Routine to find the observation matching observation #KREC1
*     in file 2.
 
 
      include '../includes/kalman_param.h'
      include 'calc_wl.h'
 
*       bounce_search           - Search function
*   ,   match                   - We found a match
*   ,   check_match             - Function to check to obs for
*                               -   a match
 
      logical bounce_search, match, check_match
 
*       krec1                   - Record number in 1st file
*   ,   krec2                   - Test for match
*   ,   last_krec2              - The previous record in file 2
*                               -   matched to file 1
 
      integer*4 krec1, krec2, last_krec2
 
***** Initialize to not found
      match_obs = 0
 
***** Let's first try the next record
      krec2 = last_krec2 + 1
 
***** Did we already go past file 2?
      if (krec2 .gt. sx_num_obs(2)) return
 
***** See if we found a match
      match = check_match(krec1,krec2)
 
***** Did we find a match?
      if (.not. match) then
 
*****     Do a "bouncing search"
          match = bounce_search(krec1,krec2)
 
*****     If we matched, say so
          if (match) write(*,300) krec1, krec2
  300     format(' Matched rec #',I5,' after search with rec #',I5)
 
      end if
 
***** If we wound up matching, assing a value to return
      if (match) match_obs = krec2
 
      end
 
CTITLE CHECK_MATCH
 
 
      logical function check_match(krec1,krec2)
 
*     J.L. Davis                   4:31 PM  FRI., 17  APR., 1987
*
*     Routine to see if two observations match
 
 
      include '../includes/kalman_param.h'
      include 'calc_wl.h'
 
*       leq                     - Lexically equal
 
      logical leq
 
*       krec1, krec2            - The record numbers in the two files
 
      integer*4 krec1, krec2
 
*       depoch                  - Difference in epochs, minutes
*   ,   tolerance               - Greatest depoch allowed, minutes
 
      real*8 depoch, tolerance
 
*       name(2,2)               - Temporary array for comparison
 
      character*8 name(2,2)
 
      data depoch / 2.0D+00 /
 
***** Default:  They match
      check_match = .true.
 
***** Determine the difference in the epochs in minutes
      depoch = 1.440D+03 * abs(sx_epoch(krec1,1) - sx_epoch(krec2,2))
 
***** If this difference is greater than the tolerance, no match
      if (depoch .gt. tolerance) check_match = .false.
 
***** What are the source names?
      name(1,1) = sx_source_names(sx_source(krec1,1),1)
      name(2,1) = sx_source_names(sx_source(krec2,2),2)
 
***** Check the source
      check_match = check_match .and. leq(name(1,1),name(2,1))
 
***** What are the site names?
      name(1,1) = sx_site_names(sx_site(1,krec1,1),1)
      name(2,1) = sx_site_names(sx_site(2,krec1,1),1)
      name(1,2) = sx_site_names(sx_site(1,krec2,2),2)
      name(2,2) = sx_site_names(sx_site(2,krec2,2),2)
 
***** Check the sites
      check_match = check_match .and.
     . ( (leq(name(1,1),name(1,2)) .and. leq(name(2,1),name(2,2))) .or.
     .   (leq(name(1,1),name(2,2)) .and. leq(name(2,1),name(1,2))) )
 
      end
 
CTITLE BOUNCE_SEARCH
 
 
      logical function bounce_search(krec1,krec2)
 
*     J.L. Davis                   5:01 PM  FRI., 17  APR., 1987
*
*     Performs a bouncing search for the observation matching the
*     KREC1th observation.  Returns true if a matching observation found.
*
 
      include '../includes/kalman_param.h'
      include 'calc_wl.h'
 
*       check_match         - Function to check the match of two
*                           -   observations
*   ,   over                - Set true when search takes us past
*                           -   EOF
*   ,   under               - Set true when search takes us to non-
*                           -   positive record numbers
 
      logical check_match, over, under
 
*       krec1               - Record number from file 1
*   ,   krec2               - Record number from file 2 to be
*                           -   returned
*   ,   move_rec            - Number of records to jump relative
*                           -   to starting record
*   ,   num_move            - The number of moves that have been made
*   ,   start               - The starting record (i.e., the value
*                           -   of KREC2 at the start
 
      integer*4 krec1, krec2, move_rec, num_move, start
 
***** Default: no match
      bounce_search = .false.
 
***** Count how many moves we make
      num_move = 0
 
***** Remember the strating record
      start = krec2
 
***** We start within bounds
      over  = .false.
      under = .false.
 
***** Search until we've found it or until we've run out of space
      do while (.not. (under .and. over)
     .             .and. .not. bounce_search)
 
*****     How many will we move?
          move_rec = (num_move + 1) / 2 + 1
 
*****     Which way?
          move_rec = (-1) ** num_move * move_rec
 
*****     Calculate the next record number
          krec2 = start + move_rec
 
*****     We've moved another record
          num_move = num_move + 1
 
*****     Are we over or under?
          if (krec2 .lt. 1)             under = .true.
          if (krec2 .gt. sx_num_obs(2)) over  = .true.
 
*****     Check for a match if we are in bounds
          if (krec2 .ge. 1 .and. krec2 .le. sx_num_obs(2)) then
 
*****         Check for a match
              bounce_search = check_match(krec1,krec2)
 
          end if
 
      end do
 
***** If we made no match, assign 0 to krec2
      if (.not. bounce_search) krec2 = 0
 
      end
 
CTITLE SET_KALOBS_FLAGS
 
      subroutine set_KalObs_flags(krec1,krec2)
 
*     J.L. Davis                   5:25 PM  FRI., 17  APR., 1987
*
*     Routine to set the KalObs flags for this observation
*
 
      include '../includes/kalman_param.h'
      include 'calc_wl.h'
 
*       bit_off             - Bit offset
*   ,   file                - File index
*   ,   krec1, krec2        - The record numbers in the EMA list.
*                           -   A zero indicates no match
*   ,   ka(2)               - Array that krec numbers are stored
*   ,   one                 - Character value for '1'
*   ,   other_file          - Index of other file
*   ,   other_qual          - FRNGE quality of other observation
*   ,   seven               - Character value for '7'
*   ,   type                - 1 = group, 2 = phase
*   ,   zero                - Character value for '0'
 
      integer*4 bit_off, file, krec1, krec2, ka(2), one, other_file,
     .    other_qual, seven, type, zero
 
      data one   / o'61' /
      data seven / o'67' /
      data zero  / o'60' /
 
***** Store krec numbers (for easy index)
      ka(1) = krec1
      ka(2) = krec2
 
***** Loop over files
      do file = 1, 2
 
*****   Check if there is a record for the first file
        if (ka(file) .ne. 0) then
 
*****     What is the other file number?
          other_file = 3 - file
 
*****     Loop over group and phase type
          do type = 1, 2
 
*****       Were we told to add this?
            if (calc_WL_cont(type)) then
 
*****         What is the bit offset?
              bit_off = (type - 1) * 6
 
*****         Bits 2/8:  No matching data
              if (ka(other_file) .eq. 0) then
                call sion(sx_ion_flag(ka(file),file),bit_off+2,1)
              else
                call sion(sx_ion_flag(ka(file),file),bit_off+2,0)
              end if
 
*****         Bits 3/9:  Matching observation has quality code of 1-7
*             Bits 6/12: Matching observation has quality code of 0
              other_qual = sx_FRNGE_code(ka(other_file),
     .                         other_file)
              if (other_qual .ge. one .and. other_qual .le. seven) then
                call sion(sx_ion_flag(ka(file),file),bit_off+3,1)
              else
                call sion(sx_ion_flag(ka(file),file),bit_off+3,0)
              end if
 
              if (other_qual .eq. zero) then
                call sion(sx_ion_flag(ka(file),file),bit_off+6,1)
              else
                call sion(sx_ion_flag(ka(file),file),bit_off+6,0)
              end if
 
*****         Bits 4/10: Ion correction available AND
*             Bits 5/11: Downweight flag for ion correction
              if (ka(other_file) .ne. 0) then
                call sion(sx_ion_flag(ka(file),file),bit_off+4,1)
                call sion(sx_ion_flag(ka(file),file),bit_off+5,0)
              else
                call sion(sx_ion_flag(ka(file),file),bit_off+4,0)
                call sion(sx_ion_flag(ka(file),file),bit_off+5,1)
              end if
 
            end if
 
          end do
 
        end if
 
      end do
 
      end
 
CTITLE SION
 
      subroutine sion(flag,bit,value)
 
*     J.L. Davis                   6:15 PM  FRI., 17  APR., 1987
*
*     Version of sbit needed to set a bit in the ion flag stored in
*     EMA
*
*     Restrictions: 0 < bit < 33
 
 
 
*       bit                 - Which bit to set
*   ,   flag                - EMA word
*   ,   temp_flag           - Non EMA word used in SBIT call
*   ,   value               - 0 or 1
 
      integer*4 bit, flag, temp_flag, value
 
***** Set temporary word
      temp_flag = flag
 
***** Call SBIT with non EMA word
      call sbit(temp_flag,bit,value)
 
***** Move back to EMA
      flag = temp_flag
 
      end
 
CTITLE EFFECTIVE_FREQS
 
      subroutine effective_freqs(ob,f)
 
*     J.L. Davis                   3:14 PM  SAT., 18  APR., 1987
*
*     Calculates the effective frequencies for the group delay, phase
*     delay, and phase delay rate and stores in the EMA array
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
      include 'calc_wl.h'
 
 
*       f                           - File number
*   ,   j                           - Loop counter
*   ,   ob                          - Observation number in EMA
 
      integer*4 f, j, ob
 
*       del_f                       - f - phs_ref_freq
*   ,   denom                       - Denominator
*   ,   freq                        - "Absolute" frequency
*   ,   num_ch                      - Real*8 num_channels
*   ,   numer                       - Numerator
*   ,   sum_del_f                   - Sum[ delta_f ]
*   ,   sum_del_f2                  - Sum[ delta_f^2 ]
*   ,   sum_inv_f                   - Sum[ 1 / f ]
*   ,   sum_del_f_ov_f              - Sum[ delta_f / f ]
*   ,   sum_f2                      - Sum[ f^2 ]
 
      real*8 del_f, denom, freq, num_ch, numer, sum_del_f, sum_del_f2,
     .    sum_inv_f, sum_del_f_ov_f, sum_f2
 
***** Initialize summing variables
      sum_del_f      = 0
      sum_del_f2     = 0
      sum_inv_f      = 0
      sum_del_f_ov_f = 0
      sum_f2         = 0
 
***** Loop over channels
      do j = 1, num_channels
 
*****     Increment variables
*                                          ! Phase ref freq has already
          del_f          = freq_seq(j)
*                                          !   been subtracted by READIN
          freq           = phs_ref_freq + del_f
 
          sum_del_f      = sum_del_f      + del_f
          sum_del_f2     = sum_del_f2     + del_f ** 2
          sum_inv_f      = sum_inv_f      + 1.0D+00 / freq
          sum_del_f_ov_f = sum_del_f_ov_f + del_f / freq
          sum_f2         = sum_f2         + freq ** 2
 
      end do
 
***** Number of channels in double precision
      num_ch = dble(num_channels)
 
***** Group delay effective frequency
      numer = num_ch * sum_del_f2 - sum_del_f ** 2
      denom = sum_del_f * sum_inv_f - num_ch * sum_del_f_ov_f
 
      sx_eff_frq(1,ob,f) = sqrt(numer / denom)
      sx_phs_ref(f) = phs_ref_freq
 
***** Phase delay effective frequency
      numer = phs_ref_freq * numer
      denom = sum_del_f2 * sum_inv_f - sum_del_f * sum_del_f_ov_f
 
      sx_eff_frq(2,ob,f) = sqrt(numer / denom)
 
***** Phase delay rate effective frequency
      sx_eff_frq(3,ob,f) = sqrt(sum_f2 / num_ch)
 
      end
 
CTITLE REPORT_NO_MATCH
 
      subroutine report_no_match(lu,rec)
 
*     J.L. Davis                   4:07 PM  SAT., 18  APR., 1987
*
*     Subroutine to report that no match was found for a record
 
 
*       lu              - Lu for output
*   ,   rec             - Record number
 
      integer*4 lu, rec
 
      write(lu,100) rec
  100 format(' No match found for record #',I5)
 
      end
 
CTITLE CHECK_NUM_OBS
 
      subroutine check_num_obs(max_obs,KalObs_name)
 
*     J.L. Davis                   4:11 PM  SAT., 18  APR., 1987
*
*     Routine to check that the number of observations for a given
*     KalObs file is not too great
*
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       error           - Error flag
*   ,   len_name        - Length of KalObs file name
*   ,   max_obs         - Maximum number of observations allowed
*   ,   TrimLen         - HP function
 
      integer*4 error, len_name, max_obs, TrimLen
 
*       KalObs_name     - KalObs file descriptor
 
      character*(*) KalObs_name
 
 
      character*13 prog
 
      data prog / 'CHECK_NUM_OBS' /
 
***** Check num_obs against max_obs
      if (num_obs .gt. max_obs) then
 
*****     Get the length of the KalObs name
          len_name = max(TrimLen(KalObs_name),1)
 
*****     Write a message
          write(*,100) num_obs,max_obs,KalObs_name(1:len_name)
  100     format(/,' <<<< Error >>>> Number of observations exceeds',
     .           /,'                 maximum allowed.',
     .          //,' Number of observations:  ',I5,
     .           /,' Maximum allowed:         ',I5,
     .           /,' KalObs file name:        ',A,/)
 
*****     Close the KalObs file
          call close_KalObs(error)
 
*****     Check for error
          call report_error('FMGR',error,'clos',KalObs_name,0,prog)
 
*****     Stop the program
          stop ' Abort:  Too many observations'
 
      end if
 
      end
 
CTITLE WRITE_TO_KALOBS
 
      subroutine write_to_KalObs(f)
 
*     J.L. Davis                   7:47 PM  SAT., 18  APR., 1987
*
*     Write results to KalObs files
*
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
      include 'calc_wl.h'
 
*       error           - Error flag
*   ,   f               - File number (1 or 2)
*   ,   ob              - Loop counter
*   ,   TrimLen         - HP function
 
      integer*4 error, f, ob, TrimLen
 
*       prog            - Routine name
 
      character*15 prog
 
      data prog / 'WRITE_TO_KALOBS' /
 
***** Write a message
      write(*,100) sx_KalObs(f)(1:TrimLen(sx_KalObs(f)))
  100 format(/,' Writing out ',A)
 
***** Open the KalObs file
      call open_KalObs(sx_KalObs(f),error)
 
***** Check for error
      call report_error('FMGR',error,'open',sx_KalObs(f),1,prog)
 
***** Read the header
      call rw_KalObs_header('R',error)
 
***** Check for error
      call report_error('FMGR',error,'read','header',1,prog)
 
***** Loop over the observations
      do ob = 1, num_obs
 
*****     Read a record
          call rw_KalObs_block('R','data',site,error,ob)
 
*****     Check for error
          call report_error('FMGR',error,'read',sx_KalObs(f),1,prog)
 
*****     Store ion info back into KalObs record
          call store_ion_info(ob,f)
 
*****     Write observation back to disk
          call rw_KalObs_block('W','data',site,error,ob)
 
*****     Check for error
          call report_error('FMGR',error,'writ',sx_KalObs(f),1,prog)
 
      end do
 
***** Now change any pertinent flags in the KalObs header
      call change_header
 
***** Rewrite header to file
      call increment_KalVer
      call rw_KalObs_header('W',error)
 
***** Check for error
      call report_error('FMGR',error,'writ','header',0,prog)
 
***** Close the KalObs file
      call close_Kalobs(error)
 
***** Check for error
      call report_error('FMGR',error,'clos',sx_KalObs(f),0,prog)
 
      end
 
CTITLE STORE_ION_INFO
 
      subroutine store_ion_info(ob,f)
 
*     J.L. Davis                   8:04 PM  SAT., 18  APR., 1987
*
*     Routine to store ion info for observation OB and file
*     F into KalObs observation record
*
 
      include '../includes/kalman_param.h'
      include '../includes/obs_data.h'
      include 'calc_wl.h'
 
*       f               - File number
*   ,   ion_rec         - The entries in SX_ION_CONT and
*                       -   SX_ION_SIGMA to use
*   ,   ob              - Observation number
*   ,   type            - Loop counter
 
      integer*4 f, ion_rec, ob, type
 
*       f_ref           - The effective frequency of the first file
*   ,   f_cont          - The effective frequency we would like
*                       -   the contribution to be at
*   ,   ion_con_calc    - The calculated ionospheric contribution
*                       -   at the frequency of the first file.
*   ,   ion_sig_calc    - The calculated ion sigma at the frequcny
*                       -   of the first file
*   ,   ratio           - See below
 
      real*8 f_ref, f_cont, ion_con_calc, ion_sig_calc, ratio
 
***** Where do we get the ionospheric delay from?  If this is file #1,
*         the ionospheric contributions and sigmas are in EMA record
*         #OB, and are referred to the frequency of file #1.  If this
*         is file #2, the ionospheric contributions and sigma are stored
*         in record  REC_FILE_1(OB), and are referred to the frequency
*         of file #1.
 
      if (f .eq. 1) then
 
*****     The information is stored in record OB
          ion_rec = ob
 
      else
 
*****     The information is stored in REC_FILE_1(OB)
          ion_rec = rec_file_1(ob)
 
      end if
 
***** Loop over contributions
      do type = 1, 3
 
*****     Do we want to store this type and does a record exist for this
*             correction?
          if (calc_WL_cont(type) .and. ion_rec .gt. 0) then
 
*****         What is the ion contribution and the sigma?
              ion_con_calc = sx_ion_cont(type,ion_rec)
              ion_sig_calc = sx_ion_sigma(type,ion_rec)
 
*****         What frequency is the contribution referred to?
              f_ref = sx_eff_frq(type,ion_rec,1)
 
*****         What frequency do we want to refer the contribution to?
              f_cont = sx_eff_frq(type,ob,f)
 
*****         The square of the ratio of the above is handy.  Note that
*                 for F=1, RATIO = 1
              ratio = (f_ref / f_cont) ** 2
 
*****         Put the contribution and sigma into the KalObs record at
*                 the correct frequency
              ion_corr(type)  = ratio * ion_con_calc
              ion_sigma(type) = ratio * ion_sig_calc
 
          end if
 
      end do
 
***** Store the ionosphere flag
      ion_flag = sx_ion_flag(ob,f)
 
      end
 
CTITLE CHANGE_HEADER
 
      subroutine change_header
 
*     J.L. Davis                   2:19 PM  SUN., 19  APR., 1987
*
*     Make the following changes to the KalObs header:
*
*     (i) Set bits 6 and/or 7 in AVAIL_BASELINE (ion available)
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include 'calc_wl.h'
 
*       bit             - Bit in AVAIL_BASELINE to set
 
      integer*4 bit, type
 
***** Loop over group and phase delay types
      do type = 1, 2
 
*****     What bit should be set for this data type?
          bit = 5 + type
 
*****     Set the bit if we calculated the ion delay
          if (calc_WL_cont(type)) call sbit(avail_baseline,bit,1)
 
      end do
 
      end
