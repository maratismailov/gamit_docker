CTITLE OPT_RUNSTR
 
      subroutine opt_runstr
 
      implicit none

*     Routine to decode the new option driven runstring for svsp3
*
 
      include '../includes/const_param.h'
      include 'svsp3.h'
 
* LOCAL VARIABLES
 
*  len_run     - Length of the runstring element
*  trimlen     - Length of string
* date_start(5), date_stop(5) - Start and stop dates
*               for generating results
*  i           - Loop counter
*  rcpar       - Gets runstring entry
*  nr          - Current runstring number.
 
 
      integer*4 len_run, i, rcpar, nr
 
      logical done    ! Set true when empty runstring.

*  runstring   - Elements of runstring
 
      character*256 runstring
 

* SVSP3 command line options
* -proc_noise <position m^2/ep> <clock m^2/ep> <range noise (m)>
* -ref_apr <X> <Y. <Z> -- Apriori coordinates for Reference (2nd) site
* -dat_apr <X> <Y. <Z> -- Apriori coordinates for first site.
* -out  <out spacing> <out type (XYZ/NEU)> 
* -debug <start epoch> <end epoch>
* -span <start epoch> <end epoch>
* -gnss <systems>  G-GPS, R-Glonass, E-Galileo, C-Beidou
* -anal_type <Type> From PC, P1, P2, P5
*

 

****  Start decoding
      nr = 3    ! 1 back from current position in runstring
      done = .false.
      do while ( .not. done ) 
         nr = nr + 1
         len_run = rcpar(nr, runstring)
         if ( len_run.eq. 0 ) then
            done = .true.
         else
*           Start decoding
            call caaefold(runstring)
            if( runstring(1:2).eq.'-P' ) then
*              -proc_noise <position m^2/ep> <clock m^2/ep> <range noise (m)>
               nr = nr + 1
               len_run = rcpar(nr, runstring)
               read(runstring,*) wn(1)
               wn(2) = wn(1)
               wn(3) = wn(1)

*              Get clock process noise
               len_run = rcpar(5, runstring)
               read(runstring,*) wn(4)
               wn(5) = 0.000d0
               wn(6) = 0.000d0

*           Get data noise
            len_run = rcpar(6, runstring)
            if( len_run.gt.0 ) then
                read(runstring,*) data_noise 
            else
                data_noise = 100.d0
            end if


            elseif ( runstring(1:2).eq.'-R' ) then
*              -ref_apr <X> <Y. <Z> -- Apriori coordinates for Reference (2nd) site
            elseif ( runstring(1:3).eq.'-DA' ) then
*              -dat_apr <X> <Y. <Z> -- Apriori coordinates for first site.
            elseif ( runstring(1:2).eq.'-O' ) then
*               -out  <out spacing> <out type (XYZ/NEU)> 
            elseif ( runstring(1:3).eq.'-DE' ) then
*               -debug <start epoch> <end epoch>
            elseif ( runstring(1:2).eq.'-S' ) then
*               -span <start epoch> <end epoch>
            elseif ( runstring(1:2).eq.'-G' ) then
*               -gnss <systems>  G-GPS, R-Glonass, E-Galileo, C-Beidou
            elseif ( runstring(1:2).eq.'-A' ) then
*               -anal_type <Type> From PC, P1, P2, P5
            else
                 write(*,'("UNKNOWN OPTION ",a)') trim(runstring)
            endif
         endif
      end do

***** Thata all
      return
      end

 
