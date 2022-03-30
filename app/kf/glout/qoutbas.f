 
      logical function qoutbas( options, s1, s2 )

      implicit none 
 
*     This function checks the options variable passed to
*     see if a baseline should be written.  If Bit 11 is
*     set (2048) then the baseline will be output for only
*     markov sites.  If Bit 11 is not set, functions always
*     returns true,
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
* PASSED VARIABLES
 
*   options     - Options passed.  If Bit 11 is set
*               - then only markov sites output.
*   s1, s2      - Two sites in the baseline.
 
      integer*4 options, s1, s2
 
* LOCAL VARIABLES
 
*   kbit        - Tests if bit is set (number is one
*               - greater than Bit Number.)
*   outs1, outs2    - Logicals to indicate if each site
*               - is markov
*   locs1, locs2    - Logicals to indicate if the station
*                 was used at this time
*   i           - Loop counter
*   eq_effected - Logical function that tests to see if
*                 station effected by an Earthquake.

      integer*4 i  

      logical kbit, outs1, outs2, eq_effected, locs1, locs2
 
****  Start, Check all the conditions that say we should
*     output this baseline.  Will be output if:
*     (1) Selected station (bak_prts command) 
*     (2) Effected by an earthquake.
*     (3) Markov site and we have said to print markov sites
*     
      qoutbas = .false.
      outs1 = eq_effected(s1) .or. kbit(bak_out_site,s1)
      outs2 = eq_effected(s2) .or. kbit(bak_out_site,s2)

*     See if should output because it is markov
      if( kbit(options,12) ) then
 
*         Test to see if each site has any markov elemnets
          do i = 1,3
              if( mar_site(i,1,s1).gt.0.d0 .or.
     .            mar_neu(i,1,s1).gt.0.d0 ) outs1 = .true.
          end do
 
*         Check second site
          do i = 1,3
              if( mar_site(i,1,s2).gt.0.d0 .or.
     .            mar_neu(i,1,s2).gt.0.d0 ) outs2 = .true.
          end do
      else
          outs1 = .true.
          outs2 = .true.
      end if

****  Now look at why it should not be output:  The only
*     reason so far:
*     (1) Not a named sited and site not used in this 
*         experiment. (Only appliies to back solution output)
      locs1 = .true.
      locs2 = .true.
      if( kbit(options,32) ) then
          locs1 = kbit(bak_out_site,s1)
          locs2 = kbit(bak_out_site,s2)
          do i = 1,cnum_sites
             if( ltog_sites(i).eq.s1 ) locs1 = .true.
             if( ltog_sites(i).eq.s2 ) locs2 = .true.
          end do
      end if

****  Now finally:  If we should output and it is local (or
*     a named site) set to output
      if( (locs1.and.outs1) .and. (locs2.and.outs2) ) then
          qoutbas = .true.
      end if

****  If the site has not been used then do not output.        
      if( .not.kbit(guse_site,s1) .or. 
     .    .not.kbit(guse_site,s2)      ) qoutbas = .false.
 
****  Thats all
      return
      end
 
