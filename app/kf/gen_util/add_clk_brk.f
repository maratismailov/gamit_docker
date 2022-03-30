CTITLE ADD_CLK_BRK
 
      subroutine add_clk_brk(brk_epoch, brk_mag, brk_site, num_clk_brk,
     .    clk_brk_epoch, clk_brk_mag, clk_brk_site, max_clk_brk )

      implicit none 
 
*     Routine to add a new clock break to the list of clock breaks.  If
*     the number of breaks will exceed the max number allowed then the
*     smallest is deleted.  (NOTE: all of the magnitudes dectected real
*     time are positive, interactively entered values are given a magnitude
*     of -1 and will not be deleted.
 
*   max_clk_brk - Maximum number of breaks allowed
*   num_clk_brk - Current number of breaks (this value accumulates
*               - the total number found.  It will often exceed
*               - the max value (so should be checked before use)
*   brk_site        - Site with the current break
*   clk_brk_site(max_clk_brk)   - List of clock break sites
 
      integer*4 max_clk_brk, num_clk_brk, brk_site,
     .    clk_brk_site(max_clk_brk)
 
*   brk_mag     - Magnitude of current break
*   clk_brk_mag(max_clk_brk)    - List of previous breaks
 
      real*4 brk_mag, clk_brk_mag(max_clk_brk)
 
*   brk_epoch   - JD of current break
*   clk_brk_epoch(max_clk_brk)  - list of epochs of clock breaks
 
      real*8 brk_epoch, clk_brk_epoch(max_clk_brk)
 
* LOCAL VARIABLES
 
*   i           - Loop counter
*   used_num        - NUmber of breaks actual in tables.  Will be either
*               - num_clk_brk or max_clk_brk
*   pos_brk     - positition in table to add current break
*   brk_with_min    - NUmber of smallest break (used for deleting)
 
      integer*4 i, used_num, pos_brk, brk_with_min
 
*   min_brk     - Smallest break in list
 
      real*4 min_brk
 
****  Increment count and see if we execced the limit
      used_num = num_clk_brk
      num_clk_brk = num_clk_brk + 1
 
*                                             ! Number exceeded get smallest
      if( num_clk_brk.gt.max_clk_brk ) then
          brk_with_min = 0
          if( brk_mag.gt.0 ) then
              min_brk = brk_mag
          else
              min_brk = 1.d20
          end if
 
****      Scan list and find smallest greater then zero
          do i = 1, max_clk_brk
              if( clk_brk_mag(i).lt.min_brk .and.
*                                                         ! New min
     .            clk_brk_mag(i).gt.0            ) then
                  min_brk = clk_brk_mag(i)
                  brk_with_min = i
              end if
          end do
 
****      If we found a min then delete this break, else just return
*         we can't find a value which is smaller than the one we will
*         add.
*                                         ! Delete
          if( brk_with_min.gt.0 ) then
              do i = brk_with_min, max_clk_brk - 1
                  clk_brk_epoch(i) = clk_brk_epoch(i+1)
                  clk_brk_mag(i)   = clk_brk_mag(i+1)
                  clk_brk_site(i)  = clk_brk_site(i+1)
              end do
*                         ! Nothing left to do, so retrun
          else
              RETURN
          end if
*                                         ! Currently used number
          used_num = max_clk_brk - 1
*                         ! Maximium number of breaks exceeded.
      end if
 
****  Now add in the new break.  Look for its position in the time
*     sorted table
      i = 0
*                             ! In case this is the latest one we
      pos_brk = used_num + 1
*                             ! have.
      do while ( i.lt.used_num )
          i = i + 1
          if( brk_epoch.lt.clk_brk_epoch(i) ) then
              pos_brk = i
              i = used_num + 1
          end if
      end do

***** See if this is just of duplicate of and existing break
      if( abs(brk_epoch-clk_brk_epoch(pos_brk-1)).lt.1.d-5 .and.
     .    brk_site.eq.clk_brk_site(pos_brk-1) ) then
          num_clk_brk = num_clk_brk - 1
          RETURN
      end if

****  Make space for new value
      do i = used_num, pos_brk, -1
          clk_brk_epoch(i+1) = clk_brk_epoch(i)
          clk_brk_mag(i+1)   = clk_brk_mag(i)
          clk_brk_site(i+1)  = clk_brk_site(i)
      end do
 
***** Now add in the new value
      clk_brk_epoch(pos_brk) = brk_epoch
      clk_brk_mag(pos_brk)   = brk_mag
      clk_brk_site(pos_brk)  = brk_site
 
***** Thats all
      return
      end
 
 
 
