CTITLE    .................................................................
 
      subroutine inc_cax(idate, intv, indx,tax_data,comjd)
c
c     Routine to increment the tic mark positions for a calender
c     axis.
c
c Variables 
c ---------
c idate -- the array containing date to incremented
c intv  -- the increment value
c indx  -- the index in tax_data for the incremet
c tax_data -- gives turn over and start values
c comjd -- the julian date (used to ensure valis date)
c
      integer*4 idate(6), intv, indx, tax_data(2,6)
 
c
      real*8 comjd, sectag
 
c
c.... Increment the date value by the amount given by intv
      sectag = 0.d0
      idate(indx) = idate(indx) + intv
c
c.... see if value has turned over
*                                                 ! yes we have gone over
      if( idate(indx).ge. tax_data(1,indx) .and. 
     .    indx.gt.1 ) then
c                                                   next unit
c....    reset the counters
         idate(indx-1) = idate(indx-1) + 1
         idate(indx)   = tax_data(2,indx)
c
c....    Make sure didnot turn over on indx-1
         if( idate(indx-1).ge.tax_data(1,indx-1) .and.
     .       indx.gt.2  ) then
            idate(indx-2) = idate(indx-2) + 1
            idate(indx-1) = tax_data(2,indx-1)
         end if
      end if
c
c.... Reconvert date to julian date to make sure we have valis date
      sectag = idate(6)
      call ymdhms_to_mjd(idate, sectag, comjd) 
 
      return
      end
 
