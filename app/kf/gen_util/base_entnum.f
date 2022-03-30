ctitle
 
      subroutine base_entnum(site, ic, nstep, entry)

      implicit none 
 
c
c     This routine will compute the entry number in the table
c     of baselines for a site pair when nstep data types
c     are being used.
c
c Variables
c ---------
c site -- the site pair.  May be in any order.
c ic -- the counter through the loop over values
c nstep -- the number of data types used
c entry -- the entry in the table.  The entry is computed in the
c     the same manner as the entries for a lower diagonal matrix
c     storage with a multiplication is nstep is greater than one.
c
      integer*4 site(2), ic, nstep, entry
 
c
c Local variables
c ---------------
c j and k -- the largerst and smallest site numbers
c base_num -- the number of the baseline
c modulo -- tells which modulus of nstep we are are on
c
      integer*4 j,k, base_num, modulo
 
c
c.... compute the baseline/observable number
      j = max(site(1),site(2))
      k = min(site(1),site(2))
      base_num = (j-1)*(j-2)/2 + k
c
c.... Now increase number to allow for multiple observable types
*                                            ! value will be 1 for delay
      modulo = ic - int((ic-1)/nstep)*nstep
c                                                            2 for rate
c.... Compute the entry
      entry = (base_num-1)*nstep + modulo
c
c.... Thats all
      return
      end
 
c.......................................................................
