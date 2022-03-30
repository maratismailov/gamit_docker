CTITLE    .................................................................
 
      subroutine check_bound(values, default, ref_val)
c
c     Routine to check the boundaries of the plot
c
c Variables
c ---------
c values  -- the value of the latest point read in
c default -- the scales to be set and checked
c ref_val -- the reference values
c
      real*4 default(2)
 
c
      real*8 values(1), ref_val
 
c
c Local variables
c ---------------
c new_value -- the value of data with the reference value removed
c
      real*4 new_value
 
c
c Scratch common
c --------------
c
      common new_value
 
c
c Functions
c ---------
c min -- find mimimum values -- FTN77 routine
c max -- find maximum values -- FTN77 routine
c
c.... Compute the value of the data with reference value removed
      new_value = values(1) - ref_val
c
c.... Now check scales
      default(1) = min(default(1),new_value)
      default(2) = max(default(2),new_value)
      
c
      return
      end
 
