CTITLE EMASIZE
 
 
      Integer*4 function emasize( start, end )
 
      implicit none 
 
*     Utility to find the size of an ema block delimited by the variables
*     start and end.  This is done by assigning -999 to end and then
*     searching for this value in memory. (1/65535 chance that we will
*     find the wrong value!)
 
*   end     - Last word of block
*   start   - First word of block
 
      integer*4 end, start
 
*   addressof     - Returns address of variable
 
* MOD TAH 190511: Changed to I*8 for 64-bit memory
      integer*8 addressof
 
 
***** Initialize end
 
      emasize = (addressof(end) - addressof(end))/4 + 1
 
***** Thats all
      return
      end
 
 
 
