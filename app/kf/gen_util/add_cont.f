CTITLE ADD_CONT
 
      subroutine add_cont( mask, avail, type, theoretical, cont)

      implicit none 
 
 
*     Routine to add a baseline dependent contribution to the
*     theoretical delays.  If the contribution is not available
*     the contributions mask is modified so that this contribution
*     will not be added in future calls.
*     All four theoretical values are modified (Group, Phase, SB
*     delays and rate)
 
*   avail   - bit mapped word which indicates the availabilty
*           - of particular contributions
*   mask    - mask which indicates which contributions are to
*           - be applied.
*   type    - number which gives which bit to be checked for
*           - this contribution
 
      integer*4 avail, mask, type
 
*   kbit    - Function to check if bit is set
 
      logical kbit
 
*   cont(2) - the contribution to delay and rate (assumed
*           - non-dispersive)
 
      real*4 cont(2)
 
*   theoretical(4) - the four theoretical values (group delay,
*           - phase delay, single band delay and phase delay rate)
 
      real*8 theoretical(4)
 
***** Check to see if contributions is to be applied
 
*                                 ! we are to apply contribution
      if( kbit(mask,type) ) then
 
*         Check to see if available
*                                     ! Yes, it is
          if( kbit(avail,type) ) then
              theoretical(1) = theoretical(1) + cont(1)
              theoretical(2) = theoretical(2) + cont(1)
              theoretical(3) = theoretical(3) + cont(1)
              theoretical(4) = theoretical(4) + cont(2)
*                                     ! No no avaible, set mask so
          else
*                                     ! we try again
              call sbit(mask,type,0 )
          end if
 
      end if
 
***** Thats all
      return
      end
 
