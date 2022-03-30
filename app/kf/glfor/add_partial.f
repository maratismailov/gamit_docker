CTITLE ADD_PARTIAL
 
      subroutine add_partial( indx_pnt, part_pnt, parn, a_part, deriv)

      implicit none  
 
*     Routine to add the pointers and partials to the partials matrix
 
 
*   indx_pnt        - Index to number of partials for this local
*                   - parameter
*   i               - Loop counter
*   part_pnt(2,1)   - Pointers to which partials are to be used
*   parn            - parameter number of global parameter to
*                   - be tested
*   slot            - NeIndex to numbext index in the partials array (a_part)
 
      integer*4 indx_pnt, i, part_pnt(2,1), parn, slot
 
*   a_part(1)       - one row from the partials matrix
*   deriv           - Derivative to be put in the array
 
 
      real*8 a_part(1), deriv
 
 
 
***** See if we estimating the global parameter
 
*                                  ! Not estimating anyway
      if( parn.eq.0 ) RETURN
 
*                                  ! No partial
      if( deriv.eq.0.d0 ) RETURN
 
****  Start be checking if we have any partials yet
 
*                                 ! We have none so add this one as the first
      if ( indx_pnt.eq.0 ) then
          indx_pnt = 1
          part_pnt(1,1) = parn
          part_pnt(2,1) = 1
          a_part(1)     = deriv
*                                 ! At least one partial already, see if
*                                 ! new one is contiguous with old
      else
 
*         Find the slot in a_part which we have to use
          slot = 0
          do i = 1, indx_pnt
              slot = slot + part_pnt(2,i)
          end do
 
*         Put partial in next slot
          slot = slot + 1
 
*         Now see if we have to add another pair of indices
          if( part_pnt(1,indx_pnt)+part_pnt(2,indx_pnt)
*                                                   ! Yes, it is contiguous
     .                              .eq.parn ) then
*                                                             ! add new one
              part_pnt(2,indx_pnt) = part_pnt(2,indx_pnt) + 1
          else
*                                                   ! new pair being added
              indx_pnt = indx_pnt + 1
*                                                   ! Parmeter number
              part_pnt(1,indx_pnt) = parn
*                                                   ! Just one at the moment
              part_pnt(2,indx_pnt) = 1
          end if
 
*         Save the partial
          a_part(slot) = deriv
 
      end if
 
***** Thats all
      return
      end
 
