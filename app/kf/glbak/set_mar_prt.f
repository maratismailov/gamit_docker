CTITLE SET_MAR_PRT
 
      subroutine set_mar_prt(loc_ind_mar)

      implicit none 
 
*     Routine to set the sites to be printed to the bakfile.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED VARIABLES
 
*   loc_ind_mar(3,gnum_sites)   - Parameter numbers of the
*                   - sites to be printed.
 
      integer*4 loc_ind_mar(3,gnum_sites)
 
* LOCAL Variables
 
*   i,j     - Look counters
 
 
      integer*4 i,j
 
*   eq_effected     - Logical function which returns
*                   - true if site effected by an earthquake
*                   - at this time.
*   kbit            - Function to check if bit is set
 
      logical eq_effected, kbit
 
***** Start: If the user has not said to clear the list then
*     add all the markov sites to the list
      if( .not.clear_bo_site ) then
          do i = 1, gnum_sites
              do j = 1,3
                  if( mar_site(j,1,i).ne.0 .or.
     .                (mar_neu(1,1,i)+mar_neu(2,1,i)+
     .                 mar_neu(3,1,i)) .ne.0) then
                      loc_ind_mar(j,i) = parn_site(j,1,i)
                  else
                      loc_ind_mar(j,i) = 0
                  end if
              end do
*             Check for local effects of earthquakes
              if( eq_effected(i) ) then
                  do j = 1,3
                      loc_ind_mar(j,i) = parn_site(j,1,i)
                  end do
              end if
          end do
      else
 
*         Clear the array
          do i = 1, gnum_sites
              do j = 1,3
                  loc_ind_mar(j,i) = 0
              end do
          end do
      end if
 
****  Now add the sites that have been specifically selected
*     by the user.
      do i = 1, gnum_sites
          if( kbit(bak_out_site,i) ) then
              do j = 1,3
                  loc_ind_mar(j,i) = parn_site(j,1,i)
              end do
          end if
      end do
 
****  Thats all
      return
      end
 
 
 
 
