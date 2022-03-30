CTITLE sort_eq
 
      subroutine sort_eq

      implicit none 
 
*     This routine will sort the earthquake list into ascending
*     time order.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   i,j,k       - Loop counters
*   temp_i4     - I*4 variaable
 
      integer*4 i,j,k, temp_i4
 
*   temp_8      - Temporary Real*8 for switching
 
      real*8 temp_8
 
*   temp_ch     - Temporary character for switching.
 
      character*8 temp_ch
 
***** Now do sort
      do i = 1, num_eq-1
          do j = 1, num_eq-i
*                                                ! Time sort 
              if( eq_epoch(j).gt.eq_epoch(j+1) ) then
                  call switch_ch(eq_codes(j), eq_codes(j+1),
     .                           temp_ch )
                  call switch_8(eq_epoch(j), eq_epoch(j+1), temp_8)
                  call switch_8(eq_rad(j), eq_rad(j+1), temp_8)
                  call switch_8(eq_depth(j), eq_depth(j+1), temp_8)
                  call switch_8(eq_log_tau(j), eq_log_tau(j+1), temp_8)

                  do k = 1,2
                      call switch_8(eq_dur(k,j),
     .                              eq_dur(k,j+1),temp_8)
                  end do
                  do k = 1,3
                      call switch_8(eq_pos(k,j),
     .                              eq_pos(k,j+1),temp_8)

                  end do
                  do k = 1,6
                      call switch_8(eq_apr_coseismic(k,j),
     .                              eq_apr_coseismic(k,j+1),temp_8)
                      call switch_8(eq_mar_pre(k,j),
     .                              eq_mar_pre(k,j+1),temp_8)
                      call switch_8(eq_mar_post(k,j),
     .                              eq_mar_post(k,j+1),temp_8)
                      call switch_8(eq_log_sig(k,j),
     .                              eq_log_sig(k,j+1),temp_8)
                  end do
                  call switch_I4(eq_rename(j), eq_rename(j+1), temp_i4)
              end if
          end do
      end do
 
***** Thats all
      return
      end
 
