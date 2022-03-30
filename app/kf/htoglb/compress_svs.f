CTITLE COMPRESS_SVS
 
      subroutine compress_svs( obs )
 
      implicit none

*     This routine will remove from the list of satellites those
*     which have no observations.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
*   obs(num)    - Number of obsevation per satellites
 
      integer*4 obs(qnum_svs)
 
* LOCAL VARIBALES
 
*   i,j     - Loop counters
 
      integer*4 i,j, k
 
****  Scan over list of satelites finding those with no observations
 
      i = 0
      do while (i.lt.qnum_svs )
          i = i + 1
*                                 ! Remove from list
        if( obs(i).eq.0 ) then
              do j = i, qnum_svs-1
                  obs(j) = obs(j+1)
                  qsvs_names(j)  = qsvs_names(j+1)

*                 Move the other information down as well
                  qsvi_prn(j)    = qsvi_prn(j+1)
                  qsvi_block(j)  = qsvi_block(j+1)
                  qsvi_antmod(j) = qsvi_antmod(j+1)
                  qsvi_ocode(j)  = qsvi_ocode(j+1)
                  do k = 1,3
                     qsvi_antpos(k,1,j) = qsvi_antpos(k,1,j+1)
                     qsvi_antpos(k,2,j) = qsvi_antpos(k,2,j+1)
                  enddo
              end do
 
*             Decrement qnum_svsber and pointer
              qnum_svs = qnum_svs - 1
              i = i - 1
          end if
      end do
 
****  Thats all
      return
      end
 
 
 
