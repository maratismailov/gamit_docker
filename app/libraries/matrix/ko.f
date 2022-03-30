C     Program test_invert
c
c
c     Program to test matrix inversion

      implicit none
c
*    max_dim		- Maxium size of matrix
 
      integer*4 max_dim		
*
*                                     ! 500 x 500 matrix
      parameter ( max_dim = 500 )
c
*   start, stop                 - Start and stop times
*   v1, v2                      - Agruments for start and stop
*   runtime                     - Runtime
*   second                      - Function for time
      real*4 start, stop, v1, v2, runtime, second
*
c
*   Norm_eq(max_dim, max_dim)   - Matrix to be inverted
*   Sol_vec(max_dim)	      - solution vector
*   scale(max_dim)              - scaing vector
*   rms_sum                     - average error in inverse
      real*8 Norm_eq(max_dim, max_dim), Sol_vec(max_dim)	,
     .    scale(max_dim), rms_sum
*
c
*   i,j      - Loop counters
*   ipivot(max_dim)    - Pivot elements
*   nsize                       - Size of matrix
      integer*4 i,j, ipivot(max_dim), nsize
*
c
***** Get the size of matrix
c
      write(*,100)
  100 format(/,' Enter size of matrix '$)
      read(*,*) nsize
c
      do i = 1, nsize
         sol_vec(i) = i
         do j = i, nsize
             norm_eq(i,j) = i
         end do
         do j = 1,i
             norm_eq(i,j) = j
         end do
      end do
 
      start = second(v1)
      write(*,'(" Start inversion")')
      call invert_vis( norm_eq, sol_vec, scale, ipivot, nsize, max_dim,
     .                 1 )
      stop = second(v2)
      runtime = stop - start
      write(*,150) start, stop, v1, v2, runtime
  150 format(' START ',f10.2, ' STOP ',f10.2,' V1,v2 ',2f10.4,/,
     .       ' RUNTIME ',f10.3)
c
      if ( nsize.lt.6 ) then
         do i = 1, nsize
            write(*,'("Row ",i3)') i
             write(*,200) (norm_eq(i,j),j=1,nsize)
         end do
c
         write(*,200) (sol_vec(i),i=1,nsize)
  200    format(5f10.6)
c     else
c        write(*,250) (norm_eq(i,nsize),i=1,nsize)
c        write(*,250) (sol_vec(i), i=1,nsize)
c 250    format(5d15.6)
      end if
c
      rms_sum = 0.d0
      do i = 1,nsize-1
         rms_sum = rms_sum + abs(sol_vec(i))
      end do
      write(*,300) rms_sum/(nsize-1)
  300 format(' AVERAGE error in solution ',d15.6)
c
      rms_sum = 0
      do i = 1, nsize-2
         rms_sum = rms_sum + abs(norm_eq(i,nsize))
      end do
      write(*,310) rms_sum/(nsize-2)
  310 format(' AVERAGE error in matrix  ',d15.6)
      stop ' Test_matrix Terminated'
c
      END

      real*4 function second()
      real*4 dtime, targ(2), dummy

      dummy  = dtime (targ)
      second = targ(1)
      return
      end

