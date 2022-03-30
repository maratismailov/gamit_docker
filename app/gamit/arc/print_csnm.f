CTITLE PRINT_CSNM

      subroutine print_csnm(type, norm, max_n, czone, ctess, stess)

      implicit none

*     Subroutine to output C and Snm gravity coefficents from degree 2 
*     to max_n
*
* INPUT:
       character*(*) type   ! Short string with label e.g., SETIDE for
                            ! for solid Earth tides, OCTIDE for ocean
       character*(*) norm   ! String with normalization type.
       integer*4 max_n      ! Maxumim n value
       real*8 czone(*)      ! Zonal coefficients
       real*8 ctess(*)      ! Non-zonal Cnm
       real*8 stess(*)      ! Non-zonal Snm

* LOCAL
       integer*4 i, j       ! Loop counters

****  OK: output description
      write(*, 110) trim(type), trim(norm), max_n
 110  format('Gravity CSnm for ',a,' Coefficients ',a,
     .       ' Max N 2-',i3)

      do i = 2,max_n
         do j = 0,i
            if( j.eq.0 ) then  ! Zonal
               write(*,150) trim(type), 
     .                    i,j, czone(i-1)*1.d11
 150           format(a,' C/S ',I2,1x,I2,1x,2(f15.4,'E-11'))
            else    ! Non-Zonal
                write(*,150) trim(type), 
     .                    i,j, ctess((i-1)*i/2+j-1)*1.d11,
     .                         stess((i-1)*i/2+j-1)*1.d11
            endif
         end do
      enddo

***** Thats all
      return
      end
 
