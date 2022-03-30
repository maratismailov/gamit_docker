CTITLE ASSIGN_EMA_VAL4
 
      subroutine assign_ema_val4( vals, num, destin, is, dim)
 
      implicit none 

 
*     Routine to assign the real*8 values in vals, into the real*4
*     arrays destin. If IS is equal to 999999 then all array elements
*     from 1 to dim are given values.
*     NUM is assumed to be the number of values to copy AND the
*     first dimension of the DESTIN
 
*   dim     - Second dimension of destin
*   i,j     - Loop counters
*   is      - site/source number for get_cmd (may be 999999)
*   num     - number of values to copy and first dimension of
*           - destin
*   start,stop  - start and stop site/source number if IS is 999999
 
      integer*4 dim, i,j, is, num, start,stop
 
*   destin(num,dim) - destimation values in ema.
 
      real*4 destin(num,dim)
 
*   vals(1)         - values to be copied
 
      real*8 vals(1)
 
 
***** First check that IS is greater than zero, meaning name found.
 
*                         ! value OK
      if( is.gt.0 ) then
 
*         Now see if ALL option specified
*                                 ! loop over all of dim
          if( is.eq.999999 ) then
              start = 1
              stop  = dim
*                                 ! just do the one IS value
          else
              start = is
              stop  = is
          end if
 
*****     Now copy values
          do i = start, stop
              do j = 1, num
                  destin(j,i) = vals(j)
              end do
          end do
      end if
 
***** Thats all
      return
      end
 
CTITLE ASSIGN_EMA_VAL8
 
      subroutine assign_ema_val8( vals, num, destin, is, dim)
 
      implicit none 

 
*     Routine to assign the real*8 values in vals, into the real*8
*     arrays destin. If IS is equal to 999999 then all array elements
*     from 1 to dim are given values.
*     NUM is assumed to be the number of values to copy AND the
*     first dimension of the DESTIN
 
*   dim     - Second dimension of destin
*   i,j     - Loop counters
*   is      - site/source number for get_cmd (may be 999999)
*   num     - number of values to copy and first dimension of
*           - destin
*   start,stop  - start and stop site/source number if IS is 999999
 
      integer*4 dim, i,j, is, num, start,stop
 
*   destin(num,dim) - destimation values in ema.
*   vals(1)         - values to be copied
 
      real*8 destin(num,dim), vals(1)
 
 
***** First check that IS is greater than zero, meaning name found.
 
*                         ! value OK
      if( is.gt.0 ) then
 
*         Now see if ALL option specified
*                                 ! loop over all of dim
          if( is.eq.999999 ) then
              start = 1
              stop  = dim
*                                 ! just do the one IS value
          else
              start = is
              stop  = is
          end if
 
*****     Now copy values
          do i = start, stop
              do j = 1, num
                  destin(j,i) = vals(j)
              end do
          end do
      end if
 
***** Thats all
      return
      end
 
