CTITLE    .................................................................
 
      subroutine first_tic1( vmin, ref_val, tic_space, start_tic)
c
c     Routine to compute the location of the first tic mark for
c     type 1 data, i.e, not time axis
c
c Variables
c ---------
c vmin -- the minumim value on the axis
c ref_val -- the reference vale removed from the data
c tic_space -- the spacing bewteen the tics
c start_tic -- the coordinate of the first tic mark
c
      real*4 vmin, start_tic, tic_space
 
c
      real*8 ref_val
 
c
c Local variables
c ---------------
c first -- the actual value of the minumum vale
c num_zero -- the number of tics we are away from zero
c first_tic -- the position of the first tic will reference value added
c test_ref_val -- Updated value of ref_val with the leading digits removed
c                 (needed when the refernce value is much greater than
c                  tic mark spacing)
 
      real*8 first, first_tic, test_ref_val, first_tic_count
 
c
      integer*4 num_zero
 
c
c.... Compute the actual value of minimum
      num_zero = 65536
      test_ref_val = ref_val
      first_tic_count = 0
      do while ( num_zero.gt. 32767 )
          first = test_ref_val + vmin
          num_zero = abs( first/tic_space )
          if( first.lt.0 ) then
             first_tic_count = first_tic_count - num_zero*tic_space
          else
             first_tic_count = first_tic_count + (num_zero+1)*tic_space
          end if
 
          test_ref_val = test_ref_val - first_tic_count
      end do
 
      first_tic  = first_tic_count
 
C     first = ref_val + vmin
C     num_zero = abs(first/tic_space)
c
c.... Get the position of the first tic
C     if( first.lt.0 ) then
C        first_tic = -num_zero*tic_space
C     else
C        first_tic = (num_zero+1)*tic_space
C     end if
c
c.... Now remove the reference value
      start_tic = first_tic - ref_val
c
      return
      end
 
