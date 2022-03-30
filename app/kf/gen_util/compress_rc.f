ctitle
 
      subroutine compress_rc(norm_eq, sol_vec, comp_array, nsize, ndim)

      implicit none 
 
c
c     routine to move the row and columns of norm_eq to the positions
c     specified in comp_array.  If the enter in comp_array is zero then
c     the row and col are not moved.  This routine assumes that
c     comp_array is in ascending order so that we do not need to
c     worry about overwritting any useful elements of norm_eq
c
c Variables
c ---------
c norm_eq -- the normal equations which will be compressed
c sol_vec -- the corresponding solution vector which will also be
c     compressed.
c comp_array -- array indicating where each row and column must be
c     moved to.  If the entry is zero then the row and column are not
c     moved.
c nsize -- the size of the matrix to be moved
c ndim  -- the dimensioned size of the matrix (needed becuase full
c     matrix is stored)
c
      integer*4 comp_array(1)
c
      integer*4 nsize, ndim, i
c 
      real*8 norm_eq(ndim,ndim), sol_vec(ndim)
c
c.... loop over matrix moving rows and columns if needed
      do i = 1, nsize
c
*                                      ! shift this row and col.
         if( comp_array(i).ne.0 ) then
c
            call dwmov(norm_eq(1,i),1, norm_eq(1,comp_array(i)),1,
*                           ! move the column
     .         nsize)
            call dwmov(norm_eq(i,1),ndim, norm_eq(comp_array(i),1),
*                           ! move the row
     .         ndim, nsize)
c
c....       now move the solution vector
            call dwmov(sol_vec(i),1, sol_vec(comp_array(i)),1, 1)
c                       ! only single value to be moved
         end if
c
      end do
c
c.... that's all
      return
      end
 
