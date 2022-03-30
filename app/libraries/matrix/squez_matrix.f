 
      subroutine squez_matrix(norm_eq, sol_vec, nsize, ndim, nsol,
     .   option)

      implicit none
 
c
c     routine to compress/expand a full ndim*ndim matrix into the upper
c     nsize*nsize elements.  This procedure is used to minimize
c     page faults during VMA operations.
c
c Variables
c ---------
c norm_eq -- the matrix to be compressed or expanded (EMA)
c sol_vec -- the LHS of the solution (EMA)
c nsize   -- the size of the matrix elements which are actually
c     used.
c ndim    -- the dimension of the matrix
c nsol    -- indicates if solution to be computed
c option  -- character string which specifies if the matrix should
c     be compressed or expanded
c
c Local variables
c ---------------
c addr  -- a value which is used to store the address in the
c     matric to which elements will be moved
c
c
C     real*8 norm_eq(ndim,1), sol_vec(1)
      real*8 norm_eq(*), sol_vec(*)
 
c
*   i,j     - Loop counters
      integer*4 nsize, ndim, nsol, i

* MOD TAH 190520: Converted addr and addi to I*8 to support matrices
*     with more than 32767 rows and columns 
c
*         addr      - Element number in squeezed matrix
*         addi      - Element number in non-squeezed matrix
*         I8        - Numeric 1 as I*8 to force I*8 calc
      integer*8 addr, addi, I8
 
c
      character*(*) option

      data I8 / 1 /
 
c
c.... see if we need to move the matrix at all
      if( nsize.eq.ndim ) return
c
c.... see which way we need to move the matrix
      if( option.eq.'compress' ) then
c
c....    start compressing the matrix
         do i = 2, nsize
            addr = (i-I8)*nsize + 1
            addi = (i-I8)*ndim  + 1
            call dwmov(norm_eq(addi),1,norm_eq(addr),1,nsize)
         end do
c
c....    now move the sol_vec
         if( nsol.gt.0 ) then
            addr = (I*8*nsize)*nsize + 1
            call dwmov(sol_vec(1),1, norm_eq(addr),1, nsize)
         end if
c
*          ! we need to expand the matrix
      else
c
c....    loop through the matrix backwards, starting with the
c        solution vector.  In particular, note that even the
c        VIS moves are done in a backwards order to ensure
c        that no part of the matrix is over written.
         if( nsol.gt.0 ) then
*                                     ! element number in compressed
            addr = (I*8*nsize)*(nsize + 1)
c                                          matrix for last element in col.
            call dwmov(norm_eq(addr),-1, sol_vec(nsize),-1, nsize)
         end if
c
c....    now do the matrix.  We only need to go to second column because
c        first column is in the correct position.
         do i = nsize,2, -1
            addr = i*nsize
*                                       ! Row nsize, column i
            addi = (I8-1)*ndim + nsize
            call dwmov(norm_eq(addr),-1, norm_eq(addi),-1, nsize)
         end do
c
      end if
c
c     thats all
      return
      end
 
