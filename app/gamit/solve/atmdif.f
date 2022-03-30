      subroutine atmdif(numsit,difop,pntdif,rowdif)
c
c     J.L. Davis  870309
c     mod:  J.L. Davis, M.H. Murray 870729    corrected bug
c           due to including all pairs of differences.
c           Changed so that only independent differences are used.
c
c     Builds difference operator for atmospheric constraints.
c
      include '../includes/dimpar.h'
c
c     Input variables
c     ---------------
c     numsit          The number of sites
c
      integer numsit
c
c     Output variables
c     ----------------
c     difop           The difference operator matrix, in "Bock
c                     packed format."
c     pntdif          The array which points to non-zero elements
c     rowdif          The row pointer
c
c    mod: correct dimensions on pntdf (previously dimop) and
c         rowdft (previously maxdif)   M. Murray 880309
c
      integer pntdif(2*(maxsit-1)), rowdif(maxsit+1)
c
      real*8 difop(*)
c
c     Internal variables
c     ------------------
c     i, j            Loop indices
c     row             Row in operator matrix
c
      integer
     .            i
     .    ,       row
c
c.... Initialize row pointer.  The first element is the total number
c     of non-zero elements in the operator.  For NUMSIT stations,
c     there are NUMSIT * (NUMSIT - 1) / 2 differences.  For each
c     difference, there are two non-zero elements in the operator
c     matrix (i.e., +1 & -1).
c
      rowdif(1) = numsit * (numsit - 1)
      rowdif(2) = 0
c
c.... Count the row
      row = 0
c
c.... Loop over sites
c     mod:   loop over only independent pairs
      do 100 i = 2, numsit
c
c....     What is the row number?
          row = row + 1
c
c....     Store the +/- 1's in DIFOP
          difop(2 * row - 1) = +1
          difop(2 * row    ) = -1
c
c....     The indices I and J are the columns for +1 and -1
          pntdif(2 * row - 1) = i
          pntdif(2 * row    ) = 1
c
c....     There are 2 non-zero elements for this row
          rowdif(2 + row) = rowdif(1 + row) + 2
c
  100 continue
c
      end
