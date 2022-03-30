      integer function jel(i,j)
c
c     J.L. Davis 870503
c
c     Routine to compute the index of the element for a matrix stored
c     in lower triangular form, whose indices in normal (expanded)
c     form are (I,J).
c
c     Input variables
c     ---------------
c     i, j            The indices
c
      integer i, j
c
      if (i .gt. j) then
c
          jel = i*(i-1)/2 + j
c
      else
c
          jel = j*(j-1)/2 + i
c
      end if
c
      end
