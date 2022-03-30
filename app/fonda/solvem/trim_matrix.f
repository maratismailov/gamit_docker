      subroutine trim_matrix(ncoef,a,id,m)
c
c     remove m rows and m columns from the id-th row,
c     and rearrange the matrix
c
      real*8 a
      dimension a((ncoef*(ncoef+1))/2)
      integer ncoef,istart,ic,idiag,ij,id,m,irow,icol
c
c     perform trim from the id-th row
      istart = id*(id-1)/2
      do 100 irow = id+m,ncoef
         idiag=irow*(irow-1)/2
         ij = 0
         ic = 0
         do 50 icol = 1,irow
            if (icol.ge.id.and.icol.lt.id+m) goto 50
            ij = ij+1
            if (icol.ge.id+m) ic = m
            a(istart+ij) = a(idiag+ij+ic)
 50      continue
         istart = istart+ij
 100  continue
c
      return
      end
