
c--------------------------------------------------------------
      subroutine submtx(n,m,nd,md,amat,bmat,mode1,mode2)
c
c     mode2 = 1:  subtract a submatrix from glabal matrix
c     mode2 = 2:  copy a submatrix to glabal matrix
c     mode2 = 3:  add a submatrix to glabal matrix
c     amat:    glabal matrix(symmetric)
c     bmat:    submatrix
c     n, nd:   start row index and row increment
c     m, md:   start column index and column increment
c     mode1 = 1:  symmetric part
c     mode1 = 2:  asymmetric part
c
      real*8 amat,bmat
      integer n,m,nd,md,i,j,i1,i2,j1,id,mode1,mode2
      dimension amat((n+nd)*(n+nd+1)/2),bmat(nd*md)
      i2 = 0
      do 20 i = 1,nd
         i1 = i+n-1
         do 10 j = 1,md
            if (mode1.eq.1.and.j.gt.i) goto 20
            j1 = j+m-1
            i2 = i2+1
            if(i1.ge.j1) id = i1*(i1-1)/2+j1
            if(j1.gt.i1) id = j1*(j1-1)/2+i1
            if(mode2.eq.1) bmat(i2) = amat(id)
            if(mode2.eq.2) amat(id) = bmat(i2)
            if(mode2.eq.3) amat(id) = amat(id)+bmat(i2)
 10      continue
 20   continue
c
      return
      end
