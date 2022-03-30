
c--------------------------------------------------------------
      subroutine fill_apr_mtx(l,lc,n,m,nd,md,amat,bmat,mode2)
c
c     mode2 = 1:  subtract a submatrix from glabal matrix
c     mode2 = 2:  copy a submatrix to glabal matrix
c     mode2 = 3:  add a submatrix to glabal matrix
c     amat:    glabal matrix(symmetric)
c     bmat:    submatrix
c     n, nd:   start row index and row increment
c     m, md:   start column index and column increment
c
      real*8 amat,bmat
      integer l,lc
      integer n,m,nd,md,i,j,i1,j1,i2,id,mode2,i3
      dimension amat(l*lc),bmat(nd*md)
      i2 = 0
      do 20 i = 1,nd
         i1 = i+n-1
         if (i1.le.lc) j1 = i1*(i1-1)/2
         if (i1.gt.lc) j1 = lc*(lc-1)/2+(i1-lc)*lc
         if (i1.le.lc) i3 = 0
         if (i1.gt.lc) i3 = lc-i
         do 10 j = 1,md
            if (j.gt.i) goto 20
            i2 = i2+1
            id = j1+j+i3
            if(mode2.eq.1) bmat(i2) = amat(id)
            if(mode2.eq.2) amat(id) = bmat(i2)
            if(mode2.eq.3) amat(id) = amat(id)+bmat(i2)
 10      continue
 20   continue
c
      return
      end
