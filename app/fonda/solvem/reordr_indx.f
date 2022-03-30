c
      subroutine reordr_indx(irow,icol,indx_tmp)
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer ib,i,i1,i2,i3,isav,iele
      integer indx_tmp(irow)
      integer irow,icol,j
c
c     First, get indx_row for A~
      i2 = 0
      iele = indx_tmp(irow)
      do 30 ib = 1,icol
         do 20 i = 1,iele
            i1 = indx_ele(i)
            if (i1.ne.ib) goto 20
            i2 = i2+1
 20      continue
         indx_row(ib) = i2
 30   continue

c     Second, shift row and column of indx_ele
      isav = 0
      do 50 ib = 1,irow
         i1 = indx_tmp(ib)
         if (i1.le.isav) goto 50
         do 40 i = isav+1,i1
            i2 = indx_ele(i)
            i3 = (i2-1)*irow+ib
            indx_ele(i) = i3
 40      continue
         isav = i1
 50   continue
      call sort_intr(iele,indx_ele,aprm,1)
c
c     Third, get compressed indx_ele
      do 60 ib = 1,iele
         i1 = indx_ele(ib)
         if (i1.le.irow) goto 60
         i3 = (i1-1)/irow
         indx_ele(ib) = i1-i3*irow
 60   continue
c
 100  continue
      return
      end
c
