      subroutine laxb_indx(l,n,m,ibnd,a,b,c,indx_row,indx_ele,type)
c
c     We only consider asymmetric matrix now.
c     And only consider array A uses index.
c     type = 1:  get C + A * B
c     type = 2:  get C - A * B
c
      integer l,n,m,icol,in1,in2,ix,i2,in,ibnd
      integer indx_row(l),indx_ele(l*ibnd),ib,type,ia
      real*8 a(l*ibnd),c(l*m),temp
      real*8 b(n*m)

      in1 = 0
      do 80 ia = 1,l
         in2 = indx_row(ia)
         if (in2.le.in1) goto 80
         do 70 ib = 1,m
            temp = 0.0d0
            i2 = (ia-1)*m+ib
            do 60 in = in1+1,in2
               icol = indx_ele(in)
               ix = (icol-1)*m+ib
               temp = temp+a(in)*b(ix)
 60         continue
            if (type.eq.1) c(i2) = c(i2)+temp
            if (type.eq.2) c(i2) = c(i2)-temp
 70      continue
         in1 = in2
 80   continue

      return
      end

