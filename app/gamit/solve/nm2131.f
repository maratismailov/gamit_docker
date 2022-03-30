Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      subroutine nm2131(nrat,ncat,iphilv,mpntb,mbcnt,g1,g2,key)
c
c     form n21, n31 submatrix

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer ielema,iel21,nrat,ncat,mbcnt,key,i1,ir,is,it,iv,iu
     .      , icntc,irow,il,i,j,l,m,iphilv(maxobs),mpntb(maxobs)

      real*8 g1,g2,small,x1,x2,x3
                         
      small = 1.0d-11
c
      i1 = -1
      do 10 i = 1,nrat
         is = irowct(i+1)+1
         it = irowct(i+2)
         i1 = i1+i
         if(is.gt.it) go to 10
         call zero1d(1,ncat,work)
         call zero1d(1,ncat,work1)
         do 15 j = is,it
            x1 = dr(j)
            if (dabs(x1).lt.small) goto 15
            x2 = x1*g2
            x1 = x1*g1
            ir = ipntct(j)
            iv = (ir*(ir-1))/2
            do 12 l = 1,ncat
               if (ir-l.ge.0) then
                  iu = l+iv
               else
                  il = (l*l-l)/2
                  iu = ir+il
               endif
               x3 = x1*cphi(iu)
               work(l) = work(l)+x3
               if (key.eq.1) work1(l) = work1(l)+x3
               if (l2flag.ge.2) work(l) = work(l)+x2*cphik(iu)
   12       continue
   15    continue
c
c work array is one row contribution to the normal matrix
c  enter into normal matrix in corresponding column (i.e., transpose)
c
c loop over columns of normal matrix
         icntc = 0
         do 60 m = 1,mbcnt
            if(iphilv(m).eq.0) go to 60
            icntc = icntc+1
            irow = mpntb(m)
            if(irow-i.ge.0) then
               ielema = i+irow*(irow-1)/2
            else
               ielema = irow+i1
            endif
            a(ielema) = a(ielema)+work(icntc)
            if (key.eq.1) then
               iel21 = i+(m-1)*nrat
               clc(iel21) = clc(iel21)+work1(icntc)
            endif
   60    continue
c
   10 continue
      return
      end
