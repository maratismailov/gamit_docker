Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      subroutine nu23(nrw,iphilv,mpntb,mbcnt,al1,al2,g1,g2,key)

c     form u2, u3 portion of u vector

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer iphilv(maxobs),mpntb(maxobs),nrw,mbcnt,key,ir,iu
     .      , icntc,irow,i,j,m

      real*8 al1(maxobs),al2(maxobs),g1,g2,wlc,wlk,u1

      do 10 i = 1,nrw
         ir = (i*i-i)/2
         wlc = 0.0d0
         wlk = 0.0d0
         do 15 j = 1,nrw
            if(i-j) 16,17,17
   16       iu = i+(j*(j-1))/2
            go to 18
   17       iu = j+ir
   18       wlc = wlc+cphi(iu)*al1(j)
            if (l2flag.ge.2) wlk = wlk+cphik(iu)*al2(j)
   15    continue
         work1(i) = g1*wlc
         work2(i) = g2*wlk
   10 continue
c
c loop over columns of u-matrix
          icntc = 0
          do 70 m = 1,mbcnt
            if(iphilv(m).eq.0) go to 70
            icntc = icntc+1
            irow = mpntb(m)
            u1 =  work1(icntc)
            b(irow) = b(irow)+u1
            if (l2flag.ge.2) b(irow) = b(irow)+work2(icntc)
            if (key.eq.1) blc(irow) = blc(irow)+u1
   70     continue
c
      return
      end
