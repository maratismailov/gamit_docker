Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      subroutine nmu11(nrat,ncat,al1,al2,iflg,g1,g2,key)
c
c     form n11 submatrix and u1 term
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer ix,iy,iz,is,it,nrat,ncat,iflg,key,ia,ir,iu,iv
     .      , ir1,i,j,k,l,m,n

      real*8  al1(maxobs),al2(maxobs),g1,g2,small,wlc,wlk,x1,x2

      logical debug/.false./

      small = 1.0d-11
c
      do 10 i = 1,nrat
         is = irowct(i+1)+1
         it = irowct(i+2)
         if(is.gt.it) go to 10
         call zero1d(1,ncat,work1)
         call zero1d(1,ncat,work2)
cd         if(debug.and.iepoch.eq.2031) 
cd     .      print *,'NMU11 i is it ',i,is,it 
         do 25 j = is,it     
cd            if(debug.and.iepoch.eq.2031) 
cd     .        print *,'j ipntct ct ',j,ipntct(j),ct(j)
            x1 = ct(j)
            if (dabs(x1).lt.small) goto 25
            x2 = x1*g2
            x1 = x1*g1
            ir = ipntct(j)
            ir1 = ir*(ir-1)/2
            do 15 k = 1,ncat
               if (ir-k.ge.0) then
                  iu = k+ir1
               else
                  iv = k*(k-1)/2
                  iu = ir+iv
               endif
               work1(k) = work1(k)+x1*cphi(iu)
               if (l2flag.ge.2) work2(k) = work2(k)+x2*cphik(iu)
   15       continue
   25    continue
         do 30 l = i,nrat
            iy = irowct(l+1)+1
            iz = irowct(l+2)
cd            if(debug.and.iepoch.eq.2031) print *,'iy iz ',iy,iz
            if(iy.gt.iz) go to 30
            ix = l*(l-1)/2 + i
            if(debug.and.ix.eq.139655) print *,'NMU11 ix l i ',ix,l,i
            wlc = 0.0d0
            wlk = 0.0d0
            do 50 m = iy,iz
               ia = ipntct(m)
               if(l2flag.ge.-1) wlc = wlc+work1(ia)*ct(m)
               if(l2flag.eq.-2) wlc = wlc+work1(ia)*ct(m)*g1
               if (l2flag.ge.2) wlk = wlk+work2(ia)*ct(m)
   50       continue
            a(ix) = a(ix)+wlc+wlk
            if (key.eq.1) alc(ix) = alc(ix)+wlc  

   30    continue
c  why ?
         if(iflg.eq.1) go to 10
         wlc = 0.0d0
         wlk = 0.0d0
         do 40 n = 1,ncat
            wlc = wlc+work1(n)*al1(n)
            if (l2flag.ge.2) wlk = wlk+work2(n)*al2(n)
   40    continue
         b(i) = b(i)+wlc+wlk
         if (key.eq.1) blc(i) = blc(i)+wlc     
         if(debug.and.ix.eq.139655) 
     .        print *,'  wlc wlk a ',wlc,b(i),blc(i)

   10 continue

      return
      end
