      subroutine axb(l,m,n,avec,bvec,cvec,mode1,mode2)
c	
c     matrix multiplication C(l,n) = A(l,m) x B(m,n)
c     mode1 = 1: both A and B are asymmetric
c     mode1 = 2: A is symmetric
c     mode1 = 3: B is symmetric
c
c     mode2 = 0: C = A x B
c     mode2 = 1: A = A(m,l),  C = A~ x B
c     mode2 = 2: B = B(n,m),  C = A x B~
c
      implicit real*8(a-h,o-z)
 
      integer l,m,n,mode1,mode2,il,in,ii,im,i1,i2
  
      dimension avec(l*m),bvec(m*n),cvec(l*n)
  
      do 30 il = 1,l
         do 20 in = 1,n
            ii = (il-1)*n+in
            temp = 0.0d0
            do 10 im = 1,m
               i1 = (il-1)*m+im
               if (mode2.eq.1) i1 = (im-1)*l+il
               i2 = (im-1)*n+in
               if (mode2.eq.2) i2 = (in-1)*m+im
               if (mode1.eq.2.and.il.ge.im) i1 = il*(il-1)/2+im 
               if (mode1.eq.2.and.il.lt.im) i1 = im*(im-1)/2+il 
               if (mode1.eq.3.and.in.ge.im) i2 = in*(in-1)/2+im 
               if (mode1.eq.3.and.in.lt.im) i2 = im*(im-1)/2+in 
               temp = temp+avec(i1)*bvec(i2)
 10         continue
            cvec(ii) = temp
 20      continue
 30   continue

      return
      end
