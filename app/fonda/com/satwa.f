      subroutine satwa(n,at,w,atxa,temp,mode)
c	
c     matrix multiplication (A~)x(W)x(A)
c     where A is rowwise symmetric storage matrix (nxn)
c     W is a full symmetric matrix
c     temp (nx1): working vector
c
      real*8 at,w,atxa,temp,temp1
      integer n,in1,in2,in3,kn,icol,irow,ii,ir,mode,in
  
      dimension w(n*(n+1)/2),at(n*(n+1)/2),temp(n)
      dimension atxa(n*(n+1)/2)
  
c     calculate (A~)x(W)
      do 30 irow = 1,n
         in1 = irow*(irow-1)/2
         do 20 icol = 1,n
            ii = icol*(icol-1)/2
            temp1 = 0.0d0
            do 10 in = 1,n
               in3 = in*(in-1)/2
               if (in.le.irow) then
                  in2 = in1+in
               else
                  in2 = in3+irow
               endif
               if (in.le.icol) then
                  kn = ii+in
               else
                  kn = in3+icol
               endif
               temp1 = temp1+at(in2)*w(kn)
 10         continue
            temp(icol) = temp1
 20      continue

c        calculate (A~)x(W)x(A)
         do 60 ir = 1,irow
            kn = ir*(ir-1)/2
            temp1 = 0.0d0
            do 40 icol = 1,n
               ii = icol*(icol-1)/2
               if (ir.ge.icol) in2 = kn+icol
               if (ir.lt.icol) in2 = ii+ir
               temp1 = temp1+temp(icol)*at(in2)
 40         continue
            atxa(in1+ir) = temp1
 60      continue

 30   continue

      return
      end
