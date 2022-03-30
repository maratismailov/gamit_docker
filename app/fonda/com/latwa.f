      subroutine latwa(l,n,at,w,cl,temp,mode,job)
c	
c     matrix multiplication (A~)x(W)x(A)
c     where A is rowwise storage matrix (lxn for A~, nxl for A)
c        cl is a (lxl) symmetric matrix
c     mode = 1: W is a full symmetric matrix
c     mode = 0: W is a diagonal matrix
c     temp (nx1): working vector
c     job = 1:   get Cl - (A~)x(W)x(A)
c     job = 2:   get Cl + (A~)x(W)x(A)
c
      implicit real*8(a-h,o-z)
  
      integer l,n,mode,job,irow,in1,icol,ii,in,in2,kn,in3,ir,ir1
      integer ic,ic1
  
      dimension at(l*n),cl(l*(l+1)/2),temp(n)
      dimension w(n*(n*mode+1)/(mode+1))
  
      do 30 irow = 1,l
         in1 = (irow-1)*n
         do 20 icol = 1,n
            if (mode.eq.0) then
               temp(icol) = w(icol)*at(in1+icol)
               goto 20
            endif
            if (mode.eq.1) then
               ii = icol*(icol-1)/2
               temp1 = 0.0d0
               do 10 in = 1,n
                  in2 = in1+in
                  if (in.le.icol) then
                     kn = ii+in
                  else
                     in3 = in*(in-1)/2
                     kn = in3+icol
                  endif
                  temp1 = temp1+w(kn)*at(in2)
 10            continue
            endif
            temp(icol) = temp1
 20   continue

         do 40 ir = irow,l
            ir1 = (ir-1)*n
            temp1 = 0.0d0
            do 50 ic = 1,n
               ic1 = ir1+ic
               temp1 = temp1+temp(ic)*at(ic1)
 50         continue
            ir1 = ir*(ir-1)/2+irow
            if(job.eq.1) cl(ir1) = cl(ir1)-temp1
            if(job.eq.2) cl(ir1) = cl(ir1)+temp1
 40      continue

 30   continue

      return
      end
