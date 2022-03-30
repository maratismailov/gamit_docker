      subroutine latwa_indx(l,n,iband,at,w,cl,temp,
     .                      indx_row,indx_ele,job)
c	
c     matrix multiplication (A~)x(W)x(A)
c     where A is rowwise storage matrix (lxn for A~, nxl for A)
c        cl is a (lxl) symmetric matrix
c     When A has a lot of elements with zero, we can compress A
c        with 2 index arrays.
c        indx_row denotes cumulated number of non-zero elements
c                 for every row of A~
c        indx_ele denotes the index of correnponding parameter 
c                 for every non-zero element
c     We only deal with full symmetric matrix W.
c     temp (nx1): working vector
c     job = 1:   get Cl - (A~)x(W)x(A)
c     job = 2:   get Cl + (A~)x(W)x(A)
c
      implicit real*8(a-h,o-z)
  
      integer l,n,iband,job,irow,icol,in,in1,in2,in0,in3,ir1
      integer itmp,ii,kn,ic
      dimension at(l*iband),cl(l*(l+1)/2),temp(n)
      dimension w(n*(n+1)/2)
      integer indx_row(l),indx_ele(l*iband)
  
      in0 = 0
      do 30 icol = 1,l
c
         in1 = indx_row(icol)
         if (in1.le.in0) goto 30
         do 20 itmp = 1,n
            ii = itmp*(itmp-1)/2
            temp1 = 0.0d0
            do 10 in = in0+1,in1
               in2 = indx_ele(in)
               if (in2.le.itmp) then
                  kn = ii+in2
               else
                  in3 = in2*(in2-1)/2
                  kn = in3+itmp
               endif
               temp1 = temp1+w(kn)*at(in)
 10         continue
            temp(itmp) = temp1
 20      continue

         in3 = in0
         do 40 irow = icol,l
            ir1 = indx_row(irow)
            if (ir1.le.in3) goto 40
            temp1 = 0.0d0
            do 50 ic = in3+1,ir1
               in2 = indx_ele(ic)
               temp1 = temp1+temp(in2)*at(ic)
 50         continue
            in3 = ir1
            ir1 = irow*(irow-1)/2+icol
            if(job.eq.1) cl(ir1) = cl(ir1)-temp1
            if(job.eq.2) cl(ir1) = cl(ir1)+temp1
 40      continue

         in0 = in1
 30   continue

      return
      end
