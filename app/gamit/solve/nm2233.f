c
      subroutine nm2233(iphilv,mpntb,mbcnt,g1,g2,isame,key)
c
c compute n22 and n33 portions of normal matrix
c  n22 = ((1-g)**2)w + ((1+g)**2)we) : l1 biases
c  n33 = (g*g*w + we)                : l2-l1 biases
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer ielem,ielema,irowa,iphilv(maxobs),mpntb(maxobs)
     .      , mbcnt,isame,key,icntr,iel20,iel22,iroww,icntc
     .      , ielemw,i,j

      real*8 g1,g2,a1,b1

      logical debug
        
c      if(iepoch.eq.2471) then
c        debug = .true.
c      else
        debug = .false.
c      endif 

c mpntb : bias parameter pointers for normal matrix
c
      icntr = 0       
      iel20 = 0
      iel22 = 0
      do 10 i = 1,mbcnt                         
        if(debug) print *,'NM2233 i iphilv ',i,iphilv(i)
        if(iphilv(i).eq.0) go to 10
        icntr = icntr+1
        iroww = icntr*(icntr-1)/2
        ielem = mpntb(i)
        irowa = ielem*(ielem-1)/2
        if (key.eq.1) iel20 = i*(i-1)/2
        icntc = 0
        b1 = 0.0d0
          do 11 j = 1,i   
            if(debug) print *,'NM2233 j iphilv ',iphilv(j) 
            if(iphilv(j).eq.0) go to 11
            icntc = icntc+1
            ielemw = iroww+icntc
            ielema = irowa+mpntb(j)
            a1 =  g1*cphi(ielemw)*dble(isame)
            if(l2flag.ge.2) b1 = g2*cphik(ielemw)*dble(isame)
            a(ielema) = a(ielema)+a1+b1
            if (key.eq.1) iel22 = iel20+j
            if (key.eq.1) an22(iel22) = an22(iel22)+a1              
            if(debug) then 
             print *,'NM2233 i ielema,ielemw,cphi(ielemw),cphik(ielemw)'
     .         ,             i,ielema,ielemw,cphi(ielemw),cphik(ielemw)
              print *,'     iel22 an22 ',iel22,an22(iel22) 
             endif 
   11     continue
   10 continue
      return
      end
