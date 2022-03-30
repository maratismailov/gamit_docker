c
      subroutine nm32(iphilv,mbcnt,g1,g2,isame)
c
c  compute N32 portion of normal matrix
c  N32 = (-g(1-g)w + (1+g)we)
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h' 

      integer irowa,ielema,ielemw,iphilv(maxobs),mbcnt,isame
     .      , icntr,ir1,icntc,i,j

      real*8 g1,g2,a1

      icntr = 0
      do 10 i = 1,mbcnt
        if(iphilv(i).eq.0) go to 10
        icntr = icntr+1
        irowa = ipntb2(i)
        irowa = irowa*(irowa-1)/2
        ir1 = icntr*(icntr-1)/2
        icntc = 0
          do 11 j = 1,mbcnt
            if(iphilv(j).eq.0) go to 11
            icntc = icntc+1
            if(i-j) 20,30,30
   20       ielemw = icntr+(icntc*(icntc-1))/2
            go to 40
   30       ielemw = icntc+ir1
   40       ielema = irowa+ipntb1(j)
            a1 =  g1*cphi(ielemw)
            if (l2flag.ge.2) a1 = a1+g2*cphik(ielemw)
            a(ielema) = a(ielema)+a1*dble(isame)
   11     continue
   10 continue
c
      return
      end
