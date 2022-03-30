Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      Subroutine DBAR( w )
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer ir,iv,is,it,ii,iu,i,j,k

      real*8 w(maxobs*maxdbl),w1,x
c
C     Matrix multiplication-- full matrix times a sparse matrix.
C     P - Symmetric storage matrix. (NRD x NRD)
C     D - Sparce matrix. (NRD x NCD) (NCD > NRD) Located by column.
C     IROWDT - Number of live elements by columns.
C        IROWDT(1) : Number of
C        IROWDT(2) = 0
C     IPNTDT - Index of live elements of D by column.
C     NRD - NUMBER OF ROWS IN D
C     NCD - NUMBER OF COLUMNS IN D
C     W = P x D.  (NRD x NCD)  Ordered by column.
C        Unfortunately, W is neither symmetric nor sparce.
C     ----- DND ----- 871229


c      print *,'DBAR irowdt ',irowdt

c---- Calculate W = P x D
      do 100 ir = 1,nrd
         iv = (ir*(ir-1))/2
         do 10 i = 1,ncd
            is = irowdt(i+1)+1
            it = irowdt(i+2)
            if(is.gt.it) go to 10
            ii = (i-1)*nrd+ir
            w1 = 0.0d0
            do 15 j = is,it
               k = ipntdt(j)
               x = dt(j)
               if(ir-k) 16,17,17
   16          iu = ir+(k*(k-1))/2
               go to 18
   17          iu = k+iv
   18          w1 = w1+x*dpdt(iu)
   15       continue
            w(ii) = w1
   10    continue
  100 continue
c
      return
      end


