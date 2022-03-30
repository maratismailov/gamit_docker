Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine sort1i(ndat,iarray)
c
c     sort an integer array from the smallest to the bigest

      implicit none
                     
     
      integer ndat,iarray,iswtch,i,k
      dimension iarray(ndat)
c
      if(ndat.le.1) go to 100
      do 50 k=1,32000
         iswtch=0
         do 30  i=2,ndat
            if(iarray(i).ge.iarray(i-1)) go to  30
            iswtch=iarray(i)
            iarray(i)=iarray(i-1)
            iarray(i-1)=iswtch
 30      continue
         if(iswtch.eq.0) go to 100
 50   continue
c
 100  continue
c
      return
      end
