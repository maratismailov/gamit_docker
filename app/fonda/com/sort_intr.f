c
      subroutine sort_intr(ndat,iarray,array,type)
c
c     sort an integer array with its corresponding real array
c     type = 1: from the smallest to the bigest
c     type = 2: from the bigest to the smallest
c
      integer type,iarray,i,k,iswtch,ndat
      dimension iarray(ndat)
      real*8 temp,array(ndat)
c
      if(ndat.le.1) go to 100
      do 50 k=1,32000
         iswtch=0  
         do 30  i=2,ndat
            if(type.eq.1.and.iarray(i).ge.iarray(i-1)) go to 30
            if(type.eq.2.and.iarray(i).lt.iarray(i-1)) go to 30
            iswtch=iarray(i)
            iarray(i)=iarray(i-1)
            iarray(i-1)=iswtch
            temp = array(i)
            array(i) = array(i-1)
            array(i-1) = temp
 30      continue
         if(iswtch.eq.0) go to 100
 50   continue
c
 100  continue  
c
      return
      end
