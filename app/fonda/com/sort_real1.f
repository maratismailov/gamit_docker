c
      subroutine sort_real1(ndat,array,indx,type)
c
c     sort a real array with its index
c     type = 1: from the smallest to the bigest
c     type = 2: from the bigest to the smallest
c
      integer type,i,k,iswtch,ndat,indx
      real*8  temp,array(ndat)
      dimension indx(ndat)
c
      if(ndat.le.1) go to 100
      do 50 k=1,32000
         iswtch=0  
         do 30  i=2,ndat
            if(type.eq.1.and.array(i).ge.array(i-1)) go to 30
            if(type.eq.2.and.array(i).lt.array(i-1)) go to 30
            temp = array(i)
            array(i) = array(i-1)
            array(i-1) = temp
            iswtch = indx(i)
            indx(i) = indx(i-1)
            indx(i-1) = iswtch
 30      continue
         if(iswtch.eq.0) go to 100
 50   continue
c
 100  continue  
c
      return
      end
