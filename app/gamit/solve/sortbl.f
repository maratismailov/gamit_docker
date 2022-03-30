      Subroutine SORTBL(nlen,ist,baslen,ndex)
c
c     Sort out all baselines by lenghths.
c
c     ist   : shift index
c     nlen  : baseline number
c
      implicit none
         
      integer nlen,ist,ndex,jj,i1,i2,i,j
     
      real*8 baslen(nlen),temp

      dimension ndex(nlen*ist)
                               
      logical debug/.false./

      if(debug) then 
        print *,'SORTBL nlen ist baslen ndex ' ,ist,nlen
        write(*,'(10i7)') (ndex(i),i=1,10)
        write(*,'(10f7.0)') (baslen(i),i=1,10)
      endif   
      

      do 30 j = 1,32000
         jj = 0
         do 20 i = 2,nlen
            if(baslen(i).ge.baslen(i-1)) go to 20
            temp = baslen(i)
            baslen(i) = baslen(i-1)
            baslen(i-1) = temp
c
            i2 = (i-1)*ist
            do 10 i1 = i2+1,i2+ist
               jj = ndex(i1)
               ndex(i1) = ndex(i1-ist)
               ndex(i1-ist) = jj
 10         continue

 20      continue
         if(jj.eq.0) go to 40
 30   continue

 40   continue

      return
      end
