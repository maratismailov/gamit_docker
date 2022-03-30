Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine sort4i(ndat,iar1,iar2,iar3,iar4,key)
c
c     sort 4 integer arrays from the smallest to the bigest
c     key : control variable.

      implicit none

      integer ndat,iar1,iar2,iar3,iar4,key,iswtch,i,k     
      dimension iar1(ndat),iar2(ndat),iar3(ndat),iar4(ndat)



      if(ndat.le.1) go to 100
      do 50 k=1,32000
         iswtch=0
         do 60  i=2,ndat
         go to (10,20,30,40,100), key
 10         if(iar1(i).ge.iar1(i-1)) go to 60
               goto 80
 20         if(iar2(i).ge.iar2(i-1)) go to 60
               goto 80
 30         if(iar3(i).ge.iar3(i-1)) go to 60
               goto 80
 40         if(iar4(i).ge.iar4(i-1)) go to 60
 80         iswtch=iar1(i)
            iar1(i)=iar1(i-1)
            iar1(i-1)=iswtch
c
            iswtch=iar2(i)
            iar2(i)=iar2(i-1)
            iar2(i-1)=iswtch
c
            iswtch=iar3(i)
            iar3(i)=iar3(i-1)
            iar3(i-1)=iswtch
c
            iswtch=iar4(i)
            iar4(i)=iar4(i-1)
            iar4(i-1)=iswtch
c
60       continue
         if(iswtch.eq.0) go to 100
50    continue
c
100   continue
c
      return
      end
