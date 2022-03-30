      subroutine isort2(n,id1,id2,id3)
c	
c     sort out integer arrays from smallest to largest.
c     id1, id2 : reference arrays
c     id3      : index arrays
c
      integer n, id1, id2, id3, iloop, jloop
      integer in1, in2, jn1, jn2, itemp
      dimension id1(n), id2(n), id3(n)
  
      do 30 iloop = 1,n-1
         in1 = id1(iloop)
         in2 = id2(iloop)
         do 20 jloop = iloop+1,n
            jn1 = id1(jloop)
            jn2 = id2(jloop)
            if ((in1.gt.jn1).or.(in1.eq.jn1.and.in2.gt.jn2)) then
               id1(iloop) = jn1
               id1(jloop) = in1
               id2(iloop) = jn2
               id2(jloop) = in2
               itemp = id3(iloop)
               id3(iloop) = id3(jloop)
               id3(jloop) = itemp
               in1 = jn1
               in2 = jn2
            endif
 20      continue
 30   continue

      return
      end
