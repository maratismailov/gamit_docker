      subroutine isort(n,id1,id2,id3,mode)
c	
c     sort out integer arrays from smallest to largest.
c     mode : specify which array is the reference array.
c
      integer n, id1, id2, id3, mode, iloop, jloop
      integer in,jn,itemp
      dimension id1(n), id2(n), id3(n)
  
      do 30 iloop = 1,n-1
         if (mode.eq.1) in = id1(iloop)
         if (mode.eq.2) in = id2(iloop)
         if (mode.eq.3) in = id3(iloop)
         do 20 jloop = iloop+1,n
            if (mode.eq.1) jn = id1(jloop)
            if (mode.eq.2) jn = id2(jloop)
            if (mode.eq.3) jn = id3(jloop)
            if (in.gt.jn) then
               in = jn
               itemp = id1(iloop)
               id1(iloop) = id1(jloop)
               id1(jloop) = itemp
               itemp = id2(iloop)
               id2(iloop) = id2(jloop)
               id2(jloop) = itemp
               itemp = id3(iloop)
               id3(iloop) = id3(jloop)
               id3(jloop) = itemp
            endif
 20      continue
 30   continue

      return
      end
