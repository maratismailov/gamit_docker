      integer function strjst(string)
c
c     remove blank character at the begining
c                       
      integer ls,i,it,il,nblen
      character*(*) string, blk*1,chr*1
      data blk/' '/
c
      ls = nblen(string)
      if (ls .eq. 0) return
      do i = 1,ls
         chr = string(i:i)
         if (chr.ne.' ') then
            it = i
            goto 10
         endif
      enddo
 
c     remove blank character at the begining
 10   il = 0
      do i = it,ls
         chr = string(i:i)
         il = il+1
         string(il:il) = chr
         string(i:i) = blk
      enddo

      strjst = il 
      return
      end
c
c
