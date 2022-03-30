      subroutine blank(string)
c
c     blank a character string
c                       
      integer ls,i
      integer nblen
      character*(*) string, blk*1
      data blk/' '/
c
      ls = nblen(string)
      if (ls .eq. 0) return
      do i = 1,ls
         string(i:i) = blk
      enddo
  
      return
      end
c
c
