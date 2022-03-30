      integer function count_arg(string)
c
c     count number of arguments in a string
c
      integer ls,i,il,nblen
      character*(*) string
      character*1 blk,chr
c
      count_arg = 0
      ls = nblen(string)
      if (ls .eq. 0) return
c
c     blank is the separation pointer
      il = 0
      blk = ' '
      do i = 1,ls
         chr = string(i:i)
         if (blk.eq.' '.and.chr.ne.' ') il = il+1
         blk = chr
      enddo

      count_arg = il

      return
      end
c
c
