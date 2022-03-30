      integer function lift_arg(string,arg,id_arg)
c
c     count number of characters in the 'id_arg'th  argument of a string
c
      integer ls,i,il,id_arg,la,nblen
      character*(*) string,arg
      character*1  blk,chr
c
      lift_arg = 0
      ls = nblen(string)
      if (ls .eq. 0) return
      call blank(arg)
c
c     visit the string and find the id_arg th argument
      il = 0
      la = 0
      blk = ' '
      do i = 1,ls
         chr = string(i:i)
         if (blk.eq.' '.and.chr.ne.' ') il = il+1
         if (il.eq.id_arg) then
            if (chr.eq.' ') goto 10
            la = la+1
            arg(la:la) = chr
         endif
         blk = chr
      enddo

 10   lift_arg = la

      return
      end
c
c
