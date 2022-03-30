c------------------------------------------------------
      integer function remedy_space(strinx,striny,ix,ic,mode)
c
c     replace space in strinx with character ic and
c     produce no-space string striny
c     mode = 1: remedy all spaces
c     mode = 2: remedy space inside strinx
      character*(*) strinx,striny
      character ic
      integer ix,i,mode,len,nblen
c
      len = nblen(strinx)
      do 10 i = 1,ix
         if (strinx(i:i).eq.' ') then
            striny(i:i) = ic
            if(mode.eq.2.and.i.gt.len) striny(i:i) = ' '
         else
            striny(i:i) = strinx(i:i)
         endif
 10   continue
c
      remedy_space = len
      return
      end
