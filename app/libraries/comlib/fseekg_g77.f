      subroutine fseekg(iunit,ioffset,ipos,iostat)
c
c     Wrapper routine to call fseek as an intrinsic subroutine (G77)
c     S. McClusky
c     10th January 1997

      implicit none
      
      integer*4 iunit,ioffset,ipos,iostat

      iostat = 0
      call fseek(iunit,ioffset,ipos)

****  Thats all
      return
      end

