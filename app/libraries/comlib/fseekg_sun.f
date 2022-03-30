      subroutine fseekg(iunit,ioffset,ipos,iostat)
c
c     Wrapper routine to call fseek as a library function (SUN HP DEC etc.)
c     S. McClusky
c     10th January 1997

      implicit none
      
      integer*4 iunit,ioffset,ipos,iostat,fseek

      iostat = fseek(iunit,ioffset,ipos)

****  Thats all
      return
      end

