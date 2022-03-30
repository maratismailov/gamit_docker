CTITLE DECODE_MAR
  
      subroutine decode_mar(element, num, code) 

      implicit none 
c 
c     This routine will decode the markov element code which is stored
c     in packed form.  The site/source number is stored in the right- 
c     hand 8 bits, and the code is in the lefthand 8 bits.  (See &plkbd 
c     for listing of meaning of the codes)
c 
c Variables 
c --------- 
c element -- the markov element in packed form
c num -- the site/source number 
c code -- the markov element code number
c 
      integer*4 element, num, code
c 
c.... Depack the code, get the site/source number 
      num = element/256 
c 
c.... Get the code
      code = element - num*256
c 
c.... Thats all 
      return
      end 
  
