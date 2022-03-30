CTITLE PLMDZ
 
 
      real*8 function plmdz(n,m,z,iopt)

      implicit none 
 
c     this function computes the derivative of the Legendre function
c     with respect to the argument z
 
 
*           iopt      - normalizing option (0 no; 1 yes)
*           m         - order of function
*           n         - degree of function
 
      integer*4 iopt, m, n
 
*           plm       - associated Legendre function (function)
*           sgn       - signum function
*           z         - argument of Legendre function
 
 
 
      real*8 plm, z
 
      plmdz=( (n + 1) * z * plm(n,m,z,iopt)
     .       -(n-m+1) * plm(n+1,m,z,iopt) ) / (1-z**2)
 
      return
      end
 
 
 
 
