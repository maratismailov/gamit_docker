CTITLE 'JEL'
 
 
      integer*4 function jel(k,l)

      implicit none 
 
*     Function to compute element number of element k,j in a
*     lower triangle symetric matrix
 
* MOD TAH 850214 multiplied integer with Integer*4
 
 
 
      integer*4 k, l
 
 
      integer*4 temp
 
      data temp / 1 /
 
      if (k .gt. l) then
 
         jel = temp * k * (k-1) / 2 + l
 
      else
 
         jel = temp * l * (l-1) / 2 + k
 
      end if
 
      end
 
