CTITLE  PLM
 
 
      real*8 function plm(n,m,z,iopt)

      implicit none 
 
c     this function returns associated Legendre functions of
c     degree n and order m. IOPT is input to indicate whether
c     the functions should be normalized (0 no; 1 yes )
 
 
*         iopt     - normalizing option
*         m        - order of function
*         mi       - order of current recursive calculation
*         n        - degree of function
*         ni       - degree of current recursive calculation
 
      integer*4 iopt, m, mi, n, ni
 
*         z        - argument (ie cos(theta) )
*         p_new    - work variable
*         pn(3)    - calculation vector
*         pm(2)    - calculation vector
 
 
      real*8 z, p_new, pn(3), pm(2)
            
c**rwk 070911: 'iopt' is always zero in current calls and not used here,
c               so add dummy statement to prevent compiler warning
      if( iopt.ne.0 ) print *,'PLN iopt ',iopt


c     Initialize recursive calculation
      pn(3) = z
      pn(2) = 1.d0
      pn(1) = 0.d0
      mi = 0
      ni = 1
 
c.... Calculate Legendre functions of order m
 
      do while (mi .lt. m)
 
         do while (ni .le. mi+2)
 
c           Increase degree by 1
            p_new = ((2*ni + 1)*z*pn(3) - (ni + mi)*pn(2))/(ni-mi+1)
            pn(1) = pn(2)
            pn(2) = pn(3)
            pn(3) = p_new
            ni = ni + 1
 
         end do
 
c        Increase order by 1
         ni = ni - 2
         pm(1) = ((ni+mi+1)*z*pn(1) - (ni-mi+1)*pn(2))/sqrt(1-z**2)
         ni = ni + 1
         pm(2) = ((ni+mi+1)*z*pn(2) - (ni-mi+1)*pn(3))/sqrt(1-z**2)
 
c        Store Pm into Pn and repeat
         pn(3) = pm(2)
         pn(2) = pm(1)
         pn(1) = 0.d0
         mi = mi + 1
 
      end do
 
c.... Advance degree to n
 
      if (ni .gt. n) then
 
         plm = pn(2)
 
      else
        do while (ni .lt. n)
           P_new = ((2*ni+1)*z*pn(3) - (ni + mi)*pn(2))/(ni-mi+1)
           pn(2) = pn(3)
           pn(3) = p_new
           ni = ni + 1
        end do
 
        plm = pn(3)
 
      end if
 
      return
      end
 
 
 
