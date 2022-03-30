CTITLE TRAN_TO_ORG
 
      subroutine tran_to_org( tran, org, iel, step )
 
      implicit none 

 
*     Routine to copy tran to org and add unit vector to values
 
*   i       - Loop counter
*   iel     - Element for unit vector
 
      integer*4 i, iel
 
*   org(6)  - Destination vector
*   tran(6) - Translation to be copied
*   step    - Length of unit vector
 
 
      real*8 org(6), tran(6), step
 
****  Copy tran to org
      do i = 1, 6
          org(i) = tran(i)
      end do
 
*     Add unit vector (if neccessary)
      if( iel.gt.0 ) then
          org(iel) = org(iel) + step
      end if
 
***** Thats all
      return
      end
 
