CTITLE 'LEQ'
 
 
      logical function leq(a,b)
 
      implicit none

*     J.L. Davis                   4:48 PM  FRI., 17  APR., 1987
*
*     Function to test for the lexical equality of two strings using
*     the FORTRAN 77 intrinsic functions lge and lle
 
 
*             a, b                - The nput strings
 
      character*(*) a, b
 
      leq = lle(a,b) .and. lge(a,b)
 
      end
 
 
