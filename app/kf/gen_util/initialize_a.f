CTITLE INITIALIZE_ATM

      subroutine initialize_atm( lapse, ht, tbias, adry,
     .           awet, bdry, bwet )

      implicit none 

*     Routine to initialize the default values of the atmospheric
*     parameters for lapse rate, height of trop, temp bias, and
*     the mapping functions.

* VARIABLES
*     lapse  - lapse rate K/km
*     ht     - Height of tropopuase (km)
*     tbias  - bias in surface temperature (K)

*     adry(4) - 4 coefficients in dry mapping function
*     awet(4) - 4 coefficinets in continued fraction for wet
*     bdry(4) - 4 coefficients for continues fraction for dry bending
*     bwet(4) - 4 coefficients for continued fraction for wet bending

      real*4 lapse, ht, tbias, adry(4), awet(4), bdry(4), bwet(4)

*     j      - Loop counter

      integer*4 j

****  Initialize the new variables
      lapse = -5.6
      ht = 10.0
      tbias = 0.0

      adry(1) = 0.001185d0
      adry(2) = 0.001144d0
      adry(3) = -0.0090d0
      adry(4) = 0.d0

      awet(1) = 0.00035d0
      awet(2) = 0.017d0
      awet(3) = 0.d0
      awet(4) = 0.d0

      do j = 1,4 
         bdry(j) = 0.d0
         bwet(j) = 0.d0
      end do
 
****  Thats all
      return
      end 

