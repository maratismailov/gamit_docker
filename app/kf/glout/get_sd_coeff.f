CTITLE get_sd_coeff
 
      subroutine get_sd_coeff( in, bs, apr_val, apr_arg, apr_num,
     .                         coeffs )

      implicit none 
 
*     This routine will match the in'th coefficient estimate with
*     the apriori value used for the coefficient.

* in -- NUmber of argumnet to match.  WARNING: This routine must be
*       keep in line with nut_coeff_pa.f (so that index and frequency go
*       together) 
* apr_arg(6,apr_num) - Five arguments for the terms (l,l',D,F, Om).  These
*       may be opposit sign for the retrograde part of the bs
* bs(2)    - Tidal bs (-1,-2,1,2 for retrograde diurnal, ret semi,
*              prograde diurnal, pro semi. and the argument sign to be
*              used i.e., are the fundamental arguments reversed.
* apr_num    - Number of values in the apriori list.
* apr_val(2,apr_num) - Apriori values for the coefficients.
* coeffs(2) -- Coefficients for this series.  Returned as sine and 
*              cosine values.

      integer*4 in, apr_num, apr_arg(6,apr_num),bs(2)
      real*8    apr_val(2,apr_num), coeffs(2)
 
*   i,j         - Loop counters
*   x1(6,15)    - Coefficent mulipliers for Nutation arguments
*   found       - TRue if arguments matched.
 
      integer*4 i, x1(6,15)
c      logical  found
 
C
C      CONTSTANTS ARE BASED ON VALUES GIVEN IN NUTW IN CALC V5.0
C      pp182-183, with addition of space for FCN.
C
C
C                 MULTIPLE OF
C                L    L'   F    D  OMEGA FCN
      DATA X1 /  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .           0 ,  0 ,  0 ,  0 ,  0 ,  1 ,
     .           0 ,  0 ,  0 ,  0 ,  1 ,  0 ,
     .           0 ,  0 ,  2 , -2 ,  2 ,  0 ,
     .           0 ,  0 ,  2 ,  0 ,  2 ,  0 ,
     .           0 ,  0 ,  0 ,  0 ,  2 ,  0 ,
     .           0 ,  1 ,  0 ,  0 ,  0 ,  0 ,
     .           1 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .           0 ,  1 ,  2 , -2 ,  2 ,  0 ,
     .           0 ,  0 ,  2 ,  0 ,  1 ,  0 ,
     .           1 ,  0 ,  2 ,  0 ,  2 ,  0 ,
     .           0 , -1 ,  2 , -2 ,  2 ,  0 ,
     .           1 ,  0 ,  0 , -2 ,  0 ,  0 ,
     .           0 ,  0 ,  2 , -2 ,  1 ,  0 ,
     .          -1 ,  0 ,  2 ,  0 ,  2 ,  0  /
 
C
*     Clear the coeffs first (so zero if not found)
      do i = 1,2
         coeffs(i) = 0.d0
      end do

*     Loop over the apr_arg values trying to find match
      do i = 1, apr_num

*        See if direct match
         if( bs(2)*X1(1,in).eq.apr_arg(1,i) .and.
     .       bs(2)*X1(2,in).eq.apr_arg(2,i) .and.
     .       bs(2)*X1(3,in).eq.apr_arg(3,i) .and.
     .       bs(2)*X1(4,in).eq.apr_arg(4,i) .and.
     .       bs(2)*X1(5,in).eq.apr_arg(5,i) .and.
     .       bs(1)         .eq.apr_arg(6,i)       ) then
             coeffs(2) = apr_val(1,i)
             coeffs(1) = apr_val(2,i)
         end if

      end do

****  THATS ALL
      return
      end
 
