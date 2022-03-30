CTITLE get_nut_coeff
 
      subroutine get_nut_coeff( in, coeffs )

      implicit none 
 
*     This rotuine will check the args passed and find the nutation
*     series which match.  ONLY WORKS WHEN in_nut series is 
*     is used. 

      include '../includes/nut_eval.h'

* in -- NUmber of argumnet to match.  WARNING: This routine must be
*       keep in line with nut_coeff_pa.f (so that index and frequency go
*       together) 
* coeffs(4) -- Coefficients for this series. Long in and out
*              oblquity in and out)

      integer*4 in
      real*8    coeffs(4)
 
*   i,j         - Loop counters
*   x1(6,15)    - Coefficent mulipliers for Nutation arguments
*   found       - TRue if arguments matched.
 
      integer*4 i, j, x1(6,15)
      logical  found
 
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
      do i = 1,4
         coeffs(i) = 0.d0
      end do

*     Loop from 3 value since 1 is precession and 2 is FCN
      i = 0
      do while (i.lt. num_series) 
         i = i + 1
         found = .true.
         do j = 1,5
            if( X1(j,in).ne. arg_series(j,i)) found = .false.
         end do

*        if found save value
         if( found ) then
             do j = 1,4
                coeffs(j) = val_nut_coeff(j,i)
             end do

*            Set i to get out of loop
             i=num_series+1
          end if
C
      end do
 
****  THATS ALL
      return
      end
 
