
*******************************************************************
*     Parameter values for the nutation seroes evalaution
*******************************************************************
 
*   max_series           - Maximum number in series series
*   max_plan            - Maximum number of planetary terms allowed
 
      integer*4  max_series, max_plan
 
      parameter ( max_series = 400  )
      parameter ( max_plan   =   50 )

*   input_ZMOA1990  - Name of the file containing the ZMAO-1990 nutation
*                     series
*   input_plan      - NAme of the planetary nutation series file.

      character*(*) input_ZMOA1990, input_plan

C     parameter ( input_ZMOA1990 = 'ZMOA1990.00' )
      parameter ( input_ZMOA1990 = 'ZMOA1990'    )
      parameter ( input_plan     = 'ZMOA_PLAN'   )
 
 
*****************************************************************
*     Common block for nutation series evaluation
*****************************************************************
 
*       arg_series(5,max_series)  - Arguments for the
*                           - series series
*       arg_plan(5,max_plan)    - Arguments for planetary nutation
*                           - terms
*       num_series           - number of terms in series series
*       num_plan            - number of planetary nutation terms
 
      integer*4  arg_series(5,max_series), arg_plan(5,max_plan),
     .    num_series, num_plan
 
*       val_nut_coeff(4,max_series)      - Coefficients of the
*                           - apriori real Earth series (long in &
*                           - out, and obl in&out [mas]
*       val_nut_rate(4,max_series)       - Coefficients of the
*                           - apriori rate of real Earth series
*                           - (long in & out, and obl in&out [mas]
 
*       val_plan(4,max_plan)    - Planetary nutations convolved with
*                           - response of the Earth (mas)
 
      real*8  val_nut_coeff(4,max_series), val_nut_rate(4,max_series), 
     .    val_plan(4,max_plan) 
 
*   series_read              - Indicates series series has been read
*   plan_read               - Indicates planetary terms read
 
      logical  series_read, plan_read
 
*--------------------------------------------------------------------
*     Common declaration
*--------------------------------------------------------------------
 
      common / nuteval_common / arg_series, arg_plan, num_series, 
     .    num_plan, val_nut_coeff, val_nut_rate, val_plan,
     .    series_read, plan_read
 
*--------------------------------------------------------------------
