
CTITLE 'COMPUTE_STD'
 
      subroutine compute_std( jd, c)

      implicit none 
 
 
*     Routine to compute the theorectical values of the nutations
*     using all of the avaiable models.
* MOD TAH 890530 Added planetary nutations to computation
 
      include '../includes/const_param.h'
      include  '../includes/nut_eval.h'
 
*         i         - Loop counters
 
      integer*4 i
 
*       jd      - Julian date of current observation
*       c(2)    - Computed value for nutation in longtide and
*               - obliquity (mas)
*   arg             - value of the argument (rads)
*   arg_period      - Period of argument in sidereal days
*   cent            - Number of centuries to J2000
 
      real*8 jd, c(2), arg, arg_period, cent
 
***** Get the number of centuries to current date
 
*                                       ! Julian century
      cent = ( jd - dj2000 )/36525.d0
 
*     Loop over all the terms in the nutation series
      do i = 1,2
          C(i) = 0.d0
      end do
 
      do i = 1, num_series
 
          call fundamental_arg( jd, arg_series(1,i),  arg, arg_period )
 
*         Sum in the contribution to the angle
          C(1) = C(1) + ( val_nut_coeff(1,i) +
     .                    val_nut_rate (1,i)*cent )*sin(arg)+
     .                  ( val_nut_coeff(2,i) +
     .                    val_nut_rate (2,i)*cent )*cos(arg)
 
          C(2) = C(2) + ( val_nut_coeff(3,i) +
     .                    val_nut_rate (3,i)*cent )*cos(arg)+
     .                  ( val_nut_coeff(4,i) +
     .                    val_nut_rate (4,i)*cent )*sin(arg)
 
      end do
 
      do i = 1, num_plan
          call planetary_arg( jd, arg_plan(1,i), arg, arg_period)
 
*         Add contribution
          C(1) = C(1) +  val_plan(1,i) * sin(arg) +
     .                   val_plan(2,i) * cos(arg)
          C(2) = C(2) +  val_plan(3,i) * cos(arg) +
     .                   val_plan(4,i) * sin(arg)
 
      end do
 
***** Thats all
      return
      end
 
