CTITLE MEAN_POLE

      subroutine mean_pole( mjd_usr, model, mpx, mpy )

      implicit none

*     Routine to compute the mean pole postion using different models
*     MODEL Choices are:
*     IERS96  -- Linear model
*     IERS10  -- Combination of linear and cubic model
* MOD TAH 200217: Added IERS20 -- Linear model again
*     IERS20  -- IERS Conventions 7,1.4 (eqn 21) 2018/02/01 (ITRF2020)

* PASSED VARIABLES
      real*8 mjd_usr  ! MJD for calcuation (JD may be passed and will be converted)
      real*8 mpx, mpy  ! Mean X and Y pole positon (mill-arc-sec) (Same as expected
                      ! in comp_ptide): IERS20 Now called secular poles 
                      ! although variable name unchanged.
      
      character*(*) model  ! Model name to be used; IERS96 IERS10

* LOCAL VARIABLES
      real*8 mjd     ! MJD converted (days)
      real*8 dt     ! Time diffence in years

****  OK; Check that MJD if OK
      mjd = mjd_usr
      if( mjd.gt. 2 400 000.d0 ) mjd = mjd_usr - 2 400 000.5d0

      dt = (mjd-51544.0d0)/365.25d0   ! time in years since 2000.0

****  See which model
      if( model(1:6).eq.'IERS96' ) then
          dt = (mjd - 51544.0d0)/365.25d0   ! time in years
*         Linear trend from IERS Conventions 2000.
          mpx = (0.054d0+0.00083d0*dt)*1000.d0
          mpy = (0.357d0+0.00395d0*dt)*1000.d0
      else if( model(1:6).eq.'IERS10' ) then
*         This model is in two parts
          if( mjd.lt.55197.0d0 ) then   ! Data before 2010
*             Cubic trend from IERS Conventions 2010.
              mpx =  55.974d0 + 1.8243d0*dt + 0.18413d0*dt**2 
     .                                      + 0.007024d0*dt**3
              mpy = 346.346d0 + 1.7896d0*dt - 0.10729d0*dt**2
     .                                      - 0.000908d0*dt**3
          else     ! Linear trend after 2010.00
              mpx =   23.513d0 + 7.6141d0*dt
              mpy =  358.891d0 - 0.6287d0*dt
          endif
       elseif( model(1:6).eq.'IERS20' ) then
*         Linear model and now call secular pole.
*         xs = 55.0+1.677*(t-2000), ys = 320.5+3.460*(t-2000)  (21)
          mpx =  55.0d0  + 1.677d0*dt
          mpy = 320.5d0  + 3.460d0*dt
       elseif( model(1:4).eq.'ZERO' ) then
          mpx = 0.d0
          mpy = 0.d0
       else
*         Model found: Must be a program bug so fatal
          call report_stat('fatal','globk','Mean_pole',model,
     .        'Model name not IERS96/IERS10',0)
       end if

*****  Thats all 
       return 
       end
 
