CTITLE OUT_GLBP
 
      subroutine out_glbp(iout, parn, est,sol,cov, conversion,
     .   first,second,units, type)

      implicit none  
 
*     This routine will see if an estimate should be output.
*     it then calls 'write_line' to write line to iout.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
c
c Variables
c ---------
c iout -- output LU number
c parn -- the parameter number (This will tell us if parameter
c     estimated.
c est  -- the estimated value of the parameter
c sol  -- the adjustments to the parameter (ema value)
c cov  -- the covariance matrix of the adjustment (ema value)
c conversion -- any conversion factor for changing units
c first  -- first label to be output.
c second -- second label to be output
c units  -- the units of the adjustment
c type   -- type of conversion on 'est' (0 -- no conversion; 1 convert
c     to deg, min, sec)
c
      integer*4 iout, parn, type
 
c
      real*8 est, sol(1), cov(num_glb_parn,1), conversion
 
c
c
      character*(*) first, second, units
 
c
c Local variables
c ---------------
c scr_real -- scratch real*8 array.  Used as needed
c
      real*8 scr_real(3)
 
c
c
c.... See if parameter estimated i.e. does it have a parameter number
*                           ! this parameter estimated
      if( parn.gt.0 ) then
*  scr_real(1) = sol(parn)*conversion      - convert to ouput units
      scr_real(1) = sol(parn)*conversion
 
      scr_real(2) = sqrt(cov(parn,parn))*abs(conversion)
 
         call write_line(iout, parn, first,second,units, est,
     .      scr_real(1),scr_real(2), type)
      end if
c
      return
      end
 
