CTITLE WRITE_PMU_CORR
 
      subroutine write_pmu_corr( unit, cov_parm, num_par, nx, ny, nu)

      implicit none 
 
*     Routine to write out the correlation between x and y pole
*     position and UT1.
 
* PASSED VARIABLES
 
*   unit        - Output unit number
*   num_par     - Total number of parameters estimated
*   nx, ny, nu  - Parameter numbers for pole position and UT1
 
      integer*4 unit, num_par, nx, ny, nu
 
*   cov_parm(num_par, num_par)  - Full covariance matrix.
 
      real*8 cov_parm(num_par, num_par)
 
* LOCAL VARIABLES
 
*   cov_xy, cov_xu, cov_yu  - Covariances between x,y and ut1
*   sig_x, sig_y, sig_u     - Sigmas of x, y, and ut1 (unless
*                           - zero in which case set to small
*                           - number of aviod dividing by 0.0)
 
      real*8 cov_xy, cov_xu, cov_yu, sig_x, sig_y, sig_u
 
****  Set nominal values
      cov_xy = 0.0
      cov_xu = 0.0
      cov_yu = 0.0
 
      sig_x = 0.001
      sig_y = 0.001
      sig_u = 0.001
 
****  Get each of the covarinace and sigmas
 
      if( nx.gt.0 ) then
          sig_x = sqrt(cov_parm(nx,nx))
        if( sig_x.lt.1.d-6 ) sig_x = 1.d-6
          if( ny.gt.0 ) cov_xy = cov_parm(nx,ny)
          if( nu.gt.0 ) cov_xu = cov_parm(nx,nu)
      end if
 
      if( ny.gt.0 ) then
          sig_y = sqrt(cov_parm(ny,ny))
        if( sig_y.lt.1.d-6 ) sig_y = 1.d-6
          if( nu.gt.0 ) cov_yu = cov_parm(ny,nu)
      end if

      if( nu.gt.0 ) then 
          sig_u = sqrt(cov_parm(nu,nu))
          if( sig_u.lt.1.d-6 ) sig_u = 1.d-6
      end if
 
****  Now write out the line only if something estimated
      if( nx+ny+nu.gt.0 ) then
          write(unit,100) cov_xy/(sig_x*sig_y),
     .                    cov_xu/(sig_x*sig_u),
     .                    cov_yu/(sig_y*sig_u)
 100      format(5x,' Pole/UT1 correlations: XY, XU, YU  ',t45,3f11.4)
      end if
 
****  Thats all
      return
      end
 
 
