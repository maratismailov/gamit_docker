CTITLE CORRECT_OBS
 
      subroutine correct_obs( final, iel, obs_corr, ndimo )

      implicit none 
  
*     Routine to add obs_corr(iel) to final if the dimension of
*     of obs_corr (ndimo) is greater then zero
 
*   iel     - Element of obs_corr to add
*   ndimo   - Dimension of obs_corr
 
      integer*4 iel, ndimo
 
*   final   - Value to be incremented
*   obs_corr(1) - Observation correction vector
 
      real*8 final, obs_corr(1)
 
 
****  If dimension is greater than add, otherwise leave final unchanged
 
      if( ndimo.gt.0 ) then
          final = final + obs_corr(iel)
      end if
 
****  Thats all
      return
      end
 
