CTITLE OUT_TYPE_STATS
 
      subroutine out_type_stats ( iout )

      implicit none  
 
*     Routine to output the statistcs for each type of parameter
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   i       - Loop counter
*   iout    - Ouput unit number
 
      integer*4 i, iout
 
*   chi2    - Chisquared per number
*   wrms    - Wrms scatter
 
      real*4 chi2, wrms
 
*   type_names(max_chi_types)   - Names and units of the parameter
*                               - types
 
      character*16 type_names(max_chi_types)
 
      data type_names / 'X coords.    (m)'
     .,                 'Y coords.    (m)'
     .,                 'Z coords.    (m)'
     .,                 'Axis offs.   (m)'
     .,                 'RA source  (mas)'
     .,                 'Dec source (mas)'
     .,                 'Pole pos.  (mas)'
     .,                 'UT1        (mas)'
     .,                 'Nut. Angs. (mas)'
     .,                 'Tide l          '
     .,                 'Tide h          '
     .,                 'Lag angle  (deg)'
     .,                 'ETD 1        (m)'
     .,                 'ETD 2        (m)'
     .,                 'ETD 3        (m)'
     .,                 'ETD 4        (m)'
     .,                 'ETD 5        (m)'
     .,                 'ETD 6        (m)'
     .,                 'ETD 7        (m)'
     .,                 'ETD 8        (m)'
     .,                 'Gamma           '
     .,                 'N coordinate (m)'
     .,                 'E coordinate (m)'
     .,                 'U coordinate (m)' /
 
 
***** Write header
 
      write ( iout, 100)
  100 format(/' POSTFIT PARAMETER SCATTTER BY TYPE ',/,
     .        ' TYPE ',t25,'WRMS',t33,'Chi**2',t48,'#')
 
      do i = 1, max_chi_types
 
          if( sum_type_num(i).gt.0 ) then
              chi2 = sum_type_res(i)/sum_type_num(i)
              wrms = sqrt(sum_type_num(i)/sum_type_wgh(i))*sqrt(chi2)
 
              write(iout,150) type_names(i), wrms, chi2,
     .                        sum_type_num(i)
  150         format(1x,a,1x,f10.5,1x,f10.3,5x,i5)
          end if
      end do
 
***** Thats all
      return
      end
 
