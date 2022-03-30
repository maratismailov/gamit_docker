      subroutine comp_ion( ion_delay, cf_omcs, cf_obsv, lv )
 
*     Routine to compute the ionspheric delay from cf_omcs(1,2) and 
*     apply to the observations and ranges omcs.  Since solve 
*     eliminates the ion-delay, this update should not affect 
*     results from solve.  It is needed for GLONASS where all data
*     are mapped to the same frequency.

      implicit none

*  MATLAB TEST CODE
*     lif1 =  1.0/(1.0 - (fL1/fL2)^2) ;
*     lif2 = -(fL1/fL2)/(1 - (fL1/fL2)^2) ;
*     iondelay = zeros(1,4);

*     iest = lif1*l1p + lif2*l2p;

*     iondelay(1) = iest ;
*     iondelay(2) = iondelay(1)*fL1/fL2 ;
*     iondelay(3) = -iondelay(1) ;
*     iondelay(4) = -iondelay(2) ;

*     l1pf = l1p - iondelay(1);
*     l2pf = l2p - iondelay(2); 
*     l1rf = l1r - iondelay(3); 
*     l2rf = l2r - iondelay(4); 

*     p1rf = p1r - iondelay(3)*c/fL1 ;
*     p2rf = p2r + iondelay(4)*c/fL2 ;
*     fprintf('Raw phs %8.3f %8.3f Rng %8.3f %8.3f\n', l1p, l2p, l1r, l2r);
*     fprintf('Fix phs %8.3f %8.3f Rng %8.3f %8.3f\n', l1pf, l2pf, l1rf, l2rf);

*     iestf = lif1*l1pf + lif2*l2pf;
*     fprintf('Ion raw %8.3f Ion fix %8.3f\n',iest,iestf);
*     disp(lif1*l1rf + lif2*l2rf);
*     disp((p1rf-p2rf)/(1-(fL1/fL2)^2));


 
 
* INCLUDES FILES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
      real*8 ion_delay(4)   ! Computed ion delay at fl1 and fl2 for 
                            ! phase (cycles) and range (cycles). Cycle
                            ! is appropriate for frenency.  Value
                            ! computed and applied here.
      real*8 cf_omcs(max_gdata_types)  ! O-minus-C for phase and range
                            ! all in cycles (no clear what 5th value in 
                            ! this array could be) INPUT/OUTPUT
     .,     cf_obsv(max_gdata_types)  !  Observed phase and range but
                            ! range is in meters units.INPUT/OUTPUT

      integer*4 lv          ! Satellite vehicle number (used to get 
                            ! correct frequencies) INPUT


***** Compute the multiplier factors (could add to set_freqs)
*     Factors computed in set_freqs
C     lif1(lv) =   1.0d0/(1.0d0 - (fL1(lv)/fL2(lv))**2) 
C     lif2(lv) = -(fL1(lv)/fL2(lv))/(1.d0 - (fL1(lv)/fL2(lv))**2)

*     Compute values using the factors that are computed in the set_freqs,f
*     routine.
C     write(*,210) lv, cf_omcs
C210  format('APP_ION CH ',i2,' PRE OMCS ',5F10.2)
      ion_delay(1) = lif1(lv)*cf_omcs(1) + lif2(lv)*cf_omcs(2)
      ion_delay(2) = ion_delay(1)*fl1(lv)/fl2(lv)
      ion_delay(3) = -ion_delay(1)
      ion_delay(4) = -ion_delay(2)

*     Now remove ion delay from omc values
      cf_omcs(1) = cf_omcs(1) - ion_delay(1)
      cf_omcs(2) = cf_omcs(2) - ion_delay(2)
      cf_omcs(3) = cf_omcs(3) - ion_delay(3)
      cf_omcs(4) = cf_omcs(4) - ion_delay(4)

*     Now update observables (noting ranges are in meters).
      cf_obsv(1) = cf_obsv(1) - ion_delay(1)
      cf_obsv(2) = cf_obsv(2) - ion_delay(2)
      cf_obsv(3) = cf_obsv(3) - ion_delay(3)*vel_light/fl1(lv)
      cf_obsv(4) = cf_obsv(4) - ion_delay(4)*vel_light/fl2(lv)

C     write(*,220) lv, ion_delay, cf_omcs
C220  format('APP_ION CH ',i2,' ION ',4F10.2,' OMCS ',5F10.2)

***** Thats all
      return 
      end




