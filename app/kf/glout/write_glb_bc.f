CTITLE write_glb_bcomp

      subroutine write_glb_bcomp( iout, options, cm, cov_parm,
     .                            sol_parm )

      implicit none  
 
*     This routine will write out the rates of changes of baseline 
*     N,E, and U and Length.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   cm      - Component cm=1 for lengths, cm=2 for rates
*   iout    - Output LU
*   icnt    - Counter for number of values output
*   iel,jel,kel - Pointers to positions in vectors
*   i,j,k   - Loop counters
*   options - OPtions parameter for output
 
      integer*4 cm, iout, jel,kel, i,j,k,l, options
 
 
*   base(3) - Components of the baseline
*   base_fin(6) - Final baseline components or rate
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix
*           - of estimated parameter
*   covar(6,6)  - Covarianve matrix of XYZ at both sites
 
*   dellen  - change in length or rate from apriori
*   dt(2)   - Epoch differences between the referenve date for
*           - the two sites in the baseline, and the current
*           - epoch
*   length  - Length of the baseline
*   lenadj  - Baseline length after adjustment (used when apriori
*             length is zero.)
*   comp_adj(4)  - Adjustments to Length, NE and U rates of changes
*   comp_fin(4)  - Final values for Length NE and U  rates of changes
*   loc_coord(3,2)    - Local coordinates a site in the baseline of both sites
*                  in the baseline 
*   pos_geod(3,2)   - Geodetic positions of the two sites
*   pos_xyz_adj(6)  - Adjustment to the XYZ position or rates
*           - to the two sites in the basline
*   pos_neu_adj(6)  - Adjustment to the NEU position or rates
*           - to the two sites in the basline
*   scr_real(10)    - Scratch array for calculations
*   sol_parm(num_glb_parn)              - Solution vector
*           - of estimated parameters
*   temp_covar(36)  - Scrach covar need for var_comp
*   trans(4,6)    - Trans formation from XYZ to baseline length, and NEU (Formed
*                   such that first row is length, then remaing rows are NE and U
*   rot_matrix(3,3)  - Rotation matrix for each station which rotates into
*                 local NEU at the station
*   cov_lneu(4,4)    - Covariance matrix for L, N E and U
*   rhone            - Correlation between north and East
*   pos_xyz_apr(6)   - Aproiri values for the positions of the two sites
 
      real*8 base(6), base_fin(6), cov_parm(num_glb_parn,num_glb_parn),
     .    covar(6,6), dt(2), loc_coord(3,2), pos_geod(3,2), length, 
     .    comp_adj(4), comp_fin(4),lenadj,
     .    pos_xyz_adj(6), pos_neu_adj(6), pos_xyz_apr(6), 
     .    sol_parm(num_glb_parn), temp_covar(36), trans(4,6), 
     .    rot_matrix(3,3), cov_lneu(4,4), rhone

*   qoutbas   - Logical function to see if should output this baseline
*               only if it is all markov

      logical qoutbas
 
      common covar, temp_covar
 
***** Output header
      if( cm.eq.1 ) then
          write(iout,100) 'BASELINE COMPONENTS'
          write(iout,110)
  110     format('   Baseline  ',18x,' North ', 26x, 
     .           'East', 20x,'Rne', 18x,'Height',/,
     .     22x,2('     Est.        Adj.   +- ',4x),10x,
     .           '     Est.        Adj.    +- ',/,
     .     32x,2('(mm)',25x),20x,'(mm)')
      else
          write(iout,100) 'BASELINE COMPONENT RATES OF CHANGE'
          write(iout,120)
  120     format('   Baseline  ',15x,' Length ', 16x, ' North ', 16x,
     .           'East', 15x,'Rne', 3x,'Height',/,
     .     23x,3(' Est.   Adj.  +- ',3x),8x,' Est.   Adj.  +- ',/,
     .     29x,3('(mm/yr)',14x),8x,'(mm/yr)')
      end if
  100 format(/,' GLOBK: ',a)
 
****  Now loop over all baselines
      do i = 1, gnum_sites-1
 
          dt(1) = (gepoch_out-site_epoch(i))/365.25d0
 
          do j = i+1, gnum_sites

*           Check to see if we should output this baseline
            if( qoutbas( options, i,j) ) then
 
              dt(2) = (gepoch_out-site_epoch(j))/365.25d0
 
*             Get the components of the baseline
*                             ! XYZ
              do k = 1,3
 
*                 Clear adjustments
                  pos_xyz_adj(k)   = 0.d0
                  pos_xyz_adj(k+3) = 0.d0
 
*                 Get the aproiri baseline
                  base(k) =   apr_val_site(k,1,i)
     .                      + apr_val_site(k,2,i)*dt(1) 
     .                      + cont_nonsec(k,i)
                  base(k+3) = apr_val_site(k,1,j)
     .                      + apr_val_site(k,2,j)*dt(2)
     .                      + cont_nonsec(k,j)
 
*                                                        ! We will add
*                                                        ! ajustments next
                  base_fin(k)   =  apr_val_site(k,cm,i) 
                  base_fin(k+3) =  apr_val_site(k,cm,j)
                  if( cm.eq.1 ) then
                      base_fin(k)   = base_fin(k) +
     .                                apr_val_site(k,2,i)*dt(1) +
     .                                cont_nonsec(k,i)
                      base_fin(k+3) = base_fin(k+3) +
     .                                apr_val_site(k,2,j)*dt(2) +
     .                                cont_nonsec(k,j)
                  end if
 
*                 Now get the compenent of the adjustment
*                                                   ! NOTE: either value
                  if( parn_site(k,cm,i).ne.0 ) then
*                                                   ! or rate
                      jel = parn_site(k,cm,i)
                      pos_xyz_adj(k) =  sol_parm(jel)
                      base_fin(k) = base_fin(k) + pos_xyz_adj(k)
                  end if
 
*                                                   ! NOTE: either value
                  if( parn_site(k,cm,j).ne.0 ) then
*                                                   ! or rate
                      kel = parn_site(k,cm,j)
                      pos_xyz_adj(k+3) = sol_parm(kel)
                      base_fin(k+3) = base_fin(k+3) + pos_xyz_adj(k+3)
                  end if

*                 Save final position for transformation (reverse sign so that
*                 true station position is saved)
                  pos_xyz_apr(k)   =  base(k)
                  pos_xyz_apr(k+3) =  base(k+3) 
*                                 ! Looping over components
              end do
 
*****         Get the transformation for components to length
              length = sqrt((base(1)-base(4))**2 + 
     .                      (base(2)-base(5))**2 + 
     .                      (base(3)-base(6))**2   )
              if( length.gt.1.d-6) then
                  do k = 1,3
                      trans(1,k)   =  (base(k) - base(k+3))/length
                      trans(1,k+3) = -(base(k) - base(k+3))/length
                  end do
              else

*                 If the apriori length is zero, project the baseline
*                 into the adjustment direction.
                  lenadj = sqrt((base_fin(1)-base_fin(4))**2 + 
     .                      (base_fin(2)-base_fin(5))**2 + 
     .                      (base_fin(3)-base_fin(6))**2   )
                  do k = 1,3
                      trans(1,k)   =  (base_fin(k) - base_fin(k+3))
     .                                /lenadj
                      trans(1,k+3) = -(base_fin(k) - base_fin(k+3))
     .                                /lenadj
                  end do
              end if
                  

*****         Get transformation from XYZ to NEU.  Do for First site
              call rotate_geod(pos_xyz_adj(1), pos_neu_adj(1), 'XYZ',
     .                         'NEU', pos_xyz_apr(1), loc_coord(1,1), 
     .                         rot_matrix) 

*             copy rotation into transformation
              do k = 1,3
                 do l = 1,3
                    trans(k+1,l) = -rot_matrix(k,l)
                 end do
              end do

*****         Get transformation from XYZ to NEU.  Do for Second site 
              call rotate_geod(pos_xyz_adj(4), pos_neu_adj(4), 'XYZ',
     .                         'NEU', pos_xyz_apr(4), loc_coord(1,2),
     .                         rot_matrix)

*             copy rotation into transformation
              do k = 1,3
                 do l = 1,3
                    trans(k+1,l+3) =  rot_matrix(k,l)
                 end do
              end do

*****         Get adjustment  and variance
              do k = 1,4
                  call dvdot(comp_adj(k), trans(k,1) ,4,
     .                       pos_xyz_adj,1, 6)
              end do
 
*****         Get the final estimates.  Because baseline values are so large
*             we difference the geod pos estimates instead of computing through
*             the linear transformation

              if( cm.eq.1 ) then

*****             Get the final values for the length and N,E and U.  For full
*                 values we need to compute directly.
                  call loc_to_geod( loc_coord(1,1), pos_geod(1,1) )
                  call loc_to_geod( loc_coord(1,2), pos_geod(1,2) )
                  comp_fin(1) = length + comp_adj(1)
                  do k = 2,4
                      comp_fin(k) = pos_geod(k-1,2) - pos_geod(k-1,1) + 
     .                              comp_adj(k)
C                     comp_fin(k) = mod( comp_fin(k),10.d0 )
C                     comp_fin(k) = comp_fin(k) - 
C    .                                  nint(comp_fin(k)/10.d0)*10.d0
                  end do
              else
                  
*****             Get final value of rate
                  do k = 1,4
                     call dvdot(comp_fin(k), trans(k,1),4, base_fin,1,6)
                  end do
              end if
 
*             Move the covariance elements into place so that we
*             can compute sigma
              call mov_cov(parn_site(1,cm,i), parn_site(1,cm,i),3,
     .                     covar(1,1), 6, cov_parm, num_glb_parn)
              call mov_cov(parn_site(1,cm,i), parn_site(1,cm,j),3,
     .                     covar(1,4), 6, cov_parm, num_glb_parn)
              call mov_cov(parn_site(1,cm,j), parn_site(1,cm,i),3,
     .                     covar(4,1), 6, cov_parm, num_glb_parn)
              call mov_cov(parn_site(1,cm,j), parn_site(1,cm,j),3,
     .                     covar(4,4), 6, cov_parm, num_glb_parn)
 
*             Now get variance and covariances
              call var_comp(trans, covar, cov_lneu, temp_covar,4,6,1)
 
*****         Finally output the results (only final answer is not zero
*             anf the variance is positive)
*                                     !  Convert everything to mm or mm/yrs
              do k = 1,4
                 comp_fin(k) = comp_fin(k)*1000.d0
                 comp_adj(k) = comp_adj(k)*1000.d0
              end do

              if( cov_lneu(2,2)*cov_lneu(3,3).gt. 0.d0 ) then
                 rhone = cov_lneu(2,3)/sqrt(cov_lneu(2,2)*cov_lneu(3,3))
              else
                 rhone = 0.001d0
              end if

*             See if sigma is small enough to output.  Value
*             corresponds to either 1 m or 1 m/yr.
              if( cov_lneu(2,2).lt.100.d0  .and.
     .            cov_lneu(3,3).lt.100.d0  .and.
     .            cov_lneu(4,4).lt.100.d0  ) then
                 if( cm.eq.1 ) then
                     write(iout, 310) gsite_names(j), gsite_names(i), 
     .                     (comp_fin(k), comp_adj(k), 
     .                      1000.d0*sqrt(abs(cov_lneu(k,k))), k = 2,3 ),
     .                      rhone, comp_fin(4), comp_adj(4),
     .                      1000.d0*sqrt(abs(cov_lneu(4,4)))
* MOD TAH 940113: Changed format to allow full component out
 310                 format(1x,a8,'-',a8,2(2x,f15.1,1x,f7.1,1x,f6.1),
     .                       2x,f6.3,
     .                       2x,f11.1,1x,f9.1,1x,f7.1)
                 else
* MOD TAH 920904: Changed format to get one extra digit.
                     write(iout, 320) gsite_names(j), gsite_names(i), 
     .                     (comp_fin(k), comp_adj(k), 
     .                      1000.d0*sqrt(abs(cov_lneu(k,k))), k = 1,3 ),
     .                      rhone, comp_fin(4), comp_adj(4),
     .                      1000.d0*sqrt(abs(cov_lneu(4,4)))
 320                 format(1x,a8,'-',a8,3(2x, f7.2,1x,f7.2,1x,f7.2),
     .                      1x,f5.3, 2x, f7.2,1x,f7.2,1x,f7.2)
                 end if
              end if
*                     ! if we should output
            end if
*                     ! Looping over second site
          end do
*                     ! Looping over first site
      end do
 
***** Thats all
      return
      end
 
