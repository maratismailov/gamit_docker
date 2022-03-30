CTITLE WRITE_GLB_BASEL

      subroutine write_glb_basel( iout, options, cm, cov_parm,
     .                            sol_parm )

      implicit none  
 
*     Routine to write either baseline lengths or their rates of
*     change.  At some future time, we may add vector rates as well
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   cm      - Component cm=1 for lengths, cm=2 for rates
*   iout    - Output LU
*   icnt    - Counter for number of values output
*   iel,jel,kel - Pointers to positions in vectors
*   i,j,k   - Loop counters
*   options - OPtions parameter for output
 
      integer*4 cm, iout, jel,kel, i,j,k, options
 
 
*   base(3) - Components of the baseline
*   base_fin(3) - Final baseline components or rate
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix
*           - of estimated parameter
*   covar(6,6)  - Covarianve matrix of XYZ at both sites
 
*   dellen  - change in length or rate from apriori
*   dt(2)   - Epoch differences between the referenve date for
*           - the two sites in the baseline, and the current
*           - epoch
*   length  - Length of the baseline
*   lenadj  - Length or rate adjustment
*   lenfin  - Final length or rate
*   loc_coord(3)    - Local coordinates a site in the baseline
*   pos_xyz_adj(6)  - Adjustment to the XYZ position or rates
*           - to the two sites in the basline
*   pos_neu_adj(6)  - Adjustment to the NEU position or rates
*           - to the two sites in the basline
*   scr_real(10)    - Scratch array for calculations
*   sol_parm(num_glb_parn)              - Solution vector
*           - of estimated parameters
*   temp_covar(36)  - Scrach covar need for var_comp
*   trans(6)    - Trans formation from XYZ to baseline length
 
 
      real*8 base(3), base_fin(3), cov_parm(num_glb_parn,num_glb_parn),
     .    covar(6,6), dt(2), length, lenadj, lenfin,
     .    pos_xyz_adj(6), scr_real(10),
     .    sol_parm(num_glb_parn), temp_covar(36), trans(6)
 

* qoutbas   - Logical function which checks to see only markov
*       baselines should be output.  (Uses options set by user)

      logical qoutbas
 
 
      common covar, temp_covar
 
***** Output header
 
      if( cm.eq.1 ) then
          write(iout,100) 'BASELINE LENGTHS'
  100     format(/,' GLOBK: ',a,/,
     .         t6,'BASELINE',t47,'Length (m)',t64,'Adjust (m)',
     .         t77,'Sigma (m)')
      else
          write(iout,120) 'BASELINE LENGTH RATES'
  120     format(/,' GLOBK: ',a,/,
     .         t6,'BASELINE',t47,'dL/dT (m/yr)',t62,'Adjust (m/yr)',
     .         t75,'Sigma (m/yr)')
      end if
 
****  Now loop over all baselines
      do i = 1, gnum_sites-1
 
          dt(1) = (gepoch_out-site_epoch(i))/365.25d0
 
          do j = i+1, gnum_sites

*           Check to see if we should output
            if( qoutbas(options, i,j) ) then
 
              dt(2) = (gepoch_out-site_epoch(j))/365.25d0
 
*             Get the components of the baseline
*                             ! XYZ
              do k = 1,3
 
*                 Clear adjustments
                  pos_xyz_adj(k)   = 0.d0
                  pos_xyz_adj(k+3) = 0.d0
 
*                 Get the aproiri baseline
                  base(k) = - apr_val_site(k,1,i)
     .                      - apr_val_site(k,2,i)*dt(1)
     .                      + apr_val_site(k,1,j)
     .                      + apr_val_site(k,2,j)*dt(2)

                  if( cm.eq.1 ) then
                      base(k) = base(k) - cont_nonsec(k,i)
     .                                  + cont_nonsec(k,i)
                  end if 
 
*                                                        ! We will add
                  base_fin(k) = -apr_val_site(k,cm,i) +
*                                                        ! ajustments next
     .                           apr_val_site(k,cm,j)
 
*                                          ! Add velocity effect to
                  if( cm.eq.1 ) then
*                                          ! position
                      base_fin(k) = base_fin(k)               -
     .                              apr_val_site(k,2,i)*dt(1) -
     .                              cont_nonsec(k,i)          +
     .                              apr_val_site(k,2,j)*dt(2) +
     .                              cont_nonsec(k,j)
                  end if

 
*                 Now get the compenent of the adjustment
*                                                   ! NOTE: either value
                  if( parn_site(k,cm,i).ne.0 ) then
*                                                   ! or rate
                      jel = parn_site(k,cm,i)
                      pos_xyz_adj(k) = sol_parm(jel)
                      base_fin(k) = base_fin(k) - pos_xyz_adj(k)
                  end if
 
*                                                   ! NOTE: either value
                  if( parn_site(k,cm,j).ne.0 ) then
*                                                   ! or rate
                      kel = parn_site(k,cm,j)
                      pos_xyz_adj(k+3) = sol_parm(kel)
                      base_fin(k) = base_fin(k) + pos_xyz_adj(k+3)
                  end if
*                                 ! Looping over components
              end do
 
*****         Get the transformation for components to length
              length = sqrt(base(1)**2 + base(2)**2 + base(3)**2)
              if( length.gt.1.d-6 ) then
                  do k = 1,3
                      trans(k)   = -base(k)/length
                      trans(k+3) =  base(k)/length
                  end do
              else

*                 If the apriori length is zero then we project
*                 baseline into the direction of the adjustment.
                  lenfin = sqrt(base_fin(1)**2+base_fin(2)**2+
     .                          base_fin(3)**2)
                  do k = 1,3
                      trans(k)   = -base_fin(k)/lenfin
                      trans(k+3) =  base_fin(k)/lenfin
                  end do
              end if
 
*****         Get final value of length or rate
              if( cm.eq.1 ) then
                  lenfin = sqrt(base_fin(1)**2 + base_fin(2)**2 +
     .                          base_fin(3)**2 )
*                                 ! Dot velocity into length
              else
                  call dvdot(lenfin, trans(4),1, base_fin,1,3)
              end if

*****         Get adjustment  and variance
              call dvdot(lenadj, trans,1, pos_xyz_adj,1, 6)
 
*             Move the covariance elements into place so that we
*             can compute sigma
              call mov_cov(parn_site(1,cm,i), parn_site(1,cm,i),3,
     .                     covar(1,1), 6, cov_parm, num_glb_parn)
              call mov_cov(parn_site(1,cm,i), parn_site(1,cm,j),3,
     .                     covar(4,1), 6, cov_parm, num_glb_parn)
              call mov_cov(parn_site(1,cm,j), parn_site(1,cm,i),3,
     .                     covar(1,4), 6, cov_parm, num_glb_parn)
              call mov_cov(parn_site(1,cm,j), parn_site(1,cm,j),3,
     .                     covar(4,4), 6, cov_parm, num_glb_parn)
 
*             Now get variance
              call var_comp(trans, covar, scr_real, temp_covar,1,6,0)
 
*****         Finally output the results (only final answer is not zero
*             anf the variance is positive)
              if( scr_real(1).gt.0 .and. lenfin.ne.0.d0 .and.
     .            scr_real(1).lt.100.d0 ) then
                  scr_real(1) = sqrt(scr_real(1))
 
                  call write_line(iout,0, gsite_names(i),'to ',
     .                gsite_names(j),lenfin, lenadj, scr_real,0)
              end if
*                     ! If we should output this baseline
            end if
*                     ! Looping over second site
          end do
*                     ! Looping over first site
      end do
 
***** Thats all
      return
      end
 
