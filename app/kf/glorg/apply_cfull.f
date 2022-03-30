CTITLE APPLY_COND_full

      subroutine apply_cond_full(iout, parn, ndim, pmu_parts, num,
     .            cov_parm, sol_parm, site_pos, cond_var,
     .            numc, cond_part, origin_gain, cat, used,
     .            cnd_parts_bits, cnd_hgt_var, use_ratio, 
     .            type, site_names, ss_var, ss_err, ss_nsig, 
     .            ss_rel, ss_min_dh, ss_min_rms, ss_min_dne, 
     .            ss_it, options, descr)

      implicit none 
 
*     Routine to apply conditions so that the net rotation and
*     translation of the coordinate system will be be zero.
 
*     This version fully computes the results so that the baselines
*     lengths will change if system is not free to translate)

      integer*4 max_stab_site
      parameter ( max_stab_site = 4096 )
 
*     Also extended to allow for scale change
 
*   iout                - Output unit
*   ndim                - Second dimension for parn array (
*                       - 2 is postion and velocity, 1 for
*                       - postion only (solvk)
*   num                 - Number of sites in the parn array
*   numc                - NUmber of rows and columns in the
*                       - covariance matrix.
*   parn(3,ndim,num)    - Parameter numbers for the sites
*   used(*)             - Bit mapped array saying which sites
*                       - to use.
*   cnd_parts_bits - Bit mapped word which sets which partials are to be
*                    estimated.  Bits 1-3 -- rotations, bits 4-6 translations
*                    bit 7 - scale

*   ss_it - Number of stabilization iteration 
*   options  -- Print options passed to glorg.  Used to see if details
*               to be printed.
 
      integer*4 iout, ndim, num, numc, parn(3,ndim,num), used(*),
     .    cnd_parts_bits, ss_it, options 
 
*   cov_parm(numc,numc) - Covariance matrix.
*   sol_parm(numc,numc) - Solution vector
*   site_pos(3,ndim,num)- Apriori site postions
*   cond_var(7)         - Variance to be given to the constraint
*                         (normally this will be zero)
 
*   pmu_parts(3,3,num)  - Polar motion/UT1 partials.  Note: the rows
*                         are for site position componenst XYZ and
*                         columns are for x,y and UT1.
*   cond_part(7,numc)    - Condition partials.
*   origin_gain(numc,7) - Computed Kalman Gain for applying the
*                       - conditions.
*   cat(numc,7)         - Temporary storage of:
*                       -  m   T
*                       - C   A
*                       -  m+1
*   cnd_hgt_var         - Variance given to height relative to horizontal
*                         components
*   use_ratio           - Ratio of height variances for site to be used in
*                         orgin.

* MOD TAH 980417: Added iterations
*   ss_var(num)  - Variance of horizontal components of each site (from 
*     covariance matrix and used for weighting)
*   ss_err(num)  - Error in horizonatl position (used for editing)
*   ss_nsig      - Sigma ratio for site to be removed from list.
*   ss_rel       - Relative weight between constant and station dependant.
*   ss_min_dh    - MIn height sigam difference
*   ss_min_rms   - Min rms difference
*   ss_min_dne   - Min NE sigma difference (post stabilzation).
 
      real*8 cov_parm(numc,numc), sol_parm(numc),
     .    site_pos(3,ndim,num), cond_var(7),
     .    pmu_parts(3,3,num), cond_part(7,numc),
     .    origin_gain(numc,7), cat(numc,7), cnd_hgt_var, use_ratio,
     .    ss_var(num), ss_err(num), ss_nsig, ss_rel, ss_min_dh, 
     .    ss_min_rms, ss_min_dne

* Stab_var(max_stab_site) -- Variances of the position/velocity estimates
*     for the stabablization sites
      real*8 stab_var(max_stab_site)

*   type  - Type of stabilization (position or rate)
*   site_names  - Names of the sites

      character*(*) type, site_names(num)

*   descr -- Descriptor  string for run (list file name in glorg)

      character*(*) descr
 
* LOCAL VARIABLES
 
*   i,j,k,l         - Loop counters
*   ir,ic           - Row and col counters
*   iel             - Generic element in matrix
*   lel             - last elment used in cond_part
*   ipivot(7)       - Used by invert_vis
*   nc              - Number of conditions (3-rotations,
*                     6-rots and translations,
*                     7-rots, trans and scale.
*   len_des         - Length of descriptor string 
*   trimlen         - Length of string
*   nact            - Actual number of conditions to apply.  If
*                     the NAPP option in set, only the rotation
*                     is applied and the covarinace matrix not changed.
 
      integer*4 i,j,k,l,m, ir,ic, iel, lel, ipivot(7), nc, nr, ns,
     .          len_des, trimlen, nact, nes
 
*   kbit            - Check to see if bit is set.
 
 
      logical kbit
 
*   ptp(7,7)            - Normal equations for the pmu_part inversion
*   scale(7)            - Used by invert_vis
*   tmp                 - Temporary summation of covariance decrament.
*   cond_result(7)      - Actuall values for the translation, rotation
*                         and scale
*   parts(3,7)          - Partials of rotation, translation, scale
*                         with respect to site positions.
*   hgt_wgh(3,3)        - Covariance matrix with height uncertainty
*   postfit_chi         - Post fit chi**2 per degree of freedom
*   tmp_xyz(3)          - Direct calculation of change.
*   fin_xyz(3)          - Final position with adjustment
 
 
      real*8 ptp(7,7), scale(7), tmp, cond_result(7),
     .    dsol_cnd(7), parts(3,7), prefit_sum,
     .    postfit_sum, atf(7), hgt_wgh(3,3), postfit_chi,
     .    tmp_xyz(3), fin_xyz(3)

*   rad  -- Raduis of site
*   av_var_hgt, min_var_hgt, max_var_hgt -- Average, min and max height
*           variances
*   var_hgt -- Variance of height estimate.
*   var_hor -- Horizontal componentent variance
*   var_hor_all -- Sum of all horizontal variances so that we get average
*              unit weight.
*   var_tot -- Total position variance; this minus hght is the horizontal
*              varinace.

*   err_hgt -- Heigt error 
*   err_tot -- Total error in position
*   sig_postfit -- Postfit RMS scatter of positions.
*   var_postedt -- Variance to use for postfit editing (based on max
*                  of rms and ss_min_rms.
*   dh_tol      -- Tolerance to be used in deleting sites due to
*                  height sigma.
*   dsol       -- Change in solution element

      real*8 rad, var_hgt , av_var_hgt, min_var_hgt, max_var_hgt,
     .       var_hor, var_hor_all, var_tot, err_hgt, err_tot, 
     .       sig_postfit, dh_tol, var_postedt, dsol 

*   nu_var_hgt -- Number in height variances
*   jel, kel   -- Parameter pointers

      integer*4 nu_var_hgt , jel, kel
      
*   cmp_name(3) -- Names of coordinate components
      character*4 cmp_name(3)
      
      data cmp_name / 'dX', 'dY', 'dZ' /
           
 
****  Start, first computed the normal eqautions for the pmu_part
*     i.e., we want to compute the generalize inverse of pmu_part
      nc = 0
      lel = 0

      len_des = max(trimlen(descr),1)

* MOD TAH 980503: Check the number of sites; if too small turn off
*     the rotation stabbliization.
      ns = 0
      do i = 1, num
         if( kbit(used,i) ) then
             ns = ns + 1
         endif
      end do

      Write(*,100) ss_it, type, descr(1:len_des)
 100  format(' Iteration ',i3,' of ',a,' stabilization ',a)
      write(iout, 105) type, ss_it,  descr(1:len_des)
 105  format(/,100('='),/,
     .      ' Starting ',a,' stabilization iteration ',i3,
     .       1x,a)

****  Check number of sites
      if( ns.lt.3 ) then
          if( ns.eq.1 ) then
              do i = 1,5
                 call sbit(cnd_parts_bits,i,0)
              end do
              call sbit(cnd_parts_bits,7,0)
              write(iout,'(a)')
     .          ' ** Only 1 site in stabablization, Z-only used'
          else
              do i = 1,3
                 call sbit(cnd_parts_bits,i,0)
              end do
              call sbit(cnd_parts_bits,7,0)
              write(iout,'(a)')
     .          ' ** Only 2 sites in stabablization, Translation only'
          end if
      end if

*     Count number of conditions to be applied
      do i = 1,7
          if( kbit( cnd_parts_bits,i) ) nc = nc + 1
      end do
 
      do i = 1,nc
          dsol_cnd(i) = 0.d0
          do j = 1,nc
              ptp(i,j) = 0.d0
          end do
      end do
 
      nr = 0
      prefit_sum = 0
      postfit_sum = 0

****  Check the quality of the sites.  Compute the range of height variances
*     and discard any sites which have anomalously high variances
* MOD TAH 980418: Replaced mean with median
      av_var_hgt = 0
      nu_var_hgt = 0
      min_var_hgt = 1.d6 
      max_var_hgt = 0.d0 

      ns = 0

      do i = 1, num
         if( kbit(used,i) ) then

             rad = sqrt(site_pos(1,1,i)**2 + site_pos(2,1,i)**2 +
     .                  site_pos(3,1,i)**2 )
             var_hgt = 0.d0
             do j = 1,3
                jel = parn(j,1,i)
                do k = 1,3
                   kel = parn(k,1,i) 
                   if( jel.gt.0 .and. kel.gt.0 ) then
                       var_hgt = var_hgt + site_pos(j,1,i)* 
     .                       cov_parm(jel,kel)*site_pos(k,1,i)/rad**2 
                   end if
                end do
             end do
             if( var_hgt.eq.0.d0 ) then
                 call sbit(used,i,0)
             else 

*                Now see how this compares
                 if( var_hgt.lt.min_var_hgt ) min_var_hgt = var_hgt 
                 if( var_hgt.gt.max_var_hgt ) max_var_hgt = var_hgt 

*                Save value for computing median
                 ns = ns + 1
                 ss_err(ns) = var_hgt
C                av_var_hgt = av_var_hgt + var_hgt
                 nu_var_hgt = nu_var_hgt + 1
             end if 
         end if
      end do 

*     Write results
c     av_var_hgt = av_var_hgt/nu_var_hgt
*     Find the median:
      call mdian2( ss_err, ns, av_var_hgt )

* MOD TAH 980514: Compute the minimum dh sigma spread.
      dh_tol = max(sqrt(av_var_hgt)-sqrt(min_var_hgt), ss_min_dh)

      if( type(1:1).eq.'P' ) then
          write(iout,110) nu_var_hgt, sqrt(min_var_hgt)*1000, 
     .                    sqrt(max_var_hgt)*1000, 
     .                    sqrt(av_var_hgt)*1000, use_ratio*dh_tol*1000,
     .                    descr(1:len_des)
 110      format(' For ',i4,' sites in origin, min/max height sigma ',
     .           2F10.2,' mm; Median  ',F10.2,' mm, Tol ',F10.2,
     .                  ' mm ',a)
      else
          write(iout,115) nu_var_hgt, sqrt(min_var_hgt)*1000, 
     .                    sqrt(max_var_hgt)*1000, 
     .                    sqrt(av_var_hgt)*1000, use_ratio*dh_tol*1000,
     .                    descr(1:len_des)
 115      format(/,' For ',i4,' sites in origin, min/max dh/dt  sigma ',
     .           2F10.2,' mm/yr; Median  ',F10.2,' mm/yr, Tol ',
     .           F10.2,' mm/yr ',a)
      end if
      
      if( sqrt(max_var_hgt)-sqrt(av_var_hgt).gt.
     .    use_ratio*(dh_tol) ) then

*         Remove the sites that are not so good.
          do i = 1, num
             if( kbit(used,i) ) then

                 rad = sqrt(site_pos(1,1,i)**2 + site_pos(2,1,i)**2 +
     .                      site_pos(3,1,i)**2 )
                 var_hgt = 0.d0
                 do j = 1,3
                    jel = parn(j,1,i)
                    do k = 1,3
                       kel = parn(k,1,i) 
                       if( jel.gt.0 .and. kel.gt.0 ) then
                           var_hgt = var_hgt + site_pos(j,1,i)* 
     .                        cov_parm(jel,kel)*site_pos(k,1,i)/rad**2 
                       end if
                    end do
                 end do

*                Delete this site if the sigma is too large.
                 If( sqrt(var_hgt)-sqrt(av_var_hgt).gt.
     .               use_ratio*dh_tol ) then
                     if( type(1:1).eq.'P' ) then
                        write(iout,120) site_names(i),
     .                                  sqrt(var_hgt)*1000,
     .                                  use_ratio, descr(1:len_des) 
 120                    format(' Removing ',a8,' from orgin condition, '
     .                        ,'height sigma ', f10.2,' mm, Ratio Tol ',
     .                        F6.3,1x,a)
                     else
                        write(iout,125) site_names(i),
     .                                  sqrt(var_hgt)*1000,
     .                                  use_ratio, descr(1:len_des) 
 125                    format(' Removing ',a8,' from orgin condition, '
     .                        ,'dh/dt  sigma ', f10.2,
     .                         ' mm/yr, Ratio Tol ', F6.3,1x,a)
                     end if
                     
                     call sbit(used,i,0 ) 
                 end if
             end if
          end do 
      end if
                
 
*     Now start forming the normal equations
      do i = 1, num
*                                         ! We are using this site
          if( kbit(used,i) ) then
 
*             Get the partials for this site
              call form_cnd_parts( cnd_parts_bits, site_pos(1,1,i),
     .                             pmu_parts(1,1,i), parts)
 
              iel = parn(3,1,i)
 
*             Form covariance matrix so height has the assigned variance
*             relative to horizontal components
              call form_hgt_wgh( site_pos(1,1,i), cnd_hgt_var, hgt_wgh,
     .                           ss_var(i) )
 
*             Increment PTP
              do j = 1,nc
 
*                Increment solution
                 do l = 1,3
                    do m = 1,3
C	       write(*,997) j,l,m, dsol_cnd(j),
C    .  		  dsol_cnd(j)+
C    .  		     parts(m,j)*hgt_wgh(m,l)*sol_parm(iel+l-3)
C997		       format('DSOL ',3i4,2e18.7)
                       dsol = parts(m,j)*hgt_wgh(m,l)*sol_parm(iel+l-3)
                       dsol_cnd(j) = dsol_cnd(j) +
     .                         parts(m,j)*hgt_wgh(m,l)*sol_parm(iel+l-3)
c                      atf(j) =  dsol_cnd(j)
c                      dsol_cnd(j) = atf(j) + dsol
c                      atf(j) = dsol_cnd(j)
                    end do
                 end do
* MOD TAH 030126: Need to add this call for 
*                GNU Fortran (GCC 3.2.1) 3.2.1 20021119 (release).
*                dsoln_cnd incorrect with -O3 optimization without
*                call.
C MOD TAH 190226: Removed dummt call which does nothing
C                 call dumdsol(dsol_cnd, nc) 
                 do k = 1, nc
                    do l = 1,3
                       do m = 1,3
                           ptp(j,k) = ptp(j,k) +
     .                            parts(m,j)*hgt_wgh(m,l)*parts(l,k)
                       end do
                    end do
                 end do
              end do
              do l = 1,3
                 nr = nr + 1
                 do m = 1,3
                     prefit_sum = prefit_sum + sol_parm(iel+l-3)*
     .                               hgt_wgh(l,m)*sol_parm(iel+m-3)
                 end do
              end do
!             write(*,145) i, iel,prefit_sum,(sol_parm(iel+l-3), l=1,3)
!145          format('PREFIT ',2i6,E12.6, 3F12.6)
 
          end if
      end do
 
****  See if we have enough data (i.e., number of coordinates more
*     than number of conditions
      if( nr.lt.nc ) then
          write( iout, 150) nr, nc
 150      format('*** FAILURE in coorindate stabilization. Only ',
     .            i2, ' coordinates for ',i2,' conditions',/,
     .         '*** SYSTEM NOT stabilized')
          call report_stat('WARNING','GLORG','apply_cond_full',
     .                     ' ','Failure to converge',0)
          RETURN
      end if
 
****  Now invert the normal equations
      do i = 1, nc
         atf(i) = dsol_cnd(i)
      end do
* Debug: Check ptp:
C     do i = 1, nc
C        write(*,998) i,dsol_cnd(i), (ptp(i,j),j=1,nc)
C998     format('PTP ',i4,8E13.5)
C     enddo
      call invert_vis(ptp, dsol_cnd,scale, ipivot, nc,7,1)

* Debug: Check ptp: 
C     do i = 1, nc
C        write(*,999) i,dsol_cnd(i), (ptp(i,j),j=1,nc)
C999     format('PTI ',i4,8E13.5)
C     enddo
      do i = 1,nc
         postfit_sum = postfit_sum + atf(i)*dsol_cnd(i)
      end do
 
 
      postfit_chi = (prefit_sum-postfit_sum)/nr
      call write_cond_result( iout, type, cnd_parts_bits, dsol_cnd,
     .                    ptp, postfit_chi, num, site_names,used,
     .                    cond_var, ss_var, ss_it, descr(1:len_des))
 
*     Write out the solution and residuals
      sig_postfit = sqrt((prefit_sum-postfit_sum)/nr)
* MOD TAH 981003: Use the square of the rms (variance comparison,
*     old code used RMS)
      var_postedt = max((prefit_sum-postfit_sum)/nr,ss_min_rms**2)

      if( type(1:1).eq.'P' ) then
          WRite(iout,170) nr, type, ss_it, sqrt(prefit_sum/nr),
     .                sig_postfit, descr(1:len_des)
 170      format(' For ',i4,1x,a,' Iter ',i2,' Pre RMS ',f9.4,
     .           ' m; Post RMS ',f9.5,' m ',a)
      else
          WRite(iout,175) nr, type, ss_it, sqrt(prefit_sum/nr),
     .                 sig_postfit, descr(1:len_des)
 175      format(' For ',i4,1x,a,' Iter ',i2,' Pre RMS ',f9.4,
     .           ' m/yr; Post RMS ',f9.5,' m/yr ',a)
      end if
 
****  Now form the condition partials
      do i = 1,numc
        do j = 1,nc
              cond_part(j,i) = 0.d0
          end do
      end do
 
***   Now fill in the non-zero values
      do i = 1, num
          if( kbit(used,i) ) then
 
*            Form the partials
             call form_cnd_parts( cnd_parts_bits, site_pos(1,1,i),
     .                             pmu_parts(1,1,i), parts)
*             Form covariance matrix so height has the assigned variance
*             relative to horizontal components
              call form_hgt_wgh( site_pos(1,1,i), cnd_hgt_var, hgt_wgh,
     .                           ss_var(i) )
 
              do j = 1,nc
                  do k = 1,3
                      iel = parn(k,1,i)
                      if( iel.gt.0 ) then
****                      Save the last element used in cond_part
                          lel = iel
                          do l = 1,nc
                             do m = 1,3
                                cond_part(j,iel) = cond_part(j,iel) +
     .                                ptp(j,l)*parts(m,l)*hgt_wgh(m,k)
                             end do
                          end do
                     end if
                  end do
             end do

          end if
      end do

 
****  Now we have partials of conditions.  Now set the variances
*     of these linear combinations to zero.
*     First form Cov_parm*cond_part(R)
      do ir = 1,numc
          do ic = 1,nc
              cat(ir,ic) = 0.d0
 
              do i = 1,lel
                  cat(ir,ic) = cat(ir,ic) +
     .                cov_parm(ir,i)*cond_part(ic,i)
              end do
          end do

      end do
 
*                m   T
****  Now form AC   A
*                m+1
      do ir = 1,nc
          do ic = 1,nc
              ptp(ir,ic) = 0.d0
              do i = 1,lel
                  ptp(ir,ic) = ptp(ir,ic) +
     .                cond_part(ir,i)*cat(i,ic)
              end do
          end do
 
*         Now add in the constraint variance.  Normally
*         this will be zero
          ptp(ir,ir) = ptp(ir,ir) + cond_var(ir)
      end do
 
****  Now invert
      call invert_vis(ptp, origin_gain,scale, ipivot, nc,7,0)
 
****  Now finishup the origin gain
      do ir = 1, numc
          do ic = 1,nc
              origin_gain(ir,ic) = 0
              do i = 1, nc
                  origin_gain(ir,ic) = origin_gain(ir,ic) +
     .                cat(ir,i)*ptp(i,ic)
              end do
          end do
      end do
 
****  Compute the translation, rotation and scale
      do i = 1,nc
         cond_result(i) = 0.0d0
         do j = 1, lel
            cond_result(i) = cond_result(i) + cond_part(i,j)*
     .                       sol_parm(j)
         end do
      end do
****  Now finally update covarirance matrix

* MOD TAH 971105: Do direct calculation of corrections
*         Get the partials for this site
    
      if( kbit(options,19) ) then
         write(iout,220)
 220     format(/,'SDET: Details of origin stabalization',/,
     .            'DIRECT Calculation of position/velocity changes')         
         do i = 1, num
           if( kbit(used,i) ) then
             call form_cnd_parts( cnd_parts_bits, site_pos(1,1,i),
     .                            pmu_parts(1,1,i), parts)
             do j = 1,3
                ir = parn(j,1,i)
                if( ir.gt.0 ) then
                    call dwdot(tmp_xyz(j),parts(j,1),3, dsol_cnd(1),
     .                         1, nc)
                end if
                fin_xyz(j) = sol_parm(ir) - tmp_xyz(j)
                write(iout,240)  site_names(i), cmp_name(j), 
     .                       fin_xyz(j), tmp_xyz(j), 
     .                       (hgt_wgh(j,k),k=1,3)
 240            format(1x,a8,1x,a,2f10.4,' m; weights ',3f8.4)
            end do
           end if
         end do
         write(iout,250) 
 250     format('SDET: Computed position change through condition',
     .          ' application')        
         do i = 1, num
           if( kbit(used,i) ) then
             do j = 1,3
                ir = parn(j,1,i)
                if( ir.gt.0 ) then
                    call dwdot(tmp_xyz(j),origin_gain(ir,1),numc, 
     .                         cond_result(1), 1, nc)
                    fin_xyz(j) = sol_parm(ir) - tmp_xyz(j)
                end if
                write(iout,260)  site_names(i), cmp_name(j), 
     .                        fin_xyz(j), tmp_xyz(j), 
     .                        (hgt_wgh(j,k),k=1,3),ir
 260            format(1x,a8,1x,a,2f10.4,' m; weights ',3f8.4,1x,i5)
             end do
*            Compute contribution to chi**2

             prefit_sum = 0
             do j = 1,3             
                 do k = 1,3
                     prefit_sum = prefit_sum + fin_xyz(j)*
     .                               hgt_wgh(j,k)*fin_xyz(k)
                 end do
             end do
             write(iout,270) site_names(i), 
     .                       sqrt(prefit_sum/3)*1000.d0
 270         format(1x,a8,' Contribution to RMS ',F10.4, ' mm')            
           end if
         end do
      end if
      
* ENDMOD TAH
*     Check to see if we are only applying rotation with no
*     change to covariance matrix. MOD TAH 030223: 
*     MOD TAH 050818: Moved code to bottom.
      nact = nc
      do ir = 1,numc
          tmp = 0.d0

          do i = 1,nact
             tmp = tmp - origin_gain(ir,i)*cond_result(i)
          end do

*         Update solution vector
          sol_parm(ir) = sol_parm(ir) + tmp
 
*         Update covariance matrix
          do ic = 1,numc
              tmp= 0.d0
              do i = 1,nact
                  tmp = tmp + origin_gain(ir,i)*cat(ic,i)
              end do
*             Update covariance matrix 
              cov_parm(ir,ic) = cov_parm(ir,ic) - tmp

              if( ir.eq.ic .and. cov_parm(ir,ic).lt.0.d0 ) then

*                 Only report if it seems excessively negative
                  if( cov_parm(ir,ic).lt.-1.d-20 ) 
     .            write(iout,300) ir, cov_parm(ir,ic), tmp
 300              format('Negative diagonal Row ',i4,' After&before',
     .                  2d11.3)

*                 Set to absolute value so that we don't have negatives
*                 in matrix
                  cov_parm(ir,ic) = abs(cov_parm(ir,ic))
              end if
          end do
      end do


****  Now compute the new weights for the sites based on the 
*     average horizontal sigma.
      var_hor_all = 0.d0
      min_var_hgt = 1.d6
      max_var_hgt = 0.d0
* MOD TAH 031101: Initialize the variances of the sites used in the
*     stabilization. ns is number of stabization sites.  Kill the run
*     if too many stabilzation sites
      if( ns.gt. max_stab_site ) then
         call report_stat('FATAL','GLORG','Apply_cond_full',' ',
     .                    'Too many stabilization sites',ns)
      endif
      do i = 1, ns
         stab_var(i) = 0.d0
      end do

* MOD TAH 0310128: Get sigma distribution and edit sites with large sigmasa.
      ns = 0
      do i = 1, num
         if( kbit(used,i) ) then

             rad = sqrt(site_pos(1,1,i)**2 + site_pos(2,1,i)**2 +
     .                  site_pos(3,1,i)**2 )
             var_hgt = 0.d0
             var_tot = 0.d0

             err_tot = 0.d0
             err_hgt = 0.d0
             do j = 1,3
                jel = parn(j,1,i)
*               Increment the total error in the station postion
*               and the variance estimate.  (Note: for the height
*               we compute the actual error and then square it at
*               end.
                if( jel.gt.0 ) then
                   err_tot = err_tot + sol_parm(jel)**2
                   err_hgt = err_hgt + sol_parm(jel)*site_pos(j,1,i)/rad
                   var_tot = var_tot + cov_parm(jel,jel)
                end if
                do k = 1,3
                   kel = parn(k,1,i) 
                   if( jel.gt.0 .and. kel.gt.0 ) then
                       var_hgt = var_hgt + site_pos(j,1,i)* 
     .                       cov_parm(jel,kel)*site_pos(k,1,i)/rad**2 
                   end if
                end do
             end do

****         Get the horizontal sigma
             var_hor = var_tot - var_hgt
             if( var_hor.lt.min_var_hgt .and. var_hor.gt.0 ) 
     .                                    min_var_hgt = var_hor
             if( var_hor.gt.max_var_hgt ) max_var_hgt = var_hor

             ss_var(i) = var_hor
             ss_err(i) = (err_tot-err_hgt**2) + (err_hgt**2/cnd_hgt_var)
             var_hor_all = var_hor_all + var_hor
             ns = ns + 1
* MOD TAH 031101: Save the variance of the sites used in the stabilzation.
             stab_var(ns) = var_hor
          end if
      end do

* MOD TAH 031028: Added additional sigma checking on the results to delete sites that
*     do not match quality of others.
*     Check the sigmas based on the postfit variannces.  Then apply a sigma delete
*     feauture similar to that used for the heights when we started.  Algorithm here 
*     is based on horizontal sigmas
*     Find the median:
      call mdian2( stab_var, ns, av_var_hgt )

* MOD TAH 980514: Compute the minimum dh sigma spread.
      dh_tol = max(sqrt(av_var_hgt)-sqrt(min_var_hgt),ss_min_dne)
     
      if( type(1:1).eq.'P' ) then
          write(iout,410) ns, sqrt(min_var_hgt)*1000, 
     .                    sqrt(max_var_hgt)*1000, 
     .                    sqrt(av_var_hgt)*1000, 
     .                    use_ratio*dh_tol*1000,
     .                    descr(1:len_des)
 410      format(' For ',i4,' sites in origin, min/max NE sigma ',
     .           2F10.2,' mm; Median  ',F10.2,' mm, Tol ',F10.2,
     .                  ' mm ',a)
      else
          write(iout,415) ns, sqrt(min_var_hgt)*1000, 
     .                    sqrt(max_var_hgt)*1000, 
     .                    sqrt(av_var_hgt)*1000, 
     .                    use_ratio*dh_tol*1000,
     .                    descr(1:len_des)
 415      format(/,' For ',i4,' sites in origin, min/max dNE/dt sigma ',
     .           2F10.2,' mm/yr; Median  ',F10.2,' mm/yr, Tol ',
     .           F10.2,' mm/yr ',a)
      end if
      
      if( sqrt(max_var_hgt)-sqrt(av_var_hgt).gt.
     .    use_ratio*(dh_tol) ) then

          ns = 0
          nes = 0
          var_hor_all = 0

*         Remove the sites that are not so good.
          do i = 1, num
             if( kbit(used,i) ) then
                 
                 ns = ns + 1
                 var_hor = stab_var(ns)

*                Delete this site if the sigma is too large.
                 If( sqrt(var_hor)-sqrt(av_var_hgt).gt.
     .               use_ratio*dh_tol ) then
                     if( type(1:1).eq.'P' ) then
                        write(iout,420) site_names(i),
     .                                  sqrt(var_hor)*1000,
     .                                  use_ratio, descr(1:len_des) 
 420                    format(' Removing ',a8,' from orgin condition, '
     .                        ,'NE sigma ', f10.2,' mm, Ratio Tol ',
     .                        F6.3,1x,a)
                     else
                        write(iout,425) site_names(i),
     .                                  sqrt(var_hor)*1000,
     .                                  use_ratio, descr(1:len_des) 
 425                    format(' Removing ',a8,' from orgin condition, '
     .                        ,'dNE/dt  sigma ', f10.2,
     .                         ' mm/yr, Ratio Tol ', F6.3,1x,a)
                     end if
                     
                     call sbit(used,i,0 )
                 else
                     nes = nes + 1
                     var_hor_all =  var_hor_all + var_hor 
                 end if
             end if
          end do
          ns = nes
 
      end if
* END MOD TAH 031101:

****  Now normalize the variance so that the average variance is 1
      do i = 1, num
         if( kbit(used,i) ) then
             ss_var(i) = 1.d0*(1.d0-ss_rel**2) +
     .                   (ss_rel**2)*ss_var(i)/(var_hor_all/ns)

*            See if site violates the nsig condition.  ss_err is the
*            error squared.
             if( ss_err(i)/var_postedt/ss_var(i).gt.ss_nsig**2 ) then
                 call sbit(used,i,0)
                 if( type(1:1).eq.'P' ) then
                     write(iout,510) site_names(i), sqrt(ss_err(i)),
     .                      ss_var(i),
     .                      sqrt(ss_err(i)/var_postedt/ss_var(i))
 510                 format('Deleting ',a8,' Position error ',F8.4,' m,'
     .                     ,' relative variance ',F8.2,' Nsigma ',F8.2)
                 else
                     write(iout,520) site_names(i), sqrt(ss_err(i)), 
     .                      ss_var(i),
     .                      sqrt(ss_err(i)/var_postedt/ss_var(i))
 520                 format('Deleting ',a8,' Velocity error ',F8.4,
     .                      ' m/yr, relative variance ',F8.2,
     .                      ' Nsigma ',F8.2)
                end if
             end if
         end if
      end do

****  If we are just applying the rotation, go back now and add the 
*     translation back 
      if ( kbit(options,23) ) then
          write(*,'(a)') 'NAPP Set: Only applying rotation'
          do i = 1, num
*             Get the partials for this site
              call form_cnd_parts( cnd_parts_bits, site_pos(1,1,i),
     .          		   pmu_parts(1,1,i), parts)

*             Loop over the 3-translations
              do j = 1,3
                 iel = parn(j,1,i)
                 if( iel.gt.0 ) then
                     sol_parm(iel) = sol_parm(iel) + cond_result(j+3) 
                     if( nc.eq.7 ) then
                	sol_parm(iel) = sol_parm(iel) + 
     .          	       cond_result(7)*parts(j,7)
                     endif
                 end if
              end do
          end do
      end if 
                  
                
 
****  Thats all.
      return
      ENd
 
CTITLE FORM_CND_PARTS
 
      subroutine form_cnd_parts( cnd_parts_bits, site_pos, pmu_parts,
     .                           parts)

      implicit none 
 
*     This routine will form up the partials derivatives of the condition
*     variables with respect to station position adjustments.  cnd_parts_bits
*     sets which partials will be formed up o a maximum of seven (3
*     translations, 3 rotations and a scale).
 
* PASSED VARIABLES
 
*   cnd_parts_bits - Bit mapped word which sets which partials are to be
*                    estimated.  Bits 1-3 -- rotations, bits 4-6 translations
*                    bit 7 - scale
 
 
      integer*4 cnd_parts_bits
 
*   pmu_parts(3,3)  - Polar motion/UT1 partials.  Note: the rows
*                     are for site position componenst XYZ and
*                     columns are for x,y and UT1. (for this site only)
*   site_pos(3)     - Apriori site postions (XYZ m)
*   parts(3,7)      - Partials of rotation, translation, scale
*                     with respect to site positions.
 
 
      real*8 pmu_parts(3,3), site_pos(3), parts(3,7)
 
* LOCAL VARIABLES
 
*   j,k     - Loop counters
*   kbit    - Tests if bit is set
*   np      - Current partial number
 
 
      integer*4 j,k, np
      logical kbit 
 
*     Form the partials first.  Form X, then Y, then Z partials
 
      do j = 1,3
 
          np = 0
*         Rotations. See if these are to be used in the conditiion
          do k = 1,3
*            Note: we use transpose:
             if( kbit( cnd_parts_bits,k ) ) then
                 np = np + 1
                 parts(j,np) = pmu_parts(k,j)
             end if
          end do
 
*         Translations.  See if these are being used
          do k = 1,3
             if( kbit( cnd_parts_bits, k+3) ) then
                 np = np + 1
*                If the translation parameter matches the coordinate
*                component then set to 1 else set the partial to be
*                zero.
                 if( k.eq.j ) then
                     parts(j,np) = 1.d0
                 else
                     parts(j,np) = 0.d0
                 end if
             end if
          end do
 
*         scale.  See if this is to estimated
          if( kbit( cnd_parts_bits,7) ) then
              np = np + 1
              parts(j,np) = site_pos(j)*1.d-9
          END If
      end do
 
***** Thats all
      return
      end
 
 
CTITTLE WRITE_COND_RESULT
 
      subroutine write_cond_result( iout, syst, cnd_parts_bits, results,
     .                cov_result, postfit_chi, num, site_names, used,
     .                cond_var, ss_var, ss_it, descr )

      implicit none 
 
*     Routine to write out the condition results.

 
* PASSED VARIABLES
 
*   iout           - Output unit number
*   cnd_parts_bits - Bit mapped word which sets which partials are to be
*                    estimated.  Bits 1-3 -- rotations, bits 4-6 translations
*                    bit 7 - scale
*   num            - Number of sites
*   used(num/32+1) - Bit mapped word saying which sites have been used.
*   ss_it          - Iteration number

      integer*4 iout, cnd_parts_bits, num, used(num/32+1), ss_it
 
*   results(7)      - Results for rotation, translation, scale parameter
*                     estimates
*   cov_result(7,7) - Covaraince matrix of results estimate
*   postfit_chi        - Postfit chI**2 of the result
*   cond_var(7)     - Variances applied to the constrraints
*   ss_var(num)     - Relative variances applied to each site.
 
 
      real*8 results(7), cov_result(7,7), postfit_chi, cond_var(7),
     .       ss_var(num)

*   syst   - Type of conditition appplied (Position or Rate)
*   site_names(num)  - NAmes of sites
*   descr  - Descriptor for solution

      character*(*) syst, site_names(num), descr
 
* LOCAL VARIABLES
 
*   i       - Loop counters
*   kbit    - Tests if bit is set
*   np      - parameter number
 
 
      integer*4 i,  np, indx
      logical kbit
 
*   type(7)   - Label for parameter estimate
 
      character*18 typp(7), typr(7)

*   line      - Line with list of stations used.

      character*132 line
 
      data typp / 'X Rotation  (mas)',
     .            'Y Rotation  (mas)',
     .            'Z Rotation  (mas)',
     .            'X Translation (m)',
     .            'Y Translation (m)',
     .            'Z Translation (m)',
     .            'Scale       (ppb)' /
 
      data typr / 'X Rotate (mas/yr)',
     .            'Y Rotate (mas/yr)',
     .            'Z Rotate (mas/yr)',
     .            'X Trans    (m/yr)',
     .            'Y Trans    (m/yr)',
     .            'Z Trans    (m/yr)',
     .            'Scale    (ppb/yr)' /
 
***** Loop over all the parameters, printing thpose that were estimated
 
      write(iout, 100) syst
 100  format(/,' ',a,' system stabilization results ',/,
     .       ' --------------------------------------- ')
      np = 0
      do i = 1,7
          if( kbit( cnd_parts_bits,i) ) then
              np = np + 1
              if( syst(1:1).eq.'P' ) then
                  write(iout, 110) typp(i), results(np),
     .                sqrt(cov_result(np,np)*postfit_chi), ss_it,
     .                descr
              else
                  write(iout, 110) typr(i), results(np),
     .                sqrt(cov_result(np,np)*postfit_chi), ss_it,
     .                descr
              end if
110           format(1x, a,1x,F10.5, ' +- ', F10.5,' Iter ',i2,1x,a)
          end if
      end do
      
      write(iout,120) (sqrt(cond_var(i)),i=1,np)
 120  format(' Condition Sigmas used ',7F10.4)

      indx = 1
      write(iout,200)
200   format('Sites and relative sigmas used in stabilization')
      do i = 1,num
	 if( kbit(used,i) ) then
C	     line(indx:) = site_names(i)
C	     indx = indx + 10
             write(line(indx:),205) site_names(i), sqrt(ss_var(i))
 205         format(a8,1x,F6.2)
             indx = indx + 17
	     if( indx.gt.100 ) then
		 write(iout,210) line(1:indx-2)
210              format(a)
		 indx = 1
             end if
         end if
      end do
      if( indx.ne.1 ) write(iout,210) line(1:indx-2) 

 
***** Thats all
      return
      end
 
CTITLE FORM_HGT_WGH
 
      subroutine form_hgt_wgh( site_pos, cnd_hgt_var, hgt_wgh, ss_var )

      implicit none 
 
*     Routine to form the weight matrix for the site coordinates.
*     The height covariance matrix is first formed and then inverted to give
*     the weight matrix.
*
* PASSED VARIABLES
 
*   site_pos(3) - site position XYZ (m)
*   cnd_hgt_var - The variance to be given to the heights
*   hgt_wgh(3,3)    - XYZ covariance matrix which is inverted to give a
*               - weight matrix
*   ss_var      - Varinace to be given to site (normalized over all sites to
*                 unity.
 
      real*8 site_pos(3), cnd_hgt_var, hgt_wgh(3,3), ss_var
 
* LOCAL PARAMETERS
 
*   j,k     - Loop counters.
*   ipivot(3)   - Pivot elements during inversion
 
 
      integer*4 k, l, ipivot(3)
 
*   dum(3)      - Dummy solution vector for invert_vis
*   scale(3)        - scale vector for invert_vis
*   rot_mat(3,3)  - Rotation matrix from XYZ to NEU (initially,
*               - transposed to yield the desized matrix
*   loc_coord(3)    - Local coordinates of site
*   hgt_neu(3,3)    - Covarinace matrix in NEU.
*   temp_cov(3,3)   - Temp matrix used in transfromation from XYZ to
*               - NEU.
 
 
      real*8 dum(3), scale(3), rot_mat(3,3), loc_coord(3),
     .    hgt_neu(3,3), temp_cov(3,3)
 
*     Get the transformation from XYZ to NEU
      call xyz_to_neu( rot_mat, site_pos, loc_coord)
 
*     Transpose the matrix
      do k = 1,2
      call dvswp(rot_mat(k,k+1),3,
     .        rot_mat(k+1,k),1, 3-k)
      end do
 
 
*     Set up the diagonal NEU matrix
      do k = 1,3
         do l = 1,3
            hgt_neu(k,l) = 0.d0
         end do
*        If horizontal components then set weight such
*        that average weighted per coordinate is 3.
         if ( k.le.2 ) then
             hgt_neu(k,k) = (3.d0-1.d0/cnd_hgt_var)/2.d0*ss_var
         else
             hgt_neu(k,k) = cnd_hgt_var*ss_var
         end if
      end do
 
      call var_comp(rot_mat, hgt_neu, hgt_wgh,
     .                        temp_cov, 3,3,1 )
 
      call invert_vis(hgt_wgh, dum,scale, ipivot, 3,3,0)
 
***** Thats all
      return
      end

CTITLE DUMDSOL

      subroutine dumdsol(dsol_cnd, nc)

      implicit none 

*     Dummy routine to fix bug in G77 
*     GNU Fortran (GCC 3.2.1) 3.2.1 20021119 (release)
      integer*4 nc
      real*8 dsol_cnd(nc)

      integer*4 i

      do i = 1, nc
      end do

      return
      end






