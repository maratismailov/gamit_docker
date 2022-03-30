      subroutine res_ambRT(ep, OK)

      implicit none

*     Subroutine to use all available data to see if ambiguities can
*     be resolved.

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRT.h'         ! Common block
      include 'trackRTObs.h'      ! Data block

* PASSED
      integer*4 ep   ! Epoch counter
      logical   OK   ! Set false if we resolve ambiquities so that
                     ! epoch will be iterated (not needed but used
                     ! as code check),


* LOCAL
      integer*4 i,j,k   ! Loop counters

      real*8 ests(5), sig(5)  ! Values and sigmas of 5 combintions of L1
                      ! and L2 ambiquities used to test:
                      ! 1 -- MW WL  1 -1
                      ! 2 -- EX WL  1 f1/f2 
                      ! 3 -- LC Est lcf1 lcf2
                      ! 4 -- L1 Est
                      ! 5 -- L2 Est
     .,      resid(5) ! Differene between est and computed value for best L1/L2
     .,      comb_coef(2,5)   ! Coefficients
     .,      comb_wgh(5)      !  Combination weights
 
      real*8 cov_col(max_parm)   !  The column of the
                  !  covariance matrix for the parameters being forced.
     .,      sol_col  !  The change in the parameter
                      ! estimate needed to get the forced values.
     .,      dsol     ! New value for the solution estimate
     .,      var_force  !  Variance with which value should be forced.
     .,      equ_gn(max_parm)  ! Kalman gain
     .,      dchi     ! Chi**2 change
 
      logical use_comb(5)   ! Set true to use a combination
      logical still_resolving  ! Set true while are are still resolving
                            ! ambiguities.
      logical amb_OK   ! Set true if ambiguity can be resolved to integer

      integer*4 iNL1, idNL1, idNL2  ! L1 and dL2 values to start search
     .,       L12_min(2), L12_nxt(2) ! Integer offsets of min and nexy chi
     .,       iL1, IL2       ! Test values of ambiquities
     .,       iof    ! Offset for parameter estimates if L1/L2 estimated
     .,       np     ! Parameter number
     .,       pass   ! Pass counter for resolving ambiguities
     .,       jel    ! Function to return baseline number

      real*8 chi2    ! Chi**2 increment
     .,      chi2_min, chi2_nxt  ! Smallest and next chi**2 increments
     .,      chi2_cont(5)  ! Contributions from each chi type
     .,      dest(5)       ! Change in estimate for iL1/iL2 values.
     .,      rel_rank  ! Relative rank of best choice
     .,      bl        ! baselen (distance) to reference site
     .,      rms       ! RMS scatter estimate
     .,      elrad     ! Elevation in (rads)

      character*8 Fcode ! F-code for not fixing ambiguties
      character*2 chi2_types(5)  ! Data types (used for output)

      data comb_coef  /  1.d0, -1.d0,  
     .                   1.d0,  exf2,
     .                   lcf1,  lcf2,
     .                   1.d0,  0.d0,
     .                   0.d0,  1.d0  /

      data comb_wgh / 1.d0, 1.d0, 1.d0, 1.d0, 1.d0 /

      data chi2_types / 'MW','EX','LC','L1','L2' /



      idNL2 = 1
      idNL1 = 1

      comb_wgh(1) = wl_fact
      comb_wgh(2) = lg_fact


      still_resolving = .true.
      pass = 0
      do while ( still_resolving )
         still_resolving = .false.
         pass = pass + 1
      
         do i = 1, num_ambs

****        See if we have an estimated parameter (ie. not resolved yet)
            if( amb_parn(1,i).gt.0 ) then
*               OK: value not resolved
                if( wls_num(i).gt.0 ) then
*                  Compute the estimate and "formal" sigma
                   do j = 1,2 
                      ests(j) = WLS_sum(j,i)/WLS_num(i)
                      rms = sqrt(abs(WLS_sqr(j,i) 
     .                          -ests(j)**2*WLS_num(i))/
     .                               max(1,WLS_num(i)-1))
                      sig(j) = rms/sqrt(min(1.*WLS_num(i),1.*WL_avnum))
                      use_comb(j) = .true.
                   end do
*                  See if we need to set min MW-WL sigma
                   if( sig(1).lt. mwwl_minsig ) sig(1) = mwwl_minsig
*                  Modify the sigma on the EW-WL to account for 
*                  systematics from the ionosphere.
                   bl = baselens(jel(1,bf_ents(1,i)))
                   if( wls_obn(1,i).gt.0 ) then
                      elrad = RT_azel(2,wls_obn(1,i))*pi/180
                      sig(2) = max(sig(2),
     .                      sqrt(exwl_minsig**2+(bl*exwl_scale)**2))*
     .                      (1+exwl_elev/sin(elrad))
                   else   ! If no data assume 
                      sig(2) = max(sig(2),
     .                      sqrt(exwl_minsig**2+(bl*exwl_scale)**2))*
     .                      (1+exwl_elev*6.0)   ! ~1/sin(10)
                   endif 

*                  If just two points make sure we do not use
                   if( WLS_num(i).le. 2 ) then
                       sig(1) = 2*float_limit(2)
                       sig(2) = 2*float_limit(2)
                   end if
                   idNL2 = nint(ests(1))
                else
                   use_comb(1) = .false.
                   use_comb(2) = .false.
                   idNL2 = 0
                end if

****            Now get the estimate
                if( l1_only ) then
                    use_comb(3) = .false.
                    iof = 1
                else
                    use_comb(4) = .false.
                    use_comb(5) = .false.
                    iof = 0
                end if
                if( neam.eq.2 ) iof = 1
                do j = 1,neam  
                   np = amb_parn(j,i) 
                   ests(2+j+iof) = sol_vecm(np)
                   sig(2+j+iof)  = sqrt(cov_parmm(np,np))
*                  Make sure sigma does not get too small
*                  which will result in high chi**2.
                   if( sig(2+j+iof).lt.min_ambsig ) then
                       sig(2+j+iof) = min_ambsig
                   end if
                   use_comb(2+j+iof) = .true.
                   if( j.eq.1 ) then
                       iNL1 = nint(ests(2+j+iof))
                       idNL1 = int(sig(2+j+iof)*2)+1
                   end if
                end do

****            OK: Now we have estimates, compute the 
*               relative rank values
                chi2_min = 1.d6
                chi2_nxt = 1.d6
                do iL1 =  iNL1-idNL1, iNL1+idNL1
                   do iL2 = iL1-idNL2-2, iL1-idNL2+2
*                     Compute the chi**2 for these values
                      chi2 = 0.0
                      do j = 1,5
                         if( use_comb(j) ) then
                            dest(j) = comb_coef(1,j)*iL1  
     .                           + comb_coef(2,j)*iL2
                            chi2 = chi2 
     .                           + comb_wgh(j)
     .                            *((ests(j)-dest(j))/sig(j))**2
                         endif
                      end do
                      if( chi2.lt.chi2_min ) then
*                         Shift the current minimum to the next
*                         entry since we are about to update with
*                         a new minumum value
                          chi2_nxt  = chi2_min
                          L12_nxt(1)= l12_min(1)
                          L12_nxt(2)= L12_min(2)
*                         Now say the new minimum entries
                          chi2_min = chi2
                          l12_min(1) = iL1
                          L12_min(2) = iL2
*                         Save contributions for output
                          do j = 1,5
                             if( use_comb(j) ) then
                                chi2_cont(j) = comb_wgh(j)
     .                            *((ests(j)-dest(j))/sig(j))**2
                                resid(j) = ests(j)-dest(j)
                             else
                                chi2_cont(j) = 0.d0
                                resid(j) = 0   ! Not being used
                             end if
                             asv_res(j,i) = resid(j)
                             asv_sig(j,i) = sig(j)
                             asv_chi(j,i) = chi2_cont(j)
                             asv_used(j,i) = use_comb(j)

                         end do
                      elseif( chi2.lt.chi2_nxt ) then
                          chi2_nxt = chi2
                          l12_nxt(1) = iL1
                          L12_nxt(2) = iL2
                      endif
                   end do
                end do

*               See how the combinations look
                rel_rank = chi2_nxt/(chi2_min+0.25)

****            See if sigmas are OK based on data being used and
*               whether we can resolve the ambiguity
                Fcode = '------'
                amb_OK = .true.
                if( use_comb(1) .and. (sig(1).gt. float_limit(2) .or.
     .              WLS_num(i).lt.wl_mnnum) ) then 
                    amb_OK = .false.
                    Fcode(2:2) = 'W'
                end if
*               Now test LC, L1 and LC against sigma
                if( use_comb(3) .and. sig(3).gt. float_limit(1) ) then
                    amb_OK = .false.
                    Fcode(3:3) = 'S'
                end if
                if( use_comb(4) .and. sig(4).gt. float_limit(1) ) then
                    amb_OK = .false.
                    Fcode(4:4) = 'S'
                end if
                if( use_comb(5) .and. sig(5).gt. float_limit(1) ) then
                    amb_OK = .false.
                    Fcode(5:5) = 'S'
                end if
                if( chi2_min.gt.max_fit ) then
                    Fcode(6:6) = 'C'
                    amb_OK = .false.
                end if
*               Test relative rank
                if( rel_rank.lt. relrank_limit ) then
                    amb_OK = .false.
                    Fcode(1:1) = 'R'
                end if

****            See if we can fix this bias
c               if( .not.amb_OK) then
c                   call ambinf_rep(6, i, ep, pass, Fcode, rel_rank,
c    .                    chi2_min, chi2_nxt, L12_min, L12_nxt, 
c    .                    chi2_cont, chi2_types, use_comb, .true. )
c               end if

*               Save information for status reports
                asv_rbn(1,i) = rel_rank
                asv_rbn(2,i) = chi2_min
                asv_rbn(3,i) = chi2_nxt
                asv_dL12(1,i) = l12_min(1)
                asv_dL12(2,i) = l12_min(2)
                asv_fcode(i) = Fcode
 
                if( amb_OK ) then
*                   This ambiguity can be fixed
                    still_resolving = .true.
                    ambiq_all(1,i) = ambiq_all(1,i) + L12_min(1)
                    ambiq_all(2,i) = ambiq_all(2,i) + L12_min(2)

                    call update_wlsstat(ep, L12_min(1),L12_min(2), i)

****                Mark the bias flag as fixed
                    call sbit(bf_ents(5,i),2,1)   ! Show resolved
                    call sbit(bf_ents(5,i),3,0)   ! Show not a jump anymore

                    asv_resep(i) = ep

                    if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .              write(*,320) ep, pass, i, site_names(bf_ents(1,i)),
     .                  bf_ents(2,i), rel_rank, Fcode, chi2_min,
     .                  l12_min, chi2_nxt, l12_nxt, iNL1, idNL2
 320                format('BF_FIX Ep ',i6,1x,I2,' Amb ',i3,1x,a,
     .                  ' PRN ',i2.2,' RR ', 
     .                   F10.3,' Fcode ',a,' Min ',F10.3,1x,2I4, 
     .                  ' Nxt ',F10.3, 1x,2I4,' Rnge ',2i4)
                    call ambinf_rep(lus, i, ep, pass, Fcode, rel_rank,
     .                    chi2_min, chi2_nxt, L12_min, L12_nxt, 
     .                    chi2_cont, chi2_types, use_comb, 
     .                    resid, sig, .true. )
                    call ambinf_rep(6, i, ep, pass, Fcode, rel_rank,
     .                    chi2_min, chi2_nxt, L12_min, L12_nxt, 
     .                    chi2_cont, chi2_types, use_comb, 
     .                    resid, sig, .true. )

****                Now force the parameter values to the ones we need
                    do j = 3, 5
                       if( use_comb(j) ) then ! This combination used
                                              ! (LC,L1 or L2)
*                         Compute the value we need to force to
                          dest(j) = comb_coef(1,j)*L12_min(1) 
     .                         + comb_coef(2,j)*L12_min(2)
*                         Get the parameter number
                          if( j.le.4 ) then
                              np = amb_parn(1,i)   ! LC or L1
*                             Set negative value here which is used 
*                             in update_posrt to remove the parameter.
                              amb_parn(1,i) = -np   
                          else if( j.eq.5 ) then
                              np = amb_parn(2,i)   ! L2
                              amb_parn(2,i) = -np    ! Later to neg to remove
c                             amb_parn(2,i) = 0    ! Later to neg to remove
                          end if
                          var_force = 1.d-9

                          call force_trackRT(cov_parmm, sol_vecm, 
     .                         max_parm, num_parm, cov_col, sol_col, 
     .                         np, dest(j), var_force, dchi, equ_gn)

                       end if
                    end do
                else
                    asv_resep(i) = 0
                    if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .              write(*,420) ep, pass, i, site_names(bf_ents(1,i)),
     .                  bf_ents(2,i), rel_rank, Fcode, chi2_min,
     .                  l12_min, chi2_nxt, l12_nxt, iNL1, idNL2
 420                format('BF_FRE Ep ',i6,1x,I2,' Amb ',i3,1x,a,
     .                  ' PRN ',i2.2,' RR ', 
     .                   F10.3,' Fcode ',a,' Min ',F10.3,1x,2I4, 
     .                  ' Nxt ',F10.3, 1x,2I4,' Rnge ',2i4)

                    if( ep.ge.debug(3) .and. ep.le.debug(4) )
     .              call ambinf_rep(6, i, ep, pass, Fcode, rel_rank,
     .                    chi2_min, chi2_nxt, L12_min, L12_nxt, 
     .                    chi2_cont, chi2_types, use_comb, 
     .                    resid, sig, .false. )

C                   If we have large chi**2
                    if( chi2_min.gt.max_fit*2 ) then
                       if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .                 write(*,440) ep, pass, i, 
     .                     site_names(bf_ents(1,i)),bf_ents(2,i),
     .                     (use_comb(j), chi2_cont(j),j=1,5)
 440                   format('BF_CHI Ep ',i6,1x,I2,' Amb ',i3,1x,a,
     .                        ' PRN ',i2.2,' Chi^2 contributions ',
     .                       5(L2,1x,F12.2,1x))
                    end if


                end if          ! Bias fixed
*           If we have data still, update the asv_res and sig numbers
            else
                if( wls_num(i).gt.0 ) then
*                  Compute the estimate and "formal" sigma
                   do j = 1,2 
                      ests(j) = WLS_sum(j,i)/WLS_num(i)
                      rms = sqrt(abs(WLS_sqr(j,i) 
     .                          -ests(j)**2*WLS_num(i))/
     .                               max(1,WLS_num(i)-1))
                      sig(j) = rms/sqrt(min(1.*WLS_num(i),1.*WL_avnum))
                      use_comb(j) = .true.
                   end do
*                  See if we need to set min MW-WL sigma
                   if( sig(1).lt. mwwl_minsig ) sig(1) = mwwl_minsig
*                  Modify the sigma on the EW-WL to account for 
*                  systematics from the ionosphere.
                   bl = baselens(jel(1,bf_ents(1,i)))
                   if( wls_obn(1,i).gt.0 ) then
                      elrad = RT_azel(2,wls_obn(1,i))*pi/180
                      sig(2) = max(sig(2),
     .                      sqrt(exwl_minsig**2+(bl*exwl_scale)**2))*
     .                      (1+exwl_elev/sin(elrad))
                   else   ! If no data assume 
                      sig(2) = max(sig(2),
     .                      sqrt(exwl_minsig**2+(bl*exwl_scale)**2))*
     .                      (1+exwl_elev*6.0)   ! ~1/sin(10)
                   endif 
*                  Now just update the MW and EX-WL values
                   do j = 1,2
                      if( use_comb(j) ) then
                         chi2_cont(j) = comb_wgh(j)
     .                     *((ests(j))/sig(j))**2    ! dEstimate is zero
                         resid(j) = ests(j)
                      else
                         chi2_cont(j) = 0.d0
                         resid(j) = 0   ! Not being used
                      end if
                      asv_res(j,i) = resid(j)
                      asv_sig(j,i) = sig(j)
                      asv_chi(j,i) = chi2_cont(j)
                      asv_used(j,i) = use_comb(j)
                   end do
                end if
            end if              ! Parameter number > 0
         end do
      end do

      if( ep.ge.debug(3) .and. ep.le.debug(4) ) then
c         Output the total ambiguity value (for comparison with track)
          do i = 1,num_ambs
             write(*,500) i, ep, bf_ents(1:6,i),ambiq_all(:,i)
 500         format('BF_ENTS # ',i4,' EP ',i6,' Tab ',2I3,1x,2I7,1x,
     .              ' S ',o3,1x,I7,1x,' AL12 ',2F20.2)
          enddo
      endif 

****  See if we resolved ambiquities (remove later)
c      if( pass.gt.1 ) then
c         OK = .false.
c      else
c         OK = .true.
c      end if

****  That all 
      return
      end

CTILE FORCE_TRACK

      subroutine force_trackRT( cov_parm, sol_parm, dim_parm, num_parm,
     .    cov_col, sol_col, nc, force_value, force_var, dchi, equ_gn)

      implicit none
 
*     This routine will force parameters in a covariance matrix to have
*     specific values.  The parameters to be forced are given in the
*     NC array.  There are num_force of them and the values to be forced
*     to are given in the force_values array.  If local array mapping
*     is to be then map_force should be called before this routine
*     setup the mapping of the arrays needed.
 
*   dum_parm    - Dimensions of cov_parm
*   num_parm        - Number of parameters in cov_parm.  Cov_parm
*               - is assumed to be dimensioned to this size.
*   num_force   - Number of parameters to be focesd to specific
*               - values
*   nc(num_force)   - The list of parameters to be made equal.
*   ipivot(num_force)   - The pivot elements for the matrix
*               - inversion
 
      integer*4 dim_parm, num_parm,  nc
 
*   cov_parm(num_parm,num_parm) - Full covariance matrix to
*               - be modified.
*   sol_parm(num_parm)      - The solution vector.  The nc
*               - elements of this vector will be set equal.
*   cov_col(num_parm, num_force)    - The columns of the
*               -  covariance matrix for the parameters
*               - being forced.
*   sol_col(num_force)      - The change in the parameter
*               - estimates needed to get the forced values.
*   force_values(num_force) - The values the forced parameters
*               - should take on.
*   force_var(num_force) - Variance with which value should be forced.
*   avat(num_force,num_force)   - The multiplication of
*               - the partials matrix of (00001000,000001000)
*               - and the covaraiance matrix.
*   equ_gn(num_parm,num_force)  - Kalman gain matrix for
*               - getting solution.
*   scale(num_force)        - Scaling vector for solution. (passed
*               - to invert_vis)
*   dchi        - Change in Chi**2
 
      real*8 cov_parm(dim_parm,dim_parm), sol_parm(dim_parm),
     .    cov_col(dim_parm), sol_col,
     .    force_value, force_var, avat, equ_gn(dim_parm), dchi
 
* LOCAL PARAMETERS
 
*   i,j,k       - Loop counters
 
      integer*4 i,j,k
 
*   dsol, dcov  - Summation variables for computing corrections
*               - to the solution vector and covariance matrix.
 
 
      real*8 dsol, dcov
 
*                        T
***** Start, form the AVA  matrix, where A is of the form:
*     y = Ax (y is vector of zero observations)
*     and
*     A =  0 0 0 1 0 0 0 0 0 0....... to num_parm
*          0 0 0 0 0 0 0 1 0 0
*          0 0 0 0 0 0 0 0 1 0
*     where the above form would set parameters 4,8, and 9.
*     For num_force parameters being equated, there are num_force
*     rows in A.
*
*     The above form is much simpler to compute is we save the
*     columns of the covariance matrix for thhose parameters to
*     be forced first.
 
      do i = 1, num_parm
          cov_col(i) = cov_parm(i,nc)
      end do
 
*     Save the forced parts of the solution vector as well
      sol_col = sol_parm(nc)
 
*               T
****  Now do AVA
*
      avat = cov_col(nc) + force_var
* MOD TAH 000302: Added variance to diagonal for constrained 
*             forcing.              

****  Now invert this matrix.  (If "sort of equal" was desired we could
*     add value to diagonal now representing variance of y above). Pass
*     zero as number in solution vector, we dont want to multiply.
*     kgain below is dummy argument.

      avat = 1/avat

****  Before continuing compute the change in Chi**2 due to condition
      dchi = (force_value-sol_col)**2*avat
 
*     Now form the Kalman gain, equ_gn given by
*                T     T -1
*     equ_gn = VA  (AVA )
*
      do i = 1, num_parm
*             Do the multiply (could use VIS but stick to straight)
*             call dwmul(equ_gn(i,j), col_col(i,1), num_force,
*    .                    avat(1,j),1, num_force)
          equ_gn(i) = cov_col(i)*avat
      end do
 
****  Now get the change to the solution vector
*
*     x  = x  -  equ_gn*(force_value-sol_col)
*      n    o
*
      do i = 1,num_parm
          dsol = equ_gn(i)*(force_value-sol_col)
          sol_parm(i) = sol_parm(i) + dsol
      end do
 
*     Now update the covariance matrix
*
*     V  = V - equ_gn*cov_col
*      n    o
 
      do i = 1, num_parm
          do j = 1, num_parm
 
*             Do summation loop
              dcov = equ_gn(i)*cov_col(j)
              cov_parm(i,j) = cov_parm(i,j) - dcov
          end do
      end do
 
****  Thats all
      return
      end

CTITLE AMBINF_REP

      subroutine ambinf_rep( unit, na, ep, pass, Fcode, rel_rank,
     .     chi2_min, chi2_nxt, L12_min, L12_nxt, 
     .     chi2_cont, chi2_types, use_comb, resid, sig, fixed )

      implicit none

*     Routine to report information about ambiquities.

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRTObs.h'      ! Data definition
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 unit   ! Output unit number
     .,  na   ! Ambiguity being reported
     .,  ep   ! Epoch number when reported 
     .,  pass ! Iteration pass
     .,  L12_min(2), L12_nxt(2)  ! DL1/2 values for best and next
     .,  j    ! Loop counter
     .,  no   ! Obs number

      real*8 rel_rank  ! Relative rank
     .,  chi2_min, chi2_nxt  ! Chi2 values
     .,  chi2_cont(5)  ! Contributions from each type
     .,  mean(2), rms(2)  ! Mean and RMS of WLs
     .,  eNl1, eNL2    ! Float point error
     .,  resid(5), sig(5)  ! Estimates and sigmas for the 5 possible
                       ! observables

     .,  uresid(5), usig(5), uchi2(5) ! Used values
     .,  out_ae(2)      ! Output azimuth/elevation (when no data) 

      integer*4 nu          ! Num type used.

      integer*4 iNL1, iNL2   ! Integer values

      logical use_comb(5)  ! Set true when type used
     .,  fixed         ! Set true if ambiqguity resolved

      character*(*) Fcode  ! Code for bias fixing
     .,   chi2_types(5)    ! 2-charater codes for type 

      character*2 utypes(5) ! Used types

      character*6 restyp    ! String with resolve type AMBFIX AMBFRE


          
****  Output information about this bias
      if( fixed ) then
          restyp = 'AMBFIX'
      else
          restyp = 'AMBFRE'
      end if

      write(unit,125) restyp, site_names(bf_ents(1,na)),bf_ents(2,na),
     .    ep, rel_rank, Fcode, L12_min, chi2_min, chi2_nxt, na,
     .    ambiq_all(:,na)
 125  format(a,1x,a4,1x,'PRN',i2.2,' EP ',I6,' RelRank ',F8.2,1x,
     .    ' FC ',a,' dL12 ',2i3,' Dchi ',2F7.2,' AMB# ',i3,
     .    ' AL12 ',2F14.1)


      if( wls_num(na).gt.0 ) then
         do j = 1,2
            mean(j) = WLS_sum(j,na)/WLS_num(na)
            rms(j) = sqrt(abs(WLS_sqr(j,na)-mean(j)**2*WLS_num(na))/
     .           max(1,WLS_num(na)-1))
         end do
         call res_wlsRT( mean(1), mean(2), eNL1, eNL2, iNl1, iNL2)
         no = wls_obn(1,na)
         if( no.gt.0 ) then
             out_ae(:) = RT_azel(:,no)
         else
             out_ae(1) = 999.99
             out_ae(2) = 99.99
         end if
         write(unit,225) site_names(bf_ents(1,na)),bf_ents(2,na),
     .     ep, bf_ents(3:5,na), iNl1, iNL2, mean, rms,WLS_num(na),
     .     eNL1, eNL2,RT_azel(:,no)
 225     format('AMBWLS',1x,a4,1x,'PRN',i2.2,' EP ',I6,
     .          ' RG ',2I7,' FX ',o3,6x,'iL12 ',2i3,
     .          ' Means ',2F7.2,' RMS ',2F6.2,' # ',I6,
     .          ' eN12 ',2F8.2,' AzEl ',2F6.2,1x,2F14.1)

      end if

****  See what we used and output values, sigmas and chi cont
      nu = 0
      do j = 1,5
         if( use_comb(j) ) then
             nu = nu + 1
             utypes(nu) = chi2_types(j)
             uresid(nu)  = resid(j)
             usig(nu)   = sig(j)
             uchi2(nu)  = chi2_cont(j)
         endif
      end do
      write(unit,325) site_names(bf_ents(1,na)),bf_ents(2,na),
     .     ep, nu, (utypes(j),uresid(j), usig(j), uchi2(j),j=1,nu)
 325  format('AMBCON',1x,a4,1x,'PRN',i2.2,' EP ',I6,' NCont ',i1,2x,
     .      5(:,a2,' Res ',F6.3,1x,F6.3,' Chi2 ',F6.2,1x))

C     write(unit,320) na, ep, nu, 
C    .     (utypes(j),uresid(j), usig(j), uchi2(j),j=1,nu)
C320  format('AMBCNT ',i4,' Ep ',i6,' NCont ',i1,2x,
C    .      5(:,a2,' Est ',F6.3,' +- ',F6.3,' Chi2 ',F6.2,2x))

      end
