      subroutine update_posrt(ep, OK, iter )

      implicit none

*     TrackRT and trackRTr routine to update the parameter estimates
*     for the current block of data.  OK is returned true if there
*     was no need for process noise changes or data editing

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

      real*8 ilf1, ilf2   ! Inverse Wavelengths at L1 and L2

      parameter ( ilf1 = fL1/vel_light )
      parameter ( ilf2 = fL2/vel_light )

* PASSED
      integer*4 ep   ! Epoch counter
      integer*4 iter ! Iteration number.  Prediction only on first 
                     ! iteration
      logical OK     ! Set true if no process noise changes or editing

* LOCAL

      real*8 sol_upd(max_obs)        ! Test update to solution
     .,      scale(max_obs)          ! Scale vector used by invert_vis
     .,      acat(max_obs, max_obs)  ! ACA* matrix 
     .,      AC(max_parm,max_obs)    ! Computation of AC matrix
     .,      obvar(max_obs)          ! Apriori variances for data (compare to
                                     ! post-fit residuals).
     .,      dummy                   ! Dummy variable for invert_vis
     .,      ssqr(max_site), swgh(max_site) ! Sum of res^2/wgh and 1/wgh for
                                     ! WRMS calcuations
     .,      wlen                    ! Wavelength to convert cycles to mm
                                     ! (Set zero to exclude range data)


      real*8 max_err                 ! Chi**2 of max postfit error
      real*8 chisum                  ! chi**2 sum
      integer*4 iddmax               ! Double-diff with largest error

      integer*4 ipivot(max_obs)      ! Pivot elements used invert_vis

      integer*4 i, j, k   ! Loop counters
     .,     np   ! Parameter number 
     .,     nd   ! Data type number
     .,     is,iv, js,jv  ! Double differece site and SV numbers
     .,     na   ! Ambuquity number
     .,     kr, kc  ! Counters over rows and columns
     .,     date(5) ! Current date
     .,     num_ldd ! Number of phase double double diffs.

      real*8 sectag  ! Seconds tag 
 
      real*8 comb_coeff(4,6)   ! Combinations to make data types from phase
                     ! and range
      real*8 part_fact(6)      ! Factors for partials to convert range partial
                     ! to obeservable type factor (1/lambda for phase)

      real*8 uvar    ! Data variance to use in data editing.  (min_lvar
                     ! effects this value).
      character*2 comb_lab(6)  ! Combination labels
      character*2 ddtyp(max_obs)  ! Data type for each DD

      logical data_OK   ! Function returns true is data is OK
      logical amb_removed  ! Set true when ambiguity removed
      logical avbad_found  ! Set true when bad residual found
      logical kbit      ! Function to test bit.

      data comb_lab / 'L1', 'L2', 'LC', 'P1', 'P2', 'PC' /
      data comb_coeff / 1.d0, 0.d0, 0.d0, 0.d0,
     .                  0.d0, 1.d0, 0.d0, 0.d0,
     .                  lcf1, lcf2, 0.d0, 0.d0,
     .                  0.d0, 0.d0, 1.d0, 0.d0,
     .                  0.d0, 0.d0, 0.d0, 1.d0,
     .                  0.d0, 0.d0, pcf1, pcf2 /

      data part_fact / ilf1, ilf2, ilf1, 1.d0, 1.d0, 1.d0 /

****  Routine to update the estimates using a standard Kalman filter
*     formulation.  
      if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .write(*,100) ep, iter
 100  format('UPDATE FILTER Ep ',i6,' Iter ',i2)

*     Steps: Prediction step
      if( iter.eq.1 ) then
****     Save the current ambiguity state in case we need to
*        iterate due to data editing
         num_psav = num_parm
         do j = 1, neam
            do i = 1, num_ambs
               amb_save(j,i) = amb_parn(j,i)
            end do
         end do
      else
*        Reset of the ambiguity
         num_parm = num_psav
         do j = 1, neam
            do i = 1, num_ambs
               amb_parn(j,i) = amb_save(j,i)
            end do
         end do
      end if
         
      call pred_step(ep)
      if( ep.ge.debug(3) .and. ep.le. debug(4) )
     .call print_parRT(6, ep,'M','After Predict')

****  Count number of analysis types
      if( ep.eq.1 ) then 
         num_anal_type = 0
         do j = 1,6
             if( index(anal_types(1),comb_lab(j)).gt.0 ) 
     .            num_anal_type = num_anal_type+1
         end do
         if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .   write(*,110) num_anal_type, anal_types(1)
 110     format('There are ',i3,' data types in ',a)
      end if 

*     See if we need to add parameters for new ambiguity parameters
c      print *,'Amb_Parn ', num_ambs, num_parm
      do i = 1, num_ambs
         if( wls_obn(1,i).ne.0 .and. wls_obn(2,i).ne.0 .and.
     .       wls_obn(3,i).ne.0 .and. wls_obn(4,i).ne.0 ) then
            do j = 1, neam
              if( amb_parn(j,i).eq.-1 .and.
     .            .not. kbit(bf_ents(5,i),2) ) then
*                  parameter initialzed so add to estimates
                   num_parm = num_parm + 1
                   
                   if( num_parm.gt.max_parm ) then
                      call report_stat('FATAL','TRACKRT','update_posrt',
     .                   ' ',  'Too many parameters: Max allowed',
     .                        max_parm)
                   endif
                   amb_parn(j,i) = num_parm
                   do k = 1,num_parm
                      cov_parmm(k,num_parm) = 0
                      cov_parmm(num_parm,k) = 0
                   end do
                   sol_vecm(num_parm) = 0
                   cov_parmm(num_parm,num_parm) = 100.d0
               elseif( amb_parn(j,i).eq.-1 .and.
     .                 kbit(bf_ents(5,i),2) ) then
                   if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .             write(*,140) ep, site_names(bf_ents(1,i)),
     .                 bf_ents(2,i)
 140               format('RESAMB Epoch ',i6,1x,a,1x,' G ',I2.2,
     .                 'amb_parn -1 but bias fixed')
                   amb_parn(j,i) = 0
               end if
            end do
          end if
      end do

****  Scan data by site and add up number of data available.  If this is
*     too small then remove all data.
      do k = 2, num_site
         num_ldd = 0
         do i = 1,num_ambs
*           See if all data needed here is OK
            if( bf_ents(1,i).eq.k .and.
     .          wls_obn(1,i).ne.0 .and. wls_obn(2,i).ne.0 .and.
     .          wls_obn(3,i).ne.0 .and. wls_obn(4,i).ne.0 ) then
*              See if we will use this data
               if( data_OK(RT_errflag(wls_obn(1,i)),data_mask) .and.
     .             data_OK(RT_errflag(wls_obn(2,i)),data_mask) .and.
     .             data_OK(RT_errflag(wls_obn(3,i)),data_mask) .and.
     .             data_OK(RT_errflag(wls_obn(4,i)),data_mask) ) then
                  num_ldd = num_ldd + 1
               end if
            end if
          end do
****      OK See how many we have
          if( num_ldd.lt.min_ldd .and. num_ldd.gt.0 ) then
*             Not enough data, so remove all data at this time.
              call mjd_to_ymdhms(RT_MJD_obs(sblk),date, sectag)
              write(lus,180) ep, date, sectag, site_names(k), 
     .            num_ldd, min_ldd
 180          format('DELETE Site  ',i6,' Time ',I4,4(1x,I2.2),1x,
     .              F5.2,' Site ',a4,' Num DD ',i3,' Min ',i3)
*             Now mark as bad
              do i = 1,num_ambs
                 if( bf_ents(1,i).eq.k .and. wls_obn(1,i).ne.0 ) then
                    call sbit(RT_errflag(wls_obn(1,i)),4,1)
                 end if
              enddo
           endif
      end do

***** Now form the double differences and double difference partials. 
*     Double differences are based on ambiguity information
      num_dd = 0

*     Loop over combination types
      nd = 0
      do j = 1,6
         if( index(anal_types(1),comb_lab(j)).gt.0 ) then
*           See what double differences can be formed
            nd = nd + 1
            do i = 1, num_ambs
                amb_dd(i) = 0  ! Clear counter incase not used

*               See if all data needed here is OK
                if( wls_obn(1,i).ne.0 .and. wls_obn(2,i).ne.0 .and.
     .              wls_obn(3,i).ne.0 .and. wls_obn(4,i).ne.0 ) then
                if( data_OK(RT_errflag(wls_obn(1,i)),data_mask) .and.
     .              data_OK(RT_errflag(wls_obn(2,i)),data_mask) .and.
     .              data_OK(RT_errflag(wls_obn(3,i)),data_mask) .and.
     .              data_OK(RT_errflag(wls_obn(4,i)),data_mask) ) then

*                  All one-ways good so add
                   num_dd = num_dd + 1
                   if( num_dd.gt.max_obs ) 
     .             call report_stat('fatal','trackrt','update_posrt',
     .                  '','Too many double differences. Max ',max_obs)
                   amb_dd(i) = num_dd
                   dd_amb(num_dd) = i
                   ddtyp(num_dd) = comb_lab(j)
c                   print *,'Form_ddrt ',j,anal_types(1),comb_lab(j),
c     .                 num_dd, i 
*                  Form the OMC double from one-way raw_omc
                   call form_ddrt(i, comb_coeff(1,j), sol_obs(num_dd))
                
*                  Form the partial derivatives for sites and atmospheres
                   call form_apart(i, part_fact(j), num_dd)

*                  Now see if we need to add partial for phase ambiquities.
C                   print *,'AmbPart ',j, i, num_dd, amb_parn(:,i)
                   if( j.eq.1 ) then     ! L1 
                      np = amb_parn(1,i) 
                      if( np.gt.0 ) apart(np,num_dd) = 1.d0
                   else if ( j.eq.2 ) then  ! L2 
                      np = amb_parn(2,i) 
                      if( np.gt.0 ) apart(np,num_dd) = 1.d0
                   else if( j.eq.3 ) then   ! LC (see if one/two ambs
                      if ( neam.eq.1 ) then
                         np = amb_parn(1,i) 
                         if( np.gt.0 ) apart(np,num_dd) = 1.d0
                      else
                         np = amb_parn(1,i) 
                         if( np.gt.0 ) apart(np,num_dd) = lcf1
                         np = amb_parn(2,i) 
                         if( np.gt.0 ) apart(np,num_dd) = lcf2
                      endif
                   endif
                end if
                end if
            end do

****        OK now increment double difference data covariance.  In forming
*           this we only correlate within a data type i.e., Correlation 
*           between L1 and LC no included.  (Normally this will not happen).
            call form_covdd(ep, nd)
         end if   ! This observable
      end do

****  Out put what we have
      if( ep.ge.debug(5) .and. ep.le.debug(6) ) then
         do i = 1, num_dd
            write(*,220) ep, i, wls_obn(:,dd_amb(i)), sol_obs(i),
     .            sqrt(cov_obs(i,i)), ddtyp(i) 
 220        format('SOLOBS_DD Ep ',i6,' DD ',i4,' OWS ',4I4,
     .            ' OMC ',F12.3,' +- ',F12.3,1x,a2)
c            write(*,240) i, (apart(k,i),k=1,num_parm)
c 240        format('APART_DD  ',i4,250(F10.3))
c            write(*,260) i, (cov_obs(i,j),j=i,num_dd)
c 260        format('COVOBS_DD ',i4,250(F10.4))
         end do
      end if

*     Now do the filtering: First form the ACAT matrix.  Save the AC part of the 
*     matrix
      do i = 1, num_parm
         do j = 1, num_dd
            AC(i,j) = 0.d0
            do k = 1, num_parm
               AC(i,j) = AC(i,j) + apart(k,j)*cov_parmm(i,k)
            end do
         end do
      end do

*     Now form ACAT matrix
      do i = 1, num_dd
         do j = 1, num_dd
             acat(i,j) = 0.d0
             do k = 1, num_parm
                acat(i,j) = acat(i,j) + AC(k,i)*apart(k,j)
             end do
         end do
      end do

*     Now add to data covariance and invert
      do i = 1, num_dd
         obvar(i) = cov_obs(i,i)
         do j = 1, num_dd
            cov_obs(i,j) = cov_obs(i,j) + acat(i,j)
         end do
      end do 

      call invert_vis(cov_obs, dummy, scale, ipivot, num_dd, max_obs, 0)

*     Now finish forming the Kalman Gain matrix
      do i = 1, num_parm
         do j = 1, num_dd
            kgain(i,j) = 0.d0
            do k = 1, num_dd
               kgain(i,j) = kgain(i,j) + AC(i,k)*cov_obs(k,j)
            end do
         end do
      end do


****  Now get the difference between the o-minus-c and the current solution
*     estimate (note: the position corrections are always zero.  This is mainly
*     for the atmospheric delay corrections)
      do j = 1, num_dd
         sol_upd(j) = sol_obs(j)
         do i = 1, num_parm
            sol_upd(j) = sol_upd(j) - apart(i,j)*sol_vecm(i)
         end do
      end do

      do i = 1, num_parm
         sol_vecm(i) = sol_vecm(i)
         do j = 1, num_dd
            sol_vecm(i) = sol_vecm(i) + kgain(i,j)*sol_upd(j)
         end do
      end do

*     Now save the final answer and update the  covariance matrix
      do i = 1, num_parm
         do j = 1, num_parm
            cov_parmm(i,j) = cov_parmm(i,j)
            do k = 1, num_dd
               cov_parmm(i,j) = cov_parmm(i,j) - kgain(i,k)*AC(j,k)
            end do
         end do
      end do

****  Compute the postfit residuals
      max_err = 0    ! Largest chi**2
      iddmax = 0     ! Obs with largest error
      chisum = 0     ! Chi**2 suum 

****  Accumulate statistics by site
      do is = 1, num_site
         ssqr(is) = 0    ! Sum (res/sig)^2
         swgh(is) = 0    ! Sum (1/sig)^2
         num_dd_site(is) = 0
      end do

      do j = 1, num_dd
         sol_upd(j) = sol_obs(j)
         do i = 1, num_parm
            sol_upd(j) = sol_upd(j) - apart(i,j)*sol_vecm(i)
         end do

****     See the worst error
         chisum = chisum + sol_upd(j)**2/obvar(j)
         uvar = obvar(j)
         if( ddtyp(j)(1:1).eq.'L' ) then
             uvar = max(obvar(j),min_lvar)
         end if
         if( sol_upd(j)**2/uvar.gt.max_err ) then
             max_err = sol_upd(j)**2/uvar 
             iddmax  = j
         end if
         na = dd_amb(j)
         is = BF_ents(1,na)
         iv = BF_ents(2,na)
         js = BF_ents(1,wls_ref(2,na))
         jv = BF_ents(2,wls_ref(1,na))

****     Save the residuals for status report
         psv_res(j) = sol_upd(j)
         psv_sig(j) = sqrt(obvar(j))
         psv_ss(1,j) = is
         psv_ss(2,j) = iv
         psv_ss(3,j) = js
         psv_ss(4,j) = jv
         psv_ae(1,j) = RT_azel(1,wls_obn(1,na))
         psv_ae(2,j) = RT_azel(2,wls_obn(1,na))
         psv_amb(j)  = na
         psv_dtype(j) = ddtyp(j)

         if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .   write(*,510) ep, j, na, site_names(is), iv,
     .        site_names(js), jv, psv_res(j), psv_sig(j),
     .        ddtyp(j),RT_azel(2,wls_obn(1,na)),  iddmax 
 510     format('PF_RESID Ep ',i6,' DD ',2I4,1x,a4,' G ',I2.2,
     .          1x,a4,' G ',I2.2,F10.4,' +- ',F10.4,1x,a2,1x,
     .          ' Elev ',F6.2,i3)
*        Summ up statistics (maybe later limit sigma here to exclude
*        range values)
         wlen = 0
         if( ddtyp(j)(1:2).eq.'L1'. or. ddtyp(j)(1:2).eq.'LC' ) then
             wlen = vel_light/fL1*1000.00
         elseif( ddtyp(j)(1:2).eq.'L2' ) then
             wlen = vel_light/fL2*1000.00
         end if
         if( wlen.gt.0 ) then
            ssqr(is) = ssqr(is) + (sol_upd(j)**2/obvar(j))*wlen**2
            swgh(is) = swgh(is) + 1/obvar(j)
            num_dd_site(is) = num_dd_site(is) + 1
         end if

      end do

****  Finish RMS calculations
      do is = 1, num_site
         if( num_dd_site(is).gt. 0 ) then
            wrms_dd_site(is) = sqrt(ssqr(is)/swgh(is))
         else
            wrms_dd_site(is) = 0
         end if
      end do


****  See if we need to edit
      if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .write(*,515) ep, num_dd, sqrt(chisum/num_dd), max_err
 515  format('PF_STATS Ep ',i6,' NumDD ',i4,' NRMS ',F10.4,
     .       ' Max chi^2 ',F14.4)
      OK = .true.
      if( max_err.gt.rms_edtol**2 ) then
         na = dd_amb(iddmax)
         call sbit(RT_errflag(wls_obn(1,na)),4,1)
         OK = .false.
         call mjd_to_ymdhms(RT_MJD_obs(sblk),date, sectag)
         if( (ep.ge.debug(3) .and. ep.le.debug(4)) .or.
     .       obvar(iddmax).lt.1.d-6 ) 
     .   write(*,520) ep, iddmax, na, wls_obn(1,na), max_err
 520     format('Deleting epoch ',I6,' DD ',i4,' Amb Obs ',i4,1x,i4,
     .          ' Max_err ',F14.4)
         j = iddmax
         is = BF_ents(1,na)
         iv = BF_ents(2,na)
         js = BF_ents(1,wls_ref(2,na))
         jv = BF_ents(2,wls_ref(1,na))
         write(*,525) ep, date, sectag, site_names(is), iv,
     .        site_names(js), jv, sol_upd(j), sqrt(obvar(j)),
     .        sqrt(max_err),RT_azel(2,wls_obn(1,na)), ddtyp(j) 
         write(lus,525) ep, date, sectag, site_names(is), iv,
     .        site_names(js), jv, sol_upd(j), sqrt(obvar(j)),
     .        sqrt(max_err), RT_azel(2,wls_obn(1,na)), ddtyp(j) 
 525     format('DELETE Epoch ',i6,' Time ',I4,4(1x,I2.2),1x,
     .          F5.2,' DD ',a4,' G ',i2.2,1x,a4,' G ',I2.2,
     .         ' PF Res ',F10.4,' +- ',F10.4,' Err ',F10.4,
     .         ' Elev ',F6.2,' deg DT ',a2)

*****     Remove the constribution of this observation from the
*         accumuation of the MW and EX widelanes.
          call dec_wlstat(ep, iddmax, na)

****      See if we have seen this residual of this size before
*         in phase data; if more then a few times, add a cycle 
*         slip.
          if( ddtyp(j)(1:1).eq.'L' ) then  ! Phase residual
*             See if seen before
              avbad_found = .false.
              do k = 1, avbad_tot
                 if( is.eq.avbad_sse(1,k) .and. 
     .               iv.eq.avbad_sse(2,k) ) then

                    avbad_found = .true.

*                   OK, we have seen before
                    if( ep-avbad_sse(3,k).le.num_edtol ) then
*                       OK; prior bad residual within num_edtol
*                       epochs.  Now see if we have seen enough
                        if( avbad_num(k).ge. num_edtol ) then
*                           Add cycle slip and remove this av_bad
*                           entry
                            call reset_cslip( na, ep )

                            write(lus,680) ep, 
     .                          site_names(bf_ents(1,na)), 
     .                          bf_ents(2,na), psv_res(j), psv_sig(j),
     .                          RT_azel(2,wls_obn(1,na)), ddtyp(j) 

 680                        format('CSLIP Ep ',i6,' Site ',a,' G ',I2.2,
     .                         ' REPEAT DD Error ',F12.2,1x,F5.2,
     .                         ' cyc, Elev', F6.2,' deg, Type ',a)

                            call rm_avbad( k,avbad_sse,  avbad_tot)

                        else ! Valid but not above number yet so 
*                          save this epoch and increment number
                           avbad_num(k)   = avbad_num(k)+1
                           avbad_sse(3,k) = ep
                        endif
                    else   ! Last bad residual was too long ago
*                       Reset this counter
                        avbad_num(k)   = 1
                        avbad_sse(3,k) = ep
                   endif
                 end if
              end do
*             If this is first bad residual for this site/sv save
              if( .not. avbad_found ) then
                  avbad_tot = avbad_tot + 1
                  if( avbad_tot.gt.max_obs )
     .            call report_stat('fatal','trackRT','update_posrt',
     .                 '','Too many bad DD saved, Max ',max_obs)
     
                  avbad_num(k)   = 1
                  avbad_sse(1,k) = is
                  avbad_sse(2,k) = iv
                  avbad_sse(3,k) = ep
              endif
          end if
      end if

****  See if we should remove any old bad dd
      j = 0
      do while ( j.lt. avbad_tot .and. j.ge.0 )
         j = j + 1
         if( ep-avbad_sse(3,j).gt. num_edtol ) then
            call rm_avbad( j,avbad_sse,  avbad_tot)
            j = j - 1  ! Continue with entry number
         end if
      end do


***** Print the solution
      if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .call print_parRT(6, ep,'M','Standard Output')

****  OK, If solution is OK so far, see if we can resolve ambiguities.
      if( OK ) then
         call res_ambRT( ep, OK )
      end if

****  If the solution is OK, update the solution vector and 
*     covariance matrix
      if ( OK ) then 
         do i = 1, num_parm
            sol_vecp(i) = sol_vecm(i)
            do j = 1, num_parm
               cov_parmp(i,j) = cov_parmm(i,j)
            end do
         end do

*****    Set all all parameter numbers to positive
         do i = 1, num_ambs
            do j = 1, 2
*              -1 value is used for reference site and a real
*              ambiguity can't have this value because site
*              and atm parameter would be down here
               if( amb_parn(j,i).lt. -1 ) 
     .            amb_parn(j,i) = abs(amb_parn(j,i))
            end do
         end do

****     Now remove any ambiguity parameters that have been fixed to
*        integer values
         amb_removed = .false.
         do i = 1, num_ambs
*           if the bias is fixed and still has a parameter
            if( kbit(bf_ents(5,i),2) .and. amb_parn(1,i).gt.0 ) then
*               Need to remove and remap
                do j = 1,2  ! L1/LC and L2 (if separate)
                   np = amb_parn(j,i)
                   if( np.gt.0 ) then ! Remove parameters
                      amb_removed = .true.
                      do kr = np+1,num_parm
                         sol_vecp(kr-1) = sol_vecp(kr)
                         do kc = 1, num_parm
                            cov_parmp(kr-1,kc) = cov_parmp(kr,kc)
                         end do
                      end do
                      do kc = np+1, num_parm
                         do kr = 1, num_parm
                            cov_parmp(kr,kc-1) = cov_parmp(kr,kc) 
                         end do
                      end do
*                     Now update the parameter pointers.  Move all
*                     parameter numbers greater than np down by one
                      do kr = 1, num_ambs
                         do kc = 1,2
                            if( amb_parn(kc,kr).gt.np ) then
                               amb_parn(kc,kr) = amb_parn(kc,kr)-1
                            endif
                         end do
                      end do
*                     OK: Reset the values for this parameter
                      amb_parn(j,i) = 0
                      num_parm = num_parm - 1
                   end if
                end do
            endif
         enddo
         if( amb_removed ) then
             if( ep.ge.debug(3) .and. ep.le. debug(4) )
     .       call print_parRT(6, ep,'P','Ambiquity resolved')
         end if
      endif

****  Thats all
      return
      end

ctitle print_parRT

      subroutine print_parRT(un, ep, pm,  title)

*     Debug routine to print parameter estimates

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 ep   ! Epoch
     .,         un   ! Unit number for output

      character*(*) title   ! Title for output
      character*(*) pm      ! Either P or M of the Plus/Minus soln

* LOCAL
      integer*4 i,j   ! Counters
      integer*4 np    ! Parameter number
           
      character*1 crd_type(3)

      data crd_type / 'X', 'Y', 'Z' / 
 
      write(un,510) ep, num_parm, pm, title
 510  format('PARAMETER ESTIMATES Epoch ',i6,' for NP ',i4,
     .       ' Type ',a,1x,a)

      do i = 1, num_site
         do j = 1, 3
            np = site_parn(j,i)
            if ( np.gt.0 ) then
               if( pm(1:1).eq.'M' ) then
                  write(un,530) ep, np, site_names(i), crd_type(j), 
     ,                sol_vecm(np),
     .                sqrt(cov_parmm(np,np)),cov_parmm(np,np) 
               else
                  write(un,530) ep, np, site_names(i), crd_type(j),
     .                sol_vecp(np),
     .                sqrt(cov_parmp(np,np)),cov_parmp(np,np) 
               end if

 530           format('PARAMETER Ep ',i6,' NP ',i3,1x,a4,1x,
     .               a1,' dCrd ',1x,f10.4,' +- ',f10.4,' m   Var ',
     .               E13.3)
            endif
         end do
      end do

      do i = 1, num_site
         np = atm_parn(i)
         if ( np.gt.0 ) then
            if( pm(1:1).eq.'M' ) then
               write(un,540) ep, np, site_names(i), sol_vecm(np),
     .           sqrt(cov_parmm(np,np)), cov_parmm(np,np)
            else
               write(un,540) ep, np, site_names(i), sol_vecp(np),
     .           sqrt(cov_parmp(np,np)), cov_parmp(np,np)
            end if
 540        format('PARAMETER Ep ',i6,' NP ',i3,1x,a4,1x,'dATM    ',
     .           f10.4,' +- ',f10.4,' m   Var ',E13.3)
         endif
      end do

*     Output the ambiguity estimates
      do i = 1, num_ambs
         do j = 1, neam
            np = amb_parn(j,i)
            if( np.gt.0 ) then
               if( pm(1:1).eq.'M' ) then
                  write(un,560) ep, np, site_names(bf_ents(1,i)),
     .             bf_ents(2,i), j, sol_vecm(np),
     .             sqrt(cov_parmm(np,np)), cov_parmm(np,np)
               else
                  write(un,560) ep, np, site_names(bf_ents(1,i)),
     .             bf_ents(2,i), j, sol_vecp(np),
     .             sqrt(cov_parmp(np,np)), cov_parmp(np,np)
               end if
 560           format('PARAMETER Ep ',i6,' NP ',i3,1x,a4,1x,'G ',
     .               i2.2,1x,i1,2x, f10.4,' +- ',f10.4,' cyc Var ',
     .               E13.3)
            endif
         end do
      end do

****  Return 
      return
      end

CTITLE RM_AVBAD

      subroutine rm_avbad( ent, avbad_sse, avbad_tot )

      implicit none

*     Move avbad_sse list down by one entry
* PASSED
      integer*4 avbad_tot   ! Total number
     .,         avbad_sse(3,avbad_tot)  ! List 
     .,         ent     ! Entry to be removed

* LOCAL
      integer*4 k   ! counter

      do k = ent+1, avbad_tot
         avbad_sse(1,k-1) = avbad_sse(1,k)
         avbad_sse(2,k-1) = avbad_sse(2,k)
         avbad_sse(3,k-1) = avbad_sse(3,k)
      end do
      avbad_tot = avbad_tot -1

      return
      end





