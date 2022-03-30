CTITLE CSLIP_REPAIR

      subroutine cslip_repair

      implicit none

*     Routine to scan ambiquity list to see if there are slips
*     that can be repaired before trying to estimate.  The slips
*     are repaired with signal differences.

      include 'track_com.h' 
      include '../includes/const_param.h'

* LOCAL VARIABLES

      integer*4 ep     ! Epoch 
     .,         i      ! site number
     .,         pn     ! satellite number
     .,         ns, ms     ! site numbers in single differences
     .,         IdN1, IdN2 ! L1 L2 cycles offset
     .,         IdN12      ! N1-N2 cycles (from MWWL) 

      real*8 av_mwwl, av_exwl  ! Average change on MW-wl and EX-wl

      real*8 bmwwl, bexwl  ! Before mw-wl and exwl
     .,      amwwl, aexwl  ! After mw-wl and exwl
     .,      sig_awl(2), sig_bwl(2)  ! Sigmas of EX and MW mean WLs
     .,      mwwl, exwl    ! Mw and ex widelanes
     .,      owdata(4,2)     ! Phase and range owdata for the two sites
                           ! in the SD.
     .,      elev(2)       ! Elevation angle
     .,      dcberr(2,2)   ! DCB range contributions
     .,      bwls(2,max_wl)     ! Saved values of EX and ME WL Before
     .,      awls(2,max_wl)     ! Saved values of EX and ME WL After
     .,      damb(2)       ! Change in cycles to remove cycle slip.

      real*8 dN1, sN1, sN12     ! Change in N1 cycles and sigma and MW sig

* MOD TAH 180501: LSQ estimate of slip with covariance matrix dchi for test
      real*8 NLx_neq(2,2)  ! Normal equations for L1/L2 cycle estimates 
     .,      NLx_bvc(2)    ! B-vector for estimate


      integer*4 bnum, anum  ! Before and after number of values
     .,         bad         ! number of bad values in building list
     .,         df(2)       ! Data flag
     .,         j,k         ! Epoch conunter
     .,         ch          ! Channel number
     .,         it, tot_iter, tot_resolved

      logical OK(2)         ! True if data OK
      logical kbit
      logical csfix         ! Set true is cycle slip is repairable
* all_wls_res -- Logical function which returns true when all the widelanes
*     at a site have been resolved 
* done  -- Set true if we have resolved all ambiquities

      logical all_wls_res, done

      integer*4 PtoL        ! Function return position in list of PRN 
                            ! in SP3 files list of satelites
      integer*4 lv          ! List number for vechile.

****  Start looping through the ambiquities that are not fixed
      write(*,100) data_mask
 100  format('CSLIP: Using Data_Mask ',o12)
      do i = 1, num_ambs

         if( .not.kbit(bf_ents(5,i),2) ) then ! Not resolved
*            Get one-way data we need before ambiquity
             k = bf_ents(3,i)
             pn = bf_ents(2,i)
             ns = bf_ents(1,i)
*            Get second site from the wls_ref entry
             ms = bf_ents(1,wls_ref(2,i))
             bnum = 0
             bad = 0
             do while ( k.gt.max(1,bf_ents(3,i)-2*num_mwwl) )
                 k = k - 1
                 call get_ow(k,ns,pn, owdata(1,1), elev(1), OK(1),df(1))
                 call get_ow(k,ms,pn, owdata(1,2), elev(2), OK(2),df(2))

                 call get_dcb(ns,pn, dcberr(1,1)) 
                 call get_dcb(ms,pn, dcberr(1,2)) 
*                If both are OK, then sum into MW and EX estimates
                 if( OK(1) .and. OK(2) ) then
                     bnum = bnum + 1
                     call save_wls( owdata, dcberr, bwls(1,bnum))       
                 else
                     bad = bad + 1
                 endif
                 if( bnum.eq.num_mwwl ) k = 0  ! Enough data
                 if( bad .eq.num_mwwl ) k = 0  ! Too many bad data
             end do
             if ( bad.eq.num_mwwl .and. bnum.gt.0 ) then
                write(*,120) i, bf_ents(3,i), pn, site_names(ns), 
     .                       site_names(ms), bad, bnum
 120            format('CSCLIP_REPAIR: Amb # ',I4,' Ep ',i6,' PRN',
     .                  I3.3,' Sites ',a,1x,a,' Too many bad ',I3,
     .                       ' Good ',i3)
             endif

***          Now do the after side
             anum = 0
             bad = 0
             do k = bf_ents(3,i), num_epochs
                 call get_ow(k,ns,pn, owdata(1,1), elev(1), OK(1),df(1))
                 call get_ow(k,ms,pn, owdata(1,2), elev(2), OK(2),df(2))
                 call get_dcb(ns,pn, dcberr(1,1)) 
                 call get_dcb(ms,pn, dcberr(1,2)) 
*                If both are OK, then sum into MW and EX estimates
                 if( OK(1) .and. OK(2) ) then
                     anum = anum + 1
                     call save_wls( owdata, dcberr, awls(1,anum))
                 else
                     bad = bad + 1
                 endif
                 if( anum.eq.num_mwwl ) exit
                 if( bad .eq.num_mwwl ) exit
             end do

****         Now see if have enough data
             if( bnum.eq.num_mwwl .and.  anum.eq.num_mwwl ) then
                 call aver_wls( bexwl, bmwwl, sig_bwl, bnum, 
     .                          num_exwl, bwls)
                 call aver_wls( aexwl, amwwl, sig_awl, anum, 
     .                          num_exwl, awls)
                 write(*,320) i, bf_ents(3,i), pn, site_names(ns), 
     .                 site_names(ms), 
     .                 aexwl-bexwl, sqrt(sig_awl(1)**2+sig_bwl(1)**2),
     .                 amwwl-bmwwl, sqrt(sig_awl(2)**2+sig_bwl(2)**2)
 320             format('CSLIP Amb # ',I4,' Ep ',i6,' PRN',
     .                  I3.3,' Sites ',a,1x,a,' dEX ',F8.4,1x,F8.4,
     .                       ' dMW ',F8.4,1x,F8.4)

****             Now compute the number of cycles and sigmas
* MOD TAH 180323: Glonass processing: Compensate for remapping of frequencies.
*                (scale by fsv/fRef to get to integer and then inverse when
*                applying to data)
                 lv = PtoL( pn ) 
                 IdN12 = nint((amwwl-bmwwl)*fL1(lv)/fR1)
                 dN1   = ((aexwl-bexwl)-(fR1/fR2)*IdN12)/(1-(fR1/fR2))
                 IdN1  = nint(dN1*fL1(lv)/fR1)
                 IdN2  = IdN1 - IdN12
                 sN1 = sqrt(sig_awl(1)**2+sig_bwl(1)**2 +
     .               (fR1/fR2)**2*(sig_awl(2)**2+sig_bwl(2)**2))
                 sN12 = sqrt(sig_awl(2)**2+sig_bwl(2)**2)
                 write(*,325)  pn, lv, fL1(lv)*1.e-6, fR1*1.e-6, 
     .                         fL2(lv)*1.e-6, fR2*1.e-6, 
     .                         (amwwl-bmwwl)*fL1(lv)/fR1, dN1
 325             format('CSLIP ',I3,1x,i3,1x, ' fL1 ',2F10.4,1x,
     .                  ' fL2 ',2F10.4,' A-B MW-WL ',2F8.4)
****             See if resolved OK: Test absolute deviation 
                 csfix = .false.
                 if( abs((dN1-IdN1)).lt.cslip_tex .and. 
     .               abs(((amwwl-bmwwl)-IdN12)).lt.cslip_tmw .and. 
     .               sN1.lt.cslip_sex .and. sN12.lt.cslip_smw   ) then
                     csfix = .true.
                 endif

                 write(*,330) i, bf_ents(3,i), pn, site_names(ns), 
     .                 site_names(ms), IdN1, IdN2, 
     .                 dN1-IdN1, (amwwl-bmwwl)-IdN12, sN1, sN12,
     .                 csfix 
 330             format('CSRES Amb # ',I4,' Ep ',i6,' PRN',
     .                  I3.3,' Sites ',a,1x,a,' dN12 ',2I4,' Err ',
     .                       2(F7.3,1x),' +- ',2(F7.3,1x),' Fix ',L1)

* MOD TAH 180501: New code using covariance matrix
                 call NLx_formneq(lv, aexwl-bexwl, amwwl-bmwwl,  
     .                       sig_awl, sig_bwl, NLx_neq, NLx_bvc)

                 call NLx_est(lv, NLx_neq, NLx_bvc, damb, csfix)

****             OK we have the estimates of changes in number of
*                cycles to fix the slip relative to current ambiq_all.
*                Save that we will fix this slip (it will be removed
*                from the bf_ents later) and update data and ambiq_all
*                for this one ambiquity. 
* MOD TAH 180323: Glonass: Map the integers back to frequency
!                damb(1) = (ambiq_all(1,i-1)-ambiq_all(1,i)) + 
!    .                                                 IdN1*fR1/fL1(lv) 
!                damb(2) = (ambiq_all(2,i-1)-ambiq_all(2,i)) + 
!    .                                                 IdN2*fR2/fL2(lv)  
                 damb(1) = (ambiq_all(1,i-1)-ambiq_all(1,i)) - damb(1)
                 damb(2) = (ambiq_all(2,i-1)-ambiq_all(2,i)) - damb(2)

*                Apply this change to all the effected L1o and L2o data.
                 do k = bf_ents(3,i), bf_ents(4,i)
                    ch = -1
                    do j = 1,num_chan_se(ns,k)
                        if( ctop_cse(j,ns,k).eq.pn ) ch = j
                    end do
                    if( ch.gt.0 ) then
                        L1o_all_cse(ch,ns,k) = L1o_all_cse(ch,ns,k)
     .                                          - damb(1)
                        L2o_all_cse(ch,ns,k) = L2o_all_cse(ch,ns,k)
     .                                          - damb(2)
*                       If this is the first epoch, remove the bias flag
*                       if csfix is true (but if user added, do not remove).
                        if( k.eq.bf_ents(3,i) .and. csfix .and.
     .                     .not. kbit(data_flag_cse(ch,ns,k),7) ) then
                           write(*,340) i, k, site_names(ns), pn, damb
 340                       format('CSLIP Amb # ',I4,' Ep ',i6,1x,a,
     .                            ' PRN',I3.3,' Removed; Slip ',2F10.2)
                           call sbit(data_flag_cse(ch,ns,k),3,0)
                           call sbit(data_flag_cse(ch,ns,k),4,0)
                        elseif( k.eq.bf_ents(3,i) ) then
                           write(*,350) i, k, site_names(ns), pn, damb
 350                       format('CSLIP Amb # ',I4,' Ep ',i6,1x,a,
     .                            ' PRN',I3.3,' Updated; Slip ',2F10.2)

                        endif
                    endif
                 end do
                 ambiq_all(1,i) = ambiq_all(1,i-1)
                 ambiq_all(2,i) = ambiq_all(2,i-1)
             end if
          endif
      enddo

****  Now redo the ambiquity selection and calculation
      data_mask = 18
      num_ambs = 0
      do i = 1, num_site
         call flag_gaps(i)
      end do

      call select_arb_bf

*     Now we resolve the widelane ambiquities for each site
*     Initialize the saved values of the widelines
      do i = 1, num_ambs
         ambiq_all(:,i) = 0.d0
         wls_all(:,i) = 0.d0
         wls_ref(:,i) = 0
      end do

* MOD TAH 070109: Set the WLS_ref to -1 for ones fixed that this
*     time.  This is used so that we can fix to reference site is
*     posssible
      do j = 1,num_ambs
         if( bf_ents(1,j).gt.1 .and. bf_ents(5,j).eq.1 ) then
            wls_ref(1,j) = -1
         end if
      end do

      
      it = 0
      tot_iter = 0
      done = .false.
      do while ( it.lt. 50 .and. .not.done )
         it = it + 1
         tot_iter = tot_iter + 1
         done = .true.
         do i = 2, num_site
            if ( .not.all_wls_res(i) ) then
                tot_resolved = 0
                done = .false.
                call resolve_wl(i,it, tot_resolved)

*               Now see if we may need to arbitarily fix another bias
*               flag to continue.  (This would be caused by a complete
*               gap in the data where all satellites would be given
*               new bias flags).  Since we reset an arbitary bias flag
*               we reset the iteration counter.
                if( tot_resolved.eq.0 .and. tot_iter.gt. 20 .and.
     .              it.gt.10 ) then
C                   call report_bf(6,'Check_Arb')
                    call check_arb_bf(i,it, tot_iter)
                end if
            end if
         end do
      end do

****  Now loop over all the ambiquities find those that need to
*     be resolved.  Go through and mark all the fixed biases as
*     finally resolved
      do i = 1, num_ambs
         if( wls_ref(1,i).le.0 ) then
             call sbit(bf_ents(5,i),2,1)
         end if
      end do

      call recomp_wls

      call report_bf(6,'CSINIT')

***** That all
      return
      end

CTITLE SAVE_WLS

      subroutine save_wls( owdata, dcberr, wls )

      implicit none

*     Routine to compute and save single difference wide-lanes

      include '../includes/const_param.h'
      include 'track_com.h' 
             
      real*8 owdata(4,2)     ! Phase and range data for the two sites
                           ! in the SD.
     .,      dcberr(2,2)   ! DCB range contributions
     .,      wls(2)        ! Widelane saved in k'th entry (EX, MW)

*     Compute the EX and MW WL and save

      wls(1) =  (owdata(1,1) - (fR1/fR2)*owdata(2,1)) -
     .          (owdata(1,2) - (fR1/fR2)*owdata(2,2)) 
      wls(2) =  (owdata(1,1) - owdata(2,1) -
     .           dfsf*((owdata(3,1)+dcberr(1,1))*fR1/vel_light+
     .                 (owdata(4,1)+dcberr(2,1))*fR2/vel_light))-
     .          (owdata(1,2) - owdata(2,2) -
     .           dfsf*((owdata(3,2)+dcberr(1,2))*fR1/vel_light+
     .                 (owdata(4,2)+dcberr(2,2))*fR2/vel_light))

***** Thats all
      return 
      end

CTITLE AVER_WLS

      subroutine aver_wls( mexwl, mmwwl, sig_mwl, nmw, nex, wls)

      implicit none

*     Routine to compute averge value and sigmas of averages

      integer*4 nmw, nex      ! Number of MW and EX values

      real*8 mexwl, mmwwl     ! Mean EXand MW WL
     .,      sig_mwl(2)       ! Sigmas of mean for EX and MW WL
     .,      wls(2,nmw)       ! Estimates of EX and MW WL

      integer*4 k
      real*8 sv, ss           ! Sum of values and squares

*     Do EX
      sv = sum(wls(1,1:nex))
      ss = dot_product(wls(1,1:nex),wls(1,1:nex))
      mexwl = sv/nex 
      sig_mwl(1) = sqrt((ss-mexwl**2*nex)/nex**2)

*     Do MW
      sv = sum(wls(2,1:nmw))
      ss = dot_product(wls(2,1:nmw),wls(2,1:nmw))
      mmwwl = sv/nmw 
      sig_mwl(2) = sqrt((ss-mmwwl**2*nmw)/nmw**2)

!     write(*,120) 'EX ',nex, mexwl, sig_mwl(1), wls(1,1:nex)
!     write(*,120) 'MW ',nmw, mmwwl, sig_mwl(2), wls(2,1:nmw)
!120  format(a3,I3,15(1x,F10.3))

      return 
      end

