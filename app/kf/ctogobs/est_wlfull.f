CTILE EST_WLFULL

      subroutine est_wlfull(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, bf_type_cse, azel_cse)

      implicit none

*     Routine to form direct estimate of Ambiquiuties
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES

*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
     .,    bf_type_cse(num_chan, num_cfiles, num_ep)

*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation

      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)

*   obs(5)  - MW-WL, L2-L1, L1, L2, LC residuals in cycles
*   wgh(5)  - Weight of measurement (inversely proporational to
*             1/sin(el)**2 and phase noise on satellite

      real*8 obs(5,4), wgh(5,4), tot_wgh

*   OK      - Logical: Set OK if data is OK
      logical OK(4)


* LOCAL VARIABLES
      integer*4 ep    ! Epoch being extracted

      integer*4 ch    ! Channel number at this particular epoch
     .,      ltoc     ! Function to return channel number of specific
                      ! satellite list number
      integer*4 i,j,k,l, in,jn , ent, bl_ent
      integer*4 np(4) 
     .,   nump          ! Number of bias parameters

      real*8 ap(4), dobs_mwwl, dobs_exwl, blen
      real*8 DDCov, DDSmw, DDSex, scale(max_cfiles*max_gprn)
      real*8 sum_wgh(max_cfiles*max_gprn), sum_cos(max_cfiles*max_gprn),
     .       sum_sin(max_cfiles*max_gprn), rms_mwl, phs_mwl
      integer*4 ipivot(max_cfiles*max_gprn),
     .       sum_num(max_cfiles*max_gprn)

      real*8 neq(max_cfiles*max_gprn,max_cfiles*max_gprn), 
     .       bvec(max_cfiles*max_gprn,2)

      real*8 neqb(max_cfiles+3*max_gprn, max_cfiles+3*max_gprn),
     .       bvcb(max_cfiles+3*max_gprn),
     .       apbias(max_cfiles+3*max_gprn)
      real*8 res, rec, obc, DDSBS
* MOD TAH 200617 Added frequency ratio (only fr1 needed)
      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.

      integer*4 off, it, p, iv, iw, jv, jw

      integer*4 stat_dcb(max_cfiles)   ! Type of DCB by station based
                  ! on bit 28 and 29 of data flag

      logical kbit, done


      logical data_OK ! Function which returns true to data passes
                      ! the phs_mask on the data flag
     .,       good_bf ! Function to see if good bias flags
      logical ibfd(max_cfiles*max_gprn) ! Logical set true when an implicit bias
                          ! have been removed (to avoid doing it multiple
                          ! times when a new bias flag have been found).

***** OK: All done, Start formal estimation
c      call clrneq(num_cfiles*num_sat, neq, bvec)
      call report_stat('status','autcln','EST_WLFULL',' ','Start',0)
      if (num_cfiles.gt.50 ) then
         stop 'too many cfiles for this code'
      end if
      nump = num_cfiles*num_sat
      do  i = 1, nump
          bvec(i,1) = 0.d0
          bvec(i,2) = 0.d0
          do j = 1,nump
             neq(i,j) = 0.d0
          end do
      end do

      do i = 1,num_cfiles
         stat_dcb(i) = -1
      end do

      ap(1) =  1.0d0
      ap(2) = -1.0d0
      ap(3) = -1.0d0
      ap(4) =  1.0d0

****  Initialize the statistics 
      do i = 1, nump
         sum_wgh(i) = 0.0d0
         sum_cos(i) = 0.0d0
         sum_sin(i) = 0.0d0
         sum_num(i) = 0.0d0
      end do

****  Now build up the normal equations
      do ep = 1, num_ep
         do i = 1, nump
            ibfd(i) = .false.
         end do

         do i = 1, num_cfiles-1
            do j = 1, num_sat-1

*              Loop over all epochs and satellite combinations 
*              solving biases implictly when needed.
               ch =  ltoc( ctol_cse(1,i,ep),j, 
     .                     actual_max_chan)

*              See if bias flag on good data that we need to remove
*              bias flag.
               if( ch.gt.0 ) then
                 if( good_bf(data_flag_cse(ch,i,ep),0,phs_mask) ) then
                     np(1) = (i-1)*num_sat+j
                     call implicit_bf(np(1),num_cfiles*num_sat, 
     .                              max_cfiles*max_gprn, neq, bvec,2)
                     ibfd(np(1)) = .true.
                 endif
               endif 

*              Data at this time is in the interval of the last
*              bias flag so see which DD ambs it contributes to
               call get_obs(ep, i,j,  
     .            L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .            L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .            azel_cse, obs(1,1), wgh(1,1), OK(1))

*              If OK:then this oneway matches and has good data
               if( OK(1) ) then 

*                  Set dcb type for site
                   if( stat_dcb(i).eq.-1 ) then
                      stat_dcb(i)  = 0
                      if( kbit(data_flag_cse(ch,i,ep),28) ) then
                         stat_dcb(i)  = 1
                       endif
                       if( kbit(data_flag_cse(ch,i,ep),29) ) then
                         stat_dcb(i)  = 2
                       end if
                   end if

*                  Find the next station
                   do k = i+1, num_cfiles
                      ch =  ltoc( ctol_cse(1,k,ep),j, 
     .                            actual_max_chan)

                      if( ch.gt.0 ) then

*                        See if implicit bias solution needed for this
*                        site/sv combination
                         np(2) = (k-1)*num_sat+j
                         if( good_bf(data_flag_cse(ch,k,ep),
     .                                          0,phs_mask) .and.
     .                        .not. ibfd(np(2)) ) then
                            call implicit_bf(np(2),num_cfiles*num_sat, 
     .                              max_cfiles*max_gprn, neq, bvec,2)
                            ibfd(np(2)) = .true.
                         endif
                         OK(2) = data_OK(data_flag_cse(ch,k,ep),
     .                                 0,phs_mask)
                      else
                          OK(2) = .false.
                      endif

                      if ( OK(2) ) then
*                         Set dcb type for site
                          if( stat_dcb(k).eq.-1 ) then
                             stat_dcb(k)  = 0
                             if( kbit(data_flag_cse(ch,k,ep),28) ) then
                                stat_dcb(k)  = 1
                             endif
                             if( kbit(data_flag_cse(ch,k,ep),29) ) then
                                stat_dcb(k)  = 2
                             end if
                         end if
                         call get_obs(ep, k, j, L1r_phs_cse, 
     .                      L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .                      L1r_rng_cse, L2r_rng_cse, 
     .                      data_flag_cse, ctol_cse,azel_cse, 
     .                      obs(1,2), wgh(1,2), OK(2))

*                        Now find next satellites
                         do l = j+1,num_sat
*                           Check station 1 to see of OK
                            ch =  ltoc( ctol_cse(1,i,ep),l, 
     .                               actual_max_chan)
                            OK(3) = .false.
                            OK(4) = .false.
                            if( ch.gt.0 ) then
**                             See if implicit bias solution
*                              needed here
                               np(3) = (i-1)*num_sat+l
                               if( good_bf(data_flag_cse(ch,i,ep),
     .                                          0,phs_mask) .and.
     .                              .not. ibfd(np(3)) ) then
                                  call implicit_bf(np(3),
     .                               num_cfiles*num_sat,  
     .                               max_cfiles*max_gprn, neq, bvec,2)
                                  ibfd(np(3)) = .true.
                               endif

                               OK(3)=data_OK(data_flag_cse(ch,i,ep),
     .                                 0,phs_mask)
                            endif

                            ch =  ltoc( ctol_cse(1,k,ep),l, 
     .                               actual_max_chan)
                            if( ch.gt.0 ) then
                               np(4) = (k-1)*num_sat+l
                               if( good_bf(data_flag_cse(ch,k,ep),
     .                                          0,phs_mask) .and.
     .                             .not. ibfd(np(4)) ) then 
                                  call implicit_bf(np(4),
     .                               num_cfiles*num_sat,  
     .                               max_cfiles*max_gprn,neq, bvec,2)
                                  ibfd(np(4)) = .true.
                               endif

                               OK(4)=data_OK(data_flag_cse(ch,k,ep),
     .                                 0,phs_mask)

                            endif
*                           If both OK, continue
                            if( OK(3) .and. OK(4) ) then

                               call get_obs(ep, i, l, L1r_phs_cse, 
     .                           L2r_phs_cse, L1_cyc_cse,
     .                           L2_cyc_cse, L1r_rng_cse,  
     .                           L2r_rng_cse,data_flag_cse,ctol_cse, 
     .                           azel_cse, obs(1,3), wgh(1,3), OK(3))
                               call get_obs(ep, k, l, L1r_phs_cse, 
     .                           L2r_phs_cse, L1_cyc_cse,
     .                           L2_cyc_cse, L1r_rng_cse,  
     .                           L2r_rng_cse,data_flag_cse,ctol_cse, 
     .                           azel_cse, obs(1,4), wgh(1,4), OK(4))

*****                          We have data, so increment normal equations
                               tot_wgh = 1.d0/(1.d0/wgh(1,1)+
     .                                         1.d0/wgh(1,2)+
     .                                         1.d0/wgh(1,3)+
     .                                         1.d0/wgh(1,4))
                               ent = bl_ent(i,k)
                               blen = baselens(ent)
                               tot_wgh = tot_wgh*
     .                             1.d0  ! weight for MWWL 
C    .                             sqrt(1.d-4/max(blen/12.d6,1.d-4))
C    .                             (12.d6/max(blen,1000.d0))
                               dobs_mwwl = (obs(1,1)-obs(1,2))-
     .                                     (obs(1,3)-obs(1,4))
* TEST: Make LC residual (changed 2, to 5,)
                               dobs_exwl = (obs(5,1)-obs(5,2))-
     .                                     (obs(5,3)-obs(5,4))
                               np(1) = (i-1)*num_sat+j
                               np(2) = (k-1)*num_sat+j
                               np(3) = (i-1)*num_sat+l
                               np(4) = (k-1)*num_sat+l
*                              Increment equations
                               do in = 1,4
                                   bvec(np(in),1) = bvec(np(in),1)+
     .                                ap(in)*dobs_mwwl*tot_wgh
                                   bvec(np(in),2) = bvec(np(in),2)+
     .                                ap(in)*dobs_exwl*tot_wgh
                                   do jn = 1,4
                                       neq(np(in),np(jn)) = 
     .                                    neq(np(in),np(jn)) +
     .                                    ap(in)*ap(jn)*tot_wgh
                                   enddo
      
****                               Increment the statistics
                                   sum_wgh(np(in)) =sum_wgh(np(in))+
     .                                 tot_wgh
                                   sum_cos(np(in)) =sum_cos(np(in))+
     .                                 tot_wgh*cos(dobs_mwwl*2*pi)
                                   sum_sin(np(in)) =sum_sin(np(in))+
     .                                 tot_wgh*sin(dobs_mwwl*2*pi)
                                   sum_num(np(in)) =sum_num(np(in))+1
 
                               enddo
                            endif   ! OK(3) and OK(4)
                         enddo      ! Loop on l satellite
                      endif         ! Second station OK
                   enddo            ! Loop on k station
               endif                ! First station/satellite OK
            end do                      ! First satellite
         end do                         ! First station
      enddo                            ! Looping over epochs

C     write(*,410) num_cfiles, num_sat, num_cfiles*num_sat
C410  format('AMB: ',3i6)
C     do i = 1,num_cfiles*num_sat
C         if( neq(i,i).eq.0 ) neq(i,i) = 1.d0
C         write(*,420) (neq(i,j),j=1,num_cfiles*num_sat)
C420      format('NEQ: ',1600E18.10)
C     enddo
C     do i = 1, num_cfiles*num_sat
C         write(*,440) bvec(i,1)
C440      format('BV1: ',1600E18.10)
C     enddo
C     do i = 1, num_cfiles*num_sat
C         write(*,460) bvec(i,2)
C460      format('BV2: ',1600E18.10)
C     enddo


***** Now apply basic constraints; One reference site and satellite
      write(*,480) wl_ref_site, cf_codes(wl_ref_site)
 480  format('Constraining site ',i4,2x,a4)
      do i = 1,num_sat
         in = (wl_ref_site-1)*num_sat+i 
         neq(in,in) = neq(in,in) + 1.d10
      enddo

      write(*,485) wl_ref_svs, prn_list(wl_ref_svs)
 485  format('Constraining sv ',i4,' PRN ',i2.2)
      do i = 1,num_cfiles
        if( i.ne.wl_ref_site) then
           in = (i-1)*num_sat + wl_ref_svs
           neq(in,in) = neq(in,in) + 1.d10
        endif
      enddo

*     Now apply a loose constraint for dependent biases
      do i = 1, nump 
         neq(i,i) = neq(i,i) + 0.1d0
      end do
      
      call invert_vis(neq,bvec,scale,ipivot,nump,max_cfiles*max_gprn,2)

*
C     do i = 1,num_cfiles*num_sat
C         if( neq(i,i).eq.0 ) neq(i,i) = 1.d0
C         write(*,520) (neq(i,j),j=1,num_cfiles*num_sat)
C520      format('INV: ',1600E18.10)
C     enddo
C     do i = 1, num_cfiles*num_sat
C         write(*,540) bvec(i,1)
C540      format('EV1: ',1600E18.10)
C     enddo
C     do i = 1, num_cfiles*num_sat
C         write(*,560) bvec(i,2)
C560      format('EV2: ',1600E18.10)
C     enddo

****  Now print results and list statostics
      do i = 1, num_cfiles
         do j = 1, num_sat
            in = (i-1)*num_sat + j
            rms_mwl = sqrt(2*(sum_wgh(in)-sum_cos(in))/
     .                     (sum_wgh(in)+1.d-3))/(2*pi)
            phs_mwl = atan2(sum_sin(in),sum_cos(in))/(2*pi)
            write(*,600) in, cf_codes(i),stat_dcb(i), prn_list(j), 
     .                   bvec(in,1),sqrt(neq(in,in)), sum_num(in),
     .                   rms_mwl, phs_mwl, bvec(in,2)
 600        format('OWBIAS',i5,2x,a4,1x,i1,1x, I2.2,F9.3,' +- ',F7.3,1x,
     .              I6,1x,' MWRMS ',F7.3,1x,F7.3,1x, ' EXWL ',f7.3)
         enddo
      enddo 

****  For the MW-WL case now add estimation of biases on each satellite by
*     receiver type
      do i = 1, num_cfiles+3*num_sat
          apbias(i) = 0.d0
      enddo 

***** Iterate the solution
      done = .false.
      it = 0
      do while ( .not. done )
         it = it + 1
         if( it.gt.5 ) done = .true.  ! quick and dirty add converge test later
         do i = 1, num_cfiles+3*num_sat
             bvcb(i) = 0.d0
             do j = 1, num_cfiles+3*num_sat
                neqb(i,j) = 0.d0
             end do
         end do
*        Loop over one-ways
         do i = 1, num_cfiles
            off = stat_dcb(i)*num_sat
            do j = 1, num_sat
               fr1 = fL1u(j)/fL1(j)
               p = (i-1)*num_sat + j
               np(1) = i
               np(2) = num_cfiles+off+j
*              If we have data
               if( sum_num(p).gt.0 ) then
                   obc = bvec(p,1) - (apbias(np(1))+apbias(np(2)))
* MOD TAH 200617: Added frequency scaling
                   rec = obc-nint(obc*fr1)/fr1
                   tot_wgh = 1.d0/(neq(p,p)+1.d-4)
                   do k = 1,2
                      in = np(k)
                      bvcb(in) = bvcb(in) + rec*tot_wgh
                      do l = 1,2
                         neqb(np(k),np(l)) = neqb(np(k),np(l))+tot_wgh
                      enddo
                   end do
               end if
            end do
         end do

****     Formed the normal equations, now solve system with weak constraint
         nump = num_cfiles+3*num_sat
         do i = 1,nump
            neqb(i,i) = neqb(i,i) + 1
         enddo

         call invert_vis(neqb,bvcb,scale,ipivot,nump,
     .                   max_cfiles+3*max_gprn,1)

****     Update the aprioris, output the results and iterate
         res = 0.d0
         do i = 1, nump
            apbias(i) = apbias(i) + bvcb(i)
            if( abs(bvcb(i)).gt.res ) res = abs(bvcb(i))
         enddo
         write(*,620) it, res
 620     format('OWB',I2.2,' Max Adj ',F8.3)
         if( res.lt.1.d-4 ) done = .true.

         do i = 1, num_cfiles
            write(*,640) it, cf_codes(i), stat_dcb(i), apbias(i), 
     .                   bvcb(i), sqrt(neqb(i,i))
 640        format('OWB',I2.2,' Site ',a4,1x,i1,1x,' WLB ',2F8.3,
     .             ' +- ',F8.3,' cyc')
         enddo
         do i = 0,2
            do j = 1,num_sat
               p = num_cfiles + i*num_sat + j
               write(*,650) it, prn_list(j), i, apbias(p), 
     .                   bvcb(p), sqrt(neqb(p,p))
 650           format('OWB',I2.2,' PRN ',I2.2,1x,i1,4x,' WLB ',2F8.3,
     .             ' +- ',F8.3,' cyc')
            enddo
         enddo

      end do


****  Now form double differences
      do i = 1, namb
         np(1) = (amb_tab(1,i)-1)*num_sat + amb_tab(3,i)
         np(2) = (amb_tab(2,i)-1)*num_sat + amb_tab(3,i)
         np(3) = (amb_tab(1,i)-1)*num_sat + amb_tab(4,i)
         np(4) = (amb_tab(2,i)-1)*num_sat + amb_tab(4,i)
         DDCov = 0.d0
         DDSmw = 0.d0
         DDSex = 0.d0
         do j = 1,4
            DDSmw = DDsmw + ap(j)*bvec(np(j),1)
            DDSex = DDsex + ap(j)*bvec(np(j),2)
            do k = 1,4
               DDCov = DDCov + ap(j)*neq(np(j),np(k))*ap(k)
            end do
         enddo
         in = amb_tab(1,i)
         jn = amb_tab(2,i)
         iv = stat_dcb(in)*num_sat + num_cfiles + amb_tab(3,i)
         iw = stat_dcb(in)*num_sat + num_cfiles + amb_tab(4,i)
         jv = stat_dcb(jn)*num_sat + num_cfiles + amb_tab(3,i)
         jw = stat_dcb(jn)*num_sat + num_cfiles + amb_tab(4,i)
         DDsbs = DDSmw - ((apbias(iv)-apbias(jv))-
     .                    (apbias(iw)-apbias(jw)))
         k = bl_ent(amb_tab(1,i),amb_tab(2,i))
         write(*,820) cf_codes(amb_tab(1,i)),
     .                cf_codes(amb_tab(2,i)),
     .                prn_list(amb_tab(3,i)),
     .                prn_list(amb_tab(4,i)),
     .                -DDSmw, sqrt(DDCov), 
     .                stat_dcb(in),stat_dcb(jn), -DDsbs,
     .                -DDSex, sqrt(DDCov),
     .                baselens(k)/1000.d0

 820     format('ESBIAS',4x,'B1L2 ',a4,1x,a4,1x,
     .            i2.2,1x,i2.2,
     .          ' MW EST ',F9.3,' +- ',F7.3,
     .          ' DCB ',2i2,' MW BIAS ',F9.3,
     .          ' EX EST  ',F9.3,' +-  ',F7.3,
     .          ' LEN ',F8.2,' km')
      end do
      call report_stat('status','autcln','EST_WLFULL',' ','Done',0)

*     Thats all
      return
      end
 
CTITLE IMPLICIT_BF
 
      subroutine implicit_bf(pn, nump, dim, neq, bvec, nb)
 
*     Routine a make an implicit solution for the bias parameters
*     similar to the way solve handles bias parameters.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
 
* PASSED VARIABLES
 
* pn   -- Parameter number to be solved
* nump -- number of parameters
* dim  -- Dimension of matrices
* nb   -- Number of terms in bvec

      integer*4 pn, nump, dim, nb 

* neq, bvec -- Normal equations and bvec
      real*8 neq(dim, dim), bvec(dim,nb)


 
* LOCAL VARIABLES
 
* i,j -- Loop variabales

      integer*4 i,j, k

*     Do the implicit solutions if values already assigned
C      if( neq(pn,pn).eq.0.d0 ) RETURN
      if( neq(pn,pn).lt.1.d-09 ) then 
         if( neq(pn,pn).ne.0.d0 ) then
            write(*,110) pn, neq(pn,pn), (bvec(pn,j),j=1,nb)
 110        format('IMPLICIT ERR: Small non-zero diagonal NP ',i4,
     ,             ' Diagonal ',E14.5,' BVEC ',10E14.5)
            do k = 1, nb
                bvec(pn,k) = 0.d0
            enddo
            do j = 1,nump
               neq(j,pn) = 0.d0
               neq(pn,j) = 0.d0
            end do
         endif
         RETURN
      endif

      do j = 1, nump
          if( j.ne.pn ) then 
             do k = 1, nb 
                bvec(j,k) = bvec(j,k) - 
     .                      bvec(pn,k)*neq(pn,j)/neq(pn,pn)
             enddo
             do i = 1, nump
                if( i.ne.pn ) then
                    neq(i,j) = neq(i,j) -
     .                         neq(i,pn)*neq(pn,j)/neq(pn,pn)
                end if
             end do
          end if
      end do              

 
*     Now remove clear the row and column for parameter pn
      do k = 1, nb
        bvec(pn,k) = 0.d0
      enddo
      do j = 1,nump
         neq(j,pn) = 0.d0
         neq(pn,j) = 0.d0
      end do

****  Thats all
      return
      end
 

