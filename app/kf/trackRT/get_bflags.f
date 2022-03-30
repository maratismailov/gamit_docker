      subroutine get_bflags( ep )

      implicit none

*     Routine to add in the estimates of the phase offsets to the 
*     phase data.  New biases are added as needed.

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep            ! Counter for number of epochs of data processed.

* LOCAL
      integer*4 i, j, k, l, m   ! counter
     .,     is, js, ls, ks      ! Pointers from the sorted data
     .,     na, nb, sdr  ! Ambiquity number and single difference reference
                      ! obs number
     .,     max_ob(max_site)    ! Obs number with maximum elevation
     .,     max_sv(max_site)    ! Satellite number with highest elevation
     .,     cs        ! Current station for DD (normally 1 but can be
                      ! higher if reference site missing
     .,     ref_ob(max_site)  ! Observation number with reference satellite
     .,     ref_am(max_site)  ! Ambiquity number for reference satellite
                      ! at each site (saved to wls_ref(1,na) -- to give 
                      ! reference satellite
     .,     num_sdat(max_site) ! Number of good obs per site (needed if no
                      ! data at a site)
     .,     np        ! Parameter number
     .,     bit_stat  ! Set 0 or 1 depending if a new reference satellite ambiquities
                      ! are resolved.

      real*8 dmw, dex  ! MW and EX widelane differences for single differences
     .,      ddmw, ddex ! Double difference MW and EX WL
     .,      ddob(4)    ! Double diffence phase and range residuals
     .,      eNL1, eNL2 ! Estimates of cycles at L1 and L2 from double diffs
     .,      max_el(max_site)     ! Maximum elevation of satellite (degs)
     .,      ddres      ! Double difference residual for testing

      integer*4 iNL1, iNL2 !  Integer ambiguity

      logical data_OK, kbit
      logical remap    ! Set true if we need to remap 
      logical DDOK     ! Set true if data in DD OK
      integer*4 nsvOK  ! number of sites, except reference site, at which a
                       ! satellite is available (used for DD reference SVS)
     .,         ng     ! Number of stations expeced (reduced for each one with
                       ! out data at an epoch)

****  Scan to get the pointer to bf_ent ambiquity of each one-way
      data_mask = -1  ! Set for the moment

      do i = sblk, eblk
         na = 0
         do j = 1, num_ambs
            if( bf_ents(1,j).eq. RT_sitenum(i) .and.
     .          bf_ents(2,j).eq. RT_satNum(i) ) then
                na = j
            end if
         end do
*        See if found
         if( na.eq.0 ) then
            call add_bfent( ep, i, na )
         end if
         RT_ambpnt(i) = na
         if( data_OK(RT_errflag(i),data_mask) ) then
             bf_ents(6,na) = bf_ents(4,na) ! Save last end epoch
             bf_ents(4,na) = ep     ! End epoch
         end if
      end do

*     OK, now apply cycle ambiquities to the phase data
      do i = sblk, eblk
         na =  RT_ambpnt(i)
         RT_omcraw(1,i) = RT_omcraw(1,i) - ambiq_all(1,na)
         RT_omcraw(2,i) = RT_omcraw(2,i) - ambiq_all(2,na)
C        write(*,120) ep, site_names(RT_sitenum(i)), 
C    .          RT_satnum(i), na, ambiq_all(:,na), BF_ents(5,na)
C120     format('AMBAP ',i6,1x,a4,' G ',I2.2,' Amb ',i4,2F12.1,
C    .          ' FLAG ',o6)
      end do

*     Clear the pointers to data for each ambiquity (not all values
*     will be set).
      do i = 1, num_ambs
         do j = 1,4
            wls_obn(j,i) = 0
         end do
      end do 

*     Loop over the sites and see if WL's match
      do i = 2, num_site
         do j = sblk, eblk
            js = RT_ord(j)  ! Loop in site order  
            if( rt_sitenum(js).eq.i ) then
*               Site matches, now try to find this satellite
*               at the reference or eariler station
                sdr = 0
                do k = 1, 1   ! In Version 1.00, only check refsite (i-1) 
                   do l = sblk, eblk
                       ls = RT_ord(l)
                       if( sdr.eq.0 .and. 
     .                     RT_satNum(js).eq. RT_satNum(ls) .and.
     .                     RT_siteNum(ls).eq.k .and.
     .                     data_OK(RT_errflag(ls),data_mask) .and.
     .                     data_OK(RT_errflag(js),data_mask) ) then
*                          Found match of satellite which we had
*                          not found before with earlier site
                           sdr = ls
                           exit   ! Exit from this loop
                       end if
                   end do
                   if( sdr.ne.0 ) exit
                end do
*               If sdr is zero, then no match found for satellite,
*               we can not form single difference 
                if( sdr.eq.0 ) then
*                  Removed setting this bit so that data can be used
*                  in other differeces (setting this but will force
*                  base station to be used for all data).
                   if( i.eq.num_site ) then ! No more data so mark bad
                       call sbit(RT_errflag(js),3,1) 
                   end if
                   call sbit(RT_errflag(js),3,1) 

                   call sbit(bf_ents(5,RT_ambpnt(js)),3,0)
                   if( debug(9).gt.0 ) 
     .             write(*,140) site_names(i), RT_satnum(js), 
     .                j, js, RT_ambpnt(js), bf_ents(5,RT_ambpnt(js)) 
 140               format('NO SDWL Site ',a4,' G ',I2.2,' J/Js ',2I5,
     .                 ' AmbNum ',i4,' BF_ENTS(5) ',o4)
                   call sbit(bf_ents(5,RT_ambpnt(js)),3,0)
                else
                   call sbit(RT_errflag(js),3,0)
*                  OK, see how big a jump we have
                   ks = sdr
                   na = RT_ambpnt(js)
                   dmw = (RT_omcraw(1,js)- RT_omcraw(2,js)-
     ,                   dfsf*(RT_omcraw(3,js)*fL1+RT_omcraw(4,js)*fL2)/
     .                                                    vel_light) -
     .                   (RT_omcraw(1,ks)- RT_omcraw(2,ks)-
     ,                   dfsf*(RT_omcraw(3,ks)*fL1+RT_omcraw(4,ks)*fL2)/
     .                                                    vel_light)
                   dex = (RT_omcraw(1,js)- RT_omcraw(2,js)*(fL1/fL2))-
     .                   (RT_omcraw(1,ks)- RT_omcraw(2,ks)*(fL1/fL2))

                   do m = 1,4
                      curr_sdob(m,na) = RT_omcraw(m,js)-RT_omcraw(m,ks)
                   end do

*                  See what we need to compare to:
                   if( wls_ref(2,na).ne.0 ) then 
*                      OK, we have been using a station in the past
*                      See if it has changed and if it has can be map
*                      the old station to the new one.
                       if( wls_ref(2,na).ne. RT_ambpnt(ks) ) then
*                         Station has changed and wls_ref(2,na) points
*                         to the ambiguitity that was previous used as
*                         single difference reference.  
                          do l = 1,num_ambs
                             if( bf_ents(1,l).eq.RT_sitenum(ks) .and.
     .                           bf_ents(2,l).eq.RT_satnum(ks) ) then
*                                We have found offset values, so 
*                                apply.  (Later add more stats here)
                                 curr_sdmw(na) = curr_sdmw(na) -
     .                                          curr_sdmw(l)
                                 curr_sdex(na) = curr_sdex(na) -
     .                                          curr_sdex(l)
                                 curr_dd(:,na) = curr_dd(:,na) -
     .                                           curr_dd(:,l)
                                 do m = 1,4
                                    curr_sdob(m,na) = curr_sdob(m,na) -
     .                                                curr_sdob(m,l)
                                 end do
                             end if
                          end do
                       end if
                   end if
*                  Save the Station single difference amb pointers
                   wls_ref(2,na) = RT_ambpnt(ks)
                   wls_obn(1,na) = js
                   wls_obn(2,na) = ks

*                  See if we have a jump in either the MW or EX widelanes
* MOD TAH 120105:  See also if there is a gap larger than max_gap
                   if( abs(dmw-curr_sdmw(na)).gt. mwwl_jmp .or.
     .                 abs(dex-curr_sdex(na)).gt. exwl_jmp .or.
     .                 ep-bf_ents(6,na).gt. max_gap ) then
*                      Jump on the widelanes, mark as such for this
*                      ambiguity
*                      Report slip
                       if( abs(dmw-curr_sdmw(na)).gt. mwwl_jmp .or.
     .                     abs(dex-curr_sdex(na)).gt. exwl_jmp ) then 
                          write(*,180) ep, site_names(RT_sitenum(js)), 
     .                      RT_satnum(js), abs(dmw-curr_sdmw(na)), 
     .                      mwwl_jmp, abs(dex-curr_sdex(na)), exwl_jmp,
     .                      RT_azel(2,js), ep-bf_ents(6,na)
                          write(lus,180) ep, site_names(RT_sitenum(js)), 
     .                      RT_satnum(js), abs(dmw-curr_sdmw(na)), 
     .                      mwwl_jmp, abs(dex-curr_sdex(na)), exwl_jmp,
     .                      RT_azel(2,js), ep-bf_ents(6,na)
 180                      format('CSLIP Ep ',i6,' Site ',a,' G ',I2.2,
     .                      ' DMW-WL/Tol ',F12.2,1x,F5.2,
     .                      ' DEX-WL/Tol ', F12.2,1x,F5.2,' cyc, Elev',
     .                      F6.2,' deg. Gap ',I7,' epochs')
                       else  ! Write message for Data gap
                          write(*,185) ep, site_names(RT_sitenum(js)), 
     .                      RT_satnum(js), ep-bf_ents(6,na),
     .                      abs(dmw-curr_sdmw(na)), 
     .                      abs(dex-curr_sdex(na)), RT_azel(2,js)
                          write(lus,185) ep, site_names(RT_sitenum(js)), 
     .                      RT_satnum(js), ep-bf_ents(6,na),
     .                      abs(dmw-curr_sdmw(na)), 
     .                      abs(dex-curr_sdex(na)), RT_azel(2,js)
 185                      format('CSLIP Ep ',i6,' Site ',a,' G ',I2.2,
     .                      ' Data Gap ',i7,' epochs, DMW-WL ',F12.2,
     .                      ' DEX-WL ', F12.2,' cyc, Elev',F6.2,' deg')
                       endif

*                      Just update for this site
                       call sbit(bf_ents(5,na),3,1)
                       call sbit(bf_ents(5,na),2,0)  ! Set resolved bit to 0
                       bf_ents(3,na) = ep

*                      Reset the parameter estimates for these parameter
                       do m = 1,2
                          if( amb_parn(m,na).eq.0 ) then
                             amb_parn(m,na) = -1
                          end if
                       end do
                   else
                       call sbit(bf_ents(5,na),3,0)
                   endif

                   if( ep.ge.debug(5) .and. ep.le.debug(6) ) 
     .             write(*,220) ep, site_names(RT_siteNum(js)), 
     .                    site_names(RT_siteNum(ks)), RT_satNum(js),
     .                    RT_satNum(ks), na, dmw-curr_sdmw(na), 
     .                    dex-curr_sdex(na), curr_sdmw(na), 
     .                    curr_sdex(na), RT_azel(2,js), 
     .                    bf_ents(1,wls_ref(2,na)), 
     .                    RT_errflag(js), RT_errflag(ks),bf_ents(5,na) 
 220               format('RT_SDWL ',i6,1x,a4,1x,a4,' G ',I2.2,1x,I2.2,
     .                    1x,i4,4(1x,F15.2),1x,F8.2,1x,i2,1x,3o4)

*                  Save current settings
                   curr_sdmw(na) = dmw
                   curr_sdex(na) = dex

                end if
             end if
         end do
      end do

***** Now check the magnitude of double difference residuals to see 
*     if we could have a cycle slip that did not show in the ex-wl
*     and mw-wl
      do i = 2, num_site
*        Loop on the data from this site
         nb = 0
         do j = ssblk(i), seblk(i)
            if( RT_satnum(j).eq.ref_sv(i) ) then
                nb = RT_ambpnt(j)
                exit
            end if
         end do
         do j = ssblk(i), seblk(i)
*           Only test if this ambiguity is fixed.
            na = RT_ambpnt(j)
*           See if data is OK: 
            ddOK = data_OK(RT_errflag(j),data_mask) .and.
     .             data_OK(RT_errflag(wls_obn(2,na)),data_mask) .and.
     .             data_OK(RT_errflag(wls_obn(1,nb)),data_mask) .and.
     .             data_OK(RT_errflag(wls_obn(2,nb)),data_mask) 
            if( ddOK .and. 
     .          ref_sv(i).ne.RT_satNum(j) .and. ref_sv(i).gt.0 ) then
*               See if we can form the double difference
                do m = 1,4
                    ddob(m) = curr_sdob(m,na)-curr_sdob(m,nb)
                end do
*               Form the LC double difference
                ddres = lcf1*ddob(1)+lcf2*ddob(2)
                if( ep.ge. debug(3) .and. ep.le.debug(4) ) then
                    write(*,225) ep, site_names(i),  RT_satnum(j),
     .                   na, nb, ddres-curr_dd(1,na), 
     .                   ddres,curr_dd(1,na)
 225                format('DD RESID Ep ',i6,1x,a4,1x,' G ',I2.2,
     .                  ' AMB #s ',2i4,' DRes ',F8.2,' RES/Curr ',
     .                  2F8.2,' cyc')
                end if
*               See if value too large and ambiguity fixed.  Only
*               check after the first epoch after the value has
*               been fixed (really should update curr_dd when 
*               ambihuity is resolved: FIX LATER
                if( abs(ddres-curr_dd(1,na)).gt. dd_jmp .and.
     .              abs(ddres).gt. dd_jmp .and.
     .              kbit(bf_ents(5,na),2) .and. 
     .              ep.gt.asv_resep(na)+1  ) then
*                   Jump on the widelanes, mark as such for this
*                   ambiguity
                    write(*,230) ep, site_names(i), 
     .                   RT_satnum(j), ddres, dd_jmp, RT_azel(2,j)
                    write(lus,230) ep, site_names(i), 
     .                   RT_satnum(j), ddres-curr_dd(1,na), dd_jmp, 
     .                   RT_azel(2,j)
 230                format('CSLIP Ep ',i6,' Site ',a,' G ',I2.2,
     .                   ' DD RESID ',F12.2,1x,' Tol ',F5.2,
     .                   ' cyc, Elev ', F6.2,' deg')

*                   Just update for this site
                    call sbit(bf_ents(5,na),3,1)
                    call sbit(bf_ents(5,na),2,0)  ! Set resolved bit to 0
                    bf_ents(3,na) = ep
*                   Reset the parameter estimates for these parameter
                    do m = 1,2
                       if( amb_parn(m,na).eq.0 ) then
                          amb_parn(m,na) = -1
                       end if
                    end do
                 end if
*                Save current DD
                 curr_dd(1,na) = ddres
            end if
         end do
      end do

****  Before testing for a double difference satellite, see if
*     any stations have no data and therefore will not be included
*     in the test
      ng = 0   ! Number of sites with data
      do i = 1, num_site
         num_sdat(i) = 0
         do j = sblk,eblk
            js = RT_ord(j)
            if( RT_sitenum(js).eq.i .and.
     ,          data_OK(RT_errflag(js),data_mask) ) then
               num_sdat(i) = num_sdat(i) + 1
            elseif (RT_sitenum(js).gt.i ) then  ! Finished with this site
               exit
            endif
         enddo
         if( num_sdat(i).gt.0 ) ng = ng + 1   ! Increment number of sites
      enddo 

****  OK we have scanned all the single differences, not form double
*     differences and increment estimates of cycle offsets
*     Scan to see which satellite to use as reference
      do i = 1,num_site
         max_el(i) = 0
         max_ob(i) = 0
         max_sv(i) = 0
      end do

      cs = 1
*     Try this first to see if we can make the reference satellite
*     one with ambiguities resolved. 
      do i = ssblk(1), seblk(1)   ! Range for site 1
*        See if this is highest
         is = i ! RT_ord(i)
         na = RT_ambpnt(is)
*        Test against other stations
         do k = 2, num_site 
             if( RT_azel(2,is).gt. max_el(k) .and.
     .          data_OK(RT_errflag(is),data_mask) ) then
*               Scan the our site and make sure data are available
                nsvOK = 0
                do j = ssblk(k), seblk(k)
                   js = j ! RT_ord(j)
                   nb = RT_ambpnt(js)
                   if( RT_satnum(js).eq.RT_satnum(is) .and.
     .                 data_OK(RT_errflag(js),data_mask) .and. 
     .                 kbit(bf_ents(5,nb),2) ) then
*                      Valid measurement to check
                       nsvOK = 1
                   endif
                end do
                if( nsvOK.eq.1 ) then
                   max_el(k) = RT_azel(2,is)
                   max_ob(k) = is
                   max_sv(k) = RT_satNum(is)
                endif
            endif
         end do
      end do
      do k = 2,num_site
         if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .   write(*,260) 'REF AMBFIX ', ep, max_sv(k), ref_sv(k), 
     .       max_el(k), max_ob(k), site_names(k) 
      end do

c      print *,'REF AMBFIX ',ep, max_ob, ' G ', max_sv, ref_sv, max_el, 
c     .      ' Ref ', site_names(RT_siteNum(max_ob))
 260  format(a,' Ep ',i6,' New Ref/Old Ref G ',I2.2,' G ',I2.2,
     .       ' Elev ',F6.2,' deg, max_ob ',I4,1x,a)

*     If max_ob is zero still, it means we could not find a reference
*     satellite with ambiguities resolved at all sites, so now just
*     try to find the highest satellite that all sites can see.
      do k = 2,num_site
         if( max_ob(k).eq.0 ) then   ! Scan again with out need for 
                                     ! being bias fixed.
             do i = ssblk(1), seblk(1)   ! Range for site 1
*                See if this is highest
                 is = i ! RT_ord(i)
                 na = RT_ambpnt(is)
*                Test against other stations
                 if( RT_azel(2,is).gt. max_el(k) .and.
     .              data_OK(RT_errflag(is),data_mask) ) then
*                   Scan the our site and make sure data are available
                    nsvOK = 0
                    do j = ssblk(k), seblk(k)
                       js = j ! RT_ord(j)
                       nb = RT_ambpnt(js)
                       if( RT_satnum(js).eq.RT_satnum(is) .and.
     .                     data_OK(RT_errflag(js),data_mask)  ) then
*                          Valid measurement to check
                           nsvOK = 1
                       endif
                    end do
                    if( nsvOK.eq.1 ) then
                       max_el(k) = RT_azel(2,is)
                       max_ob(k) = is
                       max_sv(k) = RT_satNum(is)
                    endif
                 endif
             enddo
             if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .       write(*,260) 'REF AMBFRE ', ep, max_sv(k), ref_sv(k),  
     .           max_el(k), max_ob(k), site_names(k) 
         endif
      end do

****  Make sure we have reference satellite
      nsvOK = 0
      do k = 2,num_site
         if( max_ob(k).ne.0 ) nsvOK = nsvOK + 1
      end do
      if( nsvOK.eq. 0 ) then
*        We cannot find any double differences
         call report_stat('WARNING','TRACKRT','get_bflags',
     .                 ' ','No double differences Epoch ',ep)

         write(lus,300) ep
 300     format('WARNING: No double differences with reference site ',
     .          'at Epoch ',i6)
         return
      end if

***** Find the ambiquity number at each site associated with the
*     reference satellite
      do i = 1, num_site
         ref_am(i) = 0
         do j = 1, num_ambs
            if( bf_ents(1,j).eq.i .and. bf_ents(2,j).eq.ref_sv(i) ) then
                ref_am(i) = j
                exit
            endif
         end do
      end do

****  See if we should mark the old reference satellite as fixed.  We
*     can do this if the new one is fixed
      do i = 2, num_site
         if( ref_sv(i).ne.max_sv(i) ) then 
            if( ref_am(i).gt.0 .and. kbit(bf_ents(5,ref_am(i)),2) ) then
*               New reference is fixed, so mark the old reference satellite
*               as fixed
                bit_stat = 1   ! Set as fixed
            else
                bit_stat = 0
            endif

****        If there was a cycle slip on the old reference, turn off its
*           bias fixed status
            if( ref_am(i).gt.0 .and. kbit(bf_ents(5,ref_am(i)),3) ) then
                call sbit(bf_ents(5,ref_am(i)),2,0)
                if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .          write(*,320) ep, site_names(i),bf_ents(2:6,ref_am(i)),
     .                   ambiq_all(:,ref_am(i))
 320            format('RESET OLD REF Ep ',i6,1x,A4,1x,' G ',I2.2,
     .                   ' BF_ENT(3:6) ',2I7,1x,o2,1x,I6,' Amb ',
     .                   2F15.1)
            end if

****        Now mark the status
            do j = 1, num_ambs
               if( bf_ents(1,j).eq.i .and. 
     .             bf_ents(2,j).eq.max_sv(i) ) then
                  call sbit(bf_ents(5,j),2,bit_stat)
                  if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .            write(*,330) ep, site_names(i),j, bf_ents(:,j),
     .                 bit_stat
 330              format('Mark BF_REF Ep ',i6,1x,a,' na ',i3,
     .                  ' BF_ENTS ',4i6,1x,o3,1x,i6,' NBS ',i1)
               end if
            end do
****        
         end if
      end do

****  Now get the ambuigity number for the new reference satellite
      do i = 1, num_site
         ref_am(i) = 0
         do j = 1, num_ambs
            if( bf_ents(1,j).eq.i .and. bf_ents(2,j).eq.max_sv(i) ) then
                ref_am(i) = j
                exit
            endif
         end do
      end do

****  Mark the reference satellite ambiquities as fixed.
      do i = 2, num_site
         if( ref_am(i).gt.0 ) then
            call sbit(bf_ents(5,ref_am(i)),2,1)
         end if
      end do

***** Now for each site find the observation number with the reference
*     satellite
      do i = 1, num_site
         ref_ob(i) = 0
         do j = ssblk(i), seblk(i)
            js = j ! RT_ord(j)
            if( RT_sitenum(js) .eq. i  .and. 
     .          RT_satnum(js).eq. RT_satnum(max_ob(i)) ) then
                ref_ob(i) = js
                exit
           endif
         end do
*        Make sure we have entry 
         if( ref_ob(i).eq.0 .and. i.gt.1 ) then
*           The reference satellite is not at this site
*           Since reference satellite had to be common
*           to all stations -- must be no data
            print *,'No Data Epoch ',ep,' Site ',site_names(i)
         end if
      end do
 
****  OK, now form double differences and use these to update
*     number of cycles.  If we have switched reference satellites
*     map the double differences
      remap = .false.
      do i = 2,num_site
         if( max_sv(i).ne.ref_sv(i) .and. ref_sv(i).gt.0 ) then
            remap = .true.
            exit
         end if
      end do 

      if( remap ) then
         call remapdd(ep, ref_sv, max_sv )

****     Now remap the parameter estimates
         call remap_ambparn( ep, ref_sv, max_sv, ref_am ) 
      end if

****  OK: now reset the reference satellite
      do i = 2,num_site
         ref_sv(i) = max_sv(i)
      end do

****  Mark the reference satellite ambiquities as fixed.

      do i = sblk, eblk
         is = RT_ord(i)
         k = RT_sitenum(is)
         if( k.ne. RT_sitenum(max_ob(k)) .and.
     .       RT_satnum(is).ne.ref_sv(k) .and.
     .       ref_ob(k).gt.0 .and.
     .       data_OK(RT_errflag(is),data_mask) )  then
*            OK difference the single difference arrays
             na = RT_ambpnt(is)
             nb = RT_ambpnt(ref_ob(RT_sitenum(is)))
             ddmw = curr_sdmw(na)-curr_sdmw(nb)
             ddex = curr_sdex(na)-curr_sdex(nb)
*            Set the pointers for this double diffence
             wls_obn(3,na) = wls_obn(1,nb)
             wls_obn(4,na) = wls_obn(2,nb)
             wls_ref(1,na) = ref_am(RT_sitenum(is))

             if( kbit(bf_ents(5,RT_ambpnt(is)),3) )  then
*               Cycle slip detected, so reset number
*               of cycles
                call res_wlsRT( ddmw, ddex, eNL1, eNL2, iNl1, iNL2)
!               write(*,510) na, bf_ents(2,na), ddmw, ddex, eNL1, eNL2,
!    .                   iNl1, iNL2, ambiq_all(:,na)
!510            format('RES_WLSRT NA ',i4,' G',I2.2,' DWLS ',2F10.3,
!    .                 ' eN12 ', 2F10.2,' iN12 ',2I10,
!    .                 ' Curr ambiq ',2f12.1)
 
*               Reset the initial epoch counter
                na = RT_ambpnt(is)

                ambiq_all(1,na) = ambiq_all(1,na)+iNL1
                ambiq_all(2,na) = ambiq_all(2,na)+iNL2
                RT_omcraw(1,is) = RT_omcraw(1,is)-iNL1
                RT_omcraw(2,is) = RT_omcraw(2,is)-iNL2
                curr_sdmw(na) = curr_sdmw(na) - (iNL1-iNL2)
                curr_sdex(na) = curr_sdex(na) - (iNL1-iNL2*fL1/fL2) 
                curr_sdob(1,na) = curr_sdob(1,na) - iNL1
                curr_sdob(2,na) = curr_sdob(2,na) - iNL2 
                ddmw = ddmw-(iNL1-iNL2)
                ddex = ddex-(iNL1-iNL2*fL1/fL2)

*               Now update the WL estimated values
C               call update_wlsstat(ep, iNL1, iNL2, na)
                call acc_wlsstat(ep, ddmw, ddex, na,'Reset')

*               Reset the parameter estimate associated with this ambiqity.
                call reset_ambest(ep, na) 

                if( ep.ge.debug(5) .and. ep.le. debug(6) ) 
     .          write(*,520) ep, site_names(RT_sitenum(is)),
     .             RT_satNum(is),ambiq_all(:,na),
     .             eNL1-iNL1, eNL2-iNL2
     .,            curr_sdmw(na),curr_sdex(na) 
 520            format('RT_AMB Epoch ',I6,1x,a4,' G ',i2.2,
     .             ' NL1/2 ',2(1x,F13.1),1x,' dNL12 ',2F8.3,
     .             ' Curr ',2(F13.1,1x))
             else

                call acc_wlsstat(ep, ddmw, ddex, na,'Sum')

             endif

             do m = 1,4
                ddob(m) = curr_sdob(m,na)-curr_sdob(m,nb)
             end do

*            If there is a cycle slip here, reset the 
*            number of cycles
             k = RT_sitenum(is)
             if( ep.ge.debug(5) .and. ep.le. debug(6) )
     .       write(*,410) ep, site_names(k),
     .            site_names(RT_sitenum(max_ob(k))),
     .            RT_satnum(is),RT_satnum(max_ob(k)),
     .            RT_ambpnt(is),ddmw,ddex, 
     .            RT_azel(2,is), RT_azel(2,ref_ob(k)),
     .            bf_ents(1:5,RT_ambpnt(is))  
  410        format('RT_DDWL ',i6,1x,a4,1x,a4,' G ',i2.2,1x,
     .            I2.2,1x,I4,1x,2(F15.2,1x),2(1x,F6.2),
     .            ' BF_ENTS ',2I3,2i7,1x,o3)
             if( ep.ge.debug(5) .and. ep.le. debug(6) )
     .       write(*,440) ep, site_names(RT_sitenum(is)),
     .            site_names(RT_sitenum(max_ob(k))),
     .            RT_satnum(is),RT_satnum(max_ob(k)),
     .            RT_ambpnt(is), (ddob(m),m=1,4),
     .            RT_azel(2,is), RT_azel(2,ref_ob(k)),
     .            bf_ents(1:5,RT_ambpnt(is)), wls_obn(:,RT_ambpnt(is))
  440        format('RT_DDOB ',i6,1x,a4,1x,a4,' G ',i2.2,1x,
     .            I2.2,1x,I4,1x,4(F15.2,1x),2(1x,F6.2),
     .            ' BF_ENTS ',2I3,2i7,1x,o3,
     .            ' OBN ',4I4)

C            write(*,415) ep, site_names(RT_sitenum(i)),
C    .            site_names(RT_sitenum(max_ob)),
C    .            i, max_ob, curr_sdmw(na),curr_sdmw(nb),
C    .            curr_sdex(nb),curr_sdex(nb)
C415         format('RT_DSEN ',i6,1x,a4,1x,a4,' # ',2i4,
C    .            4F15.2)
          endif
      end do

****  Report status ambiqiuities.
      call rep_wlsstat(ep)

 
***** That all 
      return 
      end





      


