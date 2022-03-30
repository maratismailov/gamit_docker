      subroutine add_bfent( ep, ob, na )

      implicit none

*     Routine to add a new bf_ents entry.  The new entry is added
*     at the end of the site list and then all entries above are
*     moved up.

      include '../includes/const_param.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep    ! Counter for number of epochs of data processed.
      integer*4 ob    ! Observation number that needs entry
      integer*4 na    ! Ambiguity number added

* LOCAL
      integer*4 i,j, k   ! Loop variables


****  Make sure we don't overflow arrays
      if( num_ambs+1.gt.max_ambs ) then
          call report_stat('fatal','trackrt','get_blags','',
     .       'Too many bias flags Max_ambs', max_ambs)
      end if 

****  See where this entry fits
      if( RT_siteNum(ob).eq.num_site ) then
         na = num_ambs + 1
      else
         na = 0
         do i = 1, num_ambs
            if( bf_ents(1,i).gt.RT_siteNum(ob) ) then
               na = i
               exit
            endif
         end do
*        If we did not find a value then we have not
*        used this site yet so just add to the end
         if( na.eq.0 ) na = num_ambs + 1
      end if

****  See if we need to move entries up
      do j = num_ambs, na, -1
         do k = 1,6
            bf_ents(k,j+1) = bf_ents(k,j)
         end do
         do k = 1,2
            ambiq_all(k,j+1) = ambiq_all(k,j)
         end do
         do k = 1,3
            wls_ref(k,j+1) = wls_ref(k,j)
         end do
         do k = 1,2
            wls_sum(k,j+1) = wls_sum(k,j)
            wls_sqr(k,j+1) = wls_sqr(k,j)
         end do
         
         wls_num(j+1) = wls_num(j)
         amb_parn(1,j+1) = amb_parn(1,j)
         amb_parn(2,j+1) = amb_parn(2,j)

         curr_sdmw(j+1) = curr_sdmw(j)
         curr_sdex(j+1) = curr_sdex(j)

****     Now move the saved quanities
         asv_res(:,j+1)  = asv_res(:,j)
         asv_sig(:,j+1)  = asv_sig(:,j)
         asv_chi(:,j+1)  = asv_chi(:,j)
         asv_rbn(:,j+1)  = asv_rbn(:,j)
         asv_dl12(:,j+1) = asv_dl12(:,j)
         asv_resep(j+1)  = asv_resep(j) 
         asv_used(:,j+1) = asv_used(:,j)
         asv_fcode(j+1)  = asv_fcode(j)

      end do

*     Now update the observation pointers to ambiguities
      do j = sblk, eblk
         if( RT_ambpnt(j).ge.na ) then
             RT_ambpnt(j) = RT_ambpnt(j) + 1
         endif
      end do

*     Now add and initialize the new entry
      bf_ents(1,na) = RT_sitenum(ob)
      bf_ents(2,na) = RT_satNum(ob)
      bf_ents(3,na) = ep     ! Start epoch
! MOD TAH 120613: Initialize to max_gap so that value will also
!     be set for first entries.
      bf_ents(6,na) = ep-max_gap-1      ! Last end epoch
      bf_ents(4,na) = ep-max_gap-1      ! End epoch

      if( RT_sitenum(ob).eq.1 ) then  ! Reference site 
         bf_ents(5,na) = 2      ! Set bit to say resolved
      else
         bf_ents(5,na) = 5      ! Set to show slip since we just added
      end if
      ambiq_all(1,na) = 0 
      ambiq_all(2,na) = 0
*     Initialize the pointers to the double differences
*     that will be used to compute these ambiguities
      wls_ref(1,na) = 0      ! Same site, other satellite 
      wls_ref(2,na) = 0      ! Other site, same satellite
      wls_ref(3,na) = 0      ! Other site, other satellite
      wls_sum(1,na) = 0
      wls_sum(2,na) = 0
      wls_sqr(1,na) = 0
      wls_sqr(2,na) = 0
      wls_num(na) = 0

***** Update the number of ambiquities
      num_ambs = num_ambs + 1

      amb_parn(1,na) = -1       ! Initially set the parameter numbers
      amb_parn(2,na) = -1       ! to -1 to show not set yet

       
****  Thats all
      return
      end




