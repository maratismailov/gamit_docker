CTITLE SAVE_GLB_APR
 
      subroutine save_glb_apr

      implicit none
 
 
*     Routine to save some of the apriori information from the current
*     solution in the globk_ema common areas.  These values are later
*     used for output.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
 
*   i       - Loop counter
 
      integer*4 i, j, k, ne
 
****  Wobble, UT1 and nutation angles
      do i = 1,8

*         Only update if we do not have an apriori value at this time. 
* MOD TAH 131224: Added test for gwob_apr(i).eq.0.0 to assign value
          if( gwob_apr(i).eq. -9999.d0 .or. gwob_apr(i).eq.0.0 ) 
     .                               gwob_apr(i) = cwob_apr(i)
          gnut_ang_apr(i) = cnut_ang_apr(i)
 
      end do
 
      do i = 1,6
* MOD TAH 131224: Added test for gut1_apr(i).eq.0.0 to assign value
          if( gut1_apr(i).eq. -9999.d0 .or. gut1_apr(i).eq.0.0 ) 
     .                               gut1_apr(i) = cut1_apr(i)
      end do
 
      gtai_utc  = ctai_utc
      ggpst_utc = cgpst_utc 

****  Now see if we multiple pmu epoch.  Unlike other parameters
      if( cent_apr_ep.gt.0 ) then 

****      Get the multi-pmu aprioris from the hfile.  We then check if 
*         we need to use any of these values
          call read_mul_pmu
*         print *,'DEBUG: read_mul_pmu ',cepoch_expt, ne,
*    .          cwob_apr(1:4), gwob_apr(1:4), cent_apr_ep
*         print *,'DEBUG: cmul_pmu ',  cmul_pmu_apr(:,:,1)
*         We now need to get apriori values from the current 
*         hfiles if we do not yet have values.  (These values
*         could be replaced later if we are using an in_pmu file).
*         Loop over the aprioris needed and see if we have values
*         match.
          do i = 1,3
             do j = 1,2
                do k = 1, num_mul_pmu

*                  See if value of apr_val_mul_pmu is set to default
*                  indicating value is unknown. 
                   if( apr_val_mul_pmu(j,i,k).eq. -9999.d0 ) then
*                      Value is unknown.  See if we can get value 
*                      from this h-file.  To do this we need the
*                      apriori codes and get the epoch numbers from
*                      these codes.

*                      Loop over the apriori code records until we
*                      find what we need
                       call assign_mul_ep( j, i, k )
                   end if
                end do
             end do
          end do
* MOD TAH 140107: See if we need to set the gwob/gut1 values for this
*         day.  If not set, they will be zero.
*         print *,'DEBUG: Assign GWOB ',cepoch_expt,gepoch_expt,
*    .         cnum_mul_pmu(1,1), cmul_pmu_ep(1),gwob_apr(1:2)
          do k = 1, cnum_mul_pmu(1,1)
              if( abs(cepoch_expt- cmul_pmu_ep(6*(k-1)+1)).lt.
     .            tol_mul_pmu .and. 
     .            abs(gepoch_expt-cepoch_expt).gt.
     .            tol_mul_pmu ) then
                 gwob_apr(1) =  cmul_pmu_apr(1,1,k)
                 gwob_apr(2) =  cmul_pmu_apr(1,2,k)
                 gwob_apr(3) =  cmul_pmu_apr(2,1,k)
                 gwob_apr(4) =  cmul_pmu_apr(2,2,k)
                 gut1_apr(1) =  cmul_pmu_apr(1,3,k)
                 gut1_apr(2) =  cmul_pmu_apr(2,3,k)
             end if
          end do
*         print *,'DEBUG: Done  GWOB ',k, cepoch_expt,
*    .         cmul_pmu_ep(6*(k-1)+1), gwob_apr(1:4)
 
      else 

*         Now multiday entries in file.  Save the current values
*         in the correct slots
          call get_mul_pmu_ent ( cepoch_expt, ne, 'N' )
*         print *,'DEBUG: get_mul_pmu_ent ',cepoch_expt, ne,
*    .          cwob_apr(1:4), gwob_apr(1:4)
          if( ne.gt.0 ) then
              if( apr_val_mul_pmu(1,1,ne).eq.-9999.d0 ) 
     .                       apr_val_mul_pmu(1,1,ne) = cwob_apr(1)
              gwob_apr(1) = apr_val_mul_pmu(1,1,ne)
              if( apr_val_mul_pmu(1,2,ne).eq.-9999.d0 ) 
     .                       apr_val_mul_pmu(1,2,ne) = cwob_apr(2)
              gwob_apr(2) = apr_val_mul_pmu(1,2,ne)
              if( apr_val_mul_pmu(2,1,ne).eq.-9999.d0 ) 
     .                       apr_val_mul_pmu(2,1,ne) = cwob_apr(3)
              gwob_apr(3) = apr_val_mul_pmu(2,1,ne)
              if( apr_val_mul_pmu(2,2,ne).eq.-9999.d0 ) 
     .                       apr_val_mul_pmu(2,2,ne) = cwob_apr(4)
              gwob_apr(4) = apr_val_mul_pmu(2,2,ne)

              if( apr_val_mul_pmu(1,3,ne).eq.-9999.d0 ) 
     .                       apr_val_mul_pmu(1,3,ne) = cut1_apr(1)
              gut1_apr(1) = apr_val_mul_pmu(1,3,ne)
              if( apr_val_mul_pmu(2,3,ne).eq.-9999.d0 ) 
     .                       apr_val_mul_pmu(2,3,ne) = cut1_apr(2)
              gut1_apr(2) = apr_val_mul_pmu(2,3,ne)

          end if
      end if


****  Thats all
      return
      end 

CTITLE ASSIGN_MUL_EP

      subroutine assign_mul_ep( or, cm, en )

*     Routine to assign the apriori values for a polar motion/UT1 
*     component.  The main task here is to see if any of the apriori
*     in this hfile can be used as apriori for the global estimates.
*                                               
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 

* PASSED VARIABLES
*-----------------
* or -- offset, rate index
* cm -- Type index X, Y, UT1
* en -- Epoch number in the estimated multi-pmu array.

      integer*4 or, cm, en

* LOCAL VARIABLES
*----------------- 
*  i -- Loop counter 
*  in -- Pointer to first epoch for each types of component.
*        (assumes all have the same number of epochs)

      integer*4 i, in

***** Loop over the values read from the hfile and assign the values
      do i = 1, cnum_mul_pmu(or,cm)
*        See if we match epoch
         in = ((or-1) + (cm-1)*2)*cnum_mul_pmu(or,cm) + i
         if( abs(gmul_pmu_ep(en)- cmul_pmu_ep(in)).lt.
     .       tol_mul_pmu) then

*            Epoch match OK, so save values
             apr_val_mul_pmu(or,cm,en) = cmul_pmu_apr(or,cm,i)
*            print *,'DEBUG: assign_mul_ep ',or,cm,en, 
*    .                      apr_val_mul_pmu(or,cm,en) 
         end if
      end do
 
****  Thats all
      return
      end
 
CTITLE READ_MUL_PMU

      subroutine read_mul_pmu

*     This routine reads the values of the apriori's for the multiday
*     polar motion values and saves them the globk_cntl common.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* LOCAL VARIABLES
*----------------
* i,j,k  - Loop counters
* ncod   - Number of entries we need to read from code record
*          (normally 128, until last record)
* na(2,3) - Counter for number of apriori multi-pmu values we have found
* rec_needed -- Record number needed from the apriori records
* curr_apr_rec -- Record number of current record in memory 
* apr_codes(128) -- Code record read from hfile
* apr_vals(64)   -- Apriori values read from hfile
* ierr    -- IOSTAT Error
* type    -- Type of parameter
* indx    -- Sub-type of paramter
* or      -- Offset or rate type (split from indx)
* m2      -- Epoch number split from list
* ent_needed -- Computed value for entry needed in list of values read.

      integer*4 i,j,k, ncod, na(2,3), rec_needed, curr_apr_rec,
     .          apr_codes(128), ierr, len, type, indx, ent_needed,
     .          or, m2

      real*8 apr_vals(64)
*
*     See if this hfile has any multiday values.  If it does not then
*     set the number of values to zero and return

      if( cent_apr_ep.eq.0 ) then
*        No values are available in this hfile.
         do i = 1,3
            do j = 1, 2
               cnum_mul_pmu(j,i) = 0
            end do
         end do
         RETURN
      end if

****  Make sure we have apriori records
      if( cnum_apr_types.eq.0 ) RETURN

*     Initialize the counters
      do i = 1,3
         do j = 1, 2
            na(j,i) = 0
         end do
      end do
      curr_apr_rec = 0

*     Some values are avaiable.  Read the apriori codes from the hfile and
*     save those that are needed

      do i = 1, cnum_apr_types
          call readd(cglb_dcb,ierr,apr_codes,128,len,
     .               crec_apr_types+i-1)

          if( i.eq.cnum_apr_types ) then
              ncod = cnum_apr_codes - (i-1)*128  
          else
              ncod = 128
          end if

          do j = 1, ncod
             call decode_code( apr_codes(j), type, indx)
             if( type.ge.56 .and. type.le.58 ) then

*                OK, found a code for multi-pmu values.  Now get the
*                apriori value from the apriori record.  See if the 
*                entry we need is in memory
                 rec_needed = crec_apr_vals + (i-1)*2 + (j-1)/64
                 ent_needed = mod(j-1,64) + 1
                 if( curr_apr_rec.ne.rec_needed) then
                     call readd(cglb_dcb,ierr,apr_vals,128,len,
     .                          rec_needed )
                     curr_apr_rec = rec_needed
                 end if 

*                We have record, so save value
                 call decode_code(indx, or, m2)
                 na(or,type-55) = na(or,type-55)+1
                 k = na(or,type-55)
                 cmul_pmu_apr(or,type-55,k) = apr_vals(ent_needed)
            end if
         end do
      end do 

****  Thats all
      return
      end


                  
                     


