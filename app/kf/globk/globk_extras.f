CTITLE GET_MUL_PMU_ENT

      subroutine get_mul_pmu_ent ( epoch, mul_pmu_ent, option )

      implicit none 

*     Routine to match an epoch 'epoch' with the available list of
*     epochs.  If there is not match mul_pmu_ent is returned as -1
*
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
* PASSED VARIABLES
*  epoch  -- Epoch of data pmu estimate
*  mul_pmu_ent -- Entry in curr_mul_pmu_ep corresponding to this epoch
*            (with in the time tolerence set by user)
*  option -- Y if station and satellite list is to be saved.

      integer*4 mul_pmu_ent
      real*8  epoch
      character*(*) option  

* LOCAL VARIABLES
* i,j  -- Loop counter
* ns   -- Global site or svs number

      integer*4 i,j, ns

***** Loop over list
      mul_pmu_ent = -1

      do i = 1, num_mul_pmu
         if( abs(epoch-gmul_pmu_ep(i)).lt. tol_mul_pmu ) then 
            mul_pmu_ent = i

****        Now count the sites that are contributing
            if( option.eq.'Y' ) then
               do j = 1, cnum_sites
                  ns = ltog_sites(j)
                  if( ns.gt.0 ) then
                     call sbit(mul_pmu_site(1,i),ns,1)
                     if( apr_neu(1,1,ns)+apr_neu(2,1,ns).lt.0.05 .and.
     .                   apr_neu(1,1,ns)+apr_neu(2,1,ns).gt.0.0 ) then
                         call sbit(mul_pmu_fidu(1,i),ns,1)
                     end if
                  end if
               end do
               do j = 1, cnum_svs
                  ns = ltog_svs(j)
                  if( ns.gt.0 ) then
                      call sbit(mul_pmu_svs(i),ns,1)
                  endif
               end do
            end if
         end if
      end do 

****  Thats all
      return
      end

CTITLE MAP_WOB_INDX

      subroutine map_wob_indx( type, sng1, mul1, mul2 )

      implicit none 

*     Routine to map the two types of polar motion/UT1 indices
*     back and forth.  The correspondence is:
*     SINGLE  MULTIPLE
*          1    1,1     xpole
*          2    1,2     ypole
*          3    2,1     xpole rate
*          4    2,2     ypole rate 
*     type sets the direction i.e, 'STOM' Mapps single to multiple

* PASSED VARIABLES
* type  -- Type of mapping STOM or MTOS
* sng1  -- Single index
* mul1, mul2 -- Pair of multiple indices

      integer*4 sng1, mul1, mul2
      character*4 type

****  Treat the two cases
      if( type.eq.'STOM' ) then
          if( sng1.le.2 ) then
              mul1 = 1
              mul2 = sng1
          else if ( sng1.le.4 ) then
              mul1 = 2
              mul2 = sng1 - 2
          else
              write(*,150) sng1, type
 150          format('** WARNING ** Invalid entry ',i4,' Passed to ',
     .               'MAP_WOB_INDX with type ',a)
          end if
      else if( type.eq.'MTOS' ) then 
          sng1 = (mul1-1)*2 + mul2
      else
          write(*,170) type
 170      format('** WARNING ** Invalid type ',a,
     .           ' passed to MAP_WOB_INDX') 
      endif

****  Thats all 
      return
      end

CTITLE SETUP_MUL_PMU 
               
      subroutine setup_mul_pmu 

      implicit none 

*     Routine to set up the multi-day pmu epochs and initialize the
*     apriori values so that new ones will be read from the binary
*     h-file if needed.

 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'


* i,j,k  -- Loop counters
* gpsweek -- GPS week number
* mblk    -- Multiday PMU estimate block number

      integer*4 i,j,k, gpsweek, mblk

      logical kbit

* gps_start_jd -- MJD of start of GPS week 
* start_jd     -- JD of start of current week
* cepoch_mod   -- Current experiment epoch matched to pmu epochs

      real*8 gps_start_jd, start_jd, cepoch_mod 
      real*8 gepoch_mid   ! Mid_point of expreriment.  Use to get the
                          ! the starting time of the mul_pmu parameters
      data gps_start_jd / 2444243.5d0 /

***** Also initialize the standard PMU values
      if( num_mul_pmu.eq.0 ) then

*         if( abs(cepoch_expt-gmul_pmu_ep(1)).gt.tol_mul_pmu ) then
* MOD TAH 131231: Changed test to use exper versus global epoch when
*         multi-day PMU is not used.
          if( abs(cepoch_expt-gepoch_expt).gt.tol_mul_pmu ) then
             do i = 1,4
                gwob_apr(i) = -9999.d0
                gut1_apr(i) = -9999.d0
              end do
          else 
             gmul_pmu_ep(1) = cepoch_expt
          end if
          RETURN
      end if

****  See what week we are in 
      call sbit(mul_pmu_opt,32,0)
      call sbit(mul_pmu_opt,31,0)


      gpsweek = int((cepoch_expt-gps_start_jd-1)/7)
* MOD TAH 030220: Use the mid-point to get the starting week
      gepoch_mid = (gepoch_start+gepoch_end)/2
      gpsweek = int((gepoch_mid-gps_start_jd-1)/7)
*     See if
      if( kbit(mul_pmu_opt,2) ) then    ! Actual epoch given 
          start_jd = start_mul_pmu 
      else
          start_jd = gpsweek*7.d0 + gps_start_jd +1 + start_mul_pmu
*         Save the start epoch as the full value and use this value
*         for later test
          start_mul_pmu = start_jd
          call sbit(mul_pmu_opt,2,1)
      end if

*     We do two check here to see if need to re-initialize the
*     EOP values.  One for mul_pmu estimation and the other for
*     single day estimation.  We also need to check to see if
*     this experiment needs a new block of multiday PMU tables
*     be established.

      mblk =  nint((cepoch_expt-start_mul_pmu)/spacing_mul_pmu)
      cepoch_mod = start_mul_pmu + mblk*spacing_mul_pmu

*     First see if the current experiment seems to fall with in
*     bounds.  There will be problems with this code for input
*     multiday PMU estimates which do not overlap.
      if( cepoch_expt.ge. gmul_pmu_ep(1)-tol_mul_pmu .and. 
     .    cepoch_expt.le. gmul_pmu_ep(num_mul_pmu)+tol_mul_pmu ) RETURN
*     if( start_jd .eq. gmul_pmu_ep(1) .and. num_mul_pmu.gt.0 ) RETURN
*     if( abs(cepoch_expt-gmul_pmu_ep(1)).lt.tol_mul_pmu .and. 
*    .        num_mul_pmu.eq.0 ) RETURN 

*     Get the block number for the estimates.  This is based on the
*     cepoch_expt.  (Again may be problem for non-alligned estimates)
      mblk = int((cepoch_mod-start_mul_pmu)/
     .           (num_mul_pmu*spacing_mul_pmu))
      if( start_mul_pmu-cepoch_expt.gt.tol_mul_pmu ) then
          mblk = mblk - 1
      end if

* MOD TAH 020322: Keep the same block unless the user has indicated 
*     that they should be propagated
      if( .not.kbit(mul_pmu_opt,6) ) then
          mblk = 0
      endif
      
      write(*,800) start_mul_pmu,  start_jd, gepoch_mid, 
     .             gepoch_start, gepoch_end, mblk  
 800  format('MUL_PMU: Start ',2F12.3,' GEPOCH Mid Start End ',3F12.3,
     .       ' MBLK ',i4)

      do i = 1, num_mul_pmu
         gmul_pmu_ep(i) = start_jd + 
     .                    (mblk*num_mul_pmu+i-1)*spacing_mul_pmu 
         do j = 1, 3
            do k = 1, 2
               apr_val_mul_pmu(k,j,i) = -9999.d0
            end do
         end do

****     Clear out the counters used to save the number of sites
*        and satellites
         do j = 1, max_glb_site_wrds
            mul_pmu_site(j,i) = 0
            mul_pmu_fidu(j,i) = 0
         end do
         mul_pmu_svs(i) = 0

      end do 

*     Indicate the mul_pmu epochs have been updated and that the process
*     noise needs to propagated (for the first experiment, bit 31 will
*     turned off in glfor).
      call sbit(mul_pmu_opt,32,1)
      call sbit(mul_pmu_opt,31,1)

***** Thats all
      return
      end

 
CTITLE PMU_MIDPT

      subroutine pmu_midpt

      implicit none 

*     Routine to convert single day polar motion/UT1 into single 
*     entry multiday polar motion/UT1 when the output epoch of the
*     soltion does not match the PMU epoch (happens when -M is set
*     in glsave or midp option is selected in org or prt_opts.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* LOCAL VARIABLES
      integer*4 i, j
      logical kbit

     
***** Polar motion/UT1 values
* MOD TAH 070823: If we are outputing at mid-point convert the
*     dailty polar motion/UT1 to multiday with the correct 
*     epoch (this needs to be the last value)
      if( num_mul_pmu.eq. 0 .and.  
     .    abs(gepoch_out-gepoch_expt).gt.1.d-3 ) then
          do i = 1,2
             parn_mul_pmu(1,i,1) = parn_wob(i)
             apr_val_mul_pmu(1,i,1) = apr_val_wob(i) + gwob_apr(i)

             parn_mul_pmu(2,i,1) = parn_wob(i+2)
             apr_val_mul_pmu(2,i,1) = apr_val_wob(i+2) + gwob_apr(i+2)
          end do
          parn_mul_pmu(1,3,1) = parn_ut1(1)
          apr_val_mul_pmu(1,3,1) = apr_val_ut1(1) + gut1_apr(1)
          parn_mul_pmu(2,3,1) = parn_ut1(2)
          apr_val_mul_pmu(2,3,1) = apr_val_ut1(2) + gut1_apr(2)
*         Set the number
          num_mul_pmu = 1
          gmul_pmu_ep(1) = gepoch_expt  ! Must be referenced to original time
          do i = 1,4
             parn_wob(i) = 0  ! Stops regular output
          end do 
          do i = 1,2
             parn_ut1(i) = 0
          end do
          do j = 1,gnum_sites
             if( kbit(guse_site,j) ) call sbit(mul_pmu_site(1,1),j,1)
          enddo
          do j = 1, gnum_svs
             call sbit(mul_pmu_svs(1),j,1)
          enddo

      end if

****  Thats all
      return
      end


  
     
