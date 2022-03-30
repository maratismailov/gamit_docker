CTITLE GLB_O_MINUS_C

      subroutine glb_o_minus_c( sol_obs, looks_ok )

      implicit none  
 
*     This routine will computed o minus c for the parameter estimates
*     read from the global file.  (Note: total values are stored in the
*     global files)
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
*   indx        - Index for the parameter (eg site/source number)
*   i           - Loop counter
*   type        - Parameter type from code.
*   orb_el      - Orbital element number (1-9)  (not used)
*   sv_num      - Number of the satellites      (not used)
*   cnt         - Countering indicating that we have printed message
*                 saying we have converted polar motion
 
      integer*4 indx, i, type, cnt
*     integer*4  orb_el, sv_num

*   Looks_OK    - Indicates prefit residuals are OK

      logical looks_OK
 
*   sol_obs(1)  - Parameter estimates from solution.  These values
*               - will be replaced with O_minus_C.
 
      real*8 sol_obs(cnum_parn)
      
* m1, m2  -- Split of index for multiday PMU
      integer*4 m1, m2, num_md_pmu(3)
      
* true_md_pmu -- Set true if more than 1 multiday PMU entry
      logical true_md_pmu      

***** Loop over the parameter estimates in sol_obs, and based on the
*     parameter code, subtract any apriori values which we have

      looks_ok = .true.
      
*     Check for multiday PMU and see if there is only one.  In which case
*     treat the problem like regular PMU estimate
      true_md_pmu = .false.
      do i = 1, 3
          num_md_pmu(i) = 0
      end do
      do i = 1, cnum_parn

         call decode_code( gpar_codes(i), type, indx )

         if( type.eq.56 .or. type.eq.57 .or. type.eq.58 ) then
             call decode_code( indx, m1, m2)
*            if( this is offset type count up number of times
*            we see parameter            
             if( m1.eq. 1 ) then
                num_md_pmu(type-55) = num_md_pmu(type-55)+1
             end if            
         end if
      end do

      true_md_pmu = .false.
      do i = 1, 3
          if( num_md_pmu(i).gt.1 ) true_md_pmu = .true.
      end do

* MOD TAH 070824: Return single entry multiday PMU for case when mid-point
*     epoch has been used
* MOD TAH 0808225: Changed test to .eq.0 from .ne.0 
      if( .not. true_md_pmu .and. num_md_pmu(1).eq.0 ) then
          write(*,*) 'Retaining single entry multiday PMU'
          true_md_pmu = .true.
      endif

*     Now if not true multiday PMU then convert to regular type
      if( .not. true_md_pmu ) then
          cnt = 0
          do i = 1, cnum_parn
             call decode_code( gpar_codes(i), type, indx )
             if( type.eq.56 .or. type.eq.57 ) then
* MOD TAH 140106: Removed message (not that useful anymore).
!                if( cnt.eq.0 ) 
!    .           write(*,*) 'Converting single multi-day PMU to regular'
                 cnt = cnt + 1
                 call decode_code( indx, m1, m2)
                 gpar_codes(i) = 13 + ((type-55) + (m1-1)*2)*256
 
             else if( type.eq.58 ) then
                 call decode_code( indx, m1, m2)
                 gpar_codes(i) = 14 + m1*256

             end if
          end do
       end if                                           
      
 
      do i = 1, cnum_parn
 
         call decode_code( gpar_codes(i), type, indx )
         call compute_OC( sol_obs(i), type, indx )

      end do
 
****  Thats all
      return
      end

CTITLE GET_TOL_PREF

      subroutine get_tol_pref( type, indx, max_prefit_diff, tol )

      implicit none 

*     Routine to set the tolerance on the prefit residual based on the
*     type of parameter.

* PASSED VARIABLES

*  type  -- Type of parameter (see glb_hdr_def.h)
*  indx  -- Index to type of parameter (glb_hdr_def.h) 
*  max_prefit_diff -- User specficied max differenvce allowed
*  tol   -- Tolerance assigned to parameter

      integer*4 type, indx
      real*4 max_prefit_diff, tol 


****  Set the default to be 100 times specificed by the user (this
*     is so that small site coordinate prefits can be used).
      tol = max_prefit_diff*100
*     Site coordinates 
      if( type.ge. 7 .and.type.le. 9 ) tol = max_prefit_diff 
*     EOP parameters (10 times position uncertainity)
      if( type.ge.13 .and.type.le.14 ) tol = max_prefit_diff*320
      if( type.ge.56 .and.type.le.58 ) tol = max_prefit_diff*320
*     Satellite orbits
      if( type.eq.51 ) tol = max_prefit_diff*1000

****  That all for the moment
      return
      end

