CTITLE GLB_INC_PREFIT
 
      subroutine glb_inc_prefit( cov_obs, sol_obs , dchi, save_chi )

      implicit none  
 
*     Routine to increment the prefit residual Chi**2 for the
*     global solution.
*
*     The chi**2 is computed rigously using:
*        2          T
*     Chi  = sol_obs  cov_obs sol_obs
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/globk_markov.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
*   i,j     - Loop counters
 
 
      integer*4 i
 
*   cov_obs(cnum_parn,cnum_parn)    - Covariance matrix
*                   - of data
*   sol_obs(cnum_parn)              - Prefit residuals after
*                   - update for current estimate of parameters
*   loc_sum         - Sum of prefits for this experiment
*   temp_sum        - Temporary summation variable
*   dchi            - Change in chi**2
 
      real*8 cov_obs(cnum_parn,cnum_parn), sol_obs(cnum_parn), loc_sum,
     .    temp_sum, dchi

*   outline   - Line to be output by report_stat
*   progname  - Program name

      character*256 outline, progname
      
*   save_chi  - Set to Y if chi**2 increments to be saved (not saved in
*     back solution)

       character*(*) save_chi

*   lenprog  - Length of program name
*   rcpar    - Returns runstring

      integer*4 lenprog, rcpar    
 
      logical first_call, kbit 

* Local variables for checking state of system
      integer*4 ir_minv, ir_maxv, ir_maxs,
     .          ir_max_nosym, ic_max_nosym,
     .          ir_max_rho, ic_max_rho, num_nopos, j

      real*8 minv, maxv, max_rho, rho, maxs, max_nosym, avs,
     .       nosym, loc_chi, big_chi
 
      save first_call, progname
 
 
      data  first_call  / .true. /, big_chi / 99999.9d0 /
 
***** If this is first call clear the summation variables
 
      if( first_call ) then
          sum_chi     = 0.d0
          sum_chi_num = 0.d0
          first_call = .false.
          lenprog = rcpar(0, progname )
      end if
 
****  Now compute increment to chi**2
 
      loc_sum = 0.0d0
      do i = 1, cnum_used
          call DWDOT( temp_sum, cov_obs(1,i),1, sol_obs,1,
     .                cnum_used )
          loc_sum = loc_sum + sol_obs(i)*temp_sum
      end do

* MOD TAH 000812: Check to make sure we do not get a NaN.
      if( cnum_used.gt.0 ) then
          dchi = loc_sum/cnum_used 
      else
          dchi = big_chi
      end if

* MOD TAH 980519: Save the chi**2 increment so that it can be
*     written to srt_file.
      dchi_save = dchi
      if(  dchi.ge.0 .and. dchi.lt.max_chi_inc .and.
     .     dchi.ne. big_chi ) then
          if( save_chi(1:1).eq.'Y' ) then
              sum_chi     = sum_chi     + loc_sum
              sum_chi_num = sum_chi_num + cnum_used
          end if
      else
          write(log_unit,90) glb_inp_file
  90      format('** ',a40,' will not be used')         
          call report_stat('warning',progname,'glfor',glb_inp_file,
     .                 'Not used chi**2 increment too large',0)
      end if

      if( cnum_used.gt.0 ) then
          loc_chi = loc_sum/cnum_used
      else
          loc_chi = big_chi
      end if

      write( log_unit,100 ) sdata_base, glb_inp_file, cnum_used, 
     .                      loc_chi           
 100  format(' For ',a10,' Glbf ',a40,' Chi**2 NP ',i5,
     .       ' is ',f10.3)
      if( log_unit.ne.6 ) then 
           write( *,100 ) sdata_base, glb_inp_file, cnum_used, 
     .                    loc_sum/cnum_used
      end if

      if( abs(loc_sum/cnum_used).lt.1000.d0 ) then
          write(outline,110) cnum_used, loc_sum/cnum_used
 110      format('Chi**2 Increment for ',i5,' dof ',F7.3)
      else
          write(outline,120) cnum_used, loc_sum/cnum_used
 120      format('Chi**2 Increment for ',i5,' dof ',d11.3)
          call report_stat('warning',progname,'glfor',glb_inp_file,
     .                 outline,0)
      end if
      call report_stat('status',progname,'glfor',glb_inp_file,
     .                 outline,0)

****  Now check the matrix more chi**2 increment os negative
      if( loc_sum.lt.0 ) then

*         Check if cov_obs has negative diagonal.
          minv = 0.d0
          maxv = 0.d0
          max_rho = 0.d0
          max_nosym = 0.d0
          num_nopos = 0
          do i = 1, cnum_used
             if( cov_obs(i,i).lt. 0.d0 ) then
                 write(log_unit,210) i, cov_obs(i,i)
 210             format('** Negative diagonal in cov_obs for row ',i5,
     .                  ' value ',d16.8)
             end if
*            Get min and max variance
             if( cov_obs(i,i).lt. minv ) then
                 minv = cov_obs(i,i)
                 ir_minv = i
             end if
             if( cov_obs(i,i).gt. maxv ) then
                 maxv = cov_obs(i,i)
                 ir_maxv = i
             end if

*            Check the off-diagonal terms
             do j = i+1, cnum_used
                nosym = (cov_obs(j,i) - cov_obs(i,j))/cov_obs(i,j)
                if ( abs(nosym).gt.abs(max_nosym) ) then
                     ir_max_nosym = i
                     ic_max_nosym = j
                     max_nosym = nosym
                end if

                rho = cov_obs(i,j)/sqrt(cov_obs(i,i)*cov_obs(j,j))
                if( abs(rho).gt.abs(max_rho) ) then
                    ir_max_rho = i
                    ic_max_rho = j
                    max_rho = rho
                end if
                if( abs(rho).gt.1.d0 ) num_nopos = num_nopos + 1
             end do
                
          end do
          if( kbit(crt_opts,16) ) then 
             write(log_unit,220)  ir_minv, minv, ir_maxv, maxv
 220         format('**Minimum Cov NEQ for row ',i5,' Value ',d16.8,/,
     .              '**Maximum Cov NEQ for row ',i5,' Value ',d16.8)
             write(log_unit,240) ir_max_nosym, ic_max_nosym, max_nosym
 240         format('**Maximum assymetry   for row ',i4,' col ',i4,
     .              ' ratio ',d16.8)
             write(log_unit,250) ir_max_rho, ic_max_rho, max_rho
 250         format('**Maximum correlation for row ',i4,' col ',i4,
     .              ' Rho   ',d16.8)
             write(log_unit,260) num_nopos
 260         format('**There are ',i4,' Non-positive definite elements')

*            Now check the maximum value of the contrubtion to chi**2
             maxs = 0.d0
             avs  = 0.d0
             do i = 1, cnum_used
                 call DWDOT( temp_sum, cov_obs(1,i),1, sol_obs,1,
     .                       cnum_used )
                 avs = avs + abs(temp_sum)
                 if( abs(temp_sum).gt.abs(maxs) ) then
                     ir_maxs = i
                     maxs = temp_sum
                 end if
             end do
             write(log_unit,340) ir_maxs, maxs, avs/cnum_used
 340         format('**Maximum Sum contribution for row ',i4,
     .              ' Value ',d16.8,' Average ',d16.8)
         end if 
      end if
   
****  That's all
      return
      end
 
