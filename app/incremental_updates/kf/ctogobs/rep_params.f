CTITLE rep_run_params
 
      subroutine rep_run_params(un)

      implicit none
 
*     This routine will report the run parameters used during a
*     ctogobs run.  It also provides a convenient means of dumping
*     the current defaults.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
c     include '../includes/gobs_header.h'
      include 'ctogobs_com.h'
      include '../includes/cfile_def.h'
 
 
* PASSED VARIABLES 

* un            - Unit number for writting parameters to

      integer*4 un
 
* LOCAL VARIABLES
 
*   trimlen     - Returns length of non-blank portion of string
*   i,j         - Loop counter
*   il          - length of string (used when length may be empty)

      integer*4 trimlen, i, j,il

*   kbit        - Checks status of bit

      logical kbit

      character*12 full_ver, hsver

***** See if we need to write this out
      if( .not.kbit(status_rep,9) )  RETURN

****  Staring writing out the current parameters
      full_ver = hsver(ctogobs_version)

      write(un,100) caprog_name, full_ver(1:trimlen(full_ver))
 100  format(50('+'),/,1x,a8,'Version ',a,' Run parameter status')
      il = max(1,trimlen(ctogobs_com_file))
      write(un,120) 'COMMAND FILE:', ctogobs_com_file(1:il)
 120  format(1x,a,t30,a)
      il = max(1,trimlen(summary_file))
      write(un,120) 'SUMMARY FILE:', summary_file(1:il)
      il = max(1,trimlen(Gobs_outfile))
      if( caprog_name(1:4).eq.'CTOG' )
     .write(un,120) 'GOBS FILE:',Gobs_outfile(1:il)
      il = max(1,trimlen(mf_name))
      write(un,120) 'POSTFIT M-FILE:',mf_name(1:il)      
      il = max(1,trimlen(rng_clk_root))
      write(un,120) 'RANGE CLOCK ROOT:',rng_clk_root(1:il)
      il = max(1,trimlen(phs_clk_root))
      write(un,120) 'PHASE CLOCK ROOT:',phs_clk_root(1:il)
      il = max(1,trimlen(phs_res_root))
      write(un,120) 'PHASE RESIDUAL ROOT:',phs_res_root(1:il)
      il = max(1,trimlen(sng_diff_root))
      write(un,120) 'SINGLE DIFFERENCE ROOT:',sng_diff_root(1:il)
      il = max(1,trimlen(dd_outfile))
      write(un,120) 'DOUBLE DIFFERENCE REPORT:',dd_outfile(1:il)
      write(un,140) cf_iy, cf_im, cf_id, cf_gnss
 140  format(/' INPUT CFILES: Date: ',i4,1x,i2,1x,i2,1x,' SYS ',a1)
      do i = 1, num_cfiles
          write(un,150) i,cf_codes(i), rcv_types(i),
     .                  cfiles(i)(1:trimlen(cfiles(i)))
 150      format(I3,'. ',a4,2x,a4,2x,a)
      end do
      write(un,160) new_cfile_series
 160  format(' OUTPUT CFILE SERIES ',a)

****  Now do the jump/editing type parameters
      write(un,200)
 200  format(/,' EDITING AND JUMP PARAMETERS')
      write(un,220) rng_jump_tol, rng_jump_min/(1.d-6*fClk)
 220  format(' Range-clock tolerances ',t30,F7.2,' n-allan sd',1x,
     .       F7.2,' usec')
      write(un,240) reset_tol/(1.d-6*fClk)
 240  format(' Millisec clock reset',t30,F7.2,' usec')
      write(un,260) rng_res_tol, rng_res_min*vel_light/fClk,
     .              rng_res_max*vel_light/fClk
 260  format(' Bad-rainge tolerances ',t30,f7.2,' n-range sd',1x,
     .       F7.2,' min m ',F7.2,' max m')
      write(un,270) rel_clk_wght
 270  format(' Relative clock weight',t30,f7.2)
 
      if( max_scan_edit.gt.0 ) then
          write(un,273) max_scan_edit
 273      format(' Site/SVS with more than ',i4,' DDscan bias',
     .           ' flags edited')
      end if
      write(un, 274) min_ow_data
 274  format(' Minimum oneway data per satellite needed ',i4,' epochs')
              
      if( num_pre_edit.gt.0 ) then
          write(un,275) num_pre_edit
 275      format(' Data pre-editing for ',i3,' cases',/,
     .           '    # SITE   PRN  Start  Stop Epoch')
          do i = 1, num_pre_edit
              if( pre_edit(1,i).gt.0 ) then
                  write(un,280) i, cf_codes(pre_edit(1,i)),
     .                            (pre_edit(j,i),j=2,4)
 280              format(1x,i3,1x,A4,2x,I3,2x,I5,1x,I5)
              else
                  write(un,280) i, ' ALL',
     .                            (pre_edit(j,i),j=2,4)
              end if
          end do
      end if

***** Report the azimuth mask used
      call report_azmask(un)

***** Report phase clock estimation and postfit editing
      if( apply_phs_clk ) then
          write(un,300)
 300      format(/,' PHASE CLOCK AND POSTFIT EDITING')
          write(un,310) pc_max_iter, pc_non_int, pc_convd, 
     .                  pc_over_shoot,
     .                  pc_full_anal, pc_max_corr
 310      format(' Phase clocks estimated: Max iterations ',i3,/,
     .        25x,'Non-Integer offsets start iteration ',i3,/,
     .        25x,'Converged when RMS change ratio < ',f5.2,/,
     .        25x,'Mean offset over-shoot ',f5.2,/,
     .        25x,'Full Analysis for ',i3,' Iterations with max. ',
     .            'correlation of ',F5.2)

          if( fdd_L2_fact.eq.0 ) then
              write(un,315)
 315          format('Phase clocks computed with L1 only data')
          end if

          if( nol1only ) write(un,'(a)') 'Default: L1-only data deleted'
          if( .not.nol1only ) write(un,'(a)') 
     .                            'L1only: All L2 data deleted'
          if( prefit_clk ) write(un,'(a)') 
     .                            'Prefit clocks used (PREFIT_CLK)'

          if( edit_postfit ) then
              write(un,320) pc_start_edit, pf_nsig, pf_maxres, 
     .                pf_max_rms
 320          format(' Postfit residuals edited: Start Iteration ',i3,
     .               ' Data > ',f5.2,' Station RMS flagged,',/,
     .               25x,'Data < ',
     .               f5.2,' (cyc) restored. SITE/SV RMS > ',f5.2,
     .               ' (cyc) removed')
          end if
      end if

****  If gaps were flagged give the size of the gaps allowed
      if( gaps_flagged ) then
          write(un, 360)
 360      format(' Minimum epochs for gaps ',/,(10(' SITE Sze ')))
          write(un,365) (cf_codes(i), gap_size(i),i=1,num_cfiles)
 365      format(100(10(1x,A4,1x,i3,1x,:),:,/))
      end if
     

****  Bias fixing paramters
      write(un,600)
 600  format(/,' BIAS DETECTION AND FIXING PARAMETERS')
      write(un,620) (phs_fit_tol(i),i=1,4)
 620  format(' One-way mean and max differences for range clock',
     .       2F9.2,', for phase clock ',2f9.2)
      write(un,640) dchi2_ratio, dchi2_min_val, dchi2_max_sep,
     .              dchi2_gap_fact
 640  format(' Double differnce: Ratio ',F7.2,' Min Chi**2 ',f7.2,
     .       ' Max separation ',F7.1,' secs, Gap factor ',F4.2)
      write(un,650) max_wl_ret, max_dd_ret, max_lg_use, tol_one_way_fix,
     .       trim_seg
 650  format(' Max. data returns: WL ',i4,' DD LC ',I4,' DD LG ',i4,
     .       '. Max gap for one-way fix ',I4,' secs, Trim Seg ',L1)
      write(un,660) dd_wl_tol, dd_lc_tol
 660  format(' Double difference slip detection:',/,
     .       ' WL: Scale ',F7.2,' Min (cyc) ',F7.2,' Max (cyc)',F7.2,/,
     .       ' LC: Scale ',F7.2,' Min (cyc) ',F7.2,' Max (cyc)',F7.2)
      write(un,680) min_dtl_bias, min_good_bias, min_dtr_end,
     .             min_good_end
 680  format(' One-way triming: Min. time between bias flags ',
     .       F7.1,' secs, min number epochs ',I4,/,
     .       '                  Min. fraction to last bias   ',
     .       F7.3,'     , min number epochs ',I4)
      write(un,690) 
 690  format(' Ionospheric jump detector (first and values)',
     .       ' which differ',/,
     .       ' Site  MaxGap (sec) dIon scale MIN dIon (cyc)',
     .       ' MAX dIon (cyc)')
      write(un,695) cf_codes(1), dt_ion_tol, ion_rw_tol(1,1),
     .             ion_rw_tol(2,1), ion_rw_tol(3,1)
 695  format(1x,a4,2x,f10.2,3x,f10.3,3x,f10.3,3x,f10.3)
      do i = 2, num_cfiles
          if( ion_rw_tol(1,i).ne.ion_rw_tol(1,i-1) .or.
     .        ion_rw_tol(2,i).ne.ion_rw_tol(2,i-1)	) then
              write(un,695) cf_codes(i), dt_ion_tol, 
     .                     ion_rw_tol(1,i), ion_rw_tol(2,i),
     .                     ion_rw_tol(3,1)
          end if
      end do

      if( pf_remove_bf ) write(un,820) 'Postfit Remove Bias flag',' '
      if( .not.pf_remove_bf ) write(un,820) 'Postfit Remove Bias flag',
     .                                 ' not '

      if( resolve_wl ) then
         write(un,720) min_wl_tol, dchi_wl_tol, msig_wl_tol, mdev_wl_tol
 720     format(' LC_AUTCLN Options: Min Overlap ',i4,
     .          ' epochs, Dchi Min ',F6.1,' Min Sig ',F6.2,
     .          ' Max Deviation ',F6.2,' cyc')
      endif 

****  Miscellaneous
      write(un,800)
 800  format(/' MISCELLANEOUS PARAMETERS')

      if( use_gamit_elc ) write(un,820) 'GAMIT Elev. cutoff',' '
 820  format(1x,a,a,'used')
      if( .not.use_gamit_elc ) write(un,820) 'GAMIT Elev. cutoff',
     .                                       ' not '

      if( use_cview_edit ) write(un,820) 'CVIEW Editing',' '
      if( .not.use_cview_edit ) write(un,820) 'CVIEW Editing',' not '

      if( usr_ignore_gaps ) write(un,840) 'DATA GAPS',' '
 840  format(1x,a,1x,a,'ignored')
      if( .not.usr_ignore_gaps ) write(un,840) 'DATA GAPS',' not '

      if( .not.do_one_bg  ) write(un,850) 'One bias or gap',' not '
      if( do_one_bg  ) write(un,850) 'One bias or gap',' '
 850  format(1x,a,1x,a,'allowed in single differences')

C     if( .not.use_gamit_elc ) write(un,860) caprog_name, 
C    .                                       min_ctog_elev*180/pi,
C    .                                       min_out_elev*180/pi
 860  format(' Data in ',a8,' Cleaned to ',f5.2,' deg, saved to ',
     .       f5.2,' deg')
      if( use_gamit_elc ) write(un,860) caprog_name,
     .                                  min_cfile_elev*180/pi,
     .                                  min_cfile_elev*180/pi

      if( apply_phs_clk ) write(un,870)
 870  format(' Residuals updated with phase clock estimates')
      if( app_ion ) write(un,'(a)') 
     .  'ION Delay applied OMC values in c-files'

      if( remap_glonass .and. cf_gnss.eq.'R')
     .   write(un,'(" GLONASS ambiguities scaled by frequency")')
      if( .not.remap_glonass .and. cf_gnss.eq.'R')
     .   write(un,'(" GLONASS ambiguities kept as integer")')



***** Thats all for the moment
      write(un,900)
 900  format(1x)
      return 
      END
 
CTITLE REPORT_AZMASK

      subroutine report_azmask(un)

      implicit none

*     Routine to report to the summary files the azimuth mask
*     used at a station.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES:
* un  -- Unit number

      integer*4 un

* LOCAL VARIABLES
* max_out   -- Maxium number of fields in az_mask at any station
* i,j       -- Loop counters

      integer*4 max_out, i,j


****  Loop over the stations to find out maxium width of field
*     needed (and if anything needs writing).
      max_out = 0
      do i = 1, num_cfiles
         if( num_azmask(i).gt.max_out ) max_out = num_azmask(i)
      end do

*     See if we have anyhing
      if( max_out.gt.0 ) then
          write(un,120) (i, i=1, max_out)
 120      format(/,' AZIMUTH MASK TABLE BY SET SITES (degrees)',/,
     .             ' SITE ',400(1x,'AZ',i2.2) )
          write(un,140) (i, i=1, max_out-1)
 140      format(8x,400(1x,'EL',i2.2))

*         Now loop over sites and write the values
          do i = 1, num_cfiles
             if( num_azmask(i).gt.0 ) then
                write(un,220) cf_codes(i), 
     .               (nint(az_mask(1,j,i)*180./pi),j=1,num_azmask(i))
 220            format(1x,a4,1x,400(1x,I4))
                write(un,240)  
     .               (nint(az_mask(2,j,i)*180./pi),j=1,num_azmask(i)-1)
 240            format(' Elev ',2x,400(1x,I4))
             end if
          end do
      end if

****  Thats all
      return
      end

