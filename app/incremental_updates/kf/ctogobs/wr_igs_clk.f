CTITLE WR_IGS_CLK

      subroutine wr_igs_clk( ins, params_cse, par_flag_cse, svcL1_ce, 
     .                       data_flag_cse, ctol_cse, cref_upd, type)

      implicit none

*     Routine to write the clock values out in IGS RINEX clock
*     format

* INCLUDES

      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include '../includes/mfile_def.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES

*   par_flag_cse(num_param,num_ep)         - parameter estimate flag
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag saying what 
*           data is available
*   ctol_cse(num_chan, num_cfiles, num_ep) -- Conversion from channel to PRN
*   ins  -- Summary file unit number

      integer*4 par_flag_cse(num_param,num_ep), 
     .          data_flag_cse(num_chan, num_cfiles, num_ep),
     .          ctol_cse(num_chan, num_cfiles, num_ep),  ins

*   params_cse(num_param,num_ep)    - Clock parameter estimates (L1 cycles)
*   svcL1_ce(num_sat,num_ep) - Part of satellite clock removed by model
*                     (L1 cycles)


      real*8 params_cse(num_param,num_ep)
      real*4 svcL1_ce(num_sat,num_ep)

* cref_upd -- Set true is reference clock list is updated.
      logical cref_upd

* type  -- P or R for range or phase clocks

      character*1 type  

* LOCAL VARIABLES
* ierr  -- IOSTAT error 
* i     -- Loop counter
* ep    -- Epoch number being output
* date(5) -- Generic array for date information
* ref_clk_num(max_cfiles) -- Numbers for reference clocks
* dday  -- NUmber of days since GPSW zero
* gpsw, gpsd  -- GPS Week number and day of week number
* num_out_clk -- Number of sites to be output
* out_clk_num(max_cfiles) -- List of sites to output
* mfs         -- Mfile site number (used to get postion adj)
* jsite_xyz(3) -- Site position in mm


      integer*4 ierr, i, j, k, ep, date(5), num_var,
     .          ref_clk_num(max_cfiles), ns, jndx, jel,
     .          trimlen, dday, gpsw, gpsd, num_out_clk,
     .          out_clk_num(max_cfiles), mfs, ls, nc, ch, ltoc

      logical data_OK  ! Function to see if data good

      real*8 jsite_xyz(3)

* mjd_start -- MJD of first cfile epoch
* mjd       -- MJD of current epoch
* dt        -- Time difference in seconds from start
* clk_out   -- Output value of clock (seconds) 
* sectag    -- Seconds tag
* mean_clk  -- Mean reference clock at each epoch
* mean_site(max_cfiles) -- Mean site clocks
* site_llr(3) -- Site position Lat, long, rad (km)
* Sum_dmean -- Sum for Dclk estimates (allows slope removal)
* SaveRMSdm(max_cfiles) -- Saved values of RMS about slope (needed
*     for sites that are not output to te clock file).
  
      real*8 mjd_start, mjd, dt, clk_out, sectag, sum_var,
     .       mean_clk, mean_site(max_cfiles), sum_mean, rms_clk,
     .       site_llr(3), dclk, sum_dmean, dmean, rms_dmean,
     .       saveRMSdm(max_cfiles), clk_sig

      logical assign_ref, kbit

      real*8 fin_stat(3,max_cfiles+max_gprn)
      integer*4 kp

* Fin_ref_sites -- List of final reference sites
* min_ref_rms   -- Minimum RMS at reference sites
* nrfn -- Number of OK reference sites
* ref_fact -- Divider on ref_rms between Phase and range.

      integer*4 fin_ref_sites(max_cfiles), nrfn
      real*8 min_ref_rms, ref_fact

* prog, runby -- Character strings with program and runby
      character*20 prog, runby

* posl(3) -- Position written as floating point so that can be
*     translated to interger

      character*12 posl(3)

* ot -- Output status (R - Reference, O output)
      character*2 ot

* line  -- Line for outputing prn list
      character*80 line

* Phase center model used
      character*40 phs_cen_mod

****  New variables when RMS fits are not good.  In this case
*     best 5 sites are choosen
      real*8 best_rms(5)
      integer*4 best_cfs(5)

****  Get the date and time of the first measurement of
*     the day
      cref_upd = .false.
      date(1) = cf_iy
      date(2) = cf_im
      date(3) = cf_id
      date(4) = cf_ihr
      date(5) = cf_min
      call ymdhms_to_mjd(date, cf_sec, mjd_start)

*     Get number of days since GPS Week 0, day 0
      dday = mjd_start - 44244.0d0
      gpsw = dday/7
      gpsd = (dday-gpsw*7) 

***   Open the output file
      if( trimlen(igs_clk_file).eq.3 ) then
*         Generate name with GPS week and day 
          write(igs_clk_file(4:),110) gpsw, gpsd
 110      format(i4.4,i1,'.clk')
      end if
      
      ns = trimlen(igs_clk_file)
      if( type.eq.'R'  ) then
          if( igs_clk_file(ns:ns).ne.'R' ) igs_clk_file(ns+1:) = 'R'
      else
          if( igs_clk_file(ns:ns).eq.'R' ) igs_clk_file(ns:ns) = ' '
      end if

      open(200,file=igs_clk_file, status='unknown',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',igs_clk_file,0,
     .                  'WR_IGS_CLK')

****  Write out the header information
* MOD Use 3.00 version headers
      write(200,120) 3.00,'CLOCK DATA'
 120  format(F9.2,11x,a10,30x,'RINEX VERSION / TYPE')
      call systime( date, sectag)
      call getenv('INSTITUTE',runby)
      write(prog,130) 'autcln',ctogobs_version
 130  format(a,1x,a)
      write(200,135) prog, runby, date
 135  format(a20,a20,i4,'-',i2.2,'-',i2.2,1x,i2.2,':',
     .       i2.2,4x,'PGM / RUN BY / DATE')

c      write(200,140)
c 140  format('CLK ANT Z-OFFSET(M): II/IIA 1.0230; IIR 0.0000',
c     .       '              COMMENT')
      phs_cen_mod = cf_antmod_snx(1:10)
* MOD TAH 200723: Write out correct GNSS system
      write(200,140) sv_gnss, cf_nversn/100., phs_cen_mod
  140 format(a1,' GAMIT ',F5.2,'       ',a40,'SYS / PCVS APPLIED')
c     format('G Program_Name_____ IGS05_1400                              SYS / PCVS APPLIED
      write(200,150)                                       
 150  format('     2    AS    AR                           ',
     .      '               # / TYPES OF DATA')
      write(200,160) runby
 160  format(a20,40x,'ANALYSIS CENTER')
 

****  Now see how many reference clocks we are using
      ref_fact = 1.d0 
      if ( type.eq.'P' ) ref_fact = 5.d0 
      write(ins,165) num_ref_clk, rms_ref_clk/ref_fact/fClk*1.d9,
     .               (ref_clk_code(i),i=1,num_ref_clk)
 165  format('REFERENCE CLOCKS: Num ',i3,' RMS (ns) ',
     .        F10.2,' LIST ',30(a4,1x))


      assign_ref = .false.
* MOD TAH030130: Make sure the reference sites look OK
      min_ref_rms = 1.d30
      nrfn = 0
      do i = 1,num_ref_clk
         jndx = 1 
         call get_cmd(ref_clk_code(i), cf_codes, num_cfiles,
     .                    jel, jndx)
         if( jel.gt.0 ) then
*            Check RMS
             if( save_rms(jel).lt.min_ref_rms ) then
                 min_ref_rms = save_rms(jel)
             end if
             if( save_rms(jel).lt. rms_ref_clk/ref_fact ) then
*                Add to final list
                 nrfn = nrfn + 1
                 fin_ref_sites(nrfn) = jel
             end if
         end if
      end do
*     See if any made it through
      if( nrfn.gt.0 ) then
          num_ref_clk = nrfn
           do i = 1, num_ref_clk
              ref_clk_num(i) = fin_ref_sites(i)
              ref_clk_code(i) = cf_codes( ref_clk_num(i))
           end do
      else
*          None made, it reset the rms tolerance and re-assess
*          below
C          rms_ref_clk = 1.1*min_ref_rms
C          assign_ref = .true.
C          ns = num_ref_clk
C          num_ref_clk = 0
*          Scan all site clocks and find the best five.
           do i = 1,5
              best_cfs(i) = 0
              best_rms(i) = 1e20
           end do
*          Now loop over rms values
           do i = 1, num_cfiles
              j = 0
              do while ( j.lt.5 ) 
                 j = j + 1
*                Site i is netter than the Jth entry so add
                 if( save_rms(i).lt.best_rms(j) ) then
                     do k = 5,j+1,-1
                        best_rms(k) = best_rms(k-1)
                        best_cfs(k) = best_cfs(k-1)
                     end do
                     best_rms(j) = save_rms(i)
                     best_cfs(j) = i
                     j = 10
                 end if
              end do
           end do
*          OK Assign new reference clocks
           num_ref_clk = 5
           do i = 1, num_ref_clk
              ref_clk_num(i) = best_cfs(i)
              ref_clk_code(i) = cf_codes(best_cfs(i))
           end do
           write(ins,166) num_ref_clk,(ref_clk_code(i),
     .             best_rms(i)/fClk*1.d9, i=1,num_ref_clk)
 166       format('AUTO REFERENCE :Num ',i3,
     .            5(1x,a4,' RMS (ns)',f10.2),' cycles') 
           cref_upd = .true. 
                     
      end if
      
C     do i = 1, num_ref_clk
C        jndx = 1
C        call get_cmd(ref_clk_code(i), cf_codes, num_cfiles,
C    .                    jel, jndx)
C        ref_clk_num(i) = jel
C     end do

***   If we have no reference clocks, report what happened
      if( num_ref_clk.eq.0 ) then
          write(ins,167) min_ref_rms, rms_ref_clk/ref_fact/fClk*1.d9
 167      format('NO REFERENCE FOUND: MIN RMS ',f10.2,
     .           ' NEW REF RMS ',f10.2,' ns')
                      
          do i = 1, ns
             jndx = 1
             call get_cmd(ref_clk_code(i), cf_codes, num_cfiles,
     .                    jel, jndx)
             if( jel.gt.0 ) then
                 write(*,168) 'REF', ref_clk_code(i), 
     .             save_drms(jel)/fClk*1.d9, save_rms(jel)/fClk*1.d9
 168             format('CLK ',a,' Site ',a4,' dRMS, RMS ',2f10.2,
     .                  ' ns')
             endif
          enddo
          do i = 1,num_cfiles
             if( save_rms(i).lt.rms_max_clk ) then
                 write(*,168) 'ALL', cf_codes(i), 
     .             save_drms(i)/fClk*1.d9, save_rms(i)/fClk*1.d9
             endif
          enddo
      end if

*     Get the RMS of the station clocks
      ns = trimlen(igs_clk_file)
      write(ins,170) igs_clk_file(1:ns), sv_gnss
 170  format(/,'CLOCK STATISTICS FOR ',a,1x,' SYS ',a1,/,
     .       'CLK Site  Number     Mean Offset (cyc)',
     .       ' dRate (D14)  RMS Diff(ns) RMS Fit (ns) STATUS')

      num_out_clk = 0

      do i = 1, num_cfiles
         sum_var = 0
         num_var = 0
         sum_mean = 0
         sum_dmean = 0
         do ep = 1, num_ep-1
            if( par_flag_cse(i,ep).eq.0 .and.
     .          par_flag_cse(i,ep+1).eq.0 ) then
                num_var = num_var + 1

* MOD TAH 0001718: Remove the average slope while computing the
*               mean
                sum_mean = sum_mean + 
     .                    (params_cse(i,ep)-
     .                     save_clk(2,i)*(ep*sampling-save_epc(i)))
                sum_dmean = sum_dmean + (params_cse(i,ep)-
     .                               params_cse(i,ep+1))
                sum_var = sum_var + (params_cse(i,ep)-
     .                               params_cse(i,ep+1))**2
c               write(*,998) i,ep, params_cse(i,ep), par_flag_cse(i,ep)
c998            format('CLK: ',i3,i5,1x,F20.2,1x,i6)
           end if
         end do

         if( num_var.gt.0 ) then
            mean_site(i) = sum_mean/num_var

* MOD TAH 000718: Changed test from this rms of differnces to
*           rms of fit to straight line.  (Computed in fit_igs_clk) 
C           if( sqrt(sum_var/num_var).lt.rms_ref_clk
            if( save_rms(i).lt. rms_ref_clk/ref_fact .and. 
     .          assign_ref ) then
                num_ref_clk = num_ref_clk + 1
                ref_clk_code(num_ref_clk) = cf_codes(i)
                ref_clk_num(num_ref_clk) = i
            end if
            if( save_rms(i).lt.rms_max_clk ) then
                num_out_clk = num_out_clk + 1
                out_clk_num(num_out_clk) = i
            end if
            ot = ' '
            do j = 1, num_ref_clk
               if( i.eq.ref_clk_num(j) ) ot(1:1) = 'R'
            end do
            do j = 1, num_out_clk
               if( i.eq.out_clk_num(j) ) ot(2:2) = 'O'
            end do
            dmean = sum_dmean/num_var
            rms_dmean = sqrt(abs(sum_var/num_var-dmean**2))
            saveRMSdm(i) = rms_dmean
C           write(ins,180) cf_codes(i), num_var, mean_site(i),
C    .                     (save_clk(2,i)/fClk)*1.d14,
C    .                     (sqrt(sum_var/num_var)/fClk)*1.d9,
C    .                     (save_rms(i)/fClk)*1.d9, ot
C180        format('DCK ',a4,1x,i7,1x,F20.2,1x,F10.2,1x,F10.3,1x,
C    .              f10.3,2x,a2)
         end if
      end do

****  Now write out the reference clock codes
      write(200,220) num_ref_clk
 220  format(i6,54x,'# OF CLK REF')
      do i = 1, num_ref_clk
         ns = ref_clk_num(i)
         write(200,230) cf_codes(ns), long_names(ns)(1:20)
 230     format(a4,1x,a20,35x,'ANALYSIS CLK REF')
      end do

****  Now write out the receiver clocks to be output
      write(200,240) num_out_clk, cf_antmod_snx(1:5)
 240  format(I6,4x,a5,'     ',40x,'# OF SOLN STA / TRF')
      do i = 1, num_out_clk
         ns = out_clk_num(i)
         call get_mfs(cf_codes(ns), mfs)
         if( mfs.gt.0 ) then
             call set_np(mfs)
             do j = 1, 3
                site_llr(j) = cf_apr_save(j,ns) + 
     .                        mf_adjust(site_np+j-1)
             end do
         else
            do j = 1,3
               site_llr(j) = cf_apr_save(j,ns) 
            end do
         end if

*        Now convert to XYZ in mm
         jsite_xyz(1) = site_llr(3)*1.d6*cos(site_llr(1))*
     .                                   cos(site_llr(2))
         jsite_xyz(2) = site_llr(3)*1.d6*cos(site_llr(1))*
     .                                   sin(site_llr(2))
         jsite_xyz(3) = site_llr(3)*1.d6*sin(site_llr(1))
         do j = 1,3
            write(posl(j),'(f12.0)') jsite_xyz(j)
         end do

         write(200,250) cf_codes(ns), long_names(ns)(1:20),
     .                  (posl(j)(1:11),j=1,3)
         
 250     format(a4,1x,a20,a11,1x,a11,1x,a11, 'SOLN STA NAME / NUM')
      end do

****  Now output number of satellites
      write(200,260) num_sat 
 260  format(i6,54x,'# OF SOLN SATS')
      do i = 1, num_sat, 15
         if( i+15.le.num_sat ) then
            ls = 15
         else
            ls = mod(num_sat,15)
            if( ls .eq. 0 ) ls = 15
         end if
* MOD TAH 200628: Added the correct GNSS designation to ouput
         write(line,270) (cf_gnss,prn_list(i+j-1),j=1,ls)
 270     format(15(a1,i2.2,1x))
         line(61:) = 'PRN LIST'
         write(200,'(a)' ) line(1:69)
      end do

      write(200,280)
 280  format(60x,'END OF HEADER')

* MOD TAH 030122: Initialize the statistics of the fits
      do i = 1, num_cfiles+num_sat
         do j = 1,3
            fin_stat(j,i) = 0.d0
         end do
      end do
*
*     Set the mask so that we count the number of data for each
*     clock estimate
      call set_phs_mask( phs_mask, phs_bias_mask )

      do ep = 1, num_ep

*        Write the receiver clocks
         mjd = mjd_start + (ep-1)*cf_inter/86400.d0 + 2.d-7/86400.d0
         call mjd_to_ymdhms(mjd, date, sectag)
         dt  = (mjd - mjd_start)*86400.d0

*        Get average reference clock
         num_var = 0
         mean_clk = 0
         sum_var  = 0
         do i = 1, num_ref_clk
            ns = ref_clk_num(i)
            if( par_flag_cse(ns,ep).eq.0 ) then
                num_var = num_var + 1
                dclk =  params_cse(ns,ep)-
     .                  (mean_site(ns)+
     .                   save_clk(2,ns)*(ep*sampling-save_epc(ns)))
                mean_clk = mean_clk + dclk 
                sum_var  = sum_var +  dclk**2
            end if
         end do
         if( num_var.gt.0 ) then
             mean_clk = mean_clk/num_var
             rms_clk  = sqrt(abs(sum_var/num_var-mean_clk**2))
         else
             mean_clk = 0
             rms_clk  = 100
         end if

         do i = 1, num_out_clk
            ns = out_clk_num(i)
            clk_out = (params_cse(ns,ep)-mean_clk)/fClk + 
     .                 apr_clk_poly(2,ns)*dt

*           Count number of values being used to compute clock.
            nc = 0
            do j = 1, actual_max_chan
               if( data_OK(data_flag_cse(j,ns,ep),0,phs_mask)) nc = nc+1
            end do
            if( nc.gt.0 ) then
                clk_sig = 1.d-9/sqrt(1.d0*nc)
            else
                clk_sig = 200.d-9
            end if

* MOD TAH 030122: If this is the range solution, save the clock
*           parameters back for allinging the phase clocks
            if( type(1:1).eq.'R' ) then
               params_cse(ns,ep) =  params_cse(ns,ep)-mean_clk
            end if

            if( .not.kbit(par_flag_cse(ns,ep),1) ) then
                if ( ep-1-
     .               int((ep-1)/igs_clk_samp)*igs_clk_samp.eq.0 ) then
                   write(200,320) cf_codes(ns), date, sectag,
     .                    1,  clk_out, clk_sig
 320               format('AR',1x,a4,1x,i4,4i3.2,1x,f9.6,1x,i2,
     .                2x,2(1x,E19.12)) 
                endif
                dclk = clk_out - 
     .                 (apr_clk_poly(2,ns)+save_clk(2,ns)/fClk)*dt
                fin_stat(1,ns) = fin_stat(1,ns)+1
                fin_stat(2,ns) = fin_stat(2,ns)+dclk
                fin_stat(3,ns) = fin_stat(3,ns)+dclk**2

            end if
         end do

*        Now do the satellite clocks
         do i = 1, num_sat
            kp = i+num_cfiles
            clk_out = (params_cse(kp,ep) - mean_clk +
     .                 svcL1_ce(i,ep))/fClk

* MOD TAH 030215:Now count up how many stations can see this satellites
            nc = 0
            do j = 1, num_cfiles
               ch = ltoc(ctol_cse(1,j,ep), i, actual_max_chan)
               if( ch.gt.0 ) then
                  if ( data_OK(data_flag_cse(ch,j,ep),0,phs_mask) )
     .                                                        nc = nc+1
               endif
            end do
            if( nc.gt.0 ) then
                clk_sig = 1.d-9/sqrt(1.d0*nc)
            else
                clk_sig = 200.d-9
            end if

* MOD TAH 030122: If this is the range solution, save the clock
*           parameters back for allinging the phase clocks
            if( type(1:1).eq.'R' ) then
               params_cse(i+num_cfiles,ep) =  
     .                     params_cse(i+num_cfiles,ep)-mean_clk
            end if

            if( .not.kbit(par_flag_cse(i+num_cfiles,ep),1) ) then
               if( ep-1-
     .             int((ep-1)/igs_clk_samp)*igs_clk_samp.eq.0 ) then
* MOD TAH 200628: Add correct GNSS designation
                  write(200,340) cf_gnss,prn_list(i), date, sectag,
     .                       1, clk_out, clk_sig
 340              format('AS',1x,a1,i2.2,2x,i4,4i3.2,1x,f9.6,1x,i2,
     .                   2x,2(1x,E19.12))
                endif
                kp = i+num_cfiles
                dclk = clk_out - svcL1_ce(i,ep)/fClk
     .                 -save_clk(2,kp)/fClk*dt
                fin_stat(1,kp) = fin_stat(1,kp)+1
                fin_stat(2,kp) = fin_stat(2,kp)+dclk
                fin_stat(3,kp) = fin_stat(3,kp)+dclk**2
C               if( type.eq.'P' ) then
C                  write(*,998) ep, prn_list(i), dclk*1e9,
C    .                  svcL1_ce(i,ep)/fClk*1e9,
C    .                  save_clk(2,kp)/fClk*dt*1e9
C998               format('FIN ',i4,1x,'G'i2.2,1x,3e18.6)
C               endif
            end if
         end do
      end do

****  Finish up clock statistics
      do i = 1, num_cfiles
         ot = ' '
         if( fin_stat(1,i).gt.0 ) then

            do j = 1, num_ref_clk
               if( i.eq.ref_clk_num(j) ) ot(1:1) = 'R'
            end do
            do j = 1, num_out_clk
               if( i.eq.out_clk_num(j) ) ot(2:2) = 'O'
            end do
            num_var = fin_stat(1,i)
            dmean = fin_stat(2,i)/num_var
            rms_dmean = sqrt(abs(fin_stat(3,i)/num_var-dmean**2))
            write(ins,420) cf_codes(i), num_var, mean_site(i),
     .                     (save_clk(2,i)/fClk)*1.d14,
     .                     (rms_dmean)*1.d9, 
     .                     (save_rms(i)/fClk)*1.d9, ot
 420        format('CLK ',a4,1x,i7,1x,F20.2,1x,F10.2,1x,F10.3,1x,
     .              f10.3,2x,a2)
         else
*           Write stats on clocks that were so poor that we did
*           not output their values
            write(ins,430) cf_codes(i), num_var, mean_site(i),
     .                     (save_clk(2,i)/fClk)*1.d14,
     .                      saveRMSdm(i)*1.d9/fClk, 
     .                     (save_rms(i)/fClk)*1.d9, ot
 430        format('CLK ',a4,1x,i7,1x,F20.2,1x,F10.1,1x,F10.1,1x,
     .              F10.1,2x,a2)
         endif
      enddo
      do i = 1, num_sat
         kp = i + num_cfiles
         ot = ' O'    ! All the satellite clocks are output
         if( fin_stat(1,kp).gt.0 ) then
            num_var = fin_stat(1,kp)
            dmean = fin_stat(2,kp)/num_var
            rms_dmean = sqrt(abs(fin_stat(3,kp)/num_var-dmean**2))
            write(ins,440) prn_list(i), num_var, dmean*1e9,
     .                     (save_clk(2,kp)/fClk)*1.d14,
     .                     (rms_dmean)*1.d9, 
     .                     (save_rms(kp)/fClk)*1.d9, ot
 440        format('CLK  G',i2.2,1x,i7,1x,F20.2,1x,F10.2,1x,F10.3,1x,
     .              f10.3,2x,a2)
         endif
      enddo


****  Thats all
      close(200)
      return
      end






