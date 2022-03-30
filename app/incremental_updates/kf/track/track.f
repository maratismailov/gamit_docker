       program track

      implicit none

*     GPS Kinematic analysis program based on the makexk/maked software
*     of Gang Chen Ph.D 1998.
 
*     Program to compute the position of kinematic GPS 
*     GANG CHEN   of MIT, 1996

*     The runstring of the program is
*     % track -f track.bat -log log_file
*     where track.bat is a command batch file

*     MOD 970502: skip the missing-data epoch
*
      include '../includes/xfile_def.h' 
      include 'track_com.h'
 
* Main pogram variables
* i,j, k -- Loop counters
* ne, tp -- Counter for number of data types used in float and
*     data analyses and the total number of parameters needed for
*     a smoothing run.
* num_xyz  -- Number epochs in averages static position
* tot_iter -- Total number of iterations
* remap_iter -- Remapping of BF counter

      integer*4 i,j,k,  luo,  ierr, num_xyz, uni_otl,
     .          ist, iend, dep,  last_min, ne, tp, ks, l,
     .          tot_iter, remap_iter, it

      integer*4 iter_res   ! Number of iterations in the resolve code
     ,,         inc_res    ! Incremental number of resolved ambiguities

      logical kbit

      logical made_updates ! Set true in remap_bf when ambiguity pointers
                           ! updated.  If none updated, iterated again.

* ISRef_D -- Reference time MJD
* ISRef_R -- Reference time in sampling intervals sub-day part

      integer*4 ISRef_D, ISRef_R
* ISRef_F -- Floating point start value ISRef_R

      real*8 ISRef_F 

* first_start -- Epoch of the earliest starting rinex or xfile
* ref_end     -- Last common ending time for the data
* avg_xyz(3)  -- Averaged position of kinematic stations based
*                on sigma < 1 m and deviation <1 m.
* dpos        -- Position difference to see if site moved

      real*8 first_start, ref_end, avg_xyz(3), dpos

* tot_resolved -- Total number of wide-lane ambiquities resolved in each
*     resolve_wl pass.
* ep  -- Counter over epochs in the data
* trimlen -- Length of string  
* MOD AZ 190305: additional variable for reading in Sun&Moon positions
* pjd_start --  PEP julian date of start epoch

      integer*4 tot_resolved, ep, trimlen, run_time(5)
      real*8 sectag, pjd_start

      real*8 cov_col(max_parm), sol_col, dsol(1), var_force(1), dchi, 
     .       avat, equ_gn(max_parm), scale(max_parm)
      real*8 rot_mat(3,3)  ! Used in XYZ_to_GEOD call.

      integer*4 ipivot(max_parm), np

* all_wls_res -- Logical function which returns true when all the widelanes
*     at a site have been resolved 
* done  -- Set true if we have resolved all ambiquities

      logical all_wls_res, done

* output_file -- Name of output file for position estimates
* format  -- String with header line for residuals
* MOD AZ 190305: additional char variable for site names
* site_name_otl -- site name for OTL
      character*256 output_file, format
      character*64  scratch_file 
      character*4   site_name_otl, lowerc

      integer*4 fmppurge  ! comlib function to delete file

* bf_name  -- String to write with the bias flag report

      character*8 bf_name 

      character*2 dtypes(max_dtype)

****  Initialize global variables
      call init_track

* MOD AZ 190305: For Test Purposes
      pin_holder = 0


***** Get the runstring runstring
 
      call get_svsrun

***** Get parameters from batch file

      call read_bat

****  Open the output summary file
      if( trimlen(sum_file).gt.0 ) then
         open(lus,file=sum_file,status='unknown',iostat=ierr)
         call report_error('IOSYAT',ierr,'open',sum_file,0,
     .                     'TRACK/SUM_FILE')
         if( ierr.ne.0 ) lus = 6
         call systime(run_time,sectag)
         write(lus,120) track_version, run_time
 120     format('SUMMARY FILE: Track Vers ',a,' Run ',
     .          i4,'/',i2.2,'/',i2.2,' ',i2.2,':',i2.2)
         write(lus,125) bat_file(1:trimlen(bat_file)),
     .        runday(1:max(1,trimlen(runday))),
     .        runweek(1:max(1,trimlen(runweek)))
 125     format('BATCH FILE: ',a,' RunDay ',a,' RunWeek ',a)
         write(lus,130) (runstr(i)(1:max(1,trimlen(runstr(i)))),
     .                                                 i=1,max_runstr)
 130     format(' STRING Options ',99(a,1x))

      endif

*     Read the atmospheric delay file if passed by user
      if( trimlen(atm_file).gt.0 ) then
          call read_atm_file
      end if
           
*     Now loop over the data files for head information
      do i = 1, num_site
 
*        Now read the header of rinex file
         obs_lu = 50+ i  
         if(obs_file_type(i).eq."R") 
     .        call read_rinex_head(obs_file(i),debug_start,rxver(i)) 

*        Save the header information from this file 
         call save_head(i)

*        Save the start times
         data_start(i) = xf_start_mjd
         data_end(i)   = xf_end_mjd 
      end do     ! end head processing

*     Get the earliest start time
      if( usr_start.eq.0 ) then 

          first_start = min(data_start(1),data_start(2))
          ref_start = max(data_start(1),data_start(2))
          do i = 3, num_site
             if( data_start(i).lt. first_start ) then
                 ref_start = first_start
                 first_start = data_start(i)
             else if( data_start(i).lt. ref_start ) then
                 ref_start = data_start(i)
             end if
          end do

* MOD TAH 070108: Check that start times are possible
          ISRef_D = int(ref_start)
          ISRef_F = (ref_start - ISRef_D)/(usr_interval/86400.d0)
          ISRef_R = nint(ISRef_F) 
          if( abs(ISRef_F-nint(ISRef_F)).gt. 1.d-3 ) then
              write(*,140) ref_start, 
     .                     ISRef_D + ISRef_R*usr_interval/86400.d0
 140          format('**WARNING** Resetting start time from ',
     .               F16.8,' MJD to ',F16.8,' MJD to match interval')
              ref_start = ISRef_D + ISRef_R*usr_interval/86400.d0
          end if
          
      else
          ref_start = usr_start
      end if
* MOD AZ 190305: Set up Sun/Moon Tables for interpolation later

      pjd_start = ref_start + 2400001.d0 + (19.0d0 + 32.184d0)/86400.d0
      
      call ephred(27,pjd_start,0,0)


* MOD AZ 190305: Initialize Ocean Tidal Loading BLQ File read-in unit
*      write(*,*) 'OKay till now'
* MOD TAH 200225: Make ocean tide optional
      if( use_blq ) then
         do i=1, num_site
  
            site_name_otl=lowerc(site_names(i))

            site_otl(i) = site_name_otl // '.BLQ'
         
            write(*,*) 'Reading in BLQ files: ', site_otl(i)
 
            uni_otl=170+i
        
            open(uni_otl,file=site_otl(i),status='old',iostat=ierr)

            call report_error('IOSTAT',ierr,'open',site_otl(i),1,
     .                      'otlcal/init_read_blq')

         end do
      end if
*      write(*,*) 'OKay till now'
* MOD AZ 190305: End of Initialize Ocean Tidal Loading BLQ

*     Read the DCB file 
      call read_dcbs(dcb_file, ref_start)

*     Now compute the number of seconds from start of SP3 file
      ref_sec = (ref_start - sp3_refmjd)*86400.d0

      write(*,150) ref_start, ref_sec
 150  format('Starting epoch count for MJD ',F16.8,
     .       ' Start seconds ',F20.8)

****  Now read in the rinex files
      num_prn = 0
      num_epochs = 0 
      num_ambs   = 0

****  Before proceeding check the interval.
      if( usr_interval.le.0 ) then
          write(*,160) 
 160      format('**DISASTER** Interval for data is zero. Check RINEX ',
     .           'file intervals or INTERVAL command')
          stop 'TRACK: Rinex file interval is zero'
      end if
 

      do i = 1, num_site
          call read_all_rinex(i)
      end do 
      call report_setup(6)
      if( lus.ne.6 ) call report_setup(lus)

* MOD TAH 180322: Map the RINEX data to common reference frequencies for
*     GNSS processing and set the reference frequencies

      call remap_rx 

****  Now get the number of epochs to process based on the last common
*     time
      ref_end = 0
      do i = 1, num_site-1
         do j = i+1, num_site
            ref_end = max(ref_end,min(data_end(i),data_end(j)))
         end do
      end do 

* MOD TAH 200527: Use nint in computing number of epochs to ensure we
*     get the last one.
      num_epochs = nint((ref_end-ref_start)/(usr_interval/86400.d0))+1
      write(*,220) num_epochs, ref_end
 220  format('There are ',i6,' epochs to last ending time MJD ',F16.8)
      if( usr_nepochs.ne.0 ) then
          write(*,240) usr_nepochs
 240      format('User specified ',i6,' epochs to process')
          num_epochs = usr_nepochs
      end if

****  Get the coordinates are the start of the run.  In current version
*     1.24 position mapped at start of session (may change later)
      do i = 1, num_site
         do j = 1,3
            site_apr(j,i) = site_int(j,i) + 
     .                   site_vel(j,i)*(ref_start-site_ep(i))/365.25d0
         end do
* Add 1.26: Compute geod coordinates
         call XYZ_to_GEOD( rot_mat, site_apr(1,i), site_geod(1,i))
      end do

***** MOD TAH Version 1.26 111122: If ion delay files have been given
*     read them now.
      do i = 1, num_site
         if( use_ionlos(i) ) call read_ionlos(i)
      enddo

****  Run an initial position solution.  We use the first Fixed site
*     as reference, and compute the positions of the kinematic sites
*     using pseudorange data.  This allows us to compute the elevation
*     angles at the sites which are used to weight the estimates of the
*     widelane measurements.
                            
      call setup('EXATM')

****  Generate an approximate trajectory using the Psuedo-range
*     data.
      data_mask = 2
      kine_known = .false.
      do ep = 1, num_epochs
C          call est_pos(ep, 'PC', 'None NoWrite RedE EDIT', 0, 0 )
C MOD tah 050120: Removed edit option: Throws out all data when station
C         is actually moving
          call est_pos(ep, 'PC', 'None NoWrite RedE ', 0, 0 )
      end do

*     Check out the nature of the Kinematic site.  If it does not
*     seem to be moving then set to a fixed position
      call pos_anal('PC range')

      
****  Now re-sort the channels at each epoch so that they are
*     are decreasing order of elevation angle
      call reorder
      
*     Mark the kinematic trajectory as known so that on the next
*     pass we will use the estimate of positions we alreay have.
      kine_known = .true.
 
*     OK, scan all the exta-wide lanes to see if any slips
      call scan_ddlg
      call scan_sdlg

* MOD TAH 071031: Remove any bias flags requested by user
      call del_usr_bf

****  Now that we know an approximate trajectory, setup up the
*     bias flags accounting for the elevation angle cut off.
      data_mask = 18
      do i = 1, num_site
         call flag_gaps(i)
      end do

****  Now we need to look at the resolving the wide lane 
*     ambiquities (ie. L1-L2 which should be able to get from the
*     Melboure-Wubbena (MW) wide lane and the ionospheric constraint.
*     Mark all the first site ambiquities as resolved and then
*     loop over all sites above the first and find the longest, earliest
*     bias flag.

      call select_arb_bf

*     Now we resolve the widelane ambiquities for each site
*     Initialize the saved values of the widelines
      do i = 1, num_ambs
         wls_all(1,i) = 0.d0
         wls_all(2,i) = 0.d0
         wls_ref(1,i) = 0
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

      call report_bf(6,'INITIAL')
c      if( lus.ne.6 ) call report_bf(lus,'INITIAL')
      
****  See if user have given an input ambiquiuty file (check done
*     inside routine)
      call read_bf

****  OK, Now try to search for the best integer 
      data_mask = 18
      num_amb_save = 0

****  Now loop over all the ambiquities find those that need to
*     be resolved.  Go through and mark all the fixed biases as
*     finally resolved
      do i = 1, num_ambs
         if( wls_ref(1,i).le.0 ) then
             call sbit(bf_ents(5,i),2,1)
         end if
      end do

* MOD TAH 160413: See if we can remove cycle slips
      if( num_mwwl.gt.0 ) then
         print *,'CALLING CLSIP_repair for NUM MWWL ',num_mwwl
         call cslip_repair
      endif


      ! stop 'CSLIP'

*     Now search for ones to search and resolve.
      last_min = 0
      done = .false. 

*     Initialize the saved rank to 0.0 for all ambiquitues.  If ambiquity
*     can not be reliably resolved the best ranked values will be used.
      do i = 1, num_ambs
         saved_rank(i) = 0.0
      end do


****  Now scan all the bias flags again and any which are
*     still un-resolved try them again
c      kine_known = .false.

      num_tot_resolved = 1
      it = 0
      remap_iter = 0
      do while (num_tot_resolved.gt.0 .and. it.lt.100)
         it = it + 1

*        Reset the number of resolved ambiquities so that we stop
*        when none are resolved.
         num_tot_resolved = 0

****     See if we have are the iteration to do a float solution
*        for the biases
         if( it.ge. float_iter .and. float_iter.gt.0 ) then
             data_mask = 18 
             write(*,320) it
 320         format(/,'Starting FLOAT iteration at ',i3)
             call setup('FLOAT')
             
             do ep = 1,  num_epochs, float_sample
                 call est_pos(ep, float_type, 'KALF EDIT FLOAT NoWrite', 
     .                        float_sample, 1 )
             end do

* MOD TAH 100520: See if we are applying post sigmas
             do i = 1, num_site
                do j = 1, 3
                   if( pst_site(j,i).ge.0 .and. 
     .                 pos_parn(j,i).gt.0 ) then
                      np = pos_parn(j,i)
                      var_force(1) = pst_site(j,i)
                      dsol(1) = 0.d0
                      call force_track( cov_parm, sol_vec, max_parm, 
     .                      num_parm,cov_col, sol_col, np, 1,  dsol,
     .                      var_force, dchi, avat, equ_gn, 
     .                      ipivot, scale)
                   end if
                end do
             end do

*            Save any bias remaining at end
             call float_end(float_type)
 
             write(*,330) num_dd_avg, sqrt(rms_dd_avg/num_dd_avg)*1000
 330         format('For ',i8,' Double differences: Average RMS ',
     .              F10.2,' mm')

*            Now see if we can resolve any ambiquities
             inc_res = 1
             iter_res = 0
             do while ( inc_res.gt.0 )
                iter_res = iter_res + 1 
                call resolve_float(it, iter_res, inc_res)
             end do

****         MOD TAH 1.26: Set the ion delay known if we are using
*            and ionex file
             if ( use_ionex ) ion_known = .true.

         end if
*        Report the final selected biases
         write(bf_name,'(a,i2.2)') 'FLOAT',it      
         call report_bf(6,bf_name)

* MOD TAH 100214: See if we should re-map the ambiguities that are
*        not fixed due to 'O' ambiguities.  If remapping is possible
*        num_tot_resolved will return the number remapped.
         if( num_tot_resolved.eq.0 ) then
            made_updates = .false.
            do while( .not.made_updates )
               write(*,340) bf_name
 340           format('After ',a,' remapping the ambiguities')
               remap_iter = remap_iter + 1
               call remap_bf(bf_name, remap_iter, made_updates)
            end do
         end if

*        Re-compute the wls (may have been changed by remap_bf
         call recomp_wls

         call pos_anal(bf_name)

      end do 
       
      kine_known = .true.

      call report_bf(6,'FINAL')
      if( lus.ne.6 ) call report_bf(lus,'FINAL')

****  Run solution with all data if float_sample not-equal 1
      if( float_sample.ne.1 ) then
         write(*,'(a)') 'Running full rate solution'
         call setup('FLOAT')
         do ep = 1,  num_epochs
             call est_pos(ep, float_type, 'KALF EDIT FLOAT NoWrite', 
     .                    float_sample, 1 )
         end do
      end if


****  Compute final MW and EX widelandes
      call check_wl


*     OK, Try an L1 solution.  Set the data mask to use the "bad"
*     data flag (Bit 2) and the elevation cut-off (Bit 5)
      data_mask = 18

      do j = 1, num_anal_type
          call setup('FLOAT')

          call clear_stats

*         Setup the output file name.
          do k = 1, num_site
              if( index(posit_type,'GEOD').gt.0 ) then 
                  output_file = posit_root(1:trimlen(posit_root)) // 
     .                      '.GEOD.' // site_names(k) //
     .                      '.' // anal_types(j)
                  luo = 21 + k
                  open(luo, file=output_file, status='unknown', 
     .                 iostat=ierr)
                  write(*,*) 'Creating ',
     .                       output_file(1:trimlen(output_file))
                  write(luo,520)
520               format('*YY  DOY        Seconds        Latitude',
     .                   '     Longitude      Height   SigN  SigE ',
     .                   ' SigH   RMS  #     Atm       +-       ',
     .                   'Fract DOY     Epoch  #BF NotF',/,
     .                   '*                                (deg) ',
     .                   '      (deg)          (m)     (cm)  (cm)',
     .                   '  (cm)  (mm)  DD    (mm)     (mm)')
              end if
              if( index(posit_type,'NEU').gt.0 ) then 
                  output_file = posit_root(1:trimlen(posit_root)) // 
     .                      '.NEU.' // site_names(k) //
     .                      '.' // anal_types(j)
                  luo = 21 + num_site + k
                  open(luo, file=output_file, status='unknown', 
     .                 iostat=ierr)
                  write(*,*) 'Creating ',
     .                       output_file(1:trimlen(output_file))
                  write(luo, 540)
 540              format('* YY  MM DD HR MIN     Sec          dNorth',
     .                   '      +-            dEast        +-    ',
     .                   '      dHeight      +-        RMS    #  ',
     .                   '    Atm     +-         Fract DOY     Epoch',
     .                   '  #BF NotF  Rho_UA',/,
     .                   '*                                    (m)  ',
     .                   '      (m)            (m)        (m)    ',
     .                   '     (m)           (m)      (mm)   DD  ',
     .                   '    (mm)    (mm)')
              end if
              if( index(posit_type,'XYZ').gt.0 ) then 
                  output_file = posit_root(1:trimlen(posit_root)) // 
     .                      '.XYZ.' // site_names(k) //
     .                      '.' // anal_types(j)
                  luo = 21 + 2*num_site + k
                  open(luo, file=output_file, status='unknown',
     .                   iostat=ierr)
                  write(luo, 560)
 560              format('* YY  MM DD HR MIN     Sec            dX  ',
     .               '      +-              dY         +-    ',
     .               '         dZ        +-        RMS    #  ',
     .               '    Atm     +-         Fract DOY     Epoch',
     .               '  #BF NotF  Rho_UA',/,
     .               '*                                    (m)  ',
     .               '      (m)            (m)        (m)    ',
     .               '     (m)           (m)      (mm)   DD  ',
     .               '    (mm)    (mm)')
              end if
              if( index(posit_type,'DHU').gt.0 ) then 
                  output_file = posit_root(1:trimlen(posit_root)) // 
     .                      '.DHU.' // site_names(k) //
     .                      '.' // anal_types(j)
                  luo = 21 + 3*num_site + k
                  open(luo, file=output_file, status='unknown',
     .                       iostat=ierr)
                  write(luo, 580)
 580              format('* YY  MM DD HR MIN     Sec          dNorth',
     .               '      +-            dEast        +-    ',
     .               '      dHeight      +-        RMS    #  ',
     .               '    Atm     +-         Fract DOY     Epoch',
     .               '  #BF NotF  Rho_UA',/,
     .               '*                                    (mm) ',
     .               '      (mm)           (mm)       (mm)   ',
     .               '     (mm)          (mm)     (mm)   DD  ',
     .               '    (mm)    (mm)')
              end if

*             Now set up the residual output files (one per station/satellite)
              if( trimlen(resid_root).gt.0 ) then
                  write_res = .true.
                  call get_num_dtype(anal_types(j),num_dtype,
     .                               dtypes)

*                 Only output residuals for kinematic sites
                  if( site_type(k).ne.0 ) then 
                     do i = 1, num_prn
                        luo = 120+(k-1)*max_sat+prn_used(i)
                        write(output_file,610) 
     .                        resid_root(1:trimlen(resid_root)),
     .                        site_names(k), 
     .                        prn_used(i), anal_types(j)
 610                    format(a,'.',a4,'.PRN',i3.3,'.',a)
                        open(luo,file=output_file,status='unknown',
     .                       iostat=ierr)
                        call report_error('IOSTAT',ierr,'creat',
     .                       output_file, 1,'TRACK/Residual files')
                        call gen_res_format(num_dtype,'label',format)
                        write(luo,'(a)') format(1:trimlen(format))
                     end do
                  end if  
              end if  


          end do
*         Test of backward direction solutuion
*         Test of backward direction solutuion   
          if( index(back_type,'BACK').gt.0 ) then

*            If we are also smoothing then open the
*            direct access file that we will need to save
*            the intermediate covariance matrices
             if( index(back_type,'SMOOTH').gt.0 ) then

*                Check to see if the number of ambiquities
*                in the float solution is less than then
*                number needed now e.g., LC used in float,
*                but now using L1+L2
                 ne = 0
                 tp = tot_parm

                 if( index(float_type,'L1').gt.0 .or.
     .               index(float_type,'LC').gt.0 ) ne = ne + 1
                 if( index(float_type,'L2').gt.0 ) ne = ne + 1
                 if( ne.eq.1 ) then
*                    See if more types are needed
                     ne = 0 
                     if( index(anal_types(j),'L1').gt.0 .or.
     .                   index(anal_types(j),'LC').gt.0 ) ne = ne + 1
                     if( index(anal_types(j),'L2').gt.0 ) ne = ne + 1
                     if( ne.gt.1 ) then
                         tp = non_amb_parm + (tot_parm-non_amb_parm)*2
                         write(*,620) tp, tot_parm
 620                     format('Increasing total parameters to ',i4,
     .                          ' from ',i4,' used in float analysis')
                     end if
                 end if

                 darecl = (tp+1)*tp + 3*tp
                 tac_parm = tp
                 scratch_file = 'scratch' // runstr(1) 
                 open(119,file=scratch_file,access='direct',
     .                recl=4*darecl,form='unformatted')
                 call report_error('IOSTAT',ierr,'open',
     .                'Direct Access File',1,'track')
             end if              
             do ep = num_epochs, 1, -1
                 if( index(back_type,'SMOOTH').gt.0 ) then
                     call est_pos(ep, anal_types(j), 
     .                    'KALF EDIT FLOAT SMOOTH NoWrite',  -1, 1 )
                 else
                     call est_pos(ep, anal_types(j), 
     .                    'KALF EDIT FLOAT',  -1, 1 )
                 end if
             end do
             call float_end(anal_types(j))
          end if

*         Run the forward solution
          call setup('FLOAT')

          pin_holder = 1   

*         Run the forward solution 
          do ep = 1, num_epochs, 1
              if( index(back_type,'SMOOTH').gt.0 ) then
                  call est_pos(ep, anal_types(j),
     .                     'KALF EDIT FLOAT SMOOTH', +1, 1 )
              else
                  call est_pos(ep, anal_types(j),
     .                     'KALF EDIT FLOAT', +1, 1 )
              endif
          end do 

          call float_end(anal_types(j))

****      See if we should output coorinates for non-stochastic
*         kniematic sites
          do k = 1, num_kine
             ks = kine_to_site(k)
             if( mar_site(1,ks).eq.0 .and. mar_site(2,ks).eq.0 .and.
     .           mar_site(3,ks).eq.0 ) then
                 write(*,720) site_names(ks), 
     .                       anal_types(j)(1:trimlen(anal_types(j))),
     .                       (kine_out(l,k),l=1,3)
 720             format(' MEAN COORDINATES ',a4,' Data type ',a, 
     .                  ' XYZ ',3F15.4)
             end if
          end do

****      Make summary of positions.
          do k = 1, num_kine
             ks = kine_to_site(k)
             do l = 1,3
                avg_xyz(l) = 0.d0
                num_xyz    = 0
             end do
             do ep = 1, num_epochs
                if( kbit(kine_OK(1,k),ep) ) then
                   if( num_xyz.eq.0 ) then
                       do l = 1,3
                          avg_xyz(l) = kine_xyz(l,k,ep)
                       end do
                       num_xyz = 1
                   else
                       dpos = 0
                       do l = 1,3
                          dpos = dpos + (kine_xyz(l,k,ep)-avg_xyz(l))**2
                       end do
                       if( dpos.lt.1.d0 ) then
                          do l = 1,3
                             avg_xyz(l) = (num_xyz*avg_xyz(l)+
     .                            kine_xyz(l,k,ep))/(num_xyz+1)
                          end do
                          num_xyz = num_xyz + 1
                       end if
                   endif
                end if
                if( ep.ge.debug_start .and. ep.le.debug_end ) then
                    write(*,995) ep, num_xyz, (avg_xyz(l),l=1,3),
     .                           (kine_xyz(l,k,ep),l=1,3)
 995                format('FINPOS ',i6,i6,3(F13.3),' Kine ',3(F13.3))
                endif 
             end do
*            Write out the averaged STATIC position
* MOD TAH 080626: Changed ks to k in static since only kinematic sites
             if( num_xyz.gt.1 .and. .not.static(k) ) then
                write(*,810 ) site_names(ks), (avg_xyz(l),l=1,3), 
     .              num_xyz
 810            format(5x,a4,2x,3(F12.3,1x),' STATIC Position from ',
     .               i6,' epochs')
             end if
          end do

*         Close the output files
          do k = 1, num_kine
             close(21+k)
             close(21+num_kine+k)
          end do
* MOD AZ 190305: Close all the BLQ files
          do k = 1, num_site
             close(170+k)
          end do

          close(119)

****      Finishup by writing out the phase statistics from the run
          call phase_stats(6,anal_types(j))
          if( lus.ne.6 ) call phase_stats(lus,anal_types(j))

      end do

****  Write out final wide lane values
      if( trimlen(wls_root).gt.0 ) then
         print *,'Output processed Wide-lanes with ',trim(wls_root)
         do k = 1, num_site
              call output_wls(k,'WLS')
          end do
      end if

****  See if DUMP output is requested.
      if( index(posit_type,'DUMP').gt.0 ) then
*        Dump out the one-way phase and range residuals
         do k = 1, num_site
            call dump_ows(k)
         enddo
      end if

*     Clean up
      if( index(back_type,'SMOOTH').gt.0 ) then
         ierr = fmppurge(scratch_file)
         call report_error('IOSTAT',ierr,'delet',scratch_file,0,'track')
      end if
      
****  Thats all
      end


