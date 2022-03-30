CTITLE INIT_TRACK

      subroutine init_track

      implicit none 

*     Routine to initialise the variable and set the defaults
*     for the track program

      include '../includes/const_param.h'   
      include 'track_com.h'
 
* LOCAL VARIABLES
* MOD AZ 190305: additional variables for path to tables
*                They are introduced because of the unique structure
*                of the HPC I'm using. Highly subject to change.
* trimlen, home_dir_e -- variables for retrieving home path
* i, j  -- loop counters
      integer*4 i,j,trimlen
      character*256 home_dir_e
****  OK, Initialize the variables which need it
      ref_start = 0
      usr_interval = 0
      usr_nepochs  = 0


      num_prn = 0
      num_site = 0
      num_epochs = 0 
      num_ambs   = 0
      num_antmods = 0
      num_tdep = 0
      num_mwwl = 0    ! Default not to remove cycle slips (=20 reasonable)
      num_exwl = 10   ! Number of EX WL to average to repair slip
      cslip_tex = 0.15d0
      cslip_tmw = 0.30d0
      cslip_sex = 0.20d0  ! Sigma limit for EX derived dN1 estimate
      cslip_smw = 0.40d0  ! Sigma limit for MW WL difference
      NLx_RR  = 10        ! Relative rank of cslip repair
      NLx_addChi = 2.0    ! Mimimum addative chi**2
      NLx_minChi = 20     ! Max chi**2, Min value can have.

*     Now start setting the defaults     
      max_gap    = 1
      min_good   = 20
      ref_rel_hum = 0.d0

      do i = 1, max_site
         apr_atm(i) = 0.d0
         mar_atm(i) = 0.d0
         mar_atm_hgt(i) = 0.d0
         do j = 1, 3
            apr_site(j,i) = 20.0
            pst_site(j,i) = -1.0   ! Set aposterori variance to show not used
            mar_site(j,i) = 20.0
            site_offarp(j,i) = 0.d0
            sit_L12(j,1,i) = 0.d0
            sit_L12(j,2,i) = 0.d0

         end do
         do j = 1,2
            num_sit_dph(j,i) = 0
         end do

         ant_name(i) = 'NONE'
         rcv_type(i) = 'N'
         atm_scale(i) = .false.  ! Set true to use atm dh scale

      end do

      do i = 1,3
         ref_xyz(i) = 0.d0
      end do

      do i = 1, max_kine
         static(i) = .false.
      end do
      
      do i = 1, max_sat
         data_var(1,i) = (0.003d0)**2
         data_var(2,i) = (0.003d0)**2
         data_var(3,i) = 1.0d0
         data_var(4,i) = 1.0d0 
         data_var(5,i) = 0.5d0
         num_svs_dph(1,i) = 0
         num_svs_dph(2,i) = 0
         do j = 1,3
            svs_L12(j,1,i) = 0.d0
            svs_L12(j,2,i) = 0.d0
         end do
         dcbs(i) = 0.d0
      end do 


      num_amb_samp = 20
      max_tot_search = 1000000
      relrank_limit = 20.0d0
      relrank_iter  = 5.d0

      ion_var  = 10.d-0   ! 10**2 cycles constant variance
      ion_corr = 6.8d-0  ! Corresponds to 1 ppm (cycles**2)
      ion_tau  = 300.d3  ! 300 km correlation length
      ion_wght = 0.01d0
      ion_ppm  = 1.d-6
      ion_ht   = 350.d3
      max_ion_jmp = 0.5d0  ! Cycles (1-1 slip is 0.3 cycles)
      rms_edtol = 3.d0     ! 3 times the data noise for the observable

      mwwl_jmp = 5.0d0     ! Jump in MW-WL that will add bias flag.
 
      debug_start = -1
      debug_end = -2

      num_edits = 0
      num_exclude = 0
      num_abf = 0
      num_rbf = 0

      posit_root = 'track'
      sum_file = 'track.sum'
      lus = 98
      resid_root = ' '
      wls_root = ' '
      rwl_root = ' '
      runday = ' '
      runweek = ' '
      posit_type = 'NEU'
      sp3_file = ' '

      write_res = .false. 

      atm_mtt = .true.
      atm_gmf = .false.
* MOD AZ 190305: additional variable for VMF
      atm_vmf = .false.
      use_atm_ref = .false.
      atm_file = ' '
* MOD TAH 202025: Make ocean tide optional (use_blq command to turn on)
      use_blq  = .false. 
c      noetide = .true.  ! Debug
      noetide = .false.  ! Debug
      set_msec_bf = .false.  ! Added version 1.22
      use_ionex   = .false.  ! Set true when ionex_file read OK
      use_ionlos  = .false.  ! Set true when ionlos_file read OK
      ion_known   = .false.  ! Set true when ion tec values known
      ionex_file  = ' '

      ambin_file = ' '
      search_type = 'NONE'
      back_type = 'NONE'
* MOD AZ 190305: retrieve home path for leapsec and sun&moon tables
      call getenv('HOME',home_dir_e)
      leapsec_path = home_dir_e(1:max(1,trimlen(home_dir_e))) //
     .               '/gg/tables/leap.sec'
*      write(*,*) leapsec_path

      nbody_path = home_dir_e(1:max(1,trimlen(home_dir_e))) //
     .               '/gg/tables/nbody'
 
c      write(*,*) nbody_path

      float_type = 'LC'
      float_iter = 1
      float_sample = 1
      float_limit(1) = 0.25
      float_limit(2) = 0.50

      out_sig_limit  = 1.0

      elev_cutoff = 10

      wl_fact  = 1.0
      lg_fact  = 1.0
      max_fit  = 25.0

      num_anal_type = 1
      anal_types(1) = 'LC'
      posit_type = 'NEU'

      stopgo_mode = .false.
      ante_off = .false.   ! Set true when antenna informtion
                ! read.  FATAL if antmod_file command used before
                ! ANTE_OFF 

      time_unit = 'epoch'
      tu_to_ep  = 1.0

      min_lcsig = 0.01
      wl_tau = 600.0d0
      dynamic_tol = 10.0   ! RMS for dynamic site

      lambda(1) = (vel_light/fR1)
      lambda(2) = (vel_light/fR2) 

* MOD TAH 180321: Set default GNSS
      tr_gnss = 'G'

*     Defaults for GPS
      fR1 = 154*10.23d6     
      fR2 = 120*10.23d6  

      call set_combs

****  Thats all
      return
      end

CTITLE SET_COMBS

      subroutine set_combs

      implicit none

*     Routine to set the linear combinations of frequencies to compute
*     MW-WL, LC, LG. PC and TECU values based on the reference L1 and L2
*     frequencies set in fR1 and fR2

      include 'track_com.h'   
 
*     Compute the factors
      dfsf = (fR1-fR2)/(fR1+fR2) 
      sfdf = (fR1+fR2)/(fR1-fR2) 

      lcf1 = 1.d0/(1.d0 - (fR2/fR1)**2) 
      lcf2 = -(fR2/fR1)/(1.d0 - (fR2/fR1)**2)

      lgf1 = -fR2/fR1
      lgf2 = 1.d0 

      exf1 = 1.d0     
      exf2 = -fR1/fR2 

      pcf1 =  fR1**2/(fR1**2-fR2**2) 
      pcf2 = -fR2**2/(fR1**2-fR2**2) 

      l1tecu = 40.3d0/fR1**2*1.d16 
      l2tecu = 40.3d0/fR2**2*1.d16 

****  Thats all
      return 
      end

CTITLE REPORT_SETUP

      subroutine report_setup(unit)

      implicit none 

*     Routine to report the setup for the track run

      include 'track_com.h'

* PASSED VARIABLES
* unit  -- Unit number for output

      integer*4 unit

* LOCAL VARIABLES
      integer*4 i, j, k, ep, trimlen, dates(5), datee(5)
      real*8 secs, sece

* stype(2) -- Type of positioning (S, K)
      character*10 stype(2)

      character*1 sqroot   ! Extended ascii symbol for square root

      data stype / 'STATIC', 'KINEMATIC' /


****  Tell user what is happening
      write(unit,50)
  50  format(/,'TRACK SETUP PARAMETERS',/,
     .         '----------------------')
      write(unit,60) 'Batch File ',bat_file(1:trimlen(bat_file))
      if( trimlen(ambin_file).gt.0 ) 
     .write(unit,60) 'Ambiguity File ',
     .               ambin_file(1:trimlen(ambin_file))
      if( nav_file_type.eq.'SP3' ) then
         write(unit,60) 'SP3 File ',sp3_file(1:trimlen(sp3_file))
 60      format(a,t24,': ',a)
      else
         write(unit,60) 'NAV File ',nav_file(1:trimlen(nav_file))
      endif
      write(unit,60) 'SUMMARY File',sum_file(1:max(trimlen(sum_file),1))

      do i = 1, num_antmods
        write(unit,60) 'ANTEX File  ',
     .                         antmod_file(i)(1:trimlen(antmod_file(i)))
      end do

      if( use_ionex .and. trimlen(ionex_file).gt.0 ) 
     .   write(unit,60) 'IONEX File ',trim(ionex_file)

      do i = 1, num_site
         write(unit,70) site_names(i), stype(site_type(i)+1),
     .                  obs_file(i)(1:trimlen(obs_file(i)))
  70     format('Rinex files',t24,': ',a,2x,'Type ',a,2x,a)
      end do

      write(unit,80) relrank_limit
  80  format('RELATIVE RANK LIMIT    : ',F12.2)
      write(unit,90) float_type(1:trimlen(float_type)), float_iter,
     .                float_sample, float_limit, wl_fact, lg_fact, 
     .                max_fit
  90  format('FLOAT TYPE             : ',a,' Start Iter ',i2,
     .       ' Sample ',i3,' Limits ',2F5.2,' WL/LG SCALES ',2F5.2,
     .       ' MAX Fit ',F6.1)
      if( num_mwwl.gt.0 ) then  ! Report CSLIP
         write(unit,95) num_mwwl, num_exwl, 
     .          cslip_tex, cslip_tmw, cslip_sex, cslip_smw 
 95      format('RM_CSLIP               : Num MWWL ',I3,' Num EXWL ',
     .          I3,' EX/MW abs tol ',2F6.3,' EX/MW Sig Limit ',2F6.3)
      else
         write(unit,'(a)') 'RM_CSLIP option not used'
      endif

      write(unit,100) (anal_types(i)(1:trimlen(anal_types(i)))
     .                ,i=1,num_anal_type)
 100  format('DATA ANALYSIS TYPES    : ',a)
      write(unit,110) (sqrt(data_var(i,1)),i=1,4),
     .       data_var(5,1)
 110  format('Data Noise L1,L2,P1,P2 GENERIC: ',2F8.4,2F8.2,
     .       ' m, Weight ',F5.2)
      do i = 1, max_sat
         if( data_var(1,i).ne.data_var(1,1) ) then
            write(unit,115) i,(sqrt(data_var(j,i)),j=1,4),
     .              data_var(5,i)
 115        format('Data Noise L1,L2,P1,P2 PRN ',i3.2,' : ',
     .             2F8.4,2F8.2,' m, Weight ',F5.2) 
         end if
      end do

      write(unit,'("GNSS SYSTEMS: ",a)') trim(tr_gnss)

      write(unit,120) min_lcsig, wl_tau, dynamic_tol
 120  format('Min LCSig ',F6.3,' cyc, WL TAU ',F8.1,' sec,',
     .       ' Dynamic Tol ',F6.1,' m')

      write(unit,130) out_sig_limit
 130  format('Postion output sigma limit    : ',f10.3,' m.')

      write(unit,140) elev_cutoff
 140  format('Elevation angle cutoff        : ',F10.3,' deg.')
      if( atm_mtt ) then
         write(unit,'(a,1x,a)') 'MTT Mapping function and',
     .                       'seasonal model used'
      elseif (atm_gmf) then
         write(unit,145) 'GPT Pressure/Temp and',
     .                     'GMF mapping function used.', 
     .                     ref_rel_hum
 145     format(a,1x,a,' with ',F5.2,' relative humidity')
      elseif (atm_vmf) then
         write(unit,'(a,1x,a)') 'GPT3 Pressure/Temp and',
     .                       'VMF3 mapping function used'
      else
         write(unit,'(a,1x,a)') 'Troposhere Model',
     .                       'Wrongly Assigned'
      endif

****  Output the site characteristics
      sqroot = 'R'  !   char(164)
      write(unit,200) time_unit(1:trimlen(time_unit)), 
     .             sqroot,time_unit(1:1), sqroot, time_unit(1:1)
 200  format('Process noise time unit: per sqrt(',a,')',/,
     .       'SITE STOCHASTIC PROPERTIES',/,
     .       'SITE Type         Apriori XYZ (m)        ',
     .       'XYZ Process noise (m/',a,a,') Atm Stats (m,m/',a,a,
     .       ') Hgt m/(m/s)')

      do i = 1, num_site

*        See if Kine or Fixed
         if( site_type(i).eq.0 ) then
*            Fixed Site
             write(unit,210) site_names(i), sqrt(apr_atm(i)), 
     .            sqrt(mar_atm(i)*tu_to_ep), sqrt(mar_atm_hgt(i))
 210         format(a4,'    F ', 56x,F7.3,1x,F8.5,1x,F8.5)
         elseif ( atm_scale(i) ) then
*            HEIGHT SCALE ATM not da/(dh/dt) term
             write(unit,220) site_names(i), (sqrt(apr_site(j,i)),j=1,3),
     .                   (sqrt(mar_site(j,i)*tu_to_ep),j=1,3), 
     .                    sqrt(apr_atm(i)), sqrt(mar_atm(i)*tu_to_ep)
 220         format(a4,'    K ',3(F8.3,1x),1x,3(F8.3,1x),1x,
     .                         F7.3,1x,F8.5,1x,' dH SCALE')
         else   ! Normal Atm delay estimation
*            Kinematic site
             write(unit,230) site_names(i), (sqrt(apr_site(j,i)),j=1,3),
     .                   (sqrt(mar_site(j,i)*tu_to_ep),j=1,3), 
     .                    sqrt(apr_atm(i)), sqrt(mar_atm(i)*tu_to_ep), 
     .                    sqrt(mar_atm_hgt(i))
 230         format(a4,'    K ',3(F8.3,1x),1x,3(F8.3,1x),1x,
     .                         F7.3,1x,F8.5,1x,F8.5)
         end if
      end do
      do i = 1, num_site

*        See if Kine or Fixed
         if( site_type(i).ne.0 .and. pst_site(1,i).ge.0 ) then
*            Kinematic site
             write(unit,240) site_names(i), (sqrt(pst_site(j,i)),j=1,3)
 240         format(a4,'    K ',3(F8.3,1x),'    Aposterori (m)')
         end if
      end do
      if( num_tdep.gt.0 ) then
         write(unit, 250)
 250     format('TIME Dependent process noise',/,
     .          'Site     XYZ Process noise (m/sqrt(t))   Start ',
     .          ' -> End ')
         do i = 1, num_site

*           See if Kine or Fixed
            if( site_type(i).ne.0 ) then
                do j = 1, num_tdep 
                    if( nst_tdep(j).eq.i .or. nst_tdep(j).eq.-1 ) then
                        call mjd_to_ymdhms(mjd_tdep(1,j),dates,secs)
                        call mjd_to_ymdhms(mjd_tdep(2,j),datee,sece)

                        write(unit,260) site_names(i), 
     .                      (sqrt(mar_tdep(k,j)*tu_to_ep),k=1,3), 
     .                       dates, secs, datee, sece
 260                    format(a4, 3(F8.3,1x),1x,I4,4i5,1x,F5.1,' -> ',
     .                         I4,4i5,1x,F5.1)   
                    end if
                end do
             end if
         end do
      endif 


****  See if there are time dependent values

      write(unit,270) 
 270  format('SITE   ARP N    ARP E    APR U (m)',
     .       ' dL1 N    dL1 E    dL1 U    ',
     .       ' dL2 N    dL2 E    dL2 U    ',
     .       ' Antenna               PCV  #El #Az RCV')
      do i = 1, num_site
         write(unit,275) site_names(i),(site_offarp(j,i),j=1,3),
     .        ((sit_L12(j,k,i),j=1,3),k=1,2),
     .        ant_name(i), (num_sit_dph(j,i),j=1,2), rcv_type(i) 
 275     format(a4,3(1x,F8.4),1x,3(1x,F8.4),1x,3(1x,F8.4), 
     .          4x,a20,5x,2(I3,1x),3x,a1)
      end do

      do i = 1, num_site
         if( use_ionlos(i) ) then
            write(*,280) site_names(i), trim(ionlos_file(i))
 280        format('IONLOS File at ',a,' File: ',a)
         endif
      end do


****  See if stopgo status set
      if( stopgo_mode ) then
         write(unit,290) stopgo_dvar
 290     format('STOPGO MODE: Variance reduced by ',e10.2)
      endif


*     Report DCB values  
      call report_dcbstr(unit)

      write(unit,320) max_ion_jmp, ion_ppm*1.d6, ion_wght, 
     .                ion_ht/1.d3, ion_tau/1.d3,  ion_corr
 320  format('IONOSPHERE CHARACTERIZATION:',/,
     .       'ION JUMP   ',f7.2,' cycles; ION PPM ',f6.2,' ppm; ',
     .       'ION Weight ',f6.2,/,
     .       'ION Height ',f7.1,' km, ION Correlation length ',
     .       F8.1,' km, ION correlated variance ',F9.2,' cycles**2')

      write(unit,340) max_gap, min_good, rms_edtol, mwwl_jmp
 340  format('New Bias introduced for gaps larger than ',i3,' epochs.',
     .       ' Minimum Good data ',i5,' epochs',/,
     .       'Data editing if RMS greater than ',F5.2,' data sigmas',/,
     .       'MW-WL Jump limit ',F7.2,' cycles')
      if( set_msec_bf ) write(*,'(a)') 
     .       'Bias flags added at millisecond jumps'

****  Tell user about SP3 file:
      write(unit,400) sp3_file(1:trimlen(sp3_file)), num_sat, num_sp3
 400  format('SP3 File ',a,' contains ',i3,' satellites at ',i5,
     .       ' epochs')

      if( num_sat.eq.0 ) then
         write(unit,420) 
 420     format('**DISASTER** No satellites in SP3 file, check file')
         stop 'TRACK: No satellites in SP3 file'
      endif
 

      if( use_atm_ref ) then 
          call report_atm_ref( unit )
      end if



****  Write out user edits
      if( num_edits.gt.0 ) write(unit,510) 
 510  format('DATA DELETED BY USER',/,
     .       'Site  PRN   Start (YMDHMS)             End (YMDHMS)')
      do i = 1, num_edits
         call mjd_to_ymdhms(tt_edit(1,i), dates, secs)
         call mjd_to_ymdhms(tt_edit(2,i), datee, sece)
         write(unit,530) site_names(ss_edit(1,i)), ss_edit(2,i), 
     .                   dates, secs, datee, sece
 530     format(a4,2x,i3.2,1x,2(i5,4i3,1x,f6.3,2x))
      end do

*     Write out excludes
      if( num_exclude.gt.0 ) then
         write(unit,550) (ss_exclude(i),i=1,num_exclude)
 550     format('EXCLUDED PRNS: ',32(i3.2,1x))
      endif

      if( num_abf.gt.0 ) write(unit,610) 
 610  format('BIAS FLAGS ADDED BY  USER',/,
     .       'Site  PRN   Epoch  Time  (YMDHMS) ')
      do i = 1, num_abf
         call mjd_to_ymdhms(tt_abf(i), dates, secs)
         ep = nint((tt_abf(i)-ref_start)*86400.d0/usr_interval)+1
         write(unit,630) site_names(abs(ss_abf(1,i))), ss_abf(2,i), 
     .                   ep, dates, secs
 630     format(a4,3x,i3.2,1x,I7, 1x, (i5,4i3,1x,f6.3,2x))
      end do

      if( num_rbf.gt.0 ) write(unit,650) 
 650  format('BIAS FLAGS DELETED BY  USER',/,
     .       'Site  PRN   Epoch  Time  (YMDHMS) ')
      do i = 1, num_rbf
         call mjd_to_ymdhms(tt_rbf(i), dates, secs)
         ep = nint((tt_rbf(i)-ref_start)*86400.d0/usr_interval)+1
         write(unit,660) site_names(ss_rbf(1,i)), ss_rbf(2,i), 
     .                   ep, dates, secs
 660     format(a4,3x,i3.2,1x,I7, 1x, (i5,4i3,1x,f6.3,2x))
      end do
      write(unit,'(1x)')
      return 
      end

CTITLE REPORT_DCBSTR

      subroutine report_dcbstr(un)

      implicit none 

*     Routine to report the DCB values in processes

      include 'track_com.h'

* PASSED 
      integer*4 un   ! Unit number

* LOCAL
      integer*4 i, trimlen
      integer*4 pn   ! We may not know number of satellites yet

      pn = num_prn
      if( pn.eq.0 ) pn = 32  ! Reset for we don't know (in trackRT
                     ! we have not read any data yet when routine called


      if( trimlen(dcb_file).eq. 0 ) RETURN

      if( num_prn.eq. 0 ) then   ! Print out at start
          write(un,120) pn, dcb_mjd, dcb_file(1:trimlen(dcb_file)),
     .         (i,dcbs(i),i=1,pn)
      else     ! Use actual list
          write(un,120) pn, dcb_mjd, dcb_file(1:trimlen(dcb_file)),
     .         (prn_used(i),dcbs(prn_used(i)),i=1,pn)
      endif


 120  format('DCBS (m) by satellite number for ',i3,' Sats. ',
     .     ' MJD ',F9.2,' File ',a,/,
     .     10(:,10(:,' ',i3.2,1x,F6.3,2x),/))

      return
      end 

