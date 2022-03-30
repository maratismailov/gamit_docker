CTITLE INIT_TRACKRT

      subroutine init_trackRT

      implicit none 

*     Routine to initialise the variable and set the defaults
*     for the trackRT and trackRTr programs

      include 'trackRT.h'
      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
 
* LOCAL VARIABLES
* i, j  -- loop counters

      integer*4 i,j

****  OK, Initialize the variables which need it
      ref_start = 0
      usr_interval = 0
      usr_nepochs  = 0


      num_prn = 0
      num_epochs = 0 
      num_ambs   = 0
      num_antmods = 0

*     Now start setting the defaults     
      max_gap    = 10  ! Default to 10-sec gap before new
                       ! ambiguity
      min_good   = 20
      ref_rel_hum = 0.d0

      do i = 1, max_site
         apr_atm(i) = 0.d0
         mar_atm(i) = 0.d0
         mar_atm_hgt(i) = 0.d0
         do j = 1, 3
            apr_site(j,i) = 1.d0
            mar_site(j,i) = 1.d0
            site_offarp(j,i) = 0.d0
            sit_L12(j,1,i) = 0.d0
            sit_L12(j,2,i) = 0.d0

         end do
         do j = 1,2
            num_sit_dph(j,i) = 0
         end do

         ant_name(i) = 'NONE'
         rcv_type(i) = 'N'
         ref_sv(i) = 0  ! Reference satellites for Double differnces
         reset(i) = .false.
         write_pos(i) = .true.  ! Default to writing positions out
      end do

      do i = 1,3
         ref_xyz(i) = 0.d0
      end do

      do i = 1, max_kine
         static(i) = .false.
      end do
      
      do i = 1, max_sat
         data_var(1,i) = (0.003d0*fL1/vel_light)**2
         data_var(2,i) = (0.003d0*fL2/vel_light)**2
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


      relrank_limit = 20.0d0

      rms_edtol = 5.d0     ! 5 times the data noise for the observable
      min_lvar  = 0.1d0**2 ! Min-phase sigma (not residual < 0.5 cyc
                           ! delted)
      num_edtol = 4        ! number bad residuals to see if cycle slip
                           ! should be added
      min_ldd   = 4        ! Minimum number of double differences needed
                           ! (number phase values is num_ldd+1)

      mwwl_jmp = 5.0d0     ! Jump in MW-WL that will add bias flag.
      exwl_jmp = 0.20d0    ! Default jump on EX-WL to have new estimate
      dd_jmp   = 0.50d0    ! Default jump on baias fixed DD
      exwl_minsig = 0.02d0 ! Min EX_WL mean sigma
      exwl_scale  = 1.0d-6 ! Baseline dependent additive value
      exwl_elev   = 0.4d0  ! Elevation angle rewieghting
      mwwl_minsig = 0.10d0 ! Min MW-WL mean sigma
 
      do j = 1,10
         debug(j) = 0
      end do
      debug(1) = -1
      debug(2) = -2

      num_edits = 0
      num_exclude = 0
      num_abf = 0
      num_rbf = 0

      posit_root  = 'trackRT'
      sum_file    = 'trackRT'
      antmod_file = ' '
      dcb_file    = ' '
      upd_file    = ' '
      lus = 98
      resid_root = ' '
      wls_root = ' '
      rwl_root = ' '
      runday = ' '
      runweek = ' '
      posit_type = 'DHU'
      file_updint = 1/24.d0  ! Hourly file updates
      curr_timetag = ' '
      sp3_file = ' '
      sp3_dir = '.'    ! Default current directory

      write_res = .false. 
      atm_mtt = .false.     ! Set default to use GPT/GMF
      L1_only = .false.

      use_atm_ref = .false.
      atm_file = ' '
c      noetide = .true.  ! Debug
      noetide = .false.  ! Debug
      set_msec_bf = .false.  ! Added version 1.22

      ambin_file = ' '
      back_type = 'NONE'

      float_type = 'LCPC'
      neam = 1   ! Only 1 parameter for LC (and L1)
      float_iter = 1
      float_sample = 1
      float_limit(1) = 0.25
      float_limit(2) = 0.50

      WL_avnum = 25   ! Number of values of computing
                       ! sigma of WL mean 
      WL_mnnum = 10    ! Minimum of MW-WL before allowing fix

      out_sig_limit  = 1.0
      elev_cutoff = 10

      wl_fact  = 1.0   ! Weight for MW-WL
      lg_fact  = 1.0   ! Weight for EX-WL (should be decreased for 
                       ! longbaselines.  Make length dependent?)
      min_ambsig = 0.05  ! Minium sigma for ambiguity estimates
      min_exsig  = 0.01  ! Minimum sigma for EX-WL 
      max_fit  = 25.0

      num_anal_type = 1
      anal_types(1) = 'LC'
      posit_type = 'DHU'

      stopgo_mode = .false.
      ante_off = .false.   ! Set true when antenna informtion
                ! read.  FATAL if antmod_file command used before
                ! ANTE_OFF 

      time_unit = 'second'
      tu_to_ep  = 1.0

      status_int = 0
      status_type = 'PAWR'

      updread = .false.
      needupd = .false.

      usr_start = 0.d0   ! Start time (all data)

      

****  Thats all
      return 
      end

CTITLE REPORT_SETUPRT

      subroutine report_setupRT(unit)

      implicit none 

*     Routine to report the setup for the track run
      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRT.h'

* PASSED VARIABLES
* unit  -- Unit number for output

      integer*4 unit

* LOCAL VARIABLES
      integer*4 i, j, k, ep, trimlen, dates(5), datee(5)
      integer*4 jel   ! Function to return baseline number
      real*8 secs, sece, totvar, bl
      real*8 datasig(4)  ! Data sigmas in meters (phase is cyc)

* stype(2) -- Type of positioning (S, K)
      character*10 stype(2)

      character*1 sqroot   ! Extended ascii symbol for square root

      data stype / 'STATIC', 'KINEMATIC' /


****  Tell user what is happening
      write(unit,20)
 20   format(/,'TRACKRT SETUP PARAMETERS',/,
     .         '------------------------')
      write(unit,40) trim(rt_machine), rt_port
 40   format('RT Comm ',t24,a,' Port ',I6)
      write(unit,60) 'Batch File ',bat_file(1:trimlen(bat_file))
      write(unit,60) 'SP3 Dir ',sp3_dir(1:trimlen(sp3_dir)), sp3_root
 60   format(a,t24,': ',a,1x,a)

      do i = 1, num_antmods
        write(unit,60) 'ANTEX File  ',
     .                         antmod_file(i)(1:trimlen(antmod_file(i)))
      end do

C     do i = 1, num_site
C        write(unit,70) site_names(i), stype(site_type(i)+1)
C70      format('Rinex files',t24,': ',a,2x,'Type ',a)
C     end do

 130  format('FLOAT TYPE             : ',a,' Start Iter ',i2,
     .       ' Sample ',i3,' Limits ',2F5.2,' WL/LG SCALES ',2F5.2,
     .       ' MAX Fit ',F6.1)
      write(unit,140) (anal_types(i)(1:trimlen(anal_types(i)))
     .                ,i=1,num_anal_type)
 140  format('DATA ANALYSIS TYPES    : ',a)
      datasig(1) = sqrt(data_var(1,1))*vel_light/fL1
      datasig(2) = sqrt(data_var(2,1))*vel_light/fL2
      datasig(3) = sqrt(data_var(3,1))
      datasig(4) = sqrt(data_var(4,1))

      write(unit,160) (datasig(i),i=1,4),data_var(5,1)
 160  format('Data Noise L1,L2,P1,P2 GENERIC: ',2F8.4,2F8.2,
     .       ' m, Weight ',F5.2)
      do i = 1, max_sat
         if( data_var(1,i).ne.data_var(1,1) ) then
            datasig(1) = sqrt(data_var(1,i))*vel_light/fL1
            datasig(2) = sqrt(data_var(2,i))*vel_light/fL2
            datasig(3) = sqrt(data_var(3,i))
            datasig(4) = sqrt(data_var(4,i))
            write(unit,180) i,(datasig(j),j=1,4), data_var(5,i)
 180        format('Data Noise L1,L2,P1,P2 PRN ',i2.2,' : ',
     .             2F8.4,2F8.2,' m, Weight ',F5.2) 
         end if
      end do

      write(unit,190) out_sig_limit
 190  format('Postion output sigma limit    : ',f10.3,' m.')

      write(unit,200) max_gap, min_good
 200  format('BF Set: Max Gap ',i4,' Min Good ',i4,' epochs')


      write(unit,210) elev_cutoff
 210  format('Elevation angle cutoff        : ',F10.3,' deg.')
      if( atm_mtt ) then
         write(unit,'(a,1x,a)') 'MTT Mapping function and',
     .                       'seasonal model used'
      else
         write(unit,215) 'GPT Pressure/Temp and',
     .                     'GMF mapping function used.', 
     .                     ref_rel_hum
 215     format(a,1x,a,' with ',F5.2,' relative humidity')
      endif

****  Output the site characteristics
      sqroot = 'R'  !   char(164)
      write(unit,220) time_unit(1:trimlen(time_unit)), 
     .             sqroot,time_unit(1:1), sqroot, time_unit(1:1)
 220  format('Process noise time unit: per sqrt(',a,')',/,
     .       'SITE STOCHASTIC PROPERTIES',/,
     .       'SITE Type         Apriori XYZ (m)        ',
     .       'XYZ Process noise (m/',a,a,') Atm Stats (m,m/',a,a,
     .       ') Hgt m/(m/s)')

      do i = 1, num_site

*        See if Kine or Fixed
         totvar = 0
         do j = 1,3
            totvar = totvar + apr_site(j,i)+mar_site(j,i)
         end do
         if( totvar.eq.0 ) site_type(i) = 0  ! No noise so mark fixed.

         if( site_type(i).eq.0 ) then
*            Fixed Site
             write(unit,240) site_names(i), sqrt(apr_atm(i)), 
     .            sqrt(mar_atm(i)*tu_to_ep), sqrt(mar_atm_hgt(i))
 240         format(a4,'    F ', 56x,F7.3,1x,F8.5,1x,F8.5)
         else
*            Kinematic site
             write(unit,245) site_names(i), (sqrt(apr_site(j,i)),j=1,3),
     .                   (sqrt(mar_site(j,i)*tu_to_ep),j=1,3), 
     .                    sqrt(apr_atm(i)), sqrt(mar_atm(i)*tu_to_ep), 
     .                    sqrt(mar_atm_hgt(i))
 245         format(a4,'    K ',3(F8.3,1x),1x,3(F8.3,1x),1x,
     .                         F7.3,1x,F8.5,1x,F8.5)
         end if
      end do
      write(unit,250) 
 250  format('SITE   ARP N    ARP E    APR U (m)',
     .       ' dL1 N    dL1 E    dL1 U    ',
     .       ' dL2 N    dL2 E    dL2 U    ',
     .       ' Antenna               PCV  #El #Az RCV')
      do i = 1, num_site
         write(unit,255) site_names(i),(site_offarp(j,i),j=1,3),
     .        ((sit_L12(j,k,i),j=1,3),k=1,2),
     .        ant_name(i), (num_sit_dph(j,i),j=1,2), rcv_type(i)
 255     format(a4,3(1x,F8.4),1x,3(1x,F8.4),1x,3(1x,F8.4), 
     .          4x,a20,5x,2(I3,1x),3x,a1)
      end do
      write(unit,'(a)') 'Distances to Reference site'
      do i = 2,num_site
         bl = baselens(jel(i,1))
         write(unit,260) site_names(i),site_names(1), bl
 260     format('Distance ',a4,'-',a4,1x,F12.3,' m')
      end do


 
      call report_dcbs(unit)
 
      write(unit,320) exwl_jmp,  exwl_minsig, exwl_scale*1.d6,
     .        exwl_elev
 320  format('EX-WL Jump ',F7.2,' cyc, Min Sigma ',F7.3,' cyc, ',
     .       'Baseline length scaling ',f6.2,' ppm, Elev weight ',F6.3)

      write(unit,340) mwwl_jmp, mwwl_minsig, wl_avnum, wl_mnnum
 340  format('MW-WL Jump ',F7.2,' cyc, Min Sigma ',F7.3,' cyc, ',
     .       'MW_WL Max Num ',I4,' MW-WL Min Num ',i4)

      write(unit,360) rms_edtol, sqrt(min_lvar), num_edtol, 
     .        dd_jmp, min_ldd
 360  format('Data editing if RMS greater than ',F5.2,' data sigmas.',
     .       ' Min Phase sigma ',F6.2,' cyc. Reset after ',
     .        i3,' bad data',/,
     .       'Double Difference Jump Tol ',F8.2,' cyc, Min Phase DD ',
     .        I3)

      write(unit,380) relrank_limit, float_limit, wl_fact, lg_fact,
     .     min_ambsig, max_fit
 380  format('RelRank Limit ',F8.2,' Float limits ',2F7.2,' cyc, '
     .       'MW Factor ',F6.2,' EX Factor ',F6.4,' Min AmbSig ',
     .       F6.3,' cyc, Max Chi ',F7.2)
      if( set_msec_bf ) write(*,'(a)') 
     .       'Bias flags added at millisecond jumps'

****  Tell user about SP3 file:
      write(unit,400) sp3_file(1:trimlen(sp3_file)), num_sat, num_sp3
 400  format('SP3 File ',a,' contains ',i3,' satellites at ',i5,
     .       ' epochs')


****  Write out user edits
      if( num_edits.gt.0 ) write(unit,510) 
 510  format('DATA DELETED BY USER',/,
     .       'Site  PRN   Start (YMDHMS)             End (YMDHMS)')
      do i = 1, num_edits
         call mjd_to_ymdhms(tt_edit(1,i), dates, secs)
         call mjd_to_ymdhms(tt_edit(2,i), datee, sece)
         write(unit,530) site_names(ss_edit(1,i)), ss_edit(2,i), 
     .                   dates, secs, datee, sece
 530     format(a4,3x,i2.2,1x,2(i5,4i3,1x,f6.3,2x))
      end do

*     Write out excludes
      if( num_exclude.gt.0 ) then
         write(unit,550) (ss_exclude(i),i=1,num_exclude)
 550     format('EXCLUDED PRNS: ',32(i2.2,1x))
      endif

      write(unit,'(1x)')
      return 
      end

CTITLE REPORT_DCBS

      subroutine report_dcbs(un)

      implicit none 

*     Routine to report the DCB values in processes

      include 'trackRT.h'

* PASSED 
      integer*4 un   ! Unit number

* LOCAL
      integer*4 i, trimlen
      integer*4 pn   ! We may not know number of satellites yet

      pn = num_prn
      if( pn.eq.0 ) pn = 32  ! Reset for we don't know (in trackRT
                     ! we have not read any data yet when routine called


      if( trimlen(dcb_file).eq. 0 ) RETURN

      write(un,120) pn, dcb_mjd, dcb_file(1:trimlen(dcb_file)),
     .     (i,dcbs(i),i=1,pn)
 120  format('DCBS (m) by satellite number for ',i3,' Sats. ',
     .     ' MJD ',F9.2,' File ',a,/,
     .     10(:,10('G',i2.2,1x,F6.3,2x),/))

      return
      end 

CTITLE REP_STATUSRT

      subroutine report_statusrt( un, ep ) 

      implicit none 

*     Routine to report the STATUS at epoch EP for the trackRT run

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRT.h'

* PASSED VARIABLES
* unit  -- Unit number for output

      integer*4 un     ! Output unit number
     .,         ep     ! Epoch number 

* LOCAL
      integer*4 date(5)  ! YMD HM
      integer*4 i, j 
      real*8 sec      ! seconds tag

      integer*4 nu      ! Number of used Widelane types
      real*4 uout(3,5)  ! Up to 5 output types

      logical kbit      ! Tests bit
      logical OK        ! Set true to output 

      character*2 utypes(5), chi2_types(5)
      character*1 fxdstat  ! Set as X or R if fixd or free

      data chi2_types / 'MW','EX','LC','L1','L2' /

****  OK: See what status entries we need
      if( status_int.eq.0 ) RETURN
      if( ep-int(ep/status_int)*status_int.ne.0 ) RETURN

***   OK; see what we need
      call mjd_to_ymdhms(kf_curr_mjd, date, sec)
      write(un,110) ep, date, sec, status_type
 110  format(100('-'),/,
     .     'STATUS REPORT Epoch ',i6,1x,' Date ',I4,4(1x,i2.2),1x,
     .      F7.3,' Type ',a)
      
      if( index(status_type,'P').gt.0 ) then
          call print_parRT(un, ep, 'P','STATUS REPORT')
      end if

****  Ambigity estimates
      if( index(status_type,'A').gt.0 ) then
         write(un,210) ep, num_ambs
 210     format('AMBIGUITY Report Ep ',i6,' Number of ambiquities ',i4)
         do i = 1,num_ambs 
            call stat_ambinf( un, i, ep )
         end do

      end if

****  Wide lane report (also includes LC residual) 
      if( index(status_type,'W').gt.0 ) then
***      See what we used and output values, sigmas and chi cont
          write(un,310) ep, num_ambs
 310      format('WIDELANE Report Ep ',i6,' Number of ambiquities ',i4)
          do i = 1, num_ambs
             OK = .true.
             if( WLS_num(i).le.0 ) OK = .false.
             if( index(status_type,'C').gt.0 .and.
     .           bf_ents(4,i).lt.ep-10 ) OK = .false.
             if( OK ) then
                nu = 0
                do j = 1,5
                   if( asv_used(j,i) ) then
                      nu = nu + 1
                      utypes(nu) = chi2_types(j)
                      uout(1,nu) = asv_res(j,i)
                      uout(2,nu) = asv_sig(j,i)
                      uout(3,nu) = asv_chi(j,i)
                   endif
                end do
*               See if fixed or free
                if( kbit(bf_ents(5,i),2) ) then
                   fxdstat = 'X'
                else
                   fxdstat = 'R'
                end if

                write(un,320) i, site_names(bf_ents(1,i)),bf_ents(2,i),
     .              ep, bf_ents(3,i), bf_ents(4,i), WLS_num(i), nu, 
     .              fxdstat, (utypes(j),uout(:,j),j=1,nu)
 320            format('WIDELANE ',i3,1x,a4,1x,'PRN',i2.2,' EP ',I6,
     .             ' Range ',2I7,' # ',I6,' NC ',i1,1x,a1,1x, 
     .              5(:,a2,' Res ',F6.2,1x,F6.2,' Chi2 ',F6.2,1x))
             end if
          end do 

      end if

      if( index(status_type,'R').gt.0 ) then
          write(un,410) ep, num_dd
 410      format('POSTFIT RESIDUAL Report Ep ',i6,' Number of DD ',i4)
          do i = 1, num_dd
             write(un,430) i, ep, site_names(psv_ss(1,i)), psv_ss(2,i),
     .          site_names(psv_ss(3,i)), psv_ss(4,i),psv_res(i),
     .          psv_sig(i), psv_dtype(i), psv_amb(i), psv_ae(2,i)
 430         format('POST ',i3,' Ep ',i6,1x,a4,' PRN',I2.2,' - ',
     .          a4,' PRN',I2.2,1x,' Res ',F10.4,' +- ',F7.4,1x,a2,
     .          ' AMB ',i3,' Elev ',F6.2,' deg')
          end do
      end if


****  Thats all
      return
      end 


CTITLE STAT_AMBINF

      subroutine stat_ambinf( un, na, ep )

      implicit none

*     Routine to status report information about ambiquities.

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRTObs.h'      ! Data definition
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 un   ! Output unit number
     .,  na   ! Ambiguity being reported
     .,  ep   ! Epoch number when reported 

      logical kbit, OK


      character*4 restyp    ! String with resolve type AMBFIX AMBFRE

****  Output information about this bias
      if( kbit(bf_ents(5,na),2) ) then
          restyp = 'FIXD'
      else
          restyp = 'FREE'
      end if
      OK = .true.
      if( WLS_num(na).le.0 ) OK = .false.
      if( index(status_type,'C').gt.0 .and.
     .    bf_ents(4,na).lt.ep-10 ) OK = .false.
*     See if LC was estimated (when the reference satellite is
*     changed, the old ref has widelane data but no LC estimate
*     and so we do not report here)
      if( .not. asv_used(3,na) ) OK = .false.
      if( OK ) 
     .write(un,125) restyp, na, site_names(bf_ents(1,na)),
     .    bf_ents(2,na), asv_rbn(1,na), asv_Fcode(na), asv_dL12(:,na), 
     .    asv_rbn(2,na), asv_rbn(3,na), psv_ae(2,na), asv_resep(na)
 125  format('AMB',a,1x,' # ',i3,1x, a4,1x,'PRN',i2.2' RelRank ',F8.2,
     .    1x,' FCode ',a,' dL12 ',2i4,' Dchi ',2F9.2,
     .       ' Elev ',F6.2,' deg; ResEpoch ',i6)

      end

