ctitle
 
      subroutine get_missrt(buffer,itype,indx)

      implicit none
c
c     routine to read the rest of the buffer line to get the
c     vals of the miscellaneous parameters
c
      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRT_cmds.h'
*                                        ! the common block
      include 'trackRT.h'
      include 'trackRTObs.h'

c Variables
c ---------
c buffer  -- 80 character buffer which will be decoded as indicated
c        by the itype value.
c itype   -- integer value indicating the type of parameter to be
c        read.  (This type is based on the control_name array as
c        listed in the block data kal_block_data.)
c
c Local variables
c ---------------
c label   -- integer*4 variable used in the computed goto statement
c vals  -- real*4 array to temporarly save the vals read from
c         the buffer string
c ians    -- buffer for questions which have a yes/no response
c
c

      character*(*) buffer
 
      integer*4 itype, label, trimlen 
c
      real*4 vals(20)
      real*8 vals8(20), val8  
c
*   err         - Error message returned from multiread and
*               - READ_LINE.
*   indx       - Used by multiread and READ_LINE to keep track
*               - of position in buffer as it is read.
*   date(5)     - Calender data input by user
*   ivals(5)    - Return of integer values from multiread
*   pn          - PRN number

      integer*4 err, indx, jndx, i, j, date(5), ivals(5), pn

*   sectag      - Seconds tag for date

      real*8 sectag
      real*8 dtun  ! Units for output file spacing.
      character*126 dtunstr  ! String with time step and unit
 
*   cvalue      - Dummy character string used in multiread.
*   name        - Name of site
*   mode        - type of default settings
 
      character*256 cvalue(5)
      character*4 name
      character*16 mode, cval
 
c.... compute the number of the miscellaneous parameter
      label = itype - mis_start + 1

*                    ! Ensures strings are decoded from the second
c      indx = 1
*                     ! character.
c
 
c.... go to the appropriate statements
****  See if we should subs for <day> or <week>

      goto ( 100, 200, 300, 400, 500, 600, 700, 800, 900,1000,
     .      1100 , 1200, 1300, 1400, 1500, 1600, 1700,
     .      1800 , 1900, 2000, 2010, 2020, 2030, 2040, 2050,
     .      2100 , 2200, 2300, 2400, 2500, 2600, 2700, 2800, 
     .      2900 , 3000, 3100, 3200, 3300, 3400, 3500, 3600, 
     .      3700 , 3800, 3900, 4000, 4100, 4200, 4300, 4400,
     .      9900 ) label
c
c.... file name of the apriori site and station coordinates
  100 continue
*        Option not supported yet.
C        call multiread(buffer,indx,'CH',err,vals,apr_file,1)
         write(*,'(a)') '**WARNING** APR_FILE command not implemented'
         return        

c.... SP3_DIR: Directory where SP3 files will be found
  200 continue
         j = -2
         call multiread(buffer,indx,'CH',err,vals,cvalue,j)
         sp3_dir = cvalue(1)
         if( j.eq.2 ) then
            sp3_root = cvalue(2)
         else
            sp3_root = 'igs'
         end if
         call caseunfold(sp3_root)

         return        
       
c
c.... GPS SV Clock file name of the GPS satellite broadcast
  300 continue
*        Not supported in Version 1.00
C        call multiread(buffer,indx,'CH',err,vals,sv_clk_file,1)
         return        

c.... Residual output option.  Get the name of the root component.
  400 continue
         call read_line(buffer,indx,'CH',err,vals,resid_root)

         return        
         
c.... BF_SET   : Bias flags setting.  Sets maximum gap size and minimum
*     duration of good data. 
  500 continue

         call read_line(buffer,indx,'I4',err,ivals,cvalue)
         if( ivals(1).gt.0 ) max_gap = ivals(1)
         call read_line(buffer,indx,'I4',err,ivals,cvalue)
         if( ivals(1).gt.0 ) min_good = ivals(1)
         return        
c

c.... DEBUG  : output epoch range: Values should be paired
*     and turn on different types of debug (see trackRT_comm.h)
  600 continue 
*        Get the start and stop epochs
         call multiread(buffer,indx,'I4',err,ivals,cvalue,2)
         debug(1) = ivals(1)
         debug(2) = ivals(2)
*        Now see if additional values have been passed
         j = 2
         do while ( err.eq.0 .and. j.lt.10 )
            call read_line(buffer,indx,'I4',err,ivals,cvalue)
            if( err.eq.0 ) then
               j = j + 1
               debug(j) = ivals(1)
            end if
         end do

         return        
c

c     EXWL_SET: Set parameters for using the Extra-widelane to resolve
*     ambiquities. 
  700 continue
*        Read each of the values to see which ones are passed
         call read_line(buffer,indx,'R8',err,val8,cvalue)
         if( err.eq.0 .and. val8.gt.0 ) exwl_jmp = val8
         call read_line(buffer,indx,'R8',err,val8,cvalue)
         if( err.eq.0 .and. val8.gt.0 ) exwl_minsig = val8
         call read_line(buffer,indx,'R8',err,val8,cvalue)
         if( err.eq.0 .and. val8.gt.0 ) exwl_scale = val8*1.d-6
         call read_line(buffer,indx,'R8',err,val8,cvalue)
         if( err.eq.0 .and. val8.gt.0 ) exwl_elev = val8

         return        
c

c.... START_TIME <yy mm dd hh mn sec.>
  800 continue
         call multiread(buffer,indx,'I4',err,date,cvalue,5)
         call multiread(buffer,indx,'R8',err,sectag,cvalue,1)
         call ymdhms_to_mjd(date, sectag, usr_start)
         return
c
c.... INTERVAL : interval sampling (seconds)
  900 continue
         call read_line(buffer,indx,'R8',err,usr_interval,
     .                  cvalue)
         return
c
c.... NUM_EPOCHS : numbers of epochs
 1000 continue

         call read_line(buffer,indx,'I4',err,usr_nepochs,cvalue)
         return

c.... OUT_TYPE  : Type of position to output (GEOD/NEU)
 1100 continue
         call read_line(buffer,indx,'CH',err,vals,posit_type)
         call casefold(posit_type)
         return

c.... DATA_NOISE : data noise level: Specify for generic L1, L2, P1, 
c     P2 measurements
 1200 continue
         call multiread(buffer,indx,'R4',err,vals,cvalue,4) 
         call read_line(buffer,indx,'R4',err,vals(5),cvalue)
         if( err.ne.0 ) vals(5) = 0.0
*        See if we have prn
         call read_line(buffer,indx,'I4',err,pn, cvalue)
*        Convert phase data noise to cycles from mm
         vals(1) = vals(1)*fL1/vel_light
         vals(2) = vals(2)*fL2/vel_light
         if( err.eq.0 ) then
            do i = 1, 4
               data_var(i,pn) = vals(i)**2
            end do
            data_var(5,pn) = vals(5)
         else
            do pn = 1, max_sat
               do i = 1, 4
                  data_var(i,pn) = vals(i)**2
               end do
               data_var(5,pn) = vals(5)
            end do
         endif
         return
         
c.... DATA_TYPES : data type be used in the position analysis. 
 1300 continue 
         num_anal_type = 0
         err = 0
         do while ( err.eq.0 )
             call read_line(buffer,indx,'CH',err,vals,
     .                      anal_types(num_anal_type+1)) 
             call casefold(anal_types(num_anal_type+1))
             if( err.eq.0 ) num_anal_type = num_anal_type + 1
         end do

         return

c.... EDIT_SSV  : Edit command (Site Satellite Start and Stop dates with seconds)
 1400 continue

         call read_line (buffer,indx,'CH',err,vals, name)
*        Match the site name 
         jndx = -1
         do i = 1, num_site
            if( site_names(i).eq.name ) jndx = i
         end do
         i = jndx
*        See if we match site
         if( i.gt.0 ) then
             num_edits = num_edits + 1
             if( num_edits.gt.max_edits ) then
                 write(*,*) '**DISASTER** Too many edits specified ',
     .                      'Max allowed is ', max_edits
                 stop 'TRACK: Too many edits specified'
             endif
             call read_line(buffer,indx,'I4',err,ivals,cvalue)
             ss_edit(1,num_edits) = i
             ss_edit(2,num_edits) = ivals(1)

*            Now get the times Start: Y M D H M
             call multiread(buffer,indx,'I4',err,date,cvalue,5) 
             call read_line(buffer,indx,'R8',err,sectag,cvalue)
             call ymdhms_to_mjd(date, sectag, tt_edit(1,num_edits))
*            Now get the times End time : Y M D H M
             call multiread(buffer,indx,'I4',err,date,cvalue,5) 
             call read_line(buffer,indx,'R8',err,sectag,cvalue)
             call ymdhms_to_mjd(date, sectag, tt_edit(2,num_edits))
          end if
        
          return
      
c.... AMBIN_FILE : ambiguity data file name 
 1500 continue
c         call read_line(buffer,indx,'CH',err,vals,ambin_file)
          write(*,'(a)') '**WARNING** AMBIN_FIL command not implemented'

         return        
c
c.... changable position markov noise
 1600 continue
*        Not implemented in version 1.0
C        call multiread(buffer,indx,'CH',err,vals,cvalue,1)
C        if (cvalue(1)(1:1).eq.'T'.or.cvalue(1)(1:1).eq.'t')
C    .        change_wn=.true.      

         return        
c

c....  AMB_CYCLE  : ambiguity constraints.  The values passed are:
*      num_amb_samp -- Number of samples to use for ranking ambiquities.
*      relrank_limit -- Ratio needed to fix the bias
*      max_tot_search -- Limit on size of sample to search

 1700 continue
          write(*,'(a)') '**WARNING** AMB_CYCLE command not implemented'
C         call multiread(buffer,indx,'R8',err,vals8,cvalue,3)
C         if( vals8(1).gt.0 ) num_amb_samp = vals8(1)
C         if( vals8(2).gt.0 ) relrank_limit = vals8(2)
C         if( vals8(3).gt.0 ) max_tot_search = vals8(3)
          return        

c.... GEOD_OUT : Root for the geodetic position output files
 1800 continue
         print *,'Command not implemented: Used POS_ROOT'
c        call read_line(buffer,indx,'CH',err,vals,posit_root)
         return        

c.... CUT_OFF  : cutoff angle (15 degree for default)
c .... add start cutoff angle to avoid ambiguity in low elev. 980430
 1900 continue
         call read_line(buffer,indx,'R8',err,vals8,cvalue)
         elev_cutoff=vals8(1)
         return

c.... SEARCH_TYPE : Ambiquity search type
 2000 continue
         write(*,'(a)') '**WARNING** SEARCH_TYP command not implemented'
C        call read_line(buffer,indx,'CH',err,vals,search_type)
C        call casefold(search_type)

         return        

c.... ATM_FILE  : User provided zenith delays for all sites
 2010 continue
         call read_line(buffer,indx,'CH',err,vals, atm_file)
         return        

c.... POS_ROOT : Name of the root file for the NEU outout 
 2020 continue
         call read_line(buffer,indx,'CH',err,vals, cvalue)
         if( trimlen(prt_root).gt.0 ) then
             call sub_char(cvalue(1),'?',trim(prt_root))
         end if
         posit_root = cvalue(1)
*        See if the update interval has been passed. Interval in
*        days. Allow option for d h m option
         call getword(buffer,dtunstr,indx)
         j = trimlen(dtunstr)
         if( j.gt.0 ) then
             call casefold(dtunstr)
             if( dtunstr(j:j).eq.'D' ) then
                 dtunstr(j:j) = ' '
                 dtun = 1.0d0
             elseif ( dtunstr(j:j).eq.'H' ) then
                 dtunstr(j:j) = ' '
                 dtun = 1.0d0/24.d0
             elseif ( dtunstr(j:j).eq.'M' ) then
                 dtunstr(j:j) = ' '
                 dtun = 1.0d0/1440.d0
             endif
             jndx = 1
             call read_line(dtunstr,jndx,'R8',err,vals8,cvalue)
             if( err.ne.0 .and. err.ne. -1 )
     .          call report_error('IOSTAT',err,'read',buffer,0,
     .              'Output interval')
             if( err.eq.0 ) then
                file_updint = vals8(1)*dtun
             endif

         end if
     
         return        

c.... atmosphric model  used 
 2030 continue
C        Not implemented in version 1.000
C        call multiread(buffer,indx,'CH',err,vals,cvalue,4)
C        do i=1, 4
C            read(cvalue(i),'(a3)') met_type(i)
C        enddo

*.... FLOAT_TYPE: sets parameters for the float ambiquity resolution 
 2040 continue
         call read_line(buffer, indx,'I4',err,ivals,cvalue)
         if( err.eq.0 .and. ivals(1).gt.0 )  float_iter = ivals(1)
         call read_line(buffer, indx,'I4',err,ivals,cvalue)
         if( err.eq.0 .and. ivals(1).gt.0 )  float_sample = ivals(1)

         call read_line(buffer, indx,'CH',err,ivals,cvalue)
         if( err.eq.0 ) float_type = cvalue(1)
         call casefold(float_type)

****     Based on type set number of ambiguities need
         neam = 1
         if( index(float_type,'L2').gt.0 ) neam = 2

         call read_line(buffer, indx,'R4',err,vals,cvalue)
         if( err.eq.0 .and. vals(1).gt.0 )  float_limit(1) = vals(1)
         call read_line(buffer, indx,'R4',err,vals,cvalue)
         if( err.eq.0 .and. vals(1).gt.0 )  float_limit(2) = vals(1)

*        Get the new parameters for wl_fact, lg_fact, and max_fit
         call read_line(buffer, indx,'R4',err,vals,cvalue)
         if( err.eq.0 .and. vals(1).ge.0 )  wl_fact = vals(1)
         call read_line(buffer, indx,'R4',err,vals,cvalue)
         if( err.eq.0 .and. vals(1).ge.0 )  lg_fact = vals(1)
         call read_line(buffer, indx,'R4',err,vals,cvalue)
         if( err.eq.0 .and. vals(1).gt.0 )  max_fit = vals(1)
         call read_line(buffer, indx,'R4',err,vals,cvalue)
         if( err.eq.0 .and. vals(1).gt.0 )  relrank_limit = vals(1)
  
         return
         
* MODE command allows default setup
 2050 continue
         call GetWord(buffer,mode, indx)
*        See type of mode given and set the defaults accordingly
         call casefold(mode)
         if( mode(1:3).eq.'AIR' ) then
             max_gap = 10
             min_good = 120
             search_type = 'NONE'
             float_type = 'LC'
             anal_types(1) = 'LC'
             do i = 2, max_site
                apr_atm(i) = 0.d0
                mar_atm(i) = 0.d0
                mar_atm_hgt(i) = 0.d0 
             end do
         else if( mode(1:3).eq.'SHO' ) then
             max_gap = 10
             min_good = 20
             search_type = 'NONE'
             float_type = 'L1+L2'
             anal_types(1) = 'L1+L2'
             do i = 2, max_site
                apr_atm(i) = 0.d0
                mar_atm(i) = 0.d0
                mar_atm_hgt(i) = 0.d0 
             end do
        else if ( mode(1:3).eq.'LON' ) then
             max_gap = 10
             min_good = 20
             search_type = 'NONE'
             float_type = 'LC'
             anal_types(1) = 'LC'
             do i = 2, max_site
                apr_atm(i) = 0.1d0
                mar_atm(i) = 1.0d-6
                mar_atm_hgt(i) = (0.23d-3)**2
             end do
             lg_fact = 0.1d0
         end if
         return

*.... BACK Command: sets type of back solution to run
 2100 continue   
         call  GetWord(buffer, back_type, indx)
         call casefold(back_type)
*        Make sure that BACK is on if SMOOTH selected
         if( index(back_type,'SMOOTH').gt.0 .and. 
     .       index(back_type,'BACK').eq.0 ) then
             back_type(trimlen(back_type)+1:) = 'BACK'
         end if
   
         RETURN

*.... OUT_SIG_LIMIT: Limits sigmas of out put value
 2200    continue
         call read_line(buffer,indx,'R4',err,vals,cvalue)
         if( vals(1).gt.0 ) out_sig_limit = vals(1)
         RETURN

*.... RMS_EDIT_TOL <n-sigma> <min-sigma> <Reset number>
*        Sets limit on editing tolerance as a multiplier of
*        the data noise
 2300    continue
         call read_line(buffer,indx,'R4',err,vals,cvalue)
         if( err.eq.0 .and. vals(1).gt.0 ) rms_edtol = vals(1)
         call read_line(buffer,indx,'R4',err,vals,cvalue)
         if( err.eq.0 .and. vals(1).gt.0 ) min_lvar = vals(1)**2
         call read_line(buffer,indx,'I4',err,ivals,cvalue)
         if( err.eq.0 .and. ivals(1).gt.0 ) num_edtol = ivals(1)

         RETURN

*.... EXCLUDE_SVS: Exclude satellites from processing
 2400    continue
         err = 0
         num_exclude = 0
         do while( err.eq.0 )
             num_exclude = num_exclude + 1
             call read_line(buffer,indx,'I4',err, 
     .            ss_exclude(num_exclude), cvalue)
             if( err.ne.0 ) num_exclude = num_exclude - 1
         end do

c.... WLS_ROOT : Name of the root file EX and WL widelanes
c     corrected for ambiquities 
 2500 continue
         call read_line(buffer,indx,'CH',err,vals, wls_root)
         call subhome(wls_root)
     
         return        
           
c.... RWL_ROOT : Name of the root file EX and WL widelanes
c     uncorrected for ambiquities (raw widelines)
 2600 continue
         call read_line(buffer,indx,'CH',err,vals, rwl_root)
         call subhome(rwl_root)
    
         return        
           
c.... SUM_FILE : Name of the summary file for summary of run
 2700 continue
         call read_line(buffer,indx,'CH',err,vals, sum_file)
         call subhome(sum_file)
         if( trimlen(prt_root).gt.0 ) then
             call sub_char(sum_file,'?',trim(prt_root))
         end if

         return

c..... STOPGO_MODE: Sets track so that process noise is set to zero
c      during STATIC mode 
 2800  continue
         jndx = indx
         call read_line(buffer,jndx,'R4',err,vals, cvalue)
         stopgo_mode = .true.
         if( err.eq.0 ) then
            stopgo_dvar = vals(1)
         else
            stopgo_dvar = 1.d6
         endif
         return

c..... MWWL_SET <Jump> <Max Averaging number> <Min number>
 2900  continue
          jndx = indx
          call read_line(buffer,indx,'R8',err,val8,cvalue)
          if( err.eq.0 .and. val8.gt.0 ) mwwl_jmp = val8
          call read_line(buffer,indx,'R8',err,val8,cvalue)
          if( err.eq.0 .and. ivals(1).gt.0 ) mwwl_minsig = val8
          call read_line(buffer,indx,'I4',err,ivals,cvalue)
          if( err.eq.0 .and. ivals(1).gt.0 ) wl_avnum = ivals(1)
          call read_line(buffer,indx,'I4',err,ivals,cvalue)
          if( err.eq.0 .and. ivals(1).gt.0 ) wl_mnnum = ivals(1)

          return

c,.... ANTMOD_FILE
 3000  continue
          jndx = indx
          call read_line(buffer,jndx,'CH',err,vals, antmod_file(1))
          call subhome(antmod_file(1))
*         Now open and read antmod file.  If read is OK, file name 
*         added to antmod_file list
          call read_antmodRT(antmod_file(1))
          return

c..... USR_ADDBF  : Allows bias flags to be added by user
 3100  continue
          
         call read_line (buffer,indx,'CH',err,vals, name)
*        Match the site name 
         jndx = -1
         do i = 1, num_site
            if( site_names(i).eq.name ) jndx = i
         end do
         i = jndx
*        See if we match site
         if( i.gt.0 ) then
             num_abf = num_abf + 1
             if( num_abf.gt.max_edits ) then
                 write(*,*) '**DISASTER** Too many User bias flags ',
     .                      'specified. Max allowed is ', max_edits
                 stop 'TRACK: Too many user BFs specified'
             endif
             call read_line(buffer,indx,'I4',err,ivals,cvalue)
             ss_abf(1,num_abf) = i
             ss_abf(2,num_abf) = ivals(1)

*            Now get the times Start: Y M D H M
             call multiread(buffer,indx,'I4',err,date,cvalue,5) 
             call read_line(buffer,indx,'R8',err,sectag,cvalue)
             call ymdhms_to_mjd(date, sectag, tt_abf(num_abf))
          end if
        
          return

c..... USR_DELBF : Allows user to remove bias flags added by track
 3200  continue
          
         call read_line (buffer,indx,'CH',err,vals, name)
*        Match the site name 
         jndx = -1
         do i = 1, num_site
            if( site_names(i).eq.name ) jndx = i
         end do
         i = jndx
*        See if we match site
         if( i.gt.0 ) then
             num_rbf = num_rbf + 1
             if( num_rbf.gt.max_edits ) then
                 write(*,*) '**DISASTER** Too many User bias flags ',
     .                      'specified. Max allowed is ', max_edits
                 stop 'TRACK: Too many user BFs specified'
             endif
             call read_line(buffer,indx,'I4',err,ivals,cvalue)
             ss_rbf(1,num_rbf) = i
             ss_rbf(2,num_rbf) = ivals(1)

*            Now get the times Start: Y M D H M
             call multiread(buffer,indx,'I4',err,date,cvalue,5) 
             call read_line(buffer,indx,'R8',err,sectag,cvalue)
             call ymdhms_to_mjd(date, sectag, tt_rbf(num_rbf))
          end if
        
          return

c....  REF_NEU: Defines reference coordinate of NEU output
 3300  continue
          call multiread(buffer,indx,'R8',err,vals8,cvalue,3)
          if( err.eq.0 ) then
              do i = 1,3
                 ref_xyz(i) = vals8(i)
              end do
          endif
          return

c.... USE_GPTGMT: Use GPT and GMF mapping functions
 3400 continue
          atm_mtt = .false.
*         See if the reference relative humifity passed
          call read_line(buffer,indx,'R8',err,vals8, cvalue)
          if( err.eq.0 ) ref_rel_hum = vals8(1)

          return
  
c.... TIME_UNIT: Sets time unit for process noise: All values
*     converted to noise-per-second (not epoch as in track)
 3500 continue
          if( usr_interval.eq.0 ) then
              write(*,3510) 
 3510         format('ERROR: INTERVAL command must be used before',
     .               'TIME_UNIT command')
              stop 'TRACK: Commands out of order'
          endif
          call  GetWord(buffer, cval, indx)
          call casefold(cval)
          if( cval(1:1).eq.'E' ) then
              time_unit = 'epoch'
              tu_to_ep = usr_interval
          elseif( cval(1:1).eq.'S' ) then
              time_unit = 'second'
              tu_to_ep = 1
          elseif( cval(1:1).eq.'M' ) then
              time_unit = 'minute'
              tu_to_ep = 60
          elseif( cval(1:1).eq.'H' ) then
              time_unit = 'hour'
              tu_to_ep = 3600
          elseif( cval(1:1).eq.'D' ) then
              time_unit = 'day'
              tu_to_ep = 86400
          else
              write(*,3520) buffer(1:trimlen(buffer))
 3520         format('WARNING: Unknown time unit in ',a,/,
     .               'Defaulting to epoch unit')
          endif
          return 

c.... SET_MSEC_BF: When used command will force bias flags
*         to be added at mill-sec reset points.  
 3600 continue
          set_msec_bf = .true.
          return

*..... AMB_SET: set parameters for fixing ambiquities
* AMB_SET <Rel Rank> <Num Av> <LC sigma> <WL sigma> <MW wght> <EX wght> <min_ambsig> <min_exsig> <Max Chi>
 3700  continue
*       
          call read_line(buffer, indx,'R4',err,vals,cvalue)
          if( err.eq.0 .and. vals(1).gt.0 )  relrank_limit = vals(1)

          call read_line(buffer, indx,'R4',err,vals,cvalue)
          if( err.eq.0 .and. vals(1).gt.0 )  float_limit(1) = vals(1)
          call read_line(buffer, indx,'R4',err,vals,cvalue)
          if( err.eq.0 .and. vals(1).gt.0 )  float_limit(2) = vals(1)

*         Get the new parameters for wl_fact, lg_fact, and max_fit
          call read_line(buffer, indx,'R4',err,vals,cvalue)
          if( err.eq.0 .and. vals(1).ge.0 )  wl_fact = vals(1)
          call read_line(buffer, indx,'R4',err,vals,cvalue)
          if( err.eq.0 .and. vals(1).ge.0 )  lg_fact = vals(1)
          call read_line(buffer, indx,'R4',err,vals,cvalue)
          if( err.eq.0 .and. vals(1).gt.0 )  min_ambsig = vals(1)
          call read_line(buffer, indx,'R4',err,vals,cvalue)
          if( err.eq.0 .and. vals(1).gt.0 )  max_fit = vals(1)

          return

c.... CSV_ROOT : Name of the root file for the CSV format used for
*     web output 
 3800 continue
         call read_line(buffer,indx,'CH',err,vals, cvalue)
         call subhome(cvalue)
         if( trimlen(prt_root).gt.0 ) then
             call sub_char(cvalue(1),'?',trim(prt_root))
         end if
         csv_root = cvalue(1)
         return 

c.... DCB_file: Name of DCB file
 3900 continue
         call read_line(buffer,indx,'CH',err,vals, dcb_file)
         call subhome( dcb_file)
         return 

c     STATUS <type> <# epochs>
 4000 continue
         call read_line(buffer,indx, 'CH',err, vals, status_type)
         call casefold(status_type)
         call read_line(buffer,indx,'I4',err,ivals,cvalue)
         call report_error('IOSTAT',err,'decod',buffer,0,
     .        'get_missrt/STATUS')
         if( err.eq.0 .and. ivals(1).gt.0 ) status_int = ivals(1)

         return

c     UPDATE_FILE <file name>
 4100 continue
         call read_line(buffer,indx,'CH',err,vals, upd_file)
         call subhome(upd_file)
         updread = .false.
         needupd = .true.
         return 

c     RESET <All/list of sites>
 4200 continue
*        Loop reading off site names or All
         err = 0
         do while ( err.eq.0 )
            call read_line(buffer,indx,'CH',err,vals, name)
            if( err.eq.0 ) then
                call casefold(name)
                if( name(1:3).eq.'ALL' ) then
                    do j = 1,num_site
                       RESET(j) = .true.
                    end do
                else
                    jndx = 1
                    call get_cmd(name, site_names, num_site, j, jndx)
                    if( j.gt.0 ) RESET(j) = .true.
                end if
             end if
          end do
          return

*.... DD_SET <Jump (cycle)> <Min Number> 
 4300 continue
         call read_line(buffer, indx,'R8',err,val8,cvalue)
         if( err.eq.0 .and. val8.gt.0 ) dd_jmp =val8
         call read_line(buffer,indx,'I4',err,ivals,cvalue)
         if( err.eq.0 .and. ivals(1).gt.0 ) min_ldd = ivals(1)
         return

*.... BNC_VERS : Set BNC version explicitly
 4400 continue
         call read_line(buffer,indx,'I4',err,ivals,cvalue)
         if( ivals(1).lt.2500 ) then
             BNC_Vers = 2500
         elseif( ivals(1).lt.2900 ) then
             BNC_Vers = 2600
         else    ! Default latest version
             BNC_Vers = 2900
         endif
         return

***** Unknown command
 9900 continue
          write(*,9920) buffer(1:trimlen(buffer))
 9920     format('TRACK: Unknown command: ',a)
           
   
c.... Thats all
      end
 
