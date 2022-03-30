      program tscon

      implicit none

*     Program to generate PBO time series files from reason XYZ files
*
*     Usages:
*     xyz_to_obo <dir> <prod_id> <cmd file> <list Reason XYZ files/<list>.lst>
*     
*     where <dir>  -- directory to put the time series in.  
*           <prod_id> is product id with form: jpl.final_frame.  Characters
*                  5:9 will used for time series type (normally rapid or final)
*           <list of Reason XYZ files>
*     
*     PROD_ID types valid in PBO:
*     Format <cen>.<series>_<frame>+<type>
*     <cen>     - Center 3 characters (bsl cwu pbo mit) or special product
*                 code e.g., aug for Augustine Volcano, aku for Akutan volcano
*     <series>  - Orbit series: final or rapid
*     <frame>   - Frame type: loose or frame 
*     <type>    - Optional type.  If not given series name is used.  Additional
*                 type put in the final series is supplemental run (suppl).
* MOD TAH 141123: Added +RESET to product ID to reset the referencce position
* MOD TAH 150106: Added new Measures .xyz format e.g. URL
* http://sopac-ftp.ucsd.edu/pub/timeseries/measures/rawXyzTimeSeries.Measures_Combination.20150101.tar.gz
* MOD TAH 180126: Added features to (1) fix UNR Easy errror (FIX_UNR command); (2) Reapply 
*                 adjustments based on read new scale data files derived from GIPSY x-files 
*                 (sh_update_xfscale script needs to run and kept updated) (REP_SCALE command);
*                 (3) convert the Measures NEU files. (.neu extent)
* MOD MAF 191218: Added assumption that UNR/NGL .txyz2 files without a reference frame
*                 in the file name are the latest IGS14 files.
* ftp://cddis.gsfc.nasa.gov/pub/GPS_Explorer/latest/WNAM_Raw_TrendNeuTimeSeries_comb_20171216.tar.gz
* where also GLB and sopac/jpl
*
* UNIT number assignment
*   99 -- command file if specified
*  100 -- Input fime series file
*  101 -- EQ_file read
*  102 -- Aprirori coordinate file 
*  103 -- Scale replacement file
*  104 -- Used to read station.info and then antmod.dat

*  200 -- Output time series file 
*  201 -- Not used
*  202 -- Output frame transformation parameters
    
*     Stsndard PROD_ID
*     pbo.rapid_frame
*     pbo.final_frame
*     pbo.final_frame+suppl

      include 'tsfit.h'
      include 'tscon.h'

      integer*4 len_run, nr, ierr, ns, rcpar, trimlen, i
     .,      unc   ! Unit number for command file (0 is not given or error)
     .,      fl   ! File name length
     .,      indx ! String pointer


      logical list  ! Set true is reading list of files
     .,       cont  ! Set true if we should continue loopsing
     .,       done  ! Set true no more files given
     .,       use_ref  ! Set true if a pbo ts files is read and
                    ! has a reference coordinate in it.

      integer*4 indx_RESET  ! Index for +RESET

      character*256 runname ! Name of file read from runstring
      character*4 cdum

      real*8 emjd  ! Evaluation MJD for stablizing coordinates
      character *80 line    ! DEBUG 

****  OK, Read the runsting (for the moment just by position)
      len_run = rcpar(1,tsdir)
      ln_tsdir = len_run
      if( len_run.eq.0 ) then
        call proper_runstring('tscon.hlp','tscon',1)
      endif
      len_run = rcpar(2,prod_id)

* MOD TAH 051129: Extract type from prod_id: Format pbo.final_frame+suppl
      ts_ref_type = prod_id(5:9)
* MOD TAH 141123: See if +RESET used 
      indx_reset = index(prod_id,"+RESET")
      if( indx_reset.gt.0 ) then
*         remove option from string
          prod_id(indx_reset:) = ' '
          len_run = indx_reset-1
      end if
          
      if( len_run.gt.15 ) then
          ts_ref_type = prod_id(17:21)
          prod_id = prod_id(1:15)
      endif

      ln_prod_id = trimlen(prod_id)

      ts_refresh = .true.
      sigma_scale = 1.d0

***** Get command file name
      nr = 3
      len_run = rcpar(3,comfile)
      if( len_run.gt.0 ) then
!         See if scale factor passed
          if( comfile(1:2).eq.'-s' .or. comfile(1:2).eq.'-S' ) then
              indx = 4
              call read_line(comfile, indx, 'R8',ierr, sigma_scale,
     .                       cdum)
              call report_error('IOSTAT',ierr,'decoding',comfile,1,
     .            'tscon/sigma_scale')
              nr = nr + 1
!             Get command file name
              len_run = rcpar(nr,comfile)
              unc = 99
              open(unc,file=comfile,iostat=ierr, status='old')
              call report_error('IOSTAT',ierr,'open',comfile,0,
     .             'tscon')
              if( ierr.ne.0 ) unc = 0

          else
!             Get command file name
              unc = 99
              open(unc,file=comfile,iostat=ierr, status='old')
              call report_error('IOSTAT',ierr,'open',comfile,0,
     .             'tscon')
              if( ierr.ne.0 ) unc = 0
          endif
      else
          unc = 0
      endif
      
***** OK Loop over input files
      num_ent = 0
      num_site = 0
      num_code = 0 
      len_run = 1
      first_mjd = 1.d20
      last_mjd = 0.d20
      save_ref_xyz = 0
      tsprog = 'tscon'
      jpl_xyz   = .false.
      mea_xyz   = .false.
      mea_neu   = .false.
      rep_scale = .false.
      fix_unre  = .false.
      reference_frame = ''
      repscl_file = ''
      out_org = ' '
      worg = .false.

****  See if file with list or if list passed
      done = .false.
      if( indx_reset.gt.0 ) then
          write(*,'(a)') "TSCON: Reference position reset"
      end if

***** Read the commnd file first (to see we need to fix 
*     UNR East problem):
* MOD TAH 180306: Added 'P' pre-read option to get just the
*     FIX_UNRE and REP_scale commands
      if( unc.gt.0 ) then
          write(*,120) trim(comfile)
120       format('TSCON: Pre-read Command file ',a)
          call read_cmd_file(unc,'P')
      endif

      do while ( .not.done )
         nr = nr + 1
         len_run = rcpar(nr,in_file)
         if( len_run.gt.0 ) then 
***          See if this is a .lst file with a list of files
             fl = trimlen(in_file) 
*            See if JPL format option
             if( in_file(1:4).eq.'JPL ' ) then
                jpl_xyz = .true.
             else if ( in_file(1:4).eq.'MEA ' ) then
                mea_xyz = .true.
             else if( in_file(fl-3:fl).eq.'.lst' ) then
*                This is a list file, so open this
*                file and read names
                 open(102,file=in_file,status='old',iostat=ierr)
                 call report_error('IOSTAT',ierr,'open',in_file,
     .                              0,'tsfit/main')
                 do while ( ierr.eq.0 )
                     read(102,'(a)', iostat=ierr) in_file
                     if( ierr.eq.0 .and. trimlen(in_file).gt.0 .and.
     .                   in_file(1:1).ne.'#' .and. 
     .                   in_file(1:1).ne.'*' ) then
*                        Try to read in files
                         call trimlead(in_file)
                         call read_in_xyz
                     endif
                 enddo
                 close(102)
             else
                 call read_in_xyz
             endif
         else
             done = .true.
         endif 
      end do

****  Now that we have all the stations read in command file
      if( unc.gt.0 ) then
* MOD TAH 180228: Moved read command above data read in case
*         went to fix UNR problems
* MOD TAH 180306: Read the rest of the commands in the command file
          rewind(unc)
          call read_cmd_file(unc,'A')
*         Now loop over all days
* MOD TAH 130326: Only stabialize if num_stab > 0
          if ( nin_stab .gt. 3 ) then

* MOD TAH 170720: See if what to output the frame origin results
             if( worg ) then
*               Try to open output file
                open(202,file=out_org, iostat=ierr, status='unknown')
                call report_error('IOSTAT',ierr,'open',out_org,0)
                if( ierr.ne.0 ) then
                    worg = .false.
                else
                    call write_org_head(202,ierr)
                    if( ierr.ne.0 ) worg = .false.
                endif

             end if
             do i = 0,nint(last_mjd-first_mjd)
                 emjd = first_mjd + i
                 call stab_sys(emjd)
             end do
          else
             write(*,220) nin_stab
 220         format('There are only ',i2,' stab sites;',
     .              ' No frame change applied')
          endif

      end if

* MOD TAH 170721: If prod_id is NONE stop at this point
      if( prod_id(1:4).eq. 'NONE' ) then
          write(*,310) 
 310      format('# TSCON: Product ID NONE so no output timeseries')
          stop 'TSCON: Product ID NONE'
      endif

      call systime( date_rel, sec_rel) 

****  OK, now loop over all the codes
      do ns = 1, num_code

*        Generate TS file name
         ts_file = tsdir(1:ln_tsdir) // '/' // in_code(ns) // '.' // 
     .             prod_id(1:ln_prod_id) // '.pos'

*****    All these series are new time series
         open(200,file=ts_file,iostat=ierr,status='old')

         num_ts = 0
         use_ref = .false.
         if( ierr.eq.0 ) then
             do i = 1,3
                ref_xyz(i) = 0.0
             end do
* MOD TAH 121214: Replaced read_ts with ts_util read_tspos
!             call read_ts(200)
             call read_tspos(200)

             ! If ref_XYZ looks OK, then use these values
* MOD TAH 141123: Use_ref unless prod_id has +RESET
             if( sqrt(ref_xyz(1)**2+ref_xyz(2)**2+ref_xyz(3)**2).gt.
     .           6.d6 .and. indx_reset.eq.0 ) then
                use_ref = .true.
             end if
         else   ! See if have saved ref_xyz values (.pos in and no existing files)
             if( sqrt(save_ref_xyz(1,ns)**2+save_ref_xyz(2,ns)**2+
     .           save_ref_xyz(3,ns)**2).gt.
     .           6.d6 .and. indx_reset.eq.0 ) then
                use_ref = .true.
                ref_xyz = save_ref_xyz(:,ns)
             endif
         endif

*****    Even when an existing series is read, reset the reference
*        points.

         new_ts = .true.

******   OK, merge new entries with old
         call merge_ts(ns)

****     See if we should read old series
         if( .not.use_ref ) then 
             call gen_refdata(ns)
         end if

* MOD TAH 080325: Before merging correct input series for any
*        East jumps due to latitude shifts

         call remove_ejmp(ns)

*****    Now close the input file and re-open to write
         close(200)
         open(200,file=ts_file,iostat=ierr,status='unknown')
         call write_ts(200, ns)
      end do

***** Thats all
      end

CITTLE READ_IN_XYZ

      subroutine read_in_xyz

*     Read in different formated files that are XYZ and also dNEU format.
*     Add back the offsets that have been removed in Measures Clean NEU
*     format.

      implicit none

      include 'tsfit.h'
      include 'tscon.h'
      include '../includes/const_param.h'

      integer*4 mmx_off      ! Maximum number of offsets allowed in Measures
                             ! NEU CleanTrended files.
      parameter ( mmx_off = 128 ) 


      character*16 gsite_full

      integer*4 date(5), ierr, jerr, i, j, k, ns, ne, yr, doy
      integer*4 trimlen, indx, fl

      integer*4 nent_start ! Entry number for the start of the current
                  ! station.  (Data layout allows sites to be mixed for
                  ! example when reading multiple org files.  This does
                  ! not happen when reading time series files.

      real*8  gmjd, pos_xyz_fin(3),xyz_std(6), unc_geod(3),llu_std(3),
     .        pos_neu_fin(3), neu_std(6), sec, decyrs, decmjd

      real*8 xyz_cor(3)   ! XYZ correlations (XY,YZ, XZ in txyz2 format)
* See http://geodesy.unr.edu/gps_timeseries/README_txyz2.txt
* MOD TAH 140919: Changed this to XY XZ YZ since most other groups including
*     PBO use this sequence.  UNR TXYZ2 files are swapped
* MOD TAH 140919: Added jpl xyzzeries format described at:
* ftp://sideshow.jpl.nasa.gov/pub/JPL_GPS_Timeseries/repro2011b/raw/position/xyzseries/00_README.format
* MOD TAH 180228: Added Measures NEU formats (both versions) and offset 
*     information that is added back to data.

      real*8 ucov_neu(3,3), ucov_xyz(3,3) ! "Unit" NEU convarinace matrix
                   ! transformed to XYZ to scale sigma
      real*8 tcov(3,3)  ! Scratch space.
      real*8 temp       ! Used swap correlations

      real*8 neuscale  ! Scale for unit NEU matrix
      real*8 unc_llu(3)  ! Computed Lat/long/height
      real*8 loc_coord(3) ! Co-lat, long, height (rad,rad. m)
      real*8 rot_mat(3,3) ! Rotation from XYZ to NEU (and transpose)
      real*8 j2000sec     ! Seconds since J2000 (JPL xyzseries format).
      real*8 dNEU(3), sNEU(3)  ! NEU dposition and sigmas from Measures NEU files/
      real*8 dXYZ(3)           ! Corresonding dXYZ to dNEU  
      real*8 min_dt       ! Closest time while cheching scale estimates

      logical neuset ! Set true when we know how to scale XYZ sigmas

      logical xyz2, xyzu, xyzj, xyzs, xyzm, xyzn
                           ! Set true if XYZ coordinates are from Blewitt format (xyz2) or
                           ! usgs xyz files (xyzu) 
                           ! JPL in name (xyzj) 
                           ! JPL xyzseries extent (xyzs)
                           ! New Measures XYZ format (xyzn) 
                           ! Measures format with no correlations (xyzm)
      logical gdcat        ! GIPSYX GDCAT format (Added 200804)
      logical xyzr, xyzrcorr ! TrendXYZ files from Measures.
      logical neum         ! Measures .neu dNEU files RAW files
      logical neuc         ! Measures ,neu dNEU files Cleaned with different format
                           ! and different reference coordinate format.
      logical lineread     ! Set true if line with data has been read while
                           ! reading header lines (neum format files).
      logical OK           ! Set true if OK to read line for file
      logical found        ! Set true when matching epoch to repscl estimates

      integer*4 nameoff    ! Offset of to get to site name in file name
 
      integer*4 ncurr      ! Current in station.info antenna information 
                           ! for fixing UNR east error.
      integer*4 J2000_sec  ! Second from J2000
      integer*4 r,c        ! Row and column 

      integer*4 zone       ! Zone for UTM coordinates
     
      character*1 hemi    ! Hemisphere for UTM coordinates

      character*512 line
      character*4 word
      character*8 dateword
      character*10 str
      character*1 NC, WC  ! Lat and Long N/S or E/W (neuc format)
      integer*4 lt_dg, lt_mn, ln_dg, ln_mn   ! Lat deg, min and Long deg, min
      integer*4 lt_sgn, ln_sgn   ! Sign for latitude and longitude
      real*8    lt_sc, ln_sc, ht ! Lat and long seconds of arc and height (m)
      integer*4 is      ! Counter to find start of name for NEU Measures files.

*     Offset information for Measures Trended files
      integer*4 mnm_off   ! Number of offsets
     .,         io        ! Index of offsets (loop n e u)
     .,         ind       ! Index for offset epoch
      real*8 mea_off(4,mmx_off)  ! Offsets: DecYrs, dN dE dU.
     .,      off          ! Read value of offset in Measures format

      data ucov_neu / 1.d0, 0.0d0,  0.d0,
     .                0.d0, 1.0d0,  0.d0,
     .                0.d0, 0.0d0, 10.d0 / 

****  Open the input file: 
*     Check the type file: 
*     .xyz are Reason/SIO format
*     .cvs are Unavco/SCEC North/East/Up delta files (differnces in flavor)
*     .pos are PBO time series format
      fl = trimlen(in_file) 
      if( in_file(fl-3:fl).eq.'.csv' ) then
         call read_in_csv
         RETURN
      elseif( in_file(fl-3:fl).eq.'.pos' ) then
         call read_in_pbo( 0 )
         RETURN
      endif
*     Default is to assume XYZ files with name format of
*     Name P067Raw.xyz
* MOD TAH 140416: Mod'd to read UNR Blewitt P474.NA12.txyz2 files 
!P474 04AUG04 2004.5914 -0.244176525880759E+07 -0.474124654176784E+07  0.348703064728680E+07 0.912811128330071E-03 0.155242828238653E-02 0.111427436367456E-02  0.792729959364514E+00 -0.784330954369549E+00 -0.701224098125203E+00   0.0083   
* MOD TAH 140808: Mod'd to read USGS XYZ files. There is no header line.  File names are p542.xyz
*date, orbit_type, decimal_date, doy, x, x-sigma, y, y-sigma, z, z-sigma, x-y correlation, x-z correlation, y-z correlation, Name and date of solution file
!20040820	F	2004.6359	233	-2471753.8159	0.0009	-4611976.8391	0.0014	3636954.5386	0.0011	0.7923	-0.7090	-0.7839	Southern_California.20040820.stacov.point-2013/04/28-16:04:39 
* MOD TAH 140829: Added JPL xyz format.  Since name format is the same as USGS; use JPL as file file name to set format
*wget ftp://fringe.jpl.nasa.gov/incoming/JPL_xyzts_igs08.tar
!2010 10 20 -2377135.346847 -4909590.111199 3296071.051902 0.000812 0.001412 0.000990 0.781593 -0.671771 -0.784449

* MOD TAH 200720: Added xyzr for TrendXYZ files.
!# X, Y, Z, and sigma values are in meters
!# Dec Yr   Yr  DayOfYr       X              Y             Z        X sig    Y sig    Z sig    CorrXY   CorrXZ   CorrYZ
!2005.7740  2005  283  -2314778.9953  -4825178.9422  3458093.6644   0.0027   0.0039   0.0031   0.6381  -0.5147  -0.6362

      xyz2 = .false.
      xyzu = .false.
      xyzj = .false.
      xyzn = .false.
      xyzs = .false.
      xyzm = .false.
      xyzr = .false.   ! Set true for Measures TrendXYZ files.
* MOD MAF 20210324: Added test of whether or not correlations are missing in "xyzr" files
      xyzrcorr = .true.  ! Set false if correlations are missing
      neum = .false.
      neuc = .false.
      gdcat = .false.
      nameoff = 7
      mnm_off = 0     ! Offsets in Measures Trended format

      nent_start = num_ent + 1   ! Entry number for start of this station.

      if( in_file(fl-5:fl).eq.'.txyz2' ) then
          xyz2 = .true.
          if( index(in_file,'NA12').gt.0 ) then
              nameoff = 11
          elseif( index(in_file,'IGS08').gt.0 ) then
              nameoff = 12
          else  ! Assume IGS14
              !stop 'Unknown UNR txyz2 frame'
              nameoff = 6
          endif
      elseif( in_file(fl-9:fl).eq.'.xyzseries' ) then
          xyzs = .true.
          nameoff = 10
      elseif ( fl.eq.8 .or. in_file(fl-8:fl-8).eq.'/' ) then  ! <dir>/aaaa.xyz type file name
          ! See if JPL in file name
          if( jpl_xyz ) then
             xyzj = .true.
          else if ( mea_xyz ) then
             xyzn = .true.
          else
             xyzu = .true.
          endif 
          nameoff = 4
* MOD TAH 200720: Added test for mea_raw files.p495_RawJumps_TrendXYZ.xyz 
* MOD MAF 20210105: Added test for MEaSUREs files, e.g. p617RawTrend.xyz or p617RawJumpsTrend.xyz
      elseif( in_file(fl-11:fl).eq.'TrendXYZ.xyz'
     .       .or. in_file(fl-8:fl).eq.'Trend.xyz' ) then
           xyzr = .true.
           if ( in_file(fl-11:fl).eq.'TrendXYZ.xyz' ) then  ! "RawJumps_TrendXYZ.xyz"
             nameoff = 22
           elseif ( in_file(fl-13:fl).eq.'JumpsTrend.xyz' ) then  ! "RawJumpsTrend.xyz"
             nameoff = 17
           else  ! "RawTrend.xyz"
             nameoff = 12
           endif
      elseif( in_file(fl-3:fl).eq.'.neu' ) then
*         Could be different formats depending on it raw or Cleaned
*         Decide based on length of file name 
*         Eaw File name  : p067r.neu
*         Clean File Nae : p067CleanTrend.neu
          is = fl
          do while ( in_file(is:is).ne.'/' .and. is.gt.0 )
             is = is - 1
          end do
          if ( fl-is.eq.9 ) then  
             neum = .true. ! Raw format
             nameoff = 5   ! Where the site name is relative to last char +3
          else
             neuc = .true.   ! Clean format
             nameoff = fl-is-4
          endif
* MOD TAH 200804: See if GIPSYX gdcat file
      elseif ( in_file(fl-5:fl).eq.'.gdcat' ) then
          gdcat = .true.
          nameoff = 6
      else

          xyzm = .true.     ! Measures format with no correlations
          nameoff = 7
      endif 

*     Get site name from file name.
      site_name = in_file(fl-(nameoff+3):fl-nameoff) // '_GPS' 

      call casefold(site_name)
      gsite_full = site_name(1:4) // '_' // prod_id(1:9)

*     See if we can match site name
      indx = 0
      call get_cmd(site_name,gsite_names,num_site,ns,indx)
      if( ns.le.0 ) then
          num_site = num_site + 1
          if( num_site.gt.max_site ) then 
              print *,'Too many sites ',num_site, max_site
              call report_stat('FATAL','TSCON','read_in_xyz','',
     .             'Too many sites',max_site)
          endif

          num_code = num_site
          ns = num_site
          gsite_names(ns) = site_name
          in_code(ns) = site_name(1:4)
          in_full(ns) = gsite_full
      end if


*     Now open and read file.
      open(100,file=in_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_file,0,'tscon')

      if( ierr.ne.0 ) RETURN

      sec = 0.d0
      neuset = .false.
      do j = 4,6
         xyz_std(j) = 0.d0
         neu_std(j) = 0.d0
      end do

*     If this is Measures NEU file read the header to get the Reference
*     XYZ values
      lineread = .false.
      if ( neum ) then 
*        Look for lines that look like:
*        # Reference_X : -2675940.067333
*        # Reference_Y : -4452982.646223
*        # Reference_Z : 3687902.754937
         do while ( .not.lineread )
             read(100,'(a)', iostat=ierr) line
             if( line(1:15).eq.'# Reference_X :') then
                 read(line(16:),*,iostat=jerr)  ref_XYZ(1) 
             elseif( line(1:15).eq.'# Reference_Y :') then
                 read(line(16:),*,iostat=jerr)  ref_XYZ(2) 
             elseif( line(1:15).eq.'# Reference_Z :') then
                 read(line(16:),*,iostat=jerr)  ref_XYZ(3) 
             elseif( line(1:1).ne.'#' ) then
                 lineread = .true.     ! First data record
             endif
         enddo
      endif
      if( neuc ) then ! Look for Reference Position Line
C# Reference position: N35 33  6.303454 W121  0 10.650662 106.99013 Datum: WGS84
         io = 0   ! Index for NEU offsets (incremented when 'componemt' line found
         do while ( .not.lineread )
            read(100,'(a)', iostat=ierr) line
            if( line(3:20).eq.'Reference position' ) then
*               Decode line with formated read
                read(line,220) NC, lt_dg, lt_mn, lt_sc, 
     .                         WC, ln_dg, ln_mn, ln_sc
 220            format(22x,a1,I2,1x,I2,1x,F9.6,1x,a1,I3,1x,I2,1x,F9.6)
                read(line(57:),*) ht
                lt_sgn = 1
                if( NC.eq.'S' ) lt_sgn = -1
                ln_sgn = 1
                if( WC.eq.'W' ) ln_sgn = -1
                loc_coord(1) = (90-lt_sgn*(lt_dg + lt_mn/60.0d0 + 
     .                                    lt_sc/3600.d0))*pi/180   ! Co-latitude (rad)
                loc_coord(2) = ln_sgn*(ln_dg + ln_mn/60.0d0 + 
     .                                    ln_sc/3600.d0)*pi/180    ! Longitude (rad)
                loc_coord(3) = ht
                call GEOD_to_XYZ(loc_coord,ref_XYZ)
             elseif( line(6:11).eq.'offset' ) then
                read(line(15:21),*,iostat=jerr) off
                call report_error('IOSTAT',jerr,'off decod',line,0,
     .                            'READ_IN_XYZ')
                read(line(50:58),*,iostat=jerr) decyrs
*#*   offset 1:-208.83  +/- 0.28  mm (2010-04-05 [2010.2589])
                if( io.eq.1 ) then   ! Reading North first so add to list
*                                      (not clear if time sorted so check
*                                      time order
                    ind = mnm_off + 1
                    do j = 1, mnm_off
                       if( decyrs.lt.mea_off(1,j) ) then
*                          Value falls before current one so move up list
*                          and add.
                           do k = mnm_off,j,-1
                              mea_off(1:2,k+1) = mea_off(1:2,k)
                           end do
                           ind = j
                           exit
                       endif
                    end do
*                   Increment numbet offsets and save values
                    mnm_off = mnm_off + 1
                    if( mnm_off.gt.mmx_off ) 
     .              call report_stat('FATAL','tscon','read_in_xyz','',
     .                  'Too many Meausures offsets ',mmx_off)
                    mea_off(1,ind) = decyrs
                    mea_off(2,ind) = off
                 else     ! Find the entry
                    ind = 1
                    do j = 1, mnm_off
                       if( mea_off(1,j).eq.decyrs ) ind = j
                    end do
                    mea_off(io+1,ind) = off
                 endif
             elseif( line(24:32).eq.'component') then
*#                    n component
                 io = io + 1
             elseif( line(1:1).ne.'#' ) then
                 lineread = .true.     ! First data record
             endif
         enddo
*        Report any offsets removed.(accumultive
         do j = 1, mnm_off
            if( j.gt.1 ) mea_off(2:4,j) = 
     .                   mea_off(2:4,j-1) + mea_off(2:4,j)
            write(*,240) site_name(1:4),mea_off(:,j)
 240        format('OffSet ',a,' DecYrs ',F10.4,' Sum dNEU ',
     .             3(F8.2,1x),' mm')
         end do 
      endif
            
*     If this is UNR see if east is be fixed
      if( fix_unre .and. xyz2 ) then
          call get_unr_ant( site_name(1:4) )
      end if

      ncurr = 1    ! Default to first entry in antenna information 
      date    = 0      ! F90 initalization of date
      date(2) = 1      ! Set Jan for month for yr doy read
      date(4) = 12     ! Set 12 hrs for formats with only date
      do while ( ierr.eq.0 )
*        See if we have read the data line already (only first time
*        for the NEUM format).
         if( .not.lineread) read(100,'(a)', iostat=ierr) line
         lineread = .false.

* MOD TAH 200720: Set OK to read if ierr == 0 and not # in
*        first column of xyzr format
         OK = .true.
         if( ierr.ne. 0 ) OK = .false.
         if( xyzr .and. line(1:1).eq.'#' ) OK = .false.

         if( OK ) then
            if( xyz2 ) then 
                read(line,*,iostat=jerr) word, dateword, decyrs, 
     .               (pos_xyz_fin(j), j=1,3), (xyz_std(j),j=1,3),
     .               (xyz_cor(j),j=1,3)
* MOD TAH 200224: Replaced decyrs call becuse number of digits is not 
*               enough and vagueness of number of days in a year (leap-years)
                call dateword_to_mjd( dateword, decmjd )
                ! call decyrs_to_mjd(decyrs, decmjd)
                gmjd = nint(decmjd) + 0.5d0  ! Gets around truncation level in inpuy
*               Swap YZ and XZ correlations so that order is XY XZ YZ
                temp = xyz_cor(3)
                xyz_cor(3) = xyz_cor(2)
                xyz_cor(2) = temp
                call mjd_to_ymdhms(gmjd, date, sec)

****            See if we are going to fix the east error in UNR time series
                if( fix_unre ) then
                   dNEU = 0.d0  ! NEU coorection set to zero initially 
 
                   if ( ncurr.lt. num_stinf ) then 
                       if( gmjd.gt.stinf_mjd(ncurr+1) ) then   ! Moved to next entry
                           ncurr = ncurr + 1
                       endif
                   endif
                   dNEU(2) = 2.54573*stinf_neu(2,1,ncurr) - 
     .                       1.54573*stinf_neu(2,2,ncurr)
*                  Get the dXYZ values
                   call rotate_geod(dNEU, dXYZ,'NEU','XYZ',
     .                     pos_xyz_fin, loc_coord,rot_mat)
*                  Appply the correction (wrong sign originally used
                   do j = 1,3
                      pos_xyz_fin(j) = pos_xyz_fin(j) - 2*dXYZ(j)/1000
                   end do
                end if 

****            See if we are adding back scale from GIPSY X-files
                found = .true.
                jerr = 0
                if( rep_scale ) then
*                   find value in repscl estimates previously read
                    found = .false.
                    min_dt = 9999
                    do j = 1,num_repscl
                       if( abs(repscl(1,j)-gmjd).lt.min_dt ) 
     .                     min_dt = abs(repscl(1,j)-gmjd)
                       if( abs(repscl(1,j)-gmjd).lt.0.1 ) then
                          found = .true.
                          ! Add back the scale that has been removed.
                          pos_xyz_fin = pos_xyz_fin*(1.d0+repscl(2,j))
                          exit
                       endif
                    end do
                    if ( .not.found ) jerr = 98   ! Error message
                end if

            elseif ( xyzs ) then
!2006.34086242  -2467567.979594  -2639809.396636   5238172.295680         0.000892         0.000868         0.001619         0.703287        -0.718097        -0.694059     200102400.00  2006  5  5 12  0  0  
                read(line,*,iostat=jerr) decyrs, 
     .               (pos_xyz_fin(j), j=1,3), (xyz_std(j),j=1,3),
     .               (xyz_cor(j),j=1,3), j2000sec, date(1:5), sec
            elseif ( xyzr ) then
* MOD TAH 200720: New Measures XYZ format (correlations XY XZ YZ are the order expected).
!# Dec Yr   Yr  DayOfYr       X              Y             Z        X sig    Y sig    Z sig    CorrXY   CorrXZ   CorrYZ
!2005.7740  2005  283  -2314778.9953  -4825178.9422  3458093.6644   0.0027   0.0039   0.0031   0.6381  -0.5147  -0.6362
c               read(line,*,iostat=jerr) decyrs, date(1), date(3), 
c    .               (pos_xyz_fin(j), j=1,3), (xyz_std(j),j=1,3),
c    .               (xyz_cor(j),j=1,3)
* MOD MAF 20210324: Read line with possible variable number of fields in case of missing correlations
                if ( trimlen(line).gt.92 ) then  ! Correlation fields are present
                  read(line,*,iostat=jerr) decyrs, date(1), date(3),
     .                 pos_xyz_fin(1:3), xyz_std(1:3), xyz_cor(1:3)
                  xyzrcorr = .true.
                else  ! Continue without known correlations
                  read(line,*,iostat=jerr) decyrs, date(1), date(3),
     .                 pos_xyz_fin(1:3), xyz_std(1:3)
                  xyz_cor(1:3) = 0.d0
                  xyzrcorr = .false.
                end if
                xyz_std(4:6) = xyz_cor(1:3)
            elseif ( xyzu ) then
                read(line,*,iostat=jerr) dateword, word, decyrs,
     .              doy, (pos_xyz_fin(j), xyz_std(j), j = 1,3),
     .              (xyz_cor(j),j=1,3)
                read(dateword,'(i4,i2,i2)') date(1:3)
                date(4) = 12
                date(5) =  0
                call decyrs_to_mjd(decyrs, decmjd)
                gmjd = nint(decmjd) + 0.5d0  ! Gets around truncation level in inpuy
****            See if we are adding back scale from GIPSY X-files
                found = .true.
                jerr = 0
                if( rep_scale ) then
*                   find value in repscl estimates previously read
                    found = .false.
                    min_dt = 9999
                    do j = 1,num_repscl
                       if( abs(repscl(1,j)-gmjd).lt.min_dt )
     .                     min_dt = abs(repscl(1,j)-gmjd)
                       if( abs(repscl(1,j)-gmjd).lt.0.1 ) then
                          found = .true.
                          ! Add back the scale that has been removed.
                          pos_xyz_fin = pos_xyz_fin*(1.d0+repscl(2,j))
                          exit
                       endif
                    end do
                    if ( .not.found ) jerr = 98   ! Error message
                end if
!2010 10 20 -2377135.346847 -4909590.111199 3296071.051902 0.000812 0.001412 0.000990 0.781593 -0.671771 -0.784449
             elseif ( xyzj ) then
                read(line,*,iostat=jerr) date(1:3), 
     .              (pos_xyz_fin(j), j=1,3), (xyz_std(j), j = 1,3),
     .              (xyz_cor(j),j=1,3)
                date(4) = 12
                date(5) =  0
*           See if new Mesuures XYZ format (MOD TAH 150106)
!p070 2007-03-22T12:00:00 -1310383.7747 -4995596.3585 3733318.6654 0.0012 0.0030 0.0022 0.6073 -0.5618 -0.8400
            elseif ( xyzn ) then
*               Read in two parts since XYZ part might not be fixed format.
                read(line,310,iostat=jerr) word, date 
 310            format(a4,1x,I4,1x,I2,1x,I2,1x,I2,1x,I2)
                read(line(25:),*,iostat=ierr )
     .              (pos_xyz_fin(j), j=1,3), (xyz_std(j), j = 1,3),
     .              (xyz_cor(j),j=1,3)
*           Measures NEU format
            elseif( neum .or. neuc ) then
*2004.0342   -56.7186   105.4902    -3.3121  0.0070  0.0060  0.0070
                 if( neum ) then
                     read(line,*,iostat=jerr) decyrs, dNEU, sNEU
*                    Check for gross errors (total > 1 m^2)
                     if( jerr.eq.0 ) then
                         if( dNEU(1)**2+dNEU(2)**2+dNEU(3)**2 > 1000 ) 
     .                                                   jerr = -99 
                     endif 
                 else
                      read(line,320,iostat=jerr) decyrs, yr, doy, dNEU, 
     .                   sNEU
 320                  format(F9.4,1x,I4,1x,I3,6(F8.2))
C2011.4315 2011 157   27.85 -122.00-1214.36    4.10    3.50    4.20 
*                     See of we have any offsets to apply.
                      ind = 0
                      do j = 1,mnm_off
                         if( decyrs.ge.mea_off(1,j) ) ind=j
                      end do
                      if( ind.gt.0 ) then
                          dNEU = dNEU + mea_off(2:4,ind)
                      endif
                      dNEU = dNEU/1000   ! Convert to meters
                      sNEU = sNEU/1000   ! Convert to meters
                 endif
         
                 call decyrs_to_mjd(decyrs, decmjd)
                 gmjd = int(decmjd) + 0.5d0  ! Gets around truncation level in input
                 call mjd_to_ymdhms(gmjd, date, sec)
*                Now convert dNEU to dXYZ; Convert sNEU to XYZ covarinace matrix,
                 call rotate_geod(dNEU, dXYZ,'NEU','XYZ',ref_XYZ,
     .                            loc_coord,rot_mat)
*                Compute Cartessian coordinates
                 do k = 1,3
                    pos_xyz_fin(k) = ref_XYZ(k) + dXYZ(k)
                 end do

                 ucov_neu(:,:) = 0.0d0 
                 do k =1,3
                    ucov_neu(k,k) = sNEU(k)**2
                 end do
                 call var_comp(rot_mat,ucov_neu, ucov_xyz,tcov, 3,3,1)
*                Now convert to XYZ sigmas and correlationa
                 do k = 1,3
                    xyz_std(k) = sqrt(ucov_xyz(k,k))
                 enddo 
                 xyz_cor(1) = ucov_xyz(1,2)/(xyz_std(1)*xyz_std(2)) ! Rho XY
                 xyz_cor(2) = ucov_xyz(1,3)/(xyz_std(1)*xyz_std(3)) ! Rho XZ
                 xyz_cor(3) = ucov_xyz(2,3)/(xyz_std(2)*xyz_std(3)) ! Rho YZ
****            See if we are adding back scale from GIPSY X-files
                found = .true.
                jerr = 0
                if( rep_scale ) then
*                   find value in repscl estimates previously read
                    found = .false.
                    min_dt = 9999
                    do j = 1,num_repscl
                       if( abs(repscl(1,j)-gmjd).lt.min_dt ) 
     .                     min_dt = abs(repscl(1,j)-gmjd)
                       if( abs(repscl(1,j)-gmjd).lt.0.1 ) then
                          found = .true.
                          ! Add back the scale that has been removed.
                          pos_xyz_fin = pos_xyz_fin*(1.d0+repscl(2,j))
                          exit
                       endif
                    end do
                    if ( .not.found ) jerr = 98   ! Error message
                end if

* MOD TAH 200804: Read the multiline gdcat file
            elseif ( gdcat ) then
                 if( index(line,'3 PARAMETERS').eq.0 ) then
                    write(*,340) trim(line)
 340                format('Expected "3 PARAMETERS" found ',a)
                 else
*                   Read lines in extract information
!1 AB43.STA.X 609983850 -2.449678808370530e+06 3.984837652881759e-03
                    do j = 1,3
                       read(100,*) r, str,   
     .                      J2000_sec, pos_xyz_fin(j), xyz_std(j)
                    enddo
!2 1 3.276850082746537e-01
                    do j = 1,3
                       read(100,'(I1,1x,I1,1x,e22.11)') r,c, temp
                       k = c*(c-1)/2 + (r-1)
                       xyz_cor(k) = temp
                    end do
                    gmjd = 51544.5d0 + nint(J2000_sec/3600.d0)/24.d0 
                    call mjd_to_ymdhms(gmjd,date,sec)
                  endif
                  

            else      ! Original measures XYZ format.
                read(line,*,iostat=jerr) decyrs, yr, doy, 
     .               (pos_xyz_fin(j), j=1,3), (xyz_std(j),j=1,3)
                date(1) = yr
                date(2) = 1
                date(3) = doy
                date(4) = 12
                date(5) =  0
            endif
            if( jerr.ne.0 ) then
               if ( jerr.eq.99 ) then ! Bad Measures data
                  write(*,410) site_name, decyrs, dNEU
 410              format('**BAD DATA** Site ',a,' Year ',F10.4,
     .                   ' dNEU ',3F10.4,' m')
               elseif ( jerr.eq.98 ) then  ! No scale value
                  write(*,420) site_name, gmjd, min_dt
 420              format('**NO SCALE** Site ',a,' MJD ',F8.2,
     .                   ' No value in scale file; min_dt ',F8.2)
               else    ! Bad line               
                  write(*,430) jerr, trim(line)
 430              format('IOSTAT error ',i4,' reading: ',a)
               endif
            endif


****        Now compute the other quantities we need
            call ymdhms_to_mjd(date, sec, gmjd)
            if( gmjd.lt.first_mjd ) first_mjd = gmjd
            if( gmjd.gt.last_mjd  ) last_mjd  = gmjd


****        Now convert XYZ to lat/long/height
            call geod_to_geod(pos_xyz_fin, unc_geod, 
     .            'XYZ', 'GEOD','WGS84','WGS84',zone,hemi)

*           Convert co-lat, long in radians to degrees
            unc_llu(1) = (pi/2-unc_geod(1))*180/pi
            unc_llu(2) = unc_geod(2)*180/pi 
            unc_llu(3) = unc_geod(3)           

****        Now get scale to convert XYZ sigmas to NEU If this is measures format
            if( xyzm ) then   ! Use approximate formula because correlations
                                    ! are not given
               if( .not. neuset ) then 
                   call XYZ_to_GEOD(rot_mat, pos_xyz_fin, loc_coord )
c....              Now transpose rot_matrix to that NEU to XYZ direction 
                   do i = 1,2
                     call dvswp(rot_mat(i,i+1),3,rot_mat(i+1,i),1, 3-i)
                   end do
                   call var_comp(rot_mat,ucov_neu,ucov_xyz,tcov,3,3,1)
*                  Now get the average scale needed
                   neuset = .true.
               end if

****           Estimate the neu sigmas based on scale factor
               neuscale = (xyz_std(1)**2/ucov_xyz(1,1)+
     .                     xyz_std(2)**2/ucov_xyz(2,2)+
     .                     xyz_std(3)**2/ucov_xyz(3,3))/3.d0
               neu_std(1) = sqrt(neuscale)
               neu_std(2) = neu_std(1)
               neu_std(3) = sqrt(10.d0*neuscale)
* MOD MAF 20210324: Added same test as above for "xyzr" files in case of missing correlations
            elseif( xyzr .and. .not. xyzrcorr ) then   ! Use approximate formula because correlations are not given
c              Reinitialize NEU unit covariance matrix
               ucov_neu(:,:) = 0.d0
               ucov_neu(1,1) = 1.d0
               ucov_neu(2,2) = 1.d0
               ucov_neu(3,3) = 10.d0
               call XYZ_to_GEOD(rot_mat, pos_xyz_fin, loc_coord )
c....          Now transpose rot_matrix to that NEU to XYZ direction 
               do i = 1,2
                 call dvswp(rot_mat(i,i+1),3,rot_mat(i+1,i),1, 3-i)
               end do
               call var_comp(rot_mat,ucov_neu,ucov_xyz,tcov,3,3,1)
*              Now get the average scale needed
               xyzrcorr = .true.

****           Estimate the neu sigmas based on scale factor
               neuscale = (xyz_std(1)**2/ucov_xyz(1,1)+
     .                     xyz_std(2)**2/ucov_xyz(2,2)+
     .                     xyz_std(3)**2/ucov_xyz(3,3))/3.d0
               neu_std(1) = sqrt(neuscale)
               neu_std(2) = neu_std(1)
               neu_std(3) = sqrt(10.d0*neuscale)
               neu_std(4:6) = 0.d0  ! Ensure NEU correlations are reset to zeros from previous data record
            else
!              We know the correlations so compute directly.
               call XYZ_to_GEOD(rot_mat, pos_xyz_fin, loc_coord )
               do i = 1,3
                   ucov_xyz(i,i) = xyz_std(i)**2
               end do
               ucov_xyz(1,2) = xyz_std(1)*xyz_std(2)*xyz_cor(1)
               ucov_xyz(1,3) = xyz_std(1)*xyz_std(3)*xyz_cor(2)
               ucov_xyz(2,3) = xyz_std(2)*xyz_std(3)*xyz_cor(3)
               ucov_xyz(2,1) = ucov_xyz(1,2)
               ucov_xyz(3,1) = ucov_xyz(1,3)
               ucov_xyz(3,2) = ucov_xyz(2,3)
               xyz_std(4) = xyz_cor(1)
               xyz_std(5) = xyz_cor(2)
               xyz_std(6) = xyz_cor(3)

               call var_comp(rot_mat,ucov_xyz, ucov_neu, tcov, 3,3,1)
               do i = 1,3
                  neu_std(i) = sqrt(ucov_neu(i,i))
               end do
               neu_std(4) = ucov_neu(1,2)/(neu_std(1)*neu_std(2))
               neu_std(5) = ucov_neu(1,3)/(neu_std(1)*neu_std(3))
               neu_std(6) = ucov_neu(2,3)/(neu_std(2)*neu_std(3))
            endif


****        Convert mm NEU to lat/long sigmas                
            llu_std(1) = neu_std(1)/Earth_rad*180/pi*1.d9
            llu_std(2) = neu_std(2)/Earth_rad*180/pi/
     .                sin(unc_geod(1))*1.d9
            llu_std(3) = neu_std(3)

****        Now get NEU values
 
            call loc_to_geod(unc_geod, pos_neu_fin)

            if( jerr.eq.0 ) then

*****         OK save this entry
              num_ent = num_ent + 1
              if( num_ent.gt.max_ent ) then
                   print *,'Too many emties.  Max allowed ',max_ent
                   call report_stat('FATAL','tscon','read_in_xyz','',
     .             'Too many enties Max ',max_ent)
              endif
              ne = num_ent   ! Short-hand version.
              in_ns(ne)  = ns
              in_cs(ne)  = ns
              in_mjd(ne) = gmjd
              in_type(ne) = ts_ref_type
              do j = 1,3
                 in_xyz(j,ne) = pos_xyz_fin(j)
                 in_neu(j,ne) = pos_neu_fin(j)
                 in_llu(j,ne) = unc_llu(j)
              end do
* MOD TAH 140829: Scale the sigmas
              xyz_std(1:3) = xyz_std(1:3)*sigma_scale
              neu_std(1:3) = neu_std(1:3)*sigma_scale
              do j = 1,6
                 in_xyz_std(j,ne) = xyz_std(j)
                 in_neu_std(j,ne) = neu_std(j)
              end do
            end if
         end if
      end do

****  Tell user were we are
      write(*,110) in_file(1:trimlen(in_file)), num_site, num_code, 
     .             num_ent, sigma_scale 
 110  format('File: ',a,' Sites ',i5,' Codes ',i5,' Entries ',i10, 
     .       ' SigScale ',F10.3)
      return
      end 

CTITLE READ_CMD_FILE

      subroutine read_cmd_file(unc, opt)  
 
      implicit none

      include 'tsfit.h'
      include 'tscon.h'

      integer*4 unc  ! Unit number for command file
      character*(*) opt   ! Option: 'P' for pre-read for FIX_UNRE and REP_SCALE commands
                          ! 'A' to read all other commnands.

* LOCAL
      integer*4 i,j, k, ih, nt    ! Loop counter
     .,    ierr, jerr, kerr      ! IOSTAT errors
     .,    trimlen
     .,    indx     ! Pointer to position in string.
     .,    date(5)  ! year, mon, day, hr, min

      logical kbit 
      logical hierach_out  ! Set true when hierarch sites listed

      real*8 sec    ! Seconds (set to zero)


      character*256 line  ! Line read from file
      character*16 cmd    ! Command extracted from line
      character*16 cval   ! Character value
      character*8 stab_out(10)     ! List of stabalizations site

****  Initialize defaults
      num_nonsec = 0
      nin_stab = 0

      do i = 1, max_stab
         do j = 1, max_hierach
            stab_len(i,j)  = 0
            stab_site(i,j) = ' '
            stab_times(1,i) = 0
         end do
      end do
      cnd_hgt_var = 1.d0  ! Since it multiplies the height sigma
      use_ratio = 3.d0
      stab_it   = 4
      stab_nsig = 4.d0
      stab_rel = 0.5d0

      stab_min_dh  = 0.005d0
      stab_min_dne = 0.0005d0
      cnd_parts_bits = 127  ! XT, YT, ZT, XR, YR, YZ

      debug = .false.
      CMest = .false.
      stinf_file  = ' '
      antmod_file = ' '


****  OK: Start reading the command file
      ierr = 0
      do while ( ierr.eq.0 )
          read(unc,'(a)',iostat=ierr) line
          if( line(1:1).eq.' ' .and. ierr.eq.0 .and.
     .        trimlen(line).gt.0 ) then

****          OK process command
*             OK: Try to read command
              indx = 0
              call GetWord(line, cmd, indx)
              call casefold(cmd)


*@ EQ_FILE <File Name>
             if(   cmd(1:2).eq. 'EQ'  .and. opt.eq.'A' ) then   ! Eq_file
*                Get name of eq_file: and read
                 call GetWord(line, eqfile, indx)
* MOD TAH 131111: Allow ~ substitution
                 call subhome(eqfile)
                 call tsread_eqfile

*             See which comamnd passed
*@ APR_FILE <File Name>
              elseif( cmd(1:3).eq. 'APR'  .and. opt.eq.'A' ) then   ! apr_file
*                 Get name of eq_file: and read
                  call GetWord(line, apr_file, indx)
* MOD TAH 131111: Allow ~ substitution
                  call subhome(apr_file)
                  call readaprf

*@ STAB_SITE <list of stablization sites> (-name to remove)
             elseif ( cmd(1:6).eq.'STAB_S' .and. opt.eq.'A' ) then   ! stab_site
*                 Extract the list of values
                  call readstab(line, indx)

*@ POS_ORG <components>
             elseif( cmd(1:3).eq.'POS' .and. opt.eq.'A' ) then
*                 Get the components
                  cnd_parts_bits = 0
                  call casefold(line)
                  if( index(line,'XR').gt.0 ) 
     .                call sbit(cnd_parts_bits,1,1)
                  if( index(line,'YR').gt.0 ) 
     .                call sbit(cnd_parts_bits,2,1)
                  if( index(line,'ZR').gt.0 ) 
     .                call sbit(cnd_parts_bits,3,1)
                  if( index(line,'XT').gt.0 ) 
     .                call sbit(cnd_parts_bits,4,1)
                  if( index(line,'YT').gt.0 ) 
     .                call sbit(cnd_parts_bits,5,1)
                  if( index(line,'ZT').gt.0 ) 
     .                call sbit(cnd_parts_bits,6,1)
                  if( index(line,'SC').gt.0 ) 
     .                call sbit(cnd_parts_bits,7,1)
* MOD TAH 170720: Added xcm, ycm, zcm for CM motions
                  if( index(line,'XC').gt.0 ) then
                      call sbit(cnd_parts_bits,8,1)
                      CMest = .true.
                  endif
                  if( index(line,'YC').gt.0 ) then
                      call sbit(cnd_parts_bits,9,1)
                      CMest = .true.
                  endif
                  if( index(line,'ZC').gt.0 ) then
                      call sbit(cnd_parts_bits,10,1)
                      CMest = .true.
                  endif

*@ STAB_ITE [# iterations] [Site Relative weight] [n-sigma]
             elseif( cmd(1:6).eq.'STAB_I' .and. opt.eq.'A' ) then
                  call read_line(line, indx,'I4',jerr,stab_it,cval)
                  call read_line(line, indx,'R8',jerr,stab_rel,cval)
                  call read_line(line, indx,'R8',jerr,stab_nsig,cval)


*@ STAB_MIN [dHsig min pos] [dNEsig min pos]
             elseif( cmd(1:6).eq.'STAB_M' .and. opt.eq.'A' ) then
                 call read_line(line, indx,'R8',jerr,stab_min_dh,cval)
                 call read_line(line, indx,'R8',jerr,stab_min_dne,cval)


*@ CND_HGTV [Height variance] [Sigma ratio]
             elseif( cmd(1:3).eq.'CND' .and. opt.eq.'A' ) then
                 call read_line(line, indx,'R8',jerr,cnd_hgt_var,cval)
                 call read_line(line, indx,'R8',jerr,use_ratio,cval)

*@ DEBUG 
             elseif( cmd(1:3).eq.'DEB' .and. opt.eq.'A' ) then
                 debug = .true.

*@ TIME_RANGE <Start date> <End dat>
             elseif( cmd(1:3).eq.'TIM'.and. opt.eq.'A' ) then
                 call multiread(line, indx,'I4',jerr,date,cval,5)
                 sec = 0.0
                 if( date(4).ne.12 )
     .           call report_stat('WARNING','TSCON','Time Range',' ',
     .                'Normally dates are are 12:00 hrs',date(4))
                 call ymdhms_to_mjd(date,sec,first_mjd)
                 call multiread(line, indx,'I4',jerr,date,cval,5)
                 call ymdhms_to_mjd(date,sec,last_mjd)
*@ SIGMA_SCALE <sigma_scale> 
             elseif( cmd(1:3).eq.'SIG'.and. opt.eq.'A' ) then
                 call read_line(line, indx,'R8',jerr,sigma_scale,cval)

*@ OUT_ORG <file name>
              elseif( cmd(1:3).eq. 'OUT'  .and. opt.eq.'A' ) then   ! out_org file
*                 Get name of eq_file: and read
                  call GetWord(line, out_org, indx)
* MOD TAH 131111: Allow ~ substitution
                  call subhome(out_org)
                  worg = .true.
   
*@ FIX_UNRE <Station.info/N> <antmod.dat>
              elseif( cmd(1:3).eq.'FIX' .and. opt.eq.'P' ) then
                  call GetWord(line, cval, indx)
*                 See if station.info included
                  call GetWord(line, stinf_file,  indx)                  
                  call GetWord(line, antmod_file, indx)  
                  fix_unre = .true.  
                  if( stinf_file(1:2).eq.'N ' .or.
     .                stinf_file(1:2).eq.'n ' )  fix_unre = .false.            

*@ REP_SCALE <scale file>
              elseif( cmd(1:3).eq.'REP' .and. opt.eq.'P' ) then
*                 Extract the file name
                  call GetWord(line, repscl_file, indx)
                  open(103,file=repscl_file,iostat=jerr,status='old')
                  call report_error('IOSTAT',jerr,'open',repscl_file)
                  if( jerr.ne.0 ) then
                     rep_scale = .false.
                     repscl_file = ' '
                  else
                     rep_scale = .true.
*                    Now read the file
                     num_repscl = 0
                     jerr = 0
                     do while ( jerr.eq.0 )
                        read(103,'(a)',iostat=jerr ) line
                        if( line(1:1).eq.' ' .and.jerr.eq.0 ) then
                           num_repscl = num_repscl + 1
                           if( num_repscl.gt.max_repscl )
     .                     call report_stat('FATAL','TSCON',
     .                       'read/repscl',repscl_file,
     .                        apr_file,'Too scale etimates',max_repscl)
                           read(line,*,iostat=kerr) repscl(:,num_repscl)
                           call report_error('IOSTAT',kerr,'decod',
     .                         line,1,repscl_file)
                        end if
                     end do 
                     close(103)
                  endif

* Command unknown
             else
* MOD TAH 180306: See if commnad OK (needed when 'P' and 'A' options
*               introdced.
                if( cmd(1:2).eq. 'EQ'    .or. cmd(1:3).eq.'APR' .or.
     .              cmd(1:6).eq.'STAB_S' .or. cmd(1:3).eq.'POS' .or.
     .              cmd(1:6).eq.'STAB_I' .or. cmd(1:6).eq.'STAB_M' .or.
     .              cmd(1:3).eq.'CND'    .or. cmd(1:3).eq.'DEB' .or. 
     .              cmd(1:3).eq.'TIM'  .  or. cmd(1:3).eq.'SIG' .or. 
     .              cmd(1:3).eq.'OUT'    .or. cmd(1:3).eq.'FIX' .or. 
     .              cmd(1:3).eq.'REP' ) then
                    indx = 1
                 else
                    write(*,600) cmd(1:trimlen(cmd)), 
     .                           line(1:trimlen(line))
 600                format('TSCON: Unknown command ',a,
     .                     ' from line ',a)  
                 endif           
             endif
         endif
      end do

****  Output the list of commands and stabalization sites
      if( opt.eq.'A' ) then
         write(*,720) comfile(1:trimlen(comfile))
 720     format('# TSCON: Command file ',a)
         write(*,740) stab_it, stab_nsig, use_ratio, cnd_parts_bits,  
     .                stab_min_dne, stab_min_dh 
 740     format('# TSCON: ',I3,' Iterations, Edit ',F4.1,' sigma',
     .          ' UseRatio ',F5.3,' ESTBITS ',o6,
     .          ' Min NE ',F7.4,' H ',F7.4,' m')
         if( fix_unre ) write(*,'("# TSCON: UNRE error removed")')
         if( rep_scale) write(*,'("# TSCON: Scale restored from ",a)') 
     .                  repscl_file   
         if( rep_scale) write(*,'("# TSCON: Scale restored ",i5,
     .                  " values")') num_repscl
         write(*,750) trim(out_org)
 750     format('# TSCON: Output Origin file ',a)

         write(*,760) num_site
 760     format('# TSCON: List of Primary stablization sites from ',i5,
     .          ' sites')

         do ih = 1, max_hierach 
            hierach_out = .false.
            j = 0
            nt = 0
            do i = 1,nin_stab
                if( stab_len(i,ih).gt.0 ) then
                   j = j + 1 
                   stab_out(j) = stab_site(i,ih)
                   if( j.eq.10 ) then
*                     If header not listed, list now for post
*                     primary sites.
                      if( ih.gt.1 .and. .not. hierach_out ) then
*                         write header out
                          write(*,770) ih
 770                      format('# TSCON: Level ',i2,
     .                           ' stablization sites')
                          hierach_out = .true.
                      endif

                      write(*,780) (nt+k,stab_out(k),k=1,10)
 780                  format(10(I5,1x,a8,1x))
                      j = 0
                      nt = nt + 10
                   endif
                endif
            end do
            if( j.gt.0 ) then
                if( ih.gt.1 .and. .not. hierach_out ) then
*                   write header out
                    write(*,770) ih
                end if
                write(*,780) (nt+k,stab_out(k),k=1,j)
            end if
         end do

*        Check for restricted use sites
         j = 0
         do i = 1, nin_stab
            if( stab_times(1,i).gt.0 ) then
*              This site has restrictions on stabilzation
               if( j.eq.0 ) then
                   write(*,820) 
 820               format('* R set stabialization sites',/,
     .                 '* Site    HFname Start MJD  End MJD')
                   j = 1
               end if
               write(*,840) stab_site(i,1), stab_hfs(i),
     .                   stab_times(1,i), stab_times(2,i)
 840           format('* ',a,2x,3x,a16,2x,F12.4,1x,F12.4)
            end if
         end do

         write(*,'(100a1)')  ('-',i=1,100)
      endif
       
****  Thats all
      return
      end

CTITLE READAPRF

      subroutine readaprf

      implicit none

* MOD TAH 141006: Updated so that EXTENDED PERIODIC terms will be ignored
*     Implemented by adding -PER to file name (case sensitive)
 
      include 'tsfit.h'
      include 'tscon.h'

      integer*4 ierr
     .,         jerr
     .,         indx, jndx 
     .,         is    ! Station number
     .,         i     ! Loop variable
     .,         trimlen
      logical no_per  ! No periodic terms to be retained
      integer*4 ind_no_per   ! Index in string of end or no-per terms


      real*8 pos_epoch  ! JD of position estimate.
     .,      vals(10)   ! Values read from line with multiread.

      character*256 line
      character*32 cpbuff
      character*8 name, cval



****  Routine to read apriori coordinate file
      write(*,120) apr_file(1:trimlen(apr_file))
 120  format('# TSCON: Reading coordinates from ',a)
* MOD TAH 141006: See if -PER added to apriori file name
      ind_no_per = index(apr_file,'-PER')
      if( ind_no_per.gt.0 ) then
           no_per = .true.
           ind_no_per = ind_no_per-1
      else
           ind_no_per = trimlen(apr_file)
           no_per = .false.
      end if
 
      open(102, file=apr_file(1:ind_no_per), iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',apr_file,0,
     .              'readaprf')

*     Read through file finding site and source names
 
      do while ( ierr.eq.0 )
 
          read(102,'(a)', iostat=ierr) line

* MOD TAH 080724: See if reference frame name is passed
          if( (line(1:16).eq.'+REFERENCE_FRAME' .or.
     .         line(1:16).eq.'+REFERENCE FRAME')  .and.
     .         ierr.eq.0                            ) then
              cpbuff = line(17:)
              call trimlead(cpbuff)
              reference_frame = cpbuff
          endif

          call casefold(line)

          if( ierr.ne.-1) then
              call report_error('IOSTAT',ierr,'read',apr_file,
     .                          0,'readaprf')
          end if
 
*         Decode if line is not a comment, and no error
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
*             See if site name
              indx = 1
              call getword(line, name, indx) 
              call casefold(name)
              jndx = 1
              call get_cmd(name, gsite_names, num_site, is,
     .                         jndx )

*             See if we need to add site name
              if( is.le.0 ) then
                  jndx = 0
                  do i = 1,num_site
                     if(name(1:4).eq.gsite_names(i)(1:4) ) 
     .                                     jndx = jndx + 1 
                  end do 
                  if( jndx.gt.0 ) then
                      num_site = num_site + 1
                      if( num_site.gt.max_site ) then 
                          call report_stat('FATAL','TSCON','readaprf',
     .                        apr_file,'Too many sites',max_site)
                      endif

                      gsite_names(num_site) = name
                      is = num_site
                  endif
              end if

*                                                 ! site name found
              if( is.gt.0 .and. is.ne.999999 ) then
                  call multiread(line,indx,'R8',jerr, vals, cval,
     .                            7)

* MOD TAH 980417: Changed to jerr on multiread so that file will continue
*            to be read.
                 if( jerr.eq.0 ) then
 
*                   Now assign values
                    do i = 1, 3
                       apr_val_site(i,1,is) = vals(i)
                       apr_val_site(i,2,is) = vals(3+i)
                    end do
 
*                   Convert epoch to JD
                    call decyrs_to_JD( vals(7), pos_epoch)
                    site_epoch(is) = pos_epoch-2 400 000.5d0 ! Make MJD
                 else 
                    write(*,200) gsite_names(is)
 200                format('** New coordinates of ',a8,
     .                 ' not used due to currupt apr_file line')
                 end if
*                                             ! site name found
             end if
 
* MOD TAH 991110: See if extended non-secular term
             indx = 1
             call GetWord(line, cval, indx)
             call casefold(cval) 
             if( cval(1:8).eq.'EXTENDED' ) then
*                OK, we have extended non-secular model entry
                  if( no_per ) then   ! Check to see if PERIODIC 
                          ! before decoding
                      if ( index(line,'PERIODIC').eq.0 ) then
                           call decode_nonsec(line, indx)
                      endif
                  else   ! Regular decoding of EXTENDED entry
                      call decode_nonsec(line, indx)
                  endif
             end if

*                         ! No file reading error
          end if
*                     ! Looping until EOF found
      end do
 
      close(102)

***** Thats all
      return
      end
  

CTITLE  STAB_SYS
  
      subroutine stab_sys(emjd)

      implicit none

****  Routine to perform translation, rotation and scale to align
*     current coordinates with a reference frame.

      include 'tsfit.h'
      include 'tscon.h'
      include '../includes/const_param.h'

      real*8 emjd  ! Evaluation MJD 

* LOCAL Variables
      integer*4 i,j,k, l,m
     .,         is, indx, ni
     .,         fin_num_stab
     .,         nc, ns          ! Number of parameters estimated
     .,         it, ait         ! Iteration counter
     .,         ipivot(max_np)      ! pivot elements
     .,         nr, ned        ! number of elements in chi**2 and 
                               ! number edited
     .,         lc, trimlen 

      integer*4 date(5) 

      integer*4 zone  ! Zone for UTM coordinates
     
      logical kbit 

      real*8 dist    ! Earthquake separation
     .,   min_nesig, min_htsig  ! Min sigmas
     .,   av_sig_ne, av_sig_ht  ! Median sigmas
     .,   dne_tol, dht_tol      ! NE and Ht tolerances
     .,   ss_nesig(max_site), ss_htsig(max_site) ! Sigmas 
     .,   dsol_cnd(max_np), atf(max_np), ptp(max_np,max_np)  ! Parameter estimator and NEW
     .,   scale(max_np) 
     .,   dobs(3), wgh_obs(3,3)
     .,   parts(3,max_np)     ! Partials 
     .,   prefit_sum, prefit_chi, postfit_chi, postfit_sum
     .,   unc_geod(3)   ! Colat, long and height (rads)
     .,   loc_coord(3), rot_mat(3,3), dneu(3) 
     .,   neu_rms(3), neu_sum(3), neu_chi(3), neu_wgh(3)    ! RMS NEU Postfit
     .,   sectag, decyr
     .,   cov_xyz(3,3), cov_neu(3,3)    ! Covariance XYZ for rotation/trans sigmas
     .,   tcov(max_np,max_np)       ! Temp matrix for var_comp

      character*1 hemi    ! Hemisphere for UTM coordinates

      character*8 name_aprf

      character*16 code   ! HFcode: Needs to match output directory
      logical found       ! Set true when stab site found



****  First find all the data that has been collected at this epoch
      num_used = 0
      num_stab = 0

*     Get list of sites
      do i = 1, num_ent
         if( abs(emjd-in_mjd(i)).lt.1.d-3 ) then
*            OK, this epoch is OK
             num_used = num_used + 1
             used_index(num_used) = i
         end if
      end do

****  Get the stab list from hierachical lists
      do j = 1, nin_stab
         found = .false.
         do k = 1, max_hierach
            if( stab_len(j,k).gt.0 ) then 
                do i = 1, num_used
                   ns = in_ns(used_index(i))
                   if( stab_site(j,k)(1:stab_len(j,k)).eq.
     .                 gsite_names(ns)(1:stab_len(j,k)) ) then
*                      site is stabilization site; see if Time Range OK
                       if( stab_times(1,j).gt.0 .and.
     .                    emjd.ge.stab_times(1,j).and. 
     .                    emjd.lt.stab_times(2,j) ) then
                          num_stab = num_stab + 1
                          stab_index(num_stab) = used_index(i)
                       else
                          num_stab = num_stab + 1
                          stab_index(num_stab) = used_index(i)
                       endif
                       found = .true.
                   endif
                   if( found ) exit
                end do
             end if
             if( found ) exit
          end do
      end do


      if( debug ) then
         write(*,'(a)') '----------------------------------------------'
         write(*,20) emjd, num_used, num_stab 
  20     format('EMJD ',F9.3,' Sites Used ',i4,' Frame ',I4)
         write(*,30) (gsite_names(in_ns(stab_index(k))),k=1,num_stab)
  30     format((10(a8,1x):,1x))
      end if
         
****  See if we have data
      if( num_used.eq.0 ) RETURN
              
***** OK, we have data
*     Now generate reference frame coordinates at this time
      do i = 1, num_stab
*        Get base name
         name_aprf = gsite_names(in_ns(stab_index(i)))
*        Now see if any renames are needed
         do j = 1,num_rn
            if( rn_codes(1,j)(1:4).eq.name_aprf(1:4) ) then
*              Rename code matches, see if time range OK
               indx = 0
               if( emjd.ge.rn_times(1,j) .and.
     .             emjd.le.rn_times(2,j) ) then
*                  Check hfile code 
                   indx = 1
                   if ( trimlen(rn_hfiles(j)).gt.0 ) then
                       code = rn_hfiles(j)
                       lc = trimlen(code)
                       indx = index(tsdir,code(1:lc))
                   endif
               end if
               if( indx.gt.0 ) then
*                  OK: Update the name
                   name_aprf = rn_codes(2,j)
*                  See if edited
                   if( name_aprf(5:8).eq.'_XCL' .or.
     .                 name_aprf(5:8).eq.'_XPS' ) then
*                      Remove this site from list
                       if( debug ) 
     .                 write(*,40) gsite_names(in_ns(stab_index(i))),
     .                     name_aprf
  40                   format('Deleting ',a,' Name ',a)
                       stab_index(i) = -1
                   end if
               endif
            endif
         end do
****     OK, now see if EQ affected if site is still used
         if( stab_index(i).gt.0 ) then
             do j = 1,num_eq
                 call eval_dist( in_xyz(1,stab_index(i)), 
     .                eq_pos(1,j), dist )
                 if( dist.le.eq_rad(j) ) then
*                    See if time is OK
                     if( emjd.ge.eq_epoch(j) ) then
*                       Site affected save, new name
                        name_aprf(7:8) = eq_codes(j)(1:2)
                     endif
                 endif
             end do
****         Now get coordinates
             indx = 1
             call get_cmd(name_aprf, gsite_names, num_site, is, indx)
             if( is.gt.0 ) then 
*                Match to name so get coordinates
                 call get_coord(is, emjd, ref_coords(1,i))
                 if( ref_coords(1,i)**2+ref_coords(2,i)**2+
     .               ref_coords(3,i)**2.lt.1.d10 ) is = -1
                 if( debug .and. is.eq.-1 ) then
                     write(*,60) name_aprf
  60                 format('Deleting ',a,' No coordinates')
                 end if

             endif
             if( is.le.0 ) then 
                 stab_index(i) = -1
                 if( debug ) 
     .           write(*,80) name_aprf
  80             format('Deleting ',a,' no name match')
             end if
         end if
      end do

****  Now do the fit to the reference sites.  First remove
*     any sites that do not match the sigma limits
      j = 0
      min_nesig = 1.d20
      min_htsig = 1.d20
      do i = 1,num_stab
*         Get the horizontal and vertical sigmas
          if( stab_index(i).gt.0 ) then
              j = j + 1
              ss_nesig(j) = sqrt(in_neu_std(1,stab_index(i))**2+
     .                           in_neu_std(2,stab_index(i))**2)
              ss_htsig(j) = in_neu_std(3,stab_index(i))
              if( min_nesig.gt.ss_nesig(j) ) min_nesig=ss_nesig(j)
              if( min_htsig.gt.ss_htsig(j) ) min_htsig=ss_htsig(j)
          endif
      end do
*     Get Medians and then edit
      call mdian2( ss_nesig, j, av_sig_ne )
      call mdian2( ss_htsig, j, av_sig_ht )

*     Now do the editing
      dht_tol = max(av_sig_ht-min_htsig, stab_min_dh )
      dne_tol = max(av_sig_ne-min_nesig, stab_min_dne)

      fin_num_stab = 0
*     Only apply editing when there are a resonable number of sites
      if( num_stab.gt.6 ) then
         do i = 1, num_stab
*           Check horizontal quality
            if ( stab_index(i).gt.0 ) then
               if( sqrt(in_neu_std(1,stab_index(i))**2+
     .                  in_neu_std(2,stab_index(i))**2)-min_nesig.gt.
     .             use_ratio*dne_tol ) then
                  stab_index(i) = -1
               end if
            endif
*           Check height quality
            if ( stab_index(i).gt.0 ) then
               if( in_neu_std(2,stab_index(i))-min_htsig.gt.
     .             use_ratio*dht_tol ) then
                  stab_index(i) = -1
               end if
            endif
            if( stab_index(i).gt.0 ) fin_num_stab = fin_num_stab + 1
         end do
      else
         fin_num_stab = num_stab
      endif

****  Make sure we still have some sites left
      if( fin_num_stab .lt. 3 ) then
          write(*,120) fin_num_stab, num_stab, emjd
 120      format('Only ',i3,' stab sites from ',i3,' at MJD ',F8.2)

      endif

****  Now start the iteration.
      nc = 0   ! Count number of parameters
      do i = 1,max_np
         if( kbit(cnd_parts_bits,i) ) nc = nc + 1
      end do

      it = 0
      do while ( it.lt.stab_it )
         it = it + 1

****     Clear the normal equations
         do i = 1,nc
             dsol_cnd(i) = 0.d0
             do j = 1,nc
                 ptp(i,j) = 0.d0  
             end do
             ! ptp(i,i) = 1.d-2  ! Put weak constraint: 10 meters
         end do

****     Now loop over sites
         prefit_sum = 0.d0
         nr = 0   ! Number of data in chi**2 estimates
         do i = 1, num_stab
            if( stab_index(i).gt.0 ) then
*               Form data weight matrix based on sigma and
*               height varinace
                call gen_obscov(stab_index(i), wgh_obs)
*               Form partials
                call gen_part(stab_index(i),cnd_parts_bits,parts)
*               Form OMC
                do j = 1,3
                   dobs(j) = in_xyz(j,stab_index(i))-ref_coords(j,i)
                end do
                if( debug )
     .          write(*,220) i,gsite_names(in_ns(stab_index(i))),dobs,
     .                 (ref_coords(j,i), j=1,3)
 220            format('PREFIT ',i3,1x,a8,3(1x,F10.4),2x,3(F14.4,1x))

***             Make sure there are no gross errors
                if( sqrt(dobs(1)**2+dobs(2)**2+dobs(3)**2).lt.1.0 ) then

*                  Increment normal equations
                   do j = 1, nc
                      do l = 1,3
                         do m = 1,3
                            dsol_cnd(j) = dsol_cnd(j) +
     .                             parts(m,j)*wgh_obs(m,l)*dobs(l)
                        end do
                      end do

                      do k = 1, nc
                         do l = 1,3
                            do m = 1,3
                                ptp(j,k) = ptp(j,k) +
     .                               parts(m,j)*wgh_obs(m,l)*parts(l,k)
                            end do
                         end do
                      end do
                   end do
                   do l = 1,3
                       nr = nr + 1
                       do m = 1,3
                           prefit_sum = prefit_sum + dobs(l)*
     .                               wgh_obs(l,m)*dobs(m)
                       end do
                   end do
                else     ! Gross error in pre-fit so remove
                   stab_index(i) = -1
                endif
            end if
         end do
****     Now solve the system
         do i = 1, nc
             atf(i) = dsol_cnd(i)
         end do

         call invert_vis(ptp, dsol_cnd,scale, ipivot, nc,max_np,1)
         postfit_sum = 0.d0
         do i = 1,nc
             postfit_sum = postfit_sum + atf(i)*dsol_cnd(i)
         end do

         if( nr.gt.9 ) then 
            prefit_chi = prefit_sum/nr
            postfit_chi = (prefit_sum-postfit_sum)/nr
         else
            prefit_chi = 1.0
            postfit_chi = 1.0
         endif

         if( debug ) then
            write(*,230) nc, dsol_cnd(1:nc)
 230        format('SOLN: ',I2,' NP Ests ',10(1x,F8.4))
         end if
         
****     Compute the postfit residuals and see if we need to edit
         ned = 0
         nr = 0
         do j = 1,3
            neu_sum(j) = 0
            neu_chi(j) = 0
            neu_wgh(j) = 0
         end do
            
         do i = 1, num_stab
            if( stab_index(i).gt.0 ) then
*               Form data weight matrix based on sigma and
*               height varinace
                call gen_obscov(stab_index(i), wgh_obs)
*               Form partials
                call gen_part(stab_index(i),cnd_parts_bits,parts)
*               Form OMC and add correction for estimates
                do j = 1,3
                   dobs(j) = in_xyz(j,stab_index(i))-ref_coords(j,i)
                   do k = 1,nc
                      dobs(j) = dobs(j) - parts(j,k)*dsol_cnd(k)
                   end do
                end do
                call rotate_geod(dobs, dneu, 'XYZ', 'NEU', 
     .                in_xyz(1,stab_index(i)), loc_coord, rot_mat)

****            Accumulate statistics on NEU residuals
                do j = 1,3
                   neu_sum(j) = neu_sum(j) + 
     .                   dneu(j)/in_neu_std(j,stab_index(i))**2
                   neu_chi(j) = neu_chi(j) +
     .                   dneu(j)**2/in_neu_std(j,stab_index(i))**2
                   neu_wgh(j) = neu_wgh(j) +
     .                   1.d0/in_neu_std(j,stab_index(i))**2
                enddo
                nr = nr + 1
                if( debug )
     .          write(*,240) i,gsite_names(in_ns(stab_index(i))),dobs,
     .              dneu
 240            format('POSFIT ',i3,1x,a8,3(1x,F10.4),2x,3(1x,F10.4))
*               Check the chi**2 value
                prefit_sum = 0.0
                do l = 1,3
                   do m = 1,3
                      prefit_sum = prefit_sum + dobs(l)*
     .                             wgh_obs(l,m)*dobs(m)
                   end do
                end do
****            See if value is too large
                if( prefit_sum/postfit_chi.gt.3*stab_nsig**2 ) then
*                   Outlier so remove
                    stab_index(i) = -1
                    ned = ned + 1
                endif
            end if
         end do
****     Finish statistics
         do j = 1,3
            if( nr.ge.3 ) then
               neu_sum(j) = neu_sum(j)/neu_wgh(j)*1000.d0
               neu_rms(j) = sqrt(neu_chi(j)/neu_wgh(j))*1000.d0
               neu_chi(j) = sqrt(neu_chi(j)/nr)
            else
               neu_sum(j) = 0.0
               neu_rms(j) = 100.0
               neu_chi(j) = 1.0
            endif 
         end do

****     See if we edited anything
         ait = it  ! Save actual iterations
         if ( ned.eq.0 ) then
             it = stab_it + 1   ! Force exit
         else
             if( debug ) write(*,260) ait, ned
 260         format('Iteration ',i3,' Deleted ',i3,' sites')
         endif
      end do

****  Output summary of results
      call jd_to_decyrs(emjd, decyr)
      call mjd_to_ymdhms(emjd, date, sectag)

* MOD TAH 150707: Added scale estimate when scale estimated.  Parameter
*     number set in gen_part
      if( nps.eq.0  ) then 
         write(*,320) decyr, (date(j),j=1,3), nr, neu_rms,neu_chi, 
     .          sqrt(prefit_chi), sqrt(postfit_chi), 
     .          neu_sum(3),neu_rms(3)/sqrt(nr*1.0), ait
 320     format('FIT ',F9.4,1x,I4,1x,I2.2,1x,I2.2,' For ',i4,
     .          ' RefSites; WRMS NEU ',
     .          3(F6.2,1x),' mm, NRMS NEU ',3(F6.2,1x),
     .          'Pre/Post Chi ',2F6.2,' Mean HGT ',F8.2,' +- ',
     .           F6.2,' mm, Iter ',i3)
      else
         write(*,325) decyr, (date(j),j=1,3), nr, neu_rms,neu_chi, 
     .          sqrt(prefit_chi), sqrt(postfit_chi), 
     .          dsol_cnd(nps)*6371.e-3, 
     .          sqrt(ptp(nps,nps)*postfit_chi)*6371.e-3, ait
 325     format('FIT ',F9.4,1x,I4,1x,I2.2,1x,I2.2,' For ',i4,
     .          ' RefSites; WRMS NEU ',
     .          3(F6.2,1x),' mm, NRMS NEU ',3(F6.2,1x),
     .          'Pre/Post Chi ',2F6.2,' Scale dH ',F8.2,' +- ',
     .           F6.2,' mm, Iter ',i3)
      endif

* MOD TAH 170720: Formally output the estimates if requested
      if( worg ) then
          call write_org_vals( 202, nc, nr, dsol_cnd, ptp, 
     .                emjd, units_scale, postfit_chi, max_np ) 
      else
*         Ned to call write_org_head to get units_scale set
*         so don't make this call.
C         call write_org_vals( 6, nc, nr, dsol_cnd, ptp, 
C    .                emjd, units_scale, postfit_chi, max_np ) 

      endif 
 
****  Estimation now finished.  Now apply corrections
      do i = 1, num_used
         if( used_index(i).gt.0 ) then

****        Generate partials
            ni = used_index(i) 
            is = in_ns(used_index(i))
            call gen_part(used_index(i),cnd_parts_bits,parts)

****        Compute adjustment to position
            do j = 1,3
               dobs(j) = 0.d0
               do k = 1,nc
                  dobs(j) = dobs(j) + parts(j,k)*dsol_cnd(k)
               end do
            end do

****        Compute noise in adjustment
            call var_comp(parts,ptp, cov_xyz, tcov, 3,max_np,1)
            do j = 1,3
               do k = 1,3
                  cov_xyz(j,k) = cov_xyz(j,k)*postfit_chi
               end do
            end do
****        Now get the NEU covariance
            call XYZ_to_GEOD(rot_mat, in_xyz(1,ni), loc_coord )
            call var_comp(rot_mat,cov_xyz, cov_neu, tcov, 3,3,1)

****        OK; update the sigmas
            do j = 1,3
                in_xyz_std(j,ni) = sqrt(in_xyz_std(j,ni)**2+
     .                                  cov_xyz(j,j))
                in_neu_std(j,ni) = sqrt(in_neu_std(j,ni)**2+
     .                                  cov_neu(j,j))
            enddo


****        Now apply to save values
            do j = 1,3
               in_xyz(j,ni) = in_xyz(j,ni)- dobs(j)
            enddo

****        Now convert XYZ to lat/long/height
            call geod_to_geod(in_xyz(1,ni), unc_geod, 
     .            'XYZ', 'GEOD','WGS84','WGS84',zone,hemi)

*           Convert co-lat, long in radians to degrees
            in_llu(1,ni) = (pi/2-unc_geod(1))*180/pi
            in_llu(2,ni) = unc_geod(2)*180/pi 
            in_llu(3,ni) = unc_geod(3)
*           Now get the dNEU
            call loc_to_geod(unc_geod, in_neu(1,ni))
         endif
      end do


****  Thats all
      return
      end

CTITLE GET_OBSCOV

      subroutine gen_obscov(ni, wgh_obs)

      implicit none

      include 'tssum.h'
      include 'tscon.h'

      integer*4 ni

      real*8 wgh_obs(3,3) 

* LOCAL
      integer*4 i,j
     .,      ipivot(3)

      real*8 cov_neu(3,3), tcov(3,3), rot_mat(3,3)
     .,      loc_coord(3), scale(3)

****  Routine to generate obs wgh matrix which is a blend
*     of data noise, constant and height downweight
      do i = 1,3
         do j = 1,3
            cov_neu(i,j) = 0
         end do
      end do
      cov_neu(1,1) = in_neu_std(1,ni)**2
      cov_neu(2,2) = in_neu_std(2,ni)**2
      cov_neu(3,3) = in_neu_std(3,ni)**2*cnd_hgt_var

****  Now convert to XYZ covariance
      call XYZ_to_GEOD(rot_mat, in_xyz(1,ni), loc_coord )
c.... Now transpose rot_matrix to that NEU to XYZ direction 
      do i = 1,2
        call dvswp(rot_mat(i,i+1),3,rot_mat(i+1,i),1, 3-i)
      end do
      call var_comp(rot_mat,cov_neu, wgh_obs, tcov, 3,3,1)
*     Now invert to get weight matrix
      call invert_vis(wgh_obs, cov_neu,scale, ipivot, 3,3,0)

****  Thats all
      return
      end


CTITLE GEN_PART

      subroutine gen_part(ni, bits, parts)

      implicit none

      include 'tssum.h'
      include 'tscon.h'

      integer*4 ni    ! Entry number of obs
     .,    bits       ! Bits set based on what is to be estimated
                      ! Bits 1-3  -- rotations, bits 4-6 translations
                      ! Bit 7     -- scale
                      ! Bits 8-10 -- CM offsets due to degree 1 load

      real*8 parts(3,max_np)  ! XYZ partials for rotation/translation and
                      ! scale
      real*8 cmpart(3,3)      ! Parts of dXYZ in CF frame to dXYZ of CM.
                      ! Computed in subroutine CM_parts

* LOCAL 

      integer*4 i,j 
     .,    np   ! Number of parameters

      real*8 ut1  ! UT1 value (set to zero at the moment)

      logical kbit

****  Get the PMU partials
      ut1 = 0.d0
      np = 0
      nps = 0  ! Parameter number for scale.
      do i = 1,3
         do j = 1,max_np
             parts(i,j) = 0.0
         end do
      end do

*     Rotations 
      do i = 1,3
          if( kbit(bits,i) ) then 
              np = np + 1
              if( i.eq.1 ) then
                 parts(1,np) = -in_xyz(3,ni)*1.d-9
                 parts(3,np) =  in_xyz(1,ni)*1.d-9
              elseif( i.eq.2 ) then
                 parts(2,np) =  in_xyz(3,ni)*1.d-9
                 parts(3,np) = -in_xyz(2,ni)*1.d-9
              elseif( i.eq.3 ) then
                 parts(1,np) =  in_xyz(2,ni)*1.d-9
                 parts(2,np) = -in_xyz(1,ni)*1.d-9
              end if
          endif
      end do
*     Translations
      do i = 4,6
         if( kbit(bits,i) ) then
             np = np + 1
             parts(i-3,np) = 1.d0
         endif
      end do
*     Scale
      if( kbit(bits,7) ) then
          np = np + 1
*         Save the scale parameter number for use in the
*         output
          nps = np
          do i = 1,3
             parts(i,np) = in_xyz(i,ni)*1.d-9
          end do
      end if

* MOD TAH 170720: Compute the CM to CF partials if needed
      if( CMEst ) then
          call CM_parts( cmpart, in_xyz(1,ni) )
*         Loop pver XYZ CM
          do i = 1, 3   ! XYZ of CM 
             if( kbit(bits,7+i) ) then  
                np = np + 1
                do j = 1,3    ! XYZ of site coordinates
                   parts(j,np) = cmpart(j,i)   ! column i (CM xyz)
                enddo
             endif
          end do 
       endif


***** Thats all
      return
      end

CTITLE GET_COORD

      subroutine get_coord(is, emjd, coords)

      implicit none 

*     Routine to return coordinates of site is at time EMJD

      include 'tssum.h'
      include 'tscon.h'

      integer*4 is
      real*8 emjd   ! MJD
     .,      coords(3)   ! Returned XYZ coordimayes

      integer*4 i
      real*8 dsol(3)    ! XYZ adjustments for non-secular terms

****  Get the linear corrinates
      do i = 1,3
         coords(i) = apr_val_site(i,1,is) + 
     .               apr_val_site(i,2,is)*(emjd-site_epoch(is))/365.25d0
      end do

****  Now evaluate the non-secular terms (0 sets that logs should be 
*     computed)
      call eval_nonsec(is, emjd, num_nonsec, param_nonsec,
     .                 apr_val_nonsec, dsol, 0 )

      do i = 1,3
         coords(i) = coords(i) + dsol(i)
      end do

***** Thats all 
      return
      end
  
CTITLE DECODE_NONSEC

      subroutine decode_nonsec( line, indx)

      implicit none

*     Routine to decode the non-secular entries read from the globk
*     apriori files.  The entries must be unique and any replicated
*     one are replaced by the latest ones read.

      include 'tsfit.h'
      include 'tscon.h'

* PASSED VARIABLES
* indx    -- Position in line

      integer*4 indx
 
* line  -- line read from apriori file

      character*(*) line

* LOCAL VARIABLES
* is  -- Site number
* it  -- Non-secular term type (1-4 for offset, periodic, exponential and
*        logarithmic)
* i, j -- Loop counters
* jerr -- Error reading line
* date(5) -- Date read from line 
* jndx    -- Pointer in string

      integer*4 is, it, i,j, jerr, date(5), cndx

* vals -- Dummy value for read_line
* sectag -- Seconds tag for dates
* vals_read(6) -- Upto to 6 arguments read from line 

      real*8 vals, sectag, vals_read(6)

* done   -- Logical set true when finished reading line
* duplicate -- Logical set true if there are duplicate entries
*     in the nonsecular terms
* different -- Set true if new extended values differ from old 
*     values

      logical done, duplicate, different

* type -- Character string with type read from line
* cval -- Dummy character value
* fatal -- Fatal message if too many terms

      character*8 type, cval, name
      character*64 fatal

c      external globk_cmd_bd
    

****  OK, Start decoding.  See if we can find station name
      call getword(line, name, indx)
      cndx = 1 
      call get_cmd(name, gsite_names, num_site, is, cndx )
      cndx = index(line,'!')
      if( cndx.eq.0 ) cndx = index(line,'#')
      if( cndx.ne.0 ) line(cndx:) = ' '

      if( is.gt.0 .and. is.ne.999 ) then
*         OK, found site name.  Now start working through the line.
*         There can be multiple entries for each site on a line.
          done = .false.
          do while ( .not.done )
             call read_line(line, indx, 'CH',jerr,vals, type)
             if( jerr.eq.0 ) then
*                OK, see if can match type
                 call casefold(type)
                 it = 0
                 if( type(1:2).eq.'OF') it = 1
                 if( type(1:2).eq.'PE') it = 2
                 if( type(1:2).eq.'EX') it = 3
                 if( type(1:2).eq.'LO') it = 4
                 if( it.gt.0 ) then

*                    Found another non-secular entry.  Increment the 
*                    number and make not too many
                     num_nonsec = num_nonsec + 1
                     if( num_nonsec.gt. max_nonsec) then
                         write(fatal,120) max_nonsec
 120                     format('Too many nonsecular terms. ',
     .                           'Max allowed is ',i6)
                         call report_stat('FATAL','GLOBK',
     .                           'decode_nonsec',' ',fatal,0)
                     end if

*                    Save the site and type for this term
                     param_nonsec(1,num_nonsec) = is
                     param_nonsec(2,num_nonsec) = it

*                    Read the date from line.
                     call multiread(line,indx,'I4',jerr,date,cval,5)
                     sectag = 0.d0
                     call ymdhms_to_mjd(date,sectag,
     .                      apr_val_nonsec(1,num_nonsec))
*                    Get the "parameter" either period or decay
*                    time.  For offset and rate value should be
*                    set to zero.
                     call read_line(line,indx,'R8', jerr, vals, cval)
                     apr_val_nonsec(2,num_nonsec) = vals
*                    Now get the  parameters for the term.
*                    The number of terms depends on the type.
                     if( it.eq. 1 .or. it.eq.2 ) then
*                        Offset and Rate change and periodic
*                        terms.   For offset and rate, the
*                        first argument is ignored, for 
*                        periodic it is period.
                         call multiread(line,indx,'R8', jerr,
     .                        vals_read,cval,6)
*                        Now convert from NEU to XYZ
                         call nonsec_convert('TOXYZ',2,vals_read,
     .                        apr_val_nonsec(3,num_nonsec),
     .                        apr_val_site(1,1,is))
                     else if( it.eq.3 .or. it.eq.4 ) then
*                        Exponential and logarithm. In this case
*                        we have just the three NE and U amplitudes
                         call multiread(line,indx,'R8', jerr,
     .                        vals_read,cval,3)
*                        Clear the other terms in values read
                         do i = 4,6
                            vals_read(i) = 0.d0
                         enddo
*                        Now convert from NEU to XYZ
                         call nonsec_convert('TOXYZ',1,vals_read,
     .                        apr_val_nonsec(3,num_nonsec),
     .                        apr_val_site(1,1,is))
                     endif

*                    Now check if the entry just added is unique
                     do i = 1, num_nonsec-1
                        duplicate = .true.
                        do j = 1, 2
                           if( param_nonsec(j,i).ne.
     .                         param_nonsec(j,num_nonsec) ) 
     .                                        duplicate = .false.
                           if( apr_val_nonsec(j,i).ne.
     .                         apr_val_nonsec(j,num_nonsec) )
     .                                        duplicate = .false.
                        end do
                        if( duplicate ) then
*                          Print a message if the new values do not
*                          match the old values
                           different = .false.
                           do j = 3,8
                              if( apr_val_nonsec(j,i).ne.
     .                            apr_val_nonsec(j,num_nonsec) )
     .                            different = .true.
                           enddo
                           
                           if( different )
     .                     write(*,180) i,gsite_names(is),type
 180                       format('**NOTE** Non-secular term ',i5,
     .                            ' at site ',a8,' Type ',a8,
     .                            ' being replaced')
                           do j = 3,8
                              apr_val_nonsec(j,i) = 
     .                            apr_val_nonsec(j,num_nonsec)
                           end do
                           num_nonsec = num_nonsec - 1
                        end if
                     end do
            
                end if
             else
*               Nothing is left on line so we are done
                done = .true.
             end if
          end do
      end if

****  Thats all
      return
      end

CTITLE NONSEC_CONVERT

      subroutine nonsec_convert(direct, num, vals_in, vals_out, xyz) 

      implicit none

*     Routine to convert the non-secular position components from
*     NEU to XYZ or visa versa.  

* PASSED VARIABLES
* num  -- Lead dimenstion of the values (ie. if num=2, then the vals
*      are paired (e.g. cos and sin terms in NE and U; if 1 then a
*      column vector of NEU)

      integer*4 num

* direct -- String containing direction of conversion (TOXYZ or
*     TONEU)

      character*(*) direct

* vals_in(num,3) -- Input values to be converted
* vals_out(num,3) -- Output values to be returned
* xyz(3)  -- Cartesian coordinates of site

      real*8 vals_in(num,3), vals_out(num,3), xyz(3)

* LOCAL VARIABLES
* i,j  -- Loop variable

      integer*4 i,j

* loc_in(3), loc_out(3) -- In and out values of terms
* loc_coord(3) -- Geodetic latitude, longtiude and height
* rot_mat(3,3) -- Rotation matrix from XYZ to NEU

      real*8 loc_in(3), loc_out(3), loc_coord(3), rot_mat(3,3)


****  OK, see which direction we are going
      if( direct(1:5).eq.'TOXYZ' ) then

         do i = 1, num
            do j = 1, 3
               loc_in(j) = vals_in(i,j)
            end do

*           Convert the vector over to the other frame 
            call rotate_geod(loc_in, loc_out, 'NEU', 'XYZ', xyz, 
     .                       loc_coord, rot_mat)

*           Save the result
            do j = 1,3
               vals_out(i,j) = loc_out(j)
            end do
         end do  
      else if( direct(1:5).eq.'TONEU' ) then

         do i = 1, num
            do j = 1, 3
               loc_in(j) = vals_in(i,j)
            end do

*           Convert the vector over to the other frame 
            call rotate_geod(loc_in, loc_out, 'XYZ', 'NEU', xyz, 
     .                       loc_coord, rot_mat)
*           Save the result
            do j = 1,3
               vals_out(i,j) = loc_out(j)
            end do
         end do  
      else 
*        Error in call
         call report_stat('FATAL','GLOBK','nonsec_convert',
     .                    direct, 'Bad direction argument',0)
      end if

****  Thats all
      return
      end

CTITLE READSTAB

      subroutine readstab(line, indx)

      implicit none

*     Routine to read list of stabalization sites 
*     - before name removes from list

      include 'tsfit.h'
      include 'tscon.h'

      integer*4 indx   ! position in string

      character*(*) line

      integer*4 jerr, kerr
     .,         is    ! Site number
     .,         cndx  ! comment index
     .,         trimlen     ! Length of string.
     .,         ls    ! Length of site name
     .,         date(5)  ! Calender date for restricted times
     .,         in, jn, ih, ps   ! Indexes for working through lists of 
                         ! sites separated by /
      real*8  sectag  ! Seconds tag.

      character*90  sitelist ! Hierarchical list of sites names 
      character*8 site
      character*8 cdum

      character*16 rhfn  ! restricted names on stab_site with R option.
      logical remove

      integer*4 nameout

      data nameout / 0 /


      cndx = index(line,'!')
      if( cndx.eq.0 ) cndx = index(line,'#')
      if( cndx.ne.0 ) line(cndx:) = ' '

****  Work through the line
      jerr = 0
      if( nameout.eq.0 ) then 
         write(*,120) num_site, (is,gsite_names(is),is=1,num_site)
 120     format('Site Names for ',i5,' sites',/,
     .          1000(10(i5,1x,a,1x),/))
         nameout = 1
      endif

      do while ( jerr.eq.0 )
         cndx = indx
         call getword(line, sitelist, cndx)   ! Get long list is passed
         call casefold(sitelist) 
         if( trimlen(sitelist).eq.0 ) then
             jerr = -1
         elseif ( sitelist.eq.'R   ' ) then
*            A restricted application is being applied.
*            See if limit on product name
             cndx = indx
             call getword(line,rhfn,cndx)
*            See if this a number
             call check_num(rhfn, kerr)
             if( kerr.ne.0 ) then ! this is sting
                stab_hfs(is) = rhfn
                indx = cndx
             else
                stab_hfs(is) = ' '
             end if
*            Get start date
             call multiread(line,indx,'I4',kerr,date,cdum,5)
             sectag = 0.0
             call ymdhms_to_mjd(date,sectag,stab_times(1,is))
             call read_line(line,indx,'I4',kerr,date(1),cdum)
             if( kerr.eq.0 ) then
                call multiread(line,indx,'I4',kerr,date(2),cdum,4)
                call ymdhms_to_mjd(date,sectag,stab_times(2,is))
             else
                stab_times(2,is) = 88069.0d0
             endif
         else

* MOD TAH 130331: Start working through the list to see what sites
*            we have.
             remove = .false.
             in = 1
             if( sitelist(1:1).eq.'-' ) then
                 remove = .true.
                 in = 2
             elseif ( sitelist(1:1).eq.'+' ) then
                 in = 2
             endif 
*            Start looping through names
             nin_stab = nin_stab + 1
             if( nin_stab.gt.max_stab ) then
                 call report_error('TOMANY site','add',line,1,max_stab)
             endif
             
             ih = 0
             jerr = 0
             do while ( jerr.eq.0 )
                ih = ih + 1
                jn = index(sitelist(in:),'/')+in-2  ! Character before /
                if( jn.eq.in-2 ) then   !  / not found so end of string
                    jn = trimlen(sitelist)
                    if ( in.ge.jn ) then    ! No name left
                        jerr = -1
                    endif
                endif
                if( jerr.eq.0 .and. sitelist(in:jn).ne.'CLEAR' ) then   ! We have a name to process
                    site = sitelist(in:jn)
                    ls = trimlen(site)
*                   See if we can find the site already, if we can
*                   don't add or remove if - used (remove): Just add
*                   it , no removal


*                   Save the site names and length
                    stab_site(nin_stab,ih) = site
                    stab_len(nin_stab,ih) = ls

C                   do is = 1, num_site
C                      if( site(1:ls).eq.gsite_names(is)(1:ls) ) then
C                         if( .not.remove ) then
C                            stab_site(is,ih) = 1
C                         else
C                            stab_site(is,ih) = 0
C                         endif
C                      endif
C                   end do
                endif
                in = jn+2   ! Start next name (skip over /)

             enddo
             is = nin_stab    ! For saving -R values

          endif    ! Stab_site command
       end do      ! Looping contents of line

***** That all
      return
      end


CITTLE READ_IN_CSV

      subroutine read_in_csv

      implicit none

*     Routine to read SCEC and PBO data files.  These file either have 1 or 
*     9 header lines.  The difference is based on the first header line.

      include 'tsfit.h'
      include 'tscon.h'
      include '../includes/const_param.h'


      character*16 gsite_full

      integer*4 date(5), ierr, jerr, i, j, ns, ne
      integer*4 trimlen, indx
      real*8  gmjd, pos_xyz_fin(3),xyz_std(6), unc_geod(3),
     .        pos_neu_fin(3), neu_std(6), sec

      real*8 ucov_neu(3,3), ucov_xyz(3,3) ! "Unit" NEU convarinace matrix
                   ! transformed to XYZ to scale sigma
      real*8 tcov(3,3)  ! Scratch space.

      real*8 unc_llu(3)  ! Computed Lat/long/height
      real*8 loc_coord(3), int_coord(3)  ! Initial coordinates from header. 
      real*8 rot_mat(3,3) ! Rotation from XYZ to NEU (and transpose)
      real*8 con(2)    ! Scale to convert lat, long from mm to rads
      real*8 dNEU(3)   ! NEU adjustments (mm)

      logical pbo_form  ! true if pbo version with sigmsa
     .,       neuset    ! Set true once rotation from NEU to XYZ known
      
 
      integer*4 zone  ! Zone for UTM coordinates
     
      character*1 hemi    ! Hemisphere for UTM coordinates

      character*512 line
      character*8 cdum
      character*10 datestr  ! String with data 2007-01-01

      data ucov_neu / 1.d0, 0.0d0,  0.d0,
     .                0.d0, 1.0d0,  0.d0,
     .                0.d0, 0.0d0, 10.d0 /

***   First open the file 
      open(100,file=in_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_file,0,'tscon')

      if( ierr.ne.0 ) RETURN

****  Read first line and see what we have
      read(100,'(a)',iostat=ierr) line
      call report_error('IOSTAT',ierr,'read',line,0,'tscon')
      if( ierr.ne.0 ) RETURN

****  See what "flavor" of file this is
      indx = index(line,'PBO Station')
      if( indx.eq.0 ) then   ! Assume SCEC format
          pbo_form = .false.
          call sub_char(line,',',' ')
          indx = 0
          call getword(line,cdum,indx)
          call getword(line,site_name,indx)
          call casefold(site_name)
          site_name(5:8) = '_GPS'
          gsite_full = site_name(1:4) // '_' // prod_id(1:9)
*         Now get the geodetic lat long and height from line
          read(line(indx:),*,iostat=jerr) int_coord
          call report_error('IOSTAT',ierr,'read',line,0,'SCEC LLU')
*         Convert to radians as needed (lat->co-lat as well)
          int_coord(1) = pi/2 - int_coord(1)*pi/180
          int_coord(2) = int_coord(2)*pi/180
          if( int_coord(2).lt.0 ) int_coord(2)=int_coord(2)+2*pi
      else
*         Read the PBO format
          pbo_form = .true.
          indx = index(line,'Frame :')
          reference_frame = line(indx+8:)
          read(100,'(a)') line     ! Format version line 
          read(100,'(a)') line     ! 4-char line
          indx = index(line,'ID')
          site_name = line(indx+3:)
          call casefold(site_name)
*         get long site name
          read(100,'(a)') line     ! Station name
          indx = index(line,'name')
          gsite_full = line(indx+5:)
*
          read(100,'(a)') line     ! Begin   Date 
          read(100,'(a)') line     ! End     Date 
          read(100,'(a)') line     ! Release Date 
*         Now get position
          read(100,'(a)') line
          indx = index(line,'position,')
          read(line(indx+9:),*,iostat=jerr) int_coord(1)
          call report_error('IOSTAT',ierr,'read',line,0,'PBO Lat')
          indx = index(line,'Latitude,')
          read(line(indx+9:),*,iostat=jerr) int_coord(2)
          call report_error('IOSTAT',ierr,'read',line,0,'PBO Long')
          indx = index(line,'Longitude,')
          read(line(indx+10:),*,iostat=jerr) int_coord(3)
          call report_error('IOSTAT',ierr,'read',line,0,'PBO Hgt')
*         Convert to co-lat, long in radias.
          int_coord(1) = pi/2 - int_coord(1)*pi/180
          int_coord(2) = int_coord(2)*pi/180
          if( int_coord(2).lt.0 ) int_coord(2)=int_coord(2)+2*pi

*         Next is date line
          read(100,'(a)') line
      endif

****  Compute the conversion factors of dN, dE in mm to radians
      con(1) = 1.d0/earth_rad
      con(2) = 1.d0/(sin(int_coord(1))*earth_rad)
      date(4) = 12
      date(5) = 0
      sec     = 0.0

*     Initialize the correlations.
      do j = 4,6
         xyz_std(j) = 0.d0
         neu_std(j) = 0.d0
      end do

      neuset = .false.

****  Now start reading the dN, dE, dU values
      do while ( ierr.eq.0 )
         read(100,'(a)',iostat=ierr) line
         if( ierr.eq.0 ) then
*            Decode the line.  Remove commas
             call sub_char(line,',',' ')
             indx = 0
             call GetWord(line,datestr,indx)
             read(datestr,'(I4,1x,I2,1x,I2)') (date(j),j=1,3)
             call ymdhms_to_mjd(date, sec, gmjd)
             if( gmjd.lt.first_mjd ) first_mjd = gmjd
             if( gmjd.gt.last_mjd  ) last_mjd  = gmjd

*            Now read the dN, dE, dU values
             call multiread(line,indx,'R8',jerr,dNEU,cdum,3)
* NOTE: SCEC format is ENU not NEU as in PBO
*            Now get the position by adding dN,dE,dU to reference
*            (Minus on first term because dlat to dco-lat.
             if( pbo_form ) then
                unc_geod(1) = int_coord(1) - dNEU(1)*con(1)*1.d-3
                unc_geod(2) = int_coord(2) + dNEU(2)*con(2)*1.d-3
                unc_geod(3) = int_coord(3) + dNEU(3)*1.d-3
             else
                unc_geod(1) = int_coord(1) - dNEU(2)*con(1)*1.d-3
                unc_geod(2) = int_coord(2) + dNEU(1)*con(2)*1.d-3
                unc_geod(3) = int_coord(3) + dNEU(3)*1.d-3
             endif

*            Convert co-lat, long in radians to degrees
             unc_llu(1) = (pi/2-unc_geod(1))*180/pi
             unc_llu(2) = unc_geod(2)*180/pi 
             unc_llu(3) = unc_geod(3)

****         Now convert to XYZ
             call geod_to_geod(unc_geod, pos_xyz_fin, 
     .            'GEOD', 'XYZ','WGS84','WGS84',zone,hemi)

*            Now see if we can sigmas (available in PBO format)
             if( pbo_form ) then
                 call multiread(line,indx,'R8',jerr,neu_std,cdum,3)
                 do i = 1,3
                    neu_std(i) = neu_std(i) * 1.d-3  ! Convert to m
                 enddo
             else
                 neu_std(1) = 1.0d-3
                 neu_std(2) = 1.0d-3
                 neu_std(3) = 3.0d-3
             endif

*            Form the NEU covariance matrix
             do i = 1,3
                do j = 1,3
                    ucov_neu(i,j) = 0.d0
                enddo
                ucov_neu(i,i) = neu_std(i)**2
             end do

*            Now transform NEU sigmas to XYZ sigmas
             if( .not. neuset ) then 
                call XYZ_to_GEOD(rot_mat, pos_xyz_fin, loc_coord )
c....           Now transpose rot_matrix to that NEU to XYZ direction 
                do i = 1,2
                  call dvswp(rot_mat(i,i+1),3,rot_mat(i+1,i),1, 3-i)
                end do
*               Now get the average scale needed
                neuset = .true.
             end if
*            Now convert sigmas to XYX
             call var_comp(rot_mat,ucov_neu, ucov_xyz, tcov, 3,3,1)

****         Compute the XYZ sigmas and add correlations
             do i = 1,3
                xyz_std(i) = sqrt(ucov_xyz(i,i))
             end do
             xyz_std(4) = ucov_xyz(1,2)/( xyz_std(1)*xyz_std(2))
             xyz_std(5) = ucov_xyz(1,3)/( xyz_std(1)*xyz_std(3))
             xyz_std(6) = ucov_xyz(2,3)/( xyz_std(2)*xyz_std(3))

****         Get final NEU values
             call loc_to_geod(unc_geod, pos_neu_fin)

             if( jerr.eq.0 ) then

*              See if we can match site name
               indx = 0
               call get_cmd(site_name,gsite_names,num_site,ns,indx)
               if( ns.le.0 ) then
                   num_site = num_site + 1
                   if( num_site.gt.max_site ) then 
                      call report_stat('FATAL','TSCON','read_in_csv',
     .                     '','Too many sites',max_site)
                   endif

                   num_code = num_site
                   ns = num_site
                   gsite_names(ns) = site_name
                   in_code(ns) = site_name(1:4)
                   in_full(ns) = gsite_full
               endif

*****          OK save this entry
               num_ent = num_ent + 1
              if( num_ent.gt.max_ent ) then
                  call report_stat('FATAL','ts_con','read_in_csv','',
     .                'Too many enties Max ',max_ent)
              endif
               ne = num_ent
               in_ns(ne)  = ns
               in_cs(ne)  = ns
               in_mjd(ne) = gmjd
               in_type(ne) = ts_ref_type
               do j = 1,3
                  in_xyz(j,ne) = pos_xyz_fin(j)
                  in_neu(j,ne) = pos_neu_fin(j)
                  in_llu(j,ne) = unc_llu(j)
               end do
* MOD TAH 140829: Scale the sigmas
              xyz_std(1:3) = xyz_std(1:3)*sigma_scale
              neu_std(1:3) = neu_std(1:3)*sigma_scale
               do j = 1,6
                  in_xyz_std(j,ne) = xyz_std(j)
                  in_neu_std(j,ne) = neu_std(j)
               end do

            end if
         end if
      end do

****  Tell user were we are
      write(*,110) in_file(1:trimlen(in_file)), num_site, num_code, 
     .             num_ent, sigma_scale
 110  format('File: ',a,' Sites ',i5,' Codes ',i5,' Entries ',i10, 
     .       ' SigScale ',F10.3)
      return
      end 

CTITLE CM_PARTS

      subroutine cm_parts( cmpart, in_xyz )

      implicit none

*     Routine to compute the dXYZ in the CF frame for dXYZ of the CM due
*     to a degree 1 load.   Eq 6, Blewitt, Lavalee, Clarke et al.,
*     "A New Global Model of Earth Deformation: Seasonal cycle detected",
*      Science, 294, 5550, pp 2342-2345, 2001.

* PASSED
       real*8 cmpart(3,3)    ! Partial derivatives (output)
       real*8 in_xyz(3)      ! Site coordinates (XYZ m).

* LOCAL
       real*8 geod_pos(3)    ! Geodetic co-lat, long and height
       real*8 rot(3,3)       ! Rotation matrix from XYZ to NEU
       real*8 ct, st, cl, sl ! Cos/Sin Theta (co-lat) and longtitude

       real*8 lp             ! l' Love number (Han, JGR 2016 uses l'/(1+k')
                             ! l' = 0.134 ; k' = 0.021

       data lp / 0.131 / 

*      Get Geodetic coordinates
       call XYZ_to_GEOD(rot, in_xyz, geod_pos )
*      Compute cosines and sines
       ct = cos(geod_pos(1))
       st = sin(geod_pos(1))
       cl = cos(geod_pos(2))
       sl = sin(geod_pos(2))

****   We code explicitly here but could be computed as l'*rot'*diag(1,1,-2)*rot
       cmpart(1,1) =  lp*(ct**2*cl**2 - 2*cl**2*st**2 + sl**2)
       cmpart(1,2) = -lp*(3*cl*st**2*sl)
       cmpart(1,3) = -lp*(3*ct*cl*st)

       cmpart(2,1) = -lp*(3*cl*st**2*sl) 
       cmpart(2,2) =  lp*(cl**2+(3*cos(geod_pos(1)*2)-1)*sl**2/2)
       cmpart(2,3) = -lp*(3*ct*st*sl)

       cmpart(3,1) = -lp*(3*ct*cl*st)
       cmpart(3,2) = -lp*(3*ct*st*sl)
       cmpart(3,3) =  lp*(-2*ct**2+st**2)

****   Thats all
       return
       end

CTITLE WRITE_ORG_HEAD

      subroutine write_org_head(unit,ierr)

      implicit none

*     Routine to write header lines for output of origin parameters

      include 'tsfit.h'
      include 'tscon.h'
      include '../includes/const_param.h'

* PASSED VARIABLES
      integer*4 unit    ! Output unit number
     .,         ierr    ! IOSTAT error

* LOCAL VARIABLES

      real*8 scl(10)    ! Scale factors for output
      character*5 labels(max_np), units(max_np)

      integer*4 i,j     ! Loopc counters
      integer*4 pos     ! Position in line
      integer*4 trimlen ! Length of string

      logical kbit      ! Test bits

      character*256 line1, line2
      
      data labels / 'XROT ', 'YROT ', 'ZROT ', 
     .              'XTRAN', 'YTRAN', 'ZTRAN',
     .              'SCALE',
     .              'XCM  ', 'YCM  ', 'ZCM  ' /
      data units /  'urad ', 'urad ', 'urad ',
     .              ' mm  ', ' mm  ', ' mm  ', 
     .              ' ppb ',
     .              ' mm  ', ' mm  ', ' mm  ' /
      data scl   /    1.0  ,    1.0 ,    1.0 ,
     .              1000.  ,  1000. ,  1000. ,
     .                1.0  ,
     .              1000.  ,  1000. ,  1000. /


      write(unit,120, iostat=ierr) comfile(1:trimlen(comfile))
 120  format('# TSCON: Command file ',a)
      call report_error('IOSTAT',ierr,'writ',out_org,1,
     .                  'write_org_head')

      write(unit,140) stab_it, stab_nsig, use_ratio, cnd_parts_bits,  
     .             stab_min_dne, stab_min_dh 
 140  format('# TSCON: ',I3,' Iterations, Edit ',F4.1,' sigma',
     .       ' UseRatio ',F5.3,' ESTBITS ',o6,
     .       ' Min NE ',F7.4,' H ',F7.4,' m')

      line1 = '#  DecYr    MJD       Nref '
      line2 = '# '

      j = 0
      do i = 1, max_np
         if( kbit(cnd_parts_bits,i) ) then
            j = j + 1
            pos = 28 + (j-1)*18
            write(line1(pos:),210) labels(i)
 210        format(2x,a,5x,'+-')
            write(line2(pos:),220) units(i), units(i)
 220        format(2x,a,5x,a)

            units_scale(j) = scl(i)
         endif
      enddo
      write(unit,'(a)') trim(line1) 
      write(unit,'(a)') trim(line2) 

      end


CTITLE WRITE_ORG_VALS

      subroutine write_org_vals( unit, nc, nr, dsol_cnd, ptp, 
     .                emjd, units_scale, postfit_chi, max_np )

      implicit none

*     Routine to write out each daily value

* PASSED VARIABLES
      integer*4 unit               ! Output unit
     .,         nc                 ! Number of estimated parameters
                                   ! (max val is max_np)
     .,         nr                 ! Number of reference sites
     .,         max_np             ! Dimensioning for PTP covariance

      real*8 dsol_cnd(nc)          ! Solution vector
     .,      ptp(max_np,max_np)    ! Covariance matrix
     .,      emjd                  ! MJD of day (12:00 GPST)
     .,      units_scale(nc)       ! Scaling for output units
     .,      postfit_chi           ! Postfit chi**2 (scaling)


* LOCAL VARIABLES
      integer*4 date(5)            ! YMD HM
     .,         j                  ! Counter

      real*8 decyr                 ! Deciminal years
     .,      sectag                ! Seconds tag

      
      call jd_to_decyrs(emjd, decyr)
      call mjd_to_ymdhms(emjd, date, sectag)

      write(unit,210) decyr, emjd, nr,
     .      (dsol_cnd(j)*units_scale(j), 
     .       sqrt(ptp(j,j)*postfit_chi)*units_scale(j), j=1,nc)
 210  format(1x,F9.4,1x,F9.2,1x,I4,10(1x,F8.3,1x,F8.3))

      return
      end

CTITLE DATEWORD_TO_MJD

      subroutine dateword_to_mjd( dateword, mjd )

      implicit none

*     Routine to convert string date of form <YY><MON><DD> e.g. 07DEC11
*     into MJD value (returns real*8) 

* PASSED
      character*(*) dateword  ! Form 07DEC11 (2007/12/11)
      real*8 mjd              ! MJD value

* LOCAL
      character*3 months(12)
      character*3 month    ! Month read from dateword
      integer*4 date(5)    ! Yr mon day hr mn
      integer*4 j, ierr    ! Loop counter and IOSTAT error
      real*8    sectag     ! Seconds tag set to zero



      data months / 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     .              'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'  /

*     Decode the string
      date = 0
      sectag = 0
      read(dateword,'(I2,A3,I2)',iostat=ierr) date(1), month, date(3)
      if( ierr.eq.0 ) then
         ! Find month
         do j = 1,12
            if( month.eq.months(j) ) then
               date(2) = j
               exit
            endif
         end do
         if( date(2).eq.0 ) then
            write(*,120) dateword, date(1:3)
 120        format('WARNING: No month in ',a,' Decoded date ',3i3)
         endif
         call ymdhms_to_mjd(date, sectag, mjd)
      else
         call report_stat('WARNING','TSCON','dateword_to_mjd', dateword,
     .                    'IOSTAT error', ierr)
         mjd = 0
      endif

****  Thats all
      end


  
