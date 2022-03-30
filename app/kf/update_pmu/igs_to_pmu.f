      program igs_to_pmu

*     Program to read an IGS formated IERS pmu file and generate
*     a globk pmu file on the day values

* MOD TAH 210704: Added additional L (Lagrange 4-point interpolator)
*     option that will generate tabular points that when interpolated
*     should reproduce offset and rate values close to the IERS IGS
*     values on the 12:00 epochs.  IL option can be used to intergrate
*     UT1 with new approch.
* MOD TAH 210826: Updated L appproch to directly compute the tablualar
*     values with extrapolation based on MATLAB WK186x.m code.
*     All matrices are 12x14 (12 tabular values for 7 IERS offset and
*     rate values. (Note:LOD needs to be changed to dUT1/dt.  L code
*     only works for 7 input values.
*
      integer*4 max_in   ! Maximum number of values allowed in input
      integer*4 max_out  ! Max. number in output
      parameter ( max_in = 36500 )
      parameter ( max_out = max_in + 5 )

      integer*4 max_nl     ! Maximum number of leap seconds
      parameter ( max_nl = 1000 )

* PROGRAM VARIABLES
      integer*4 num_in, num_out,n  ! Number of input and output records
      integer*4 len_run, rcpar   ! Length of runstring
      integer*4 ierr,jerr        ! IOSTAT error
      integer*4 idls             ! Input index for any leap seconds
      integer*4 i,j
      integer*4 trimlen          ! Length of string
      integer*4 date(5)

      integer*4  nl  	! NUmber of leap seconds
     .,          ils    ! counter in leap-sec file for leap second
                        ! to be used.

      real*8 in(max_in,14)   ! Input IGS values (jd, pole, ut1 etc) 1e-6, 1e-7
                             ! mas and ms
      real*8 save_ut         ! Copy of the input UT1 values in case there
                             ! is a leap second and UT1 integrated).
      real*8 jdo(max_out)    ! Ouput epochs
      real*8 xao(max_out), yao(max_out), uao(max_out) ! Output X Y and UT1
      real*8 xa, ya, ua      ! Averaged values
      real*8 dls, dut        ! Any offset due to leap second
      real*8 dt, off, sectag  ! Time offset between in and out values

      real*8 dut1, dlod      ! Tidal contributions to UT1 and LOD (seconds
                             ! as returned from short_period_lod

      real*8 xjd       ! PEP JD and standard JD
     .,      leap_mjd(max_nl)   ! Julian dates of leap seconds
     .,      tai_utc(max_nl)   ! Value of TAI-UTC after the date of
                               ! each leap second.

* MOD TAH 210704: Matrices needs for Lagrange interpoator tabular
*     genration.  Derived with Matlab code Wk185x.m
      real*8 A_inv(12,14) ! Matrix to take 7 values and 7 rate estimates
                        ! to create 12-tabular points.  The center two 
                        ! are best determined.(numerical values for Matlab)
                        ! L option: Weight for X and Y pole.  Rates down-
                        ! by 25 relative to offset values/
      real*8 L_inv(12,14) ! A- matrix wtth offset and rate with same sigma.
      real*8 ILinv(12,14) ! A- matrix for UT1 where LOD is integrated.  
                        ! UT1 from integration with mean value matching input.

      real*8 vrin(14)    ! Input 7-values then 7 rates for use with A_inv.
                        ! Fitr UT1-UTC need to change sign of LOD values
      real*8 vrout(12)   ! Output 12-tabular values 
  
      real*8 Lm1(4), Lm2(4)  !  Lagrange coefficents to extraplote minus 1 and 2
                        ! days based on coefficinent formulas in gen_util/lagrange_int.f
      real*8 Lp1(4), Lp2(4)  ! Lagrange coefficents to extraplote plus 1 and 2 days

      real*8 RefXYU(3)  ! Initial reference values, removed so values near
                        ! zero and add bavk at end (nned for Lagrange option).

      character*8   int_opt    ! Runstring entry to say integrate UT1.

      character*256 infile, outfile ! Input and output file names
      character*256 line
      character*256 leap_file ! Leap second file
     .,             home_dir  ! Home directory (gg/tables expected here)

* MOD TAH 210724: Matrices from Matlab: Note: Values entered by columns
* MOD TAH 210826: Updated to compute directly: 12 out for 14 input,
*                 Matrix entered by column.Aminus Matrx 12 14 Type 1 Weight 2.50e+01
      data A_inv /   0.5839406,  0.7762637,  1.1728733,  0.4795780, 
     .              -0.1328827,  0.0847648, -0.0436536,  0.0237205, 
     .              -0.0131914,  0.0077690, -0.0087440, -0.0161149, 
     .               0.6284917,  0.3393095, -0.2785989,  0.6741911, 
     .               0.6276060, -0.2252641,  0.1327660, -0.0700963, 
     .               0.0392300, -0.0230729,  0.0259762,  0.0478748, 
     .              -0.3265894, -0.1774717,  0.1603528, -0.2495661, 
     .               0.6511056,  0.6393828, -0.2322861,  0.1376881, 
     .              -0.0752169,  0.0444655, -0.0500068, -0.0921522, 
     .               0.1745495,  0.0946732, -0.0837887,  0.1449755, 
     .              -0.2371412,  0.6442903,  0.6442903, -0.2371412, 
     .               0.1449755, -0.0837887,  0.0946732,  0.1745495, 
     .              -0.0921522, -0.0500068,  0.0444655, -0.0752169, 
     .               0.1376881, -0.2322861,  0.6393828,  0.6511056, 
     .              -0.2495661,  0.1603528, -0.1774717, -0.3265894, 
     .               0.0478748,  0.0259762, -0.0230729,  0.0392300, 
     .              -0.0700963,  0.1327660, -0.2252641,  0.6276060, 
     .               0.6741911, -0.2785989,  0.3393095,  0.6284917, 
     .              -0.0161149, -0.0087440,  0.0077690, -0.0131914, 
     .               0.0237205, -0.0436536,  0.0847648, -0.1328827, 
     .               0.4795780,  1.1728733,  0.7762637,  0.5839406, 
     .              -2.6937665, -1.6037574, -0.4338962,  0.1512196, 
     .              -0.0859566,  0.0455143, -0.0243164,  0.0131074, 
     .              -0.0073020,  0.0042988, -0.0048387, -0.0089177, 
     .              -0.1269353, -0.0679933,  0.0575668, -0.1143734, 
     .               0.1113025, -0.0502413,  0.0282137, -0.0150115, 
     .               0.0083901, -0.0049357,  0.0055566,  0.0102409, 
     .               0.0653144,  0.0355393, -0.0319603,  0.0517739, 
     .              -0.1098256,  0.1091235, -0.0491917,  0.0278965, 
     .              -0.0153428,  0.0090603, -0.0101909, -0.0187790, 
     .              -0.0351054, -0.0190366,  0.0168614, -0.0289935, 
     .               0.0494848, -0.1088267,  0.1088267, -0.0494848, 
     .               0.0289935, -0.0168614,  0.0190366,  0.0351054, 
     .               0.0187790,  0.0101909, -0.0090603,  0.0153428, 
     .              -0.0278965,  0.0491917, -0.1091235,  0.1098256, 
     .              -0.0517739,  0.0319603, -0.0355393, -0.0653144, 
     .              -0.0102409, -0.0055566,  0.0049357, -0.0083901, 
     .               0.0150115, -0.0282137,  0.0502413, -0.1113025, 
     .               0.1143734, -0.0575668,  0.0679933,  0.1269353, 
     .               0.0089177,  0.0048387, -0.0042988,  0.0073020, 
     .              -0.0131074,  0.0243164, -0.0455143,  0.0859566, 
     .              -0.1512196,  0.4338962,  1.6037574,  2.6937665 /

*     Aminus Matrx 12 14 Type 3 Weight 1.00e+00
      data L_inv /   0.7193969,  0.8573152,  0.8022998,  0.4901680, 
     .               0.1674721,  0.0627166,  0.0234866,  0.0088277, 
     .               0.0033711,  0.0014299,  0.0010342,  0.0020360, 
     .               0.1753357,  0.0892058,  0.1236480,  0.3349971, 
     .               0.3311523,  0.1079241,  0.0404377,  0.0151653, 
     .               0.0057903,  0.0024559,  0.0017763,  0.0034970, 
     .               0.0657925,  0.0334275,  0.0463247,  0.1093645, 
     .               0.3297210,  0.3291927,  0.1072232,  0.0402645, 
     .               0.0153401,  0.0065053,  0.0047049,  0.0092626, 
     .               0.0246794,  0.0125362,  0.0173364,  0.0409689, 
     .               0.1073970,  0.3290192,  0.3290192,  0.1073970, 
     .               0.0409689,  0.0173364,  0.0125362,  0.0246794, 
     .               0.0092626,  0.0047049,  0.0065053,  0.0153401, 
     .               0.0402645,  0.1072232,  0.3291927,  0.3297210, 
     .               0.1093645,  0.0463247,  0.0334275,  0.0657925, 
     .               0.0034970,  0.0017763,  0.0024559,  0.0057903, 
     .               0.0151653,  0.0404377,  0.1079241,  0.3311523, 
     .               0.3349971,  0.1236480,  0.0892058,  0.1753357, 
     .               0.0020360,  0.0010342,  0.0014299,  0.0033711, 
     .               0.0088277,  0.0234866,  0.0627166,  0.1674721, 
     .               0.4901680,  0.8022998,  0.8573152,  0.7193969, 
     .              -2.6677489, -1.5855180, -0.6209688,  0.1788639, 
     .               0.0638510,  0.0247201,  0.0092880,  0.0034936, 
     .               0.0013343,  0.0005660,  0.0004093,  0.0008059, 
     .              -0.1759693, -0.0892452, -0.1212524, -0.2960066, 
     .               0.2875357,  0.1037385,  0.0396354,  0.0149123, 
     .               0.0056968,  0.0024164,  0.0017478,  0.0034409, 
     .              -0.0658295, -0.0334285, -0.0461716, -0.1069011, 
     .              -0.2906960,  0.2895334,  0.1045207,  0.0400182, 
     .               0.0152955,  0.0064900,  0.0046942,  0.0092417, 
     .              -0.0246740, -0.0125324, -0.0173213, -0.0408009, 
     .              -0.1049029, -0.2899155,  0.2899155,  0.1049029, 
     .               0.0408009,  0.0173213,  0.0125324,  0.0246740, 
     .              -0.0092417, -0.0046942, -0.0064900, -0.0152955, 
     .              -0.0400182, -0.1045207, -0.2895334,  0.2906960, 
     .               0.1069011,  0.0461716,  0.0334285,  0.0658295, 
     .              -0.0034409, -0.0017478, -0.0024164, -0.0056968, 
     .              -0.0149123, -0.0396354, -0.1037385, -0.2875357, 
     .               0.2960066,  0.1212524,  0.0892452,  0.1759693, 
     .              -0.0008059, -0.0004093, -0.0005660, -0.0013343, 
     .              -0.0034936, -0.0092880, -0.0247201, -0.0638510, 
     .              -0.1788639,  0.6209688,  1.5855180,  2.6677489 / 

*     Aminus Matrx 12 14 Type 4 Weight 1.00e-06
      data ILinv /  -0.5655476,  0.1870474,  0.2160360,  0.2171514, 
     .               0.2171930,  0.2171935,  0.2171926,  0.2171924, 
     .               0.2172038,  0.2175133,  0.2255568,  0.4343819, 
     .               0.2261840,  0.1174746,  0.1132873,  0.1131262, 
     .               0.1131197,  0.1131189,  0.1131184,  0.1131182, 
     .               0.1131241,  0.1132854,  0.1174746,  0.2262350, 
     .               0.2262484,  0.1174826,  0.1132931,  0.1131319, 
     .               0.1131260,  0.1131257,  0.1131251,  0.1131250, 
     .               0.1131309,  0.1132921,  0.1174816,  0.2262485, 
     .               0.2262498,  0.1174824,  0.1132928,  0.1131316, 
     .               0.1131257,  0.1131258,  0.1131258,  0.1131257, 
     .               0.1131316,  0.1132928,  0.1174824,  0.2262498, 
     .               0.2262485,  0.1174816,  0.1132921,  0.1131309, 
     .               0.1131250,  0.1131251,  0.1131257,  0.1131260, 
     .               0.1131319,  0.1132931,  0.1174826,  0.2262484, 
     .               0.2262350,  0.1174746,  0.1132854,  0.1131241, 
     .               0.1131182,  0.1131184,  0.1131189,  0.1131197, 
     .               0.1131262,  0.1132873,  0.1174746,  0.2261840, 
     .               0.4343819,  0.2255568,  0.2175133,  0.2172038, 
     .               0.2171924,  0.2171926,  0.2171935,  0.2171930, 
     .               0.2171514,  0.2160360,  0.1870474, -0.5655476, 
     .              -3.3014138, -1.9162949, -0.9028873,  0.0605936, 
     .               0.0977049,  0.0991339,  0.0991885,  0.0991905, 
     .               0.0991958,  0.0993372,  0.1030106,  0.1983799, 
     .              -1.4526589, -0.7541878, -0.7259124, -0.6892159, 
     .               0.2366426,  0.2723039,  0.2736764,  0.2737289, 
     .               0.2737453,  0.2741355,  0.2842729,  0.5474586, 
     .              -1.2262498, -0.6367360, -0.6139762, -0.6117288, 
     .              -0.5760353,  0.3497842,  0.3854438,  0.3868170, 
     .               0.3868902,  0.3874436,  0.4017712,  0.7737390, 
     .              -0.9999906, -0.5192542, -0.5007350, -0.4999696, 
     .              -0.4985697, -0.4629091,  0.4629091,  0.4985697, 
     .               0.4999696,  0.5007350,  0.5192542,  0.9999906, 
     .              -0.7737390, -0.4017712, -0.3874436, -0.3868902, 
     .              -0.3868170, -0.3854438, -0.3497842,  0.5760353, 
     .               0.6117288,  0.6139762,  0.6367360,  1.2262498, 
     .              -0.5474586, -0.2842729, -0.2741355, -0.2737453, 
     .              -0.2737289, -0.2736764, -0.2723039, -0.2366426, 
     .               0.6892159,  0.7259124,  0.7541878,  1.4526589, 
     .              -0.1983799, -0.1030106, -0.0993372, -0.0991958, 
     .              -0.0991905, -0.0991885, -0.0991339, -0.0977049, 
     .              -0.0605936,  0.9028873,  1.9162949,  3.3014138 / 


      data Lm1 /   4.0,   -6.0,    4.0,   -1.0 /
      data Lm2 /  10.0,  -20.0,   15.0,   -4.0 /
      data Lp1 /  -1.0,    4.0,   -6.0,    4.0 /
      data Lp2 /  -4.0,   15.0,  -20.0,   10.0 /

****  Get the names of the in and output files
      len_run = rcpar(1,infile)
      if( len_run.eq.0 ) then
          call proper_runstring('igs_to_pmu.hlp','igs_to_pmu',1)
      endif
      open(50,file=infile,iostat=ierr,status='old')
      call report_error('IOSTAT',ierr,'open',infile,1,'igs_to_pmu')

      len_run = rcpar(2,outfile)
      if( len_run.eq.0 ) then
          call proper_runstring('igs_to_pmu.hlp','igs_to_pmu',1)
      endif
      open(60,file=outfile,iostat=ierr,status='unknown')
      call report_error('IOSTAT',ierr,'open',outfile,1,'igs_to_pmu')

* MOD TAH 200517: See if UT1 integrate option passed.
      len_run = rcpar(3,int_opt)
      if( len_run.eq.0 ) int_opt = 'N'
      call casefold(int_opt)
      write(60,110) infile(1:trimlen(infile))
      write(*,110) infile(1:trimlen(infile))
 110  format('* GLOBK PMU file generated from IGS file ',a)
      if( int_opt(1:1).eq.'I' ) then
        write(60,120)
        write( *,120)
 120    format('* UT1: Integrated -LOD from first data point')
      endif
      if( index(int_opt,'L').gt.0 ) then
         write(60,125)
         write(*,125)
 125     format('* Tabular points generated by inverse of 4-point ',
     .          ' Lagrange Interpolator (200704 Added)')
      endif


****************
* MOD TAH 090101 Open and read the leap second file
      leap_file = 'leap.sec'
      open(55,file=leap_file, status='old', iostat=ierr)
      if( ierr.ne.0 ) then

*         Could not open file.  Try to generate name
          call getenv('HOME', home_dir)
          leap_file = home_dir(1:trimlen(home_dir)) //
     .                '/gg/tables/leap.sec'
          open(55,file=leap_file, status='old',iostat=ierr)
          if( ierr.ne.0 ) then
              write(*,140) leap_file(1:trimlen(leap_file)), ierr
 140          format('**ERROR** Could not open leap.sec, also tried ',
     .               a,' IOSTAT error ',i5,/,
     .               ' Cannot convert UTC to TAI')
              stop 'IGS_to_PMU: Could not open leap.sec file'
          end if
      end if

****  Now read the leap.sec file and make tai minus utc
      nl = 1
      leap_mjd(1) = 44786.00d0
      tai_utc(1) = +20.d7    ! Make IGS unit (0.1 usec)

      read(55,'(a)' )  line
      read(55,'(a)' )  line
      ierr = 0 
      do while ( ierr.eq.0 ) 
         read(55,'(a)', iostat=ierr) line
         if( ierr.eq.0 .and. trimlen(line).gt.0 ) then
             read(line,*,iostat=jerr) xjd
             if( jerr.eq.0 ) then
                 nl = nl + 1
                 leap_mjd(nl) = xjd + 0.5d0 - 2400000.5d0
                 tai_utc(nl) = tai_utc(nl-1) + 1.d7  ! IGS inits (0.1 usec)
             end if
         end if
      end do

      close(55)

****  OK Read the input file
      ierr = 0
      num_in = 0
      idls = 0
      dls = 0.d0
      do while ( ierr.eq.0 )
         read(50,'(a)',iostat=ierr) line
         if( ierr.eq.0 .and.line(1:1).ne.' ' ) then
             num_in = num_in + 1
             if( num_in.gt.max_in ) then
                 call report_stat('FATAL','IGS_to_PMU','read',' ',
     .                'Too many input values ',num_in)
             end if
             read(line,*,iostat=jerr) (in(num_in,j),j=1,14)
             call report_error('IOSTAT',jerr,'read',line,0,'igs_to_pmu')
             if( jerr.ne.0 ) num_in = num_in - 1

*            See if there is a leap second
C            if( jerr.eq.0 .and. num_in.gt.1 .and. idls.eq.0 ) then
C                dut = in(num_in,4)-in(num_in-1,4)
C                if ( abs(dut).gt.0.5d7 ) then
C                   dls = nint(dut/1.d7)*1.d7
C                   idls = num_in
C                end if
C            endif
C            if( jerr.eq.0 .and. idls.gt.0 ) then
C                in(num_in,4)=in(num_in,4)-dls
C            end if
* MOD TAH 090101: Convert input UTC to TAI by adding in leap seconds
*            so that time series is continuous
             ils = 0
             do i = 1, nl
                if( leap_mjd(i).le.in(num_in,1) ) then
                  ils = i
                endif
             enddo
*            Remove the leap second to make continous series
             if( ils.gt.0 ) then
                 in(num_in,4) = in(num_in,4) - tai_utc(ils)
             else
                 call report_stat('FATAL','IGS_to_PMU','add',' ',
     .                'leap second index 0 at record',num_in)
             endif

         endif
      end do
      write(*,160) num_in, infile(1:trimlen(infile))
 160  format('IGS_TO_PMU: ',i4,' records read from ',a)
C     if( idls.gt.0 ) then
C         write(*,140) idls, dls/1.d7
C140      format('LeapSecond detected at index ',i4,' Value ',F6.3,
C    .           ' seconds')
C     endif

* MOD TAH 031202: Regularize the UT1 and LOD before interpolation
      do i = 1,num_in
         call short_period_lod(in(i,1)+2400000.5d0,dut1,dlod)
*        Remove from read values
         in(i,4) = in(i,4) - dut1*1.d7
         in(i,5) = in(i,5) - dlod*1.d7
C        write(*,180) in(i,1), in(i,4), in(i,5), dut1*1.d7, dlod*1.d7
C180     format('UT1R ',F12.2,4F14.1)
      end do

* MOD TAH 210705: Get reference XYU values and remove from input
      RefXYU = in(1,2:4)
      in(1:num_in,2) = in(1:num_in,2) - RefXYU(1)
      in(1:num_in,3) = in(1:num_in,3) - RefXYU(2)
      in(1:num_in,4) = in(1:num_in,4) - RefXYU(3)


* MOD TAH 200517: See if we are going to integrate UT1 values.
* MOD TAH 210704: Only apply here if L option not selected
* MOD TAH 210705: Moved this code to after removing the short-period
*     changes so that derivatives should ebe smoother.
      if ( int_opt(1:2).eq.'I ' ) then
         do i = 2,num_in
            save_ut = in(i,4)
            in(i,4) = in(i-1,4) - (in(i-1,5)+in(i,5))/2
* MOD TAH 210704: Removed this code because UT1-AT is being used
*           independent of input so no leap seconds.
C           if( abs(save_ut-in(i,4)).gt..5d7 ) then ! Add back
C                                       ! leap second
C               if( in(i,4).lt.save_ut) then
C                  in(i,4) = in(i,4) + 1.d7
C               else
C                  in(i,4) = in(i,4) + 1.d7
C               endif
C           endif
         end do
      endif

*     OK Now generate the output values (values go from 3 before to
*     3 after the values in the input) (Since we only interpolate 
*     between the input values, the number of interpolated values is
*     num_in - 1, hence add 5 to get the 6 extra values
      num_out = num_in + 5
      dt = 0.5d0

* MOD TAH 210704: Non-Lagrange code (original).
      if( index(int_opt,'L').eq.0 ) then
         do i = 1, num_in -1
            jdo(i+3) = in(i,1)+dt
C        r1 = in(i,13)
C        dr = in(i+1,13)-r1
C        xf = in(i,2) + r1*dt+dr*dt^2/2
C        r1 = mit(i+1,13)
C        dr = r1- mit(i,13)
C        xb = in(i+1,2) - r1*dt+dr*dt^2/2
            xa = (in(i+1,2)-in(i+1,13)*dt+(in(i+1,13)-in(i,13))*dt**2/2+
     .            in(i,2)+in(i,13)*dt+(in(i+1,13)-in(i,13))*dt**2/2)/2
            ya = (in(i+1,3)-in(i+1,14)*dt+(in(i+1,14)-in(i,14))*dt**2/2+
     .            in(i,3)+in(i,14)*dt+(in(i+1,14)-in(i,14))*dt**2/2)/2
            ua = (in(i+1,4)+in(i+1,5)*dt-(in(i+1,5)-in(i,5))*dt**2/2+
     .            in(i,4)-in(i,5)*dt-(in(i+1,5)-in(i,5))*dt**2/2)/2
C         write(*,220) i, jdo(i+3), xa, ya, ua
C 220     format(i4,F12.2,3F12.1)
            xao(i+3) = xa*1.e-6
            yao(i+3) = ya*1.e-6
            uao(i+3) = ua*1.e-7 
         end do

****  Now extrapolate the first three and last three values
         do i = 1,3
            off = -(3-i)-dt
            jdo(i) = in(1,1)+off
            xao(i) = (in(1,2)+in(1,13)*off)*1.e-6
            yao(i) = (in(1,3)+in(1,14)*off)*1.e-6
            uao(i) = (in(1,4)-in(1,5)*off)*1.e-7
         enddo
         do i = 1,3
            off = i-dt
            n = num_in
            jdo(n+i+2) = in(n,1)+off
            xao(n+i+2) = (in(n,2)+in(n,13)*off)*1.e-6
            yao(n+i+2) = (in(n,3)+in(n,14)*off)*1.e-6
            uao(n+i+2) = (in(n,4)-in(n,5)*off)*1.e-7
         end do
      ELSE
* MOD TAH 210704: New interpolator and LOD integrator that is 
*        more consistent with 4-pt Lagrangian interpolation.
*        Calculation in two parts here.  Generate tabular
*        points around values (i.e., 7 Offset+rate inputs
*        generates 8 tables points (-0.5 to +0.5 from start and end).
*        Then we need to extrapolate 2 points additional at start
*        in end (only one needed for Largrange interpolation, 2nd
*        need because arc starts early and ends late).
* MOD TAH 210826: Updated to direct calculation but only value
*        for 7 tabular points.
         if( num_in.ne.7 ) then
            write(*,230) num_in
 230        format(/,'** ERROR ** L option only works with 7 input ',
     .               'values. ',I3,' is incorrect number')
            stop 'IGS_TO_PMU: L option only workd for 7 days'
         endif

*        Generate time tags of all data:
         do i = 1,num_in+5
            jdo(i) = in(1,1)+(i-4)+dt
         enddo

*        Set up the input arrays (7 offsets then 7 rates).
*        X-pole Values and Rates
         vrin      = (/ in(1:7,2), in(1:7,13) /)
         xao(1:12) =  matmul(A_inv,vrin)*1.e-6
*        Y-pole Values and Rates
         vrin      = (/ in(1:7,3), in(1:7,14) /)
         yao(1:12) =  matmul(A_inv,vrin)*1.e-6
*        UT1 values and rates used for interpolation.
*        Change LOD sign so that it is dUT/dt.
         vrin = (/ in(1:7,4), -in(1:7,5) /)
         if ( index(int_opt,'I').eq.0 ) then
*           Both UT1 and LOD values are used with equal weight
            uao(1:12) =  matmul(L_inv,vrin)*1e-7
         else
*           Get UT1 by integrating LOD and aligning average
*           value of UT1.
            uao(1:12) =  matmul(ILinv,vrin)*1.e-7
         endif

      ENDIF

****  Now undo any leap second corrections we have made
C     if( idls.ne.0 ) then
C         write(60,240) idls, dls*1.e-7
C240      format('* Leap second detected at index ',i4,' of ',
C    .           F6.3,' seconds')
C         do i = idls+2, num_out
C            uao(i) = uao(i) + dls*1e-7
C         end do
C     end if

* MOD TAH 210705: Add back the RefXYU value to to the output
*     tabkes )change units to match output units)
      xao(1:num_in+5) = xao(1:num_in+5) + RefXYU(1)*1.d-6
      yao(1:num_in+5) = yao(1:num_in+5) + RefXYU(2)*1.d-6
      uao(1:num_in+5) = uao(1:num_in+5) + RefXYU(3)*1.d-7

* MOD TAH 090101: Now add back the leap seconds that we removed.
      do i = 1, num_out
         ils = 0
         do j = 1, nl
            if( leap_mjd(j).le.jdo(i) ) then
               ils = j
            endif
         end do
         uao(i) = uao(i) + tai_utc(ils)*1.d-7  ! Converted back to secs.
      end do

***** Now write the out values
      do i = 1, num_out
          call jd_to_ymdhms(jdo(i), date, sectag)
*         Add back the tidal contributions
          call short_period_lod(jdo(i)+2400000.5d0,dut1,dlod)
*          write(*,150) jdo(i), uao(i)*1.d7, dlod*1.d7, 
*     .                 dut1*1.d7, dlod*1.d7
 
          uao(i) = uao(i) + dut1

          write(60,320) date, xao(i), yao(i), uao(i)
 320      format(1x,i4,4i3,F11.6,'  0.0001 ',F11.6,'  0.0001 ',
     .           F12.7,'  0.00001')
      end do
*
      close(60)

      
      end
      




