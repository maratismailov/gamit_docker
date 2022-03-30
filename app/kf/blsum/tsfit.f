      program tsfit

*     Program to fit PBO generated times series using a globk 
*     eathquake file input and other optional parameters (such as 
*     periodic signals).
* MOD TAH 141205: Added Kalman filter estimator (KFILT command)

*     Usage:
*     tsfit <command file> <summary file> <list of files/file containing list>
*     where <command file> is command file that contains parameters to
*                          be estimated and option globk eq file.  If the
*                          name NONE is used, just linear rates will be
*                          fit.
*           <summary file> output summary file that contains statistics
*                          and estimated parameters
*           <list of files/file containing list> is either a list of tssum
*                          files to process or a file containing a list
*                          of files.
*
      implicit none

      include 'tsfit.h'

* MAIN Program variables
      integer*4 len_run   ! Length of runstring
     .,         nr        ! Entry in runstinng
     .,         trimlen   ! Length of string
     .,         ierr      ! IOSTAT error
     .,         rcpar     ! Function to read runstring
     .,         ss, se    ! Start and end of solution types.
     .,         iter      ! Iteration through fit
     .,         i         ! Loop counter


      character*256 line    ! Line read from file
      character*256 inname  ! File name with .pos removed

      logical done   ! Set true when job completed
     .,       cont   ! Remains true will run should comtinue
     .,       list   ! Set true if readung from list
     .,       edits  ! Set true while there are data being edited 
                     ! based on n-sigma limits
      


***** OK: Start decoding the runstring.
      len_run = rcpar(1,cmdfile)
      if( len_run.eq.0 ) then
          call proper_runstring('tsfit.hlp','tsfit',1)
      endif

      len_run = rcpar(2,sumfile)
      if( len_run.eq.0 ) then
          call proper_runstring('tsfit.hlp','tsfit',1)
      endif

      tsprog = 'tsfit'
      kfopt  = ' '

      call read_tscmd 

      call sort_tsfeq

****  Open the output summary file
      uns = 200
      unv = 0
      una = 0
      unp = 0
      open(uns,file=sumfile,iostat=ierr,status='unknown')
      call report_error('IOSTAT',ierr,'open',sumfile,1,'tsfit')


****  Output the tsfit summary file header.  The head lists the parameters
*     to be estimated.
      call tsfit_head

****  OK: Now start working our way through the files.
      done = .false.
      list = .false.
      nr = 2
      do while ( .not.done )
          cont = .true.
          if( .not.list ) then 
              nr = nr + 1
              len_run = rcpar(nr,infile)
              if( len_run.eq.0 ) then
                  done = .true.
                  cont = .false.
              endif
              unr = 100
          else
*             Get the next file name
              unr = 101
              read(100,'(a)',iostat=ierr) infile
              if( ierr.ne.0 .or. infile(1:1).eq.'#' .or.
     .            infile(1:1).eq.'*' .or. 
     .            trimlen(infile).eq.0 ) then
                 call trimlead(infile)
                 cont = .false.
                 if( ierr.ne.0 ) done = .true.
              endif
          endif

          if( cont ) then
*             See if this is a pos file or list of files.
              open(unr,file=infile,iostat=ierr, status='old')
              call report_error('IOSTAT',ierr,'open',infile,0,
     .                          'tsfit/pos')
              if( ierr.ne. 0 ) cont = .false.
          end if
          if( cont ) then
*             Check first record
              read(unr,'(a)',iostat=ierr) line
              if ( line(1:11).eq.'PBO Station') then
                  write(*,220) infile(1:trimlen(infile))
 220              format('Processing ',a)
*                 Get the solution type
*                 Name format: P067.pbo.final_frame.pos
*                 Extract the file name (start of name
*                 or after /)
                  ts_file = infile
                  inname = infile
                  call sub_char(inname,'.pos',' ')
           
                  ss = 1
                  do i = trimlen(inname),1,-1
                     if( inname(i:i).eq.'/' ) then
                         ss = i + 1
                         exit
                     endif
                  end do
C                 ss = trimlen(infile)-23
                  if( ss.lt.1 ) ss = 1
                  se = ss + 13
                  solnstr = inname(ss:se)
*
                  call read_tspos(unr)
                  call remove_ejmp(0)
                  call edit_ts
                  edits = .true.
                  iter = 0
                  do while ( edits ) 
                      iter = iter + 1
                      call fit_tspos(edits, iter)
                  end do
* MOD TAH 141205: Once the edits are complete, run the Kalman
*                 filter code to generate final estimates
                  if( trimlen(kfopt).gt.0 ) then
                       call kf_tspos
                  endif            

                  close(unr)
              elseif ( .not.list ) then
*                 List of files
                  list = .true.
                  listfile = infile
                  write(*,240) infile(1:trimlen(infile))
 240              format('Reading list ',a)
                  rewind(100)
              else
                  write(*,260) infile(1:trimlen(infile)),
     .                         listfile(1:trimlen(listfile))
 260              format('PBO POS File ',a,' in list ',a,
     .                   ' is not a PBO tssum position file')
              endif
          endif
      end do

****  Thats all
      end



CTITLE READ_TSCMD

      subroutine read_tscmd

      implicit none

*     Rouitne to read the command file  

      include 'tsfit.h'

* LOCAL VARIABLES
      integer*4 ierr, jerr  ! IOSTAT errors
     .,         indx        ! Position in string
     .,         trimlen     ! Lenght of string
     .,         j           ! Loop counter
     .,         date(5)     ! Date for time range
     .,         num         ! counter

      real*8    val          ! Value read from string
      real*8    sectag       ! Seconds value read from date

      character*8 cmd        ! Commdand read from line.
      character*8 cdum       ! Dummy character for read_line.
      
      character*16 rval(max_ptyp) ! Array read from line

      character*256 line     ! Line read from file
     .,         cval         ! String read from string



****  Initialise
      num_per = 0
      num_pat = 0
      num_eq = 0
      num_rn = 0
      num_off = 0
      av_res = 0.d0
      use_constraints = .false.
      report_edits = .false.
      report_renames = .false.
      report_constraints = .false.
      do j = 1,3
         max_sigma(j) = 0.10d0
      end do
c      do j = 1, max_ptyp
c         rcode(j) = 0
c      end do	     
      detroot = 'ts_'
      real_sigma = .false.
      time_range(1) = 1.d0
      time_range(2) = 1.d20
      eq_out = .false.
      eq_rad_scale = 1.0d0
      mean_only = .false.
      resid_sum = .false.
      eqo_root = 'tsfit_'
      min_sigscale = 10 ! Minumum for rescaling sigmas 
      min_rsig = 30     ! Minimum number for real sigma (1.0.3) 
      max_persig = 5.0d0  ! Max sigma for periodic sigma
      max_eqsig  = 2.0d0  ! Max sigma for log/exp postseismic
      htscale    = 4.0d0  ! Scale for heights for sigmas
      max_rmchi  = 3.0d0  ! Max ratio of parameter estimate to sigma
                          ! for the parameter to be removed.
      rf_out     = .false.  ! Set false so that reference frame line
                          ! will written later.
      pbov_header = .false. ! Set true after writing pbovel header
      nsigma     = 0.0    ! No n-sigma delete by default. Value can
                          ! be passed with MEAN+4.0

****  See if NONE used as name
      if( cmdfile(1:5).eq.'NONE ' ) then
          detroot = 'NONE'
          RETURN
      endif
      if( cmdfile(1:4).eq.'MEAN' ) then
          mean_only = .true.
*         See if sigma limit past
          indx = index(cmdfile,'+')
          if( indx.gt.0 ) then
              call read_line(cmdfile,indx+1,'R8',ierr,val,cdum)
              if( ierr.eq.0 ) nsigma = val
          end if
*         MOD TAH 160223: See of additional option of RES added
          if( index(cmdfile,'RES').gt.0 ) then          
               resid_sum = .true.
          endif
          detroot = 'NONE'
          RETURN
      end if

****  Open the command file
      open(100,file=cmdfile,iostat=ierr,status='old')
      call report_error('IOSTAT',ierr,'open',cmdfile,1,'tsfit')

      
****  Loop over commands
      do while ( ierr.eq.0 )
         read(100,'(a)',iostat=ierr) line
         if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .       trimlen(line).gt.0 ) then
*            OK: Try to read command
             indx = 0
             call GetWord(line, cmd, indx)
             call casefold(cmd)

*            See which comamnd passed
*@ EQ_FILE <File Name>
             if(   cmd(1:2).eq. 'EQ'    ) then   ! Eq_file
*                Get name of eq_file: and read
                 call GetWord(line, eqfile, indx)
                 call subhome( eqfile )

* MOD TAH 170406: See if optional scale passed
                 call read_line(line,indx,'R8',jerr,val,cval)
                 if( jerr.eq.0 ) then
                     eq_rad_scale = val
                 endif           
                 call tsread_eqfile

*@ PERIODIC <Period (days)> [Apriori Sigma]
             elseif ( cmd(1:2).eq.'PE'    ) then   ! Periodic
*                Get the Period and optional sigma
                 call read_line(line,indx,'R8',jerr,val,cval)
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/Period')
                 if( jerr.eq.0 ) then
                     num_per = num_per + 1
                     per_per(num_per) = val
*                    See if sigma passed
                     call read_line(line,indx,'R8',jerr,val,cval)
                     if( jerr.eq.0 ) then
                        per_sig(num_per) = val
                     else
                        per_sig(num_per) = 0.d0
                     endif
                 endif 

*@ SITES <Pattern of using site>
             elseif ( cmd(1:2).eq.'SI'    ) then   ! Sites
*                Start strippinng off the patterns
                 jerr = 0
                 do while ( jerr.eq.0 )
                     call read_line(line,indx,'CH',jerr,val,cval)
                     if( jerr.eq.0 ) then
                         num_pat = num_pat + 1
                         site_pat(num_pat) = cval
                     endif
                enddo

*@ USE_CONSTRAINTS
             elseif ( cmd(1:2).eq.'US'    ) then   ! Use_constraints
                 use_constraints = .true.

*@ REP_EDITS
            elseif ( cmd(1:5).eq.'REP_E' ) then   ! Report Edits
                 report_edits = .true.
*                See if rename file passed
                 call GetWord(line, rename_file, indx)
                 if ( trimlen(rename_file).gt.0 ) then
                    call wild_card(rename_file,sumfile)
                    report_renames = .true.
                 endif


*@ REP_CONSTRAINTS
             elseif ( cmd(1:5).eq.'REP_C' ) then   ! Report Constraints
                 report_constraints = .true.

*@ AV_RES <days>
             elseif ( cmd(1:2).eq.'AV'    ) then   ! Average Residuals
                 call read_line(line,indx,'R8',jerr,val,cval)
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/av_res')
                 if( jerr.eq.0 ) then
                     av_res = val
                 endif

*@ REAL_SIGMA
             elseif( cmd(1:4).eq.'REAL' ) then
                 real_sigma = .true.

*@ DETROOT <det_root>
             elseif( cmd(1:4).eq.'DETR' ) then                 
                 call GetWord(line, detroot, indx)
* MOD 131111: Allow ~ and wildcard substution
                call wild_card(detroot,sumfile)
C               call subhome( detroot)

*@ VELFILE <vel file name>
             elseif( cmd(1:4).eq.'VELF' ) then                 
                 call GetWord(line, velfile, indx)
                 call wild_card(velfile,sumfile)

*@ NSIGMA <nsigma limit>
             elseif( cmd(1:4).eq.'NSIG' ) then
                 call read_line(line,indx,'R8',jerr,val,cval)
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/nsigma')
                 if( jerr.eq.0 ) then
                     nsigma = val
                 endif

*@ MAX_SIGMA <Sig N> <Sig E> <Sig U> meters
             elseif( cmd(1:5).eq.'MAX_S' ) then
*                Get start time
                 call multiread(line,indx,'R8',jerr,max_sigma,cval,3)
            

*@ TIME_RANGE <Start Date> <End Date>
             elseif( cmd(1:4).eq.'TIME' ) then
*                Get start time
                 call multiread(line,indx,'I4',jerr,date,cval,5)
                 val = 0.0
                 call ymdhms_to_mjd(date, val, time_range(1))
*                See if end time passed
                 call read_line(line,indx,'I4',jerr,date,cval)
                 if( jerr.eq.0 ) then
*                    Read rest of line
                     call multiread(line,indx,'I4',jerr,date(2),cval,4)
                     val = 0.0
                     call ymdhms_to_mjd(date, val, time_range(2))
                 endif

*@ OUT_APRF <file name> 
             elseif( cmd(1:5).eq.'OUT_A' ) then
                 call GetWord(line, outapr_file, indx)
                 call wild_card(outapr_file,sumfile)

*@ OUT_EQROOT <root to output earthquake files> [out days]
             elseif( cmd(1:5).eq.'OUT_E' ) then
                  call GetWord(line, eqo_root, indx)
                  call wild_card(eqo_root,sumfile)
                  call read_line(line,indx,'R8',jerr,outlog_days,cval)
                  if( jerr.ne.0 ) outlog_days=0
                  eq_out = .true.

*@ OUT_PBOV <file name>
             elseif( cmd(1:6).eq.'OUT_PB' ) then
                  call GetWord(line, out_pbovel, indx)
                  call wild_card(out_pbovel,sumfile)

*@ OUT_POSF <file name>
             elseif( cmd(1:6).eq.'OUT_PO' ) then
                  call GetWord(line, out_posf, indx)
                  call wild_card(out_posf,sumfile)


*@ RESROOT <root to output residual files> [out days]
             elseif( cmd(1:4).eq.'RESR' ) then
                  call read_line(line, indx, 'CH', jerr, val, cval)
                  call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/resroot')
                  if( jerr.eq.0 ) then
                     resroot = cval
                     call subhome(resroot)
                  endif

*@ RESXCL <parametes to exclude from residual calculation > [offset rate periodic break earthquake log exp drate]
             elseif( cmd(1:4).eq.'RESX' ) then
                  num = 0    ! Sets multiread to read all values on line.
                  call multiread(line,indx,'CH',jerr,val,rval,num)
		  j = 1
		  do while (j .le. num )
		    call casefold(rval(j))
		    if (    rval(j)(1:2) .eq. "OF" )   then
		       rcode(j) = 1
		    else if ( rval(j)(1:2) .eq. "RA" ) then
		       rcode(j) = 2
		    else if ( rval(j)(1:2) .eq. "PE" ) then
		       rcode(j) = 3
*                      Periodic are patameter type numbers 3 and 4 (set both)
		       j = j + 1
		       rcode(j) = 4
		    else if ( rval(j)(1:2) .eq. "BR" ) then
		       rcode(j) = 5
		    else if ( rval(j)(1:2) .eq. "EA" ) then
		       rcode(j) = 6
		    else if ( rval(j)(1:2) .eq. "LO" ) then
		       rcode(j) = 7
		    else if ( rval(j)(1:2) .eq. "EX" ) then
		        rcode(j) = 8
		    else if ( rval(j)(1:2) .eq. "DR" ) then
		        rcode(j) = 9
		    else 
                        call report_error('IOSTAT',-1,'decod',rval(j),
     .                  0,'read_tscmd/RESXCL - unknown paramater type')
                    endif
		    j = j + 1
		  enddo	
		  	 
*@ MINNUM <Minimum number of values to rescale sigma> <min Real-sigma>
             elseif (cmd(1:4).eq.'MINN' ) then
                  call read_line(line, indx, 'R8', jerr, val,  cval)
                  call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/MINSCALE RESCALE')
                  if( jerr.eq.0 ) then
                      min_sigscale = nint(val)
                  endif
*                 Real-sigma min is optional
                  call read_line(line, indx, 'R8', jerr, val,  cval)
                  if( jerr.ne.0 .and. jerr.ne.-1 )
     .            call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/MINSCALE REALSIGMA')
                  if( jerr.eq.0 ) then
                      min_rsig = nint(val)
                  endif


*@ MAX_PARSIG <Max per sigma (mm)> <Max log/exp sigma (mm)> <Ht scale> <max_rmchi>
             elseif (cmd(1:5).eq.'MAX_P' ) then
                  call read_line(line, indx, 'R8', jerr, val,  cval)
                  call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/MAX_PARSIG PERIODIC')
                  if( jerr.eq.0 .and. val.ne.-1.d0 ) then
                      max_persig = val
                  endif
*                 Max earthquake log/exp sigma is optional
                  call read_line(line, indx, 'R8', jerr, val,  cval)
                  if( jerr.ne.0 .and. jerr.ne.-1 )
     .            call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/MAX_PARSIG LOG/EXP')
                  if( jerr.eq.0 .and. val.ne.-1.d0 ) then
                      max_eqsig = val
                  endif
*                 Scale for height is optional
                  call read_line(line, indx, 'R8', jerr, val,  cval)
                  if( jerr.ne.0 .and. jerr.ne.-1 )
     .            call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/MAX_PARSIG HTSCALE')
                  if( jerr.eq.0 .and. val.ne.-1.d0 ) then
                      htscale = val
                  endif
*                 Max remove chi (optional)
                  call read_line(line, indx, 'R8', jerr, val,  cval)
                  if( jerr.ne.0 .and. jerr.ne.-1 )
     .            call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'read_tscmd/MAX_PARSIG MAXCHI')
                  if( jerr.eq.0 .and. val.ne.-1.d0 ) then
                      max_rmchi = val
                  endif
*@ MEAN_ONLY 
             elseif( cmd(1:2).eq.'ME' ) then
*                 Set logical mean_only to estimate only mean
                  mean_only = .true.
* MOD TAH 160223: See if RES has been added
                  if( index(cmd,'RES').gt.0 ) resid_sum = .true.

*@ KFILT <RW/FOGM/WN> <Min RW Variance>
             elseif( cmd(1:2).eq.'KF' ) then
* MOD TAH 141205: Added KFILT option.
*                 Use Kalman filter rather than WLS estimator. 
*                 This option used with RealSigma to get the process
*                 noise and optional correlation time. 
                  call GetWord(line, kfopt, indx)
                  if( trimlen(kfopt).eq.0 ) kfopt = 'RW'
                  call casefold(kfopt)
                  call read_line(line, indx, 'R8', jerr, min_rwvar,  
     .                           cval)
                  if( jerr.ne.0 ) min_rwvar = 0.05 ! mm^2/yr
                 
             else
                 write(*,800) line(1:trimlen(line))
 800             format('** WARNING ** Unknow command: ',a)
             endif
          end if


      enddo

****  Thats all
      close(100)

      return
      end


CTITLE TSFIT_HEAD

      subroutine tsfit_head

      implicit none

      include 'tsfit.h'
 
      integer*4 trimlen, ierr

****  Write out header for summary file (ensum format).

      write(uns,120) cmdfile(1:trimlen(cmdfile))
 120  format('* TSFIT NEU Secular trend components from command file ',
     .       a)
      if( .not.resid_sum) write(uns,160)
 160  format('  Stations    enu  S#   #      mean length    sig   ',
     .     ' wrms   nrms     slope    sig     wrms  nrms   dur   mean',
     .     /,'                                  (m)        (m)    (mm)',
     .       '          (mm/yr)  (mm/yr)   (mm)        (yrs)  (yrs)')
      if( resid_sum) write(uns,170)
 170  format('  Stations    enu  S#   #         Residual    sig   ',
     .     ' wrms   nrms     slope    sig     wrms  nrms   dur   mean',
     .     /,'                                   (mm)      (mm)   (mm)',
     .       '          (mm/yr)  (mm/yr)   (mm)        (yrs)  (yrs)')
 
      if( trimlen(velfile).gt.0 ) then
          unv = 202
          open(unv,file=velfile,iostat=ierr, status='unknown')
          call report_error('IOSTAT',ierr,'open',velfile,0,'tsfit_head')
          if( ierr.eq.0 ) then
             write(unv,220) cmdfile(1:trimlen(cmdfile))
 220         format('* TSFIT Velocity file from command file ',a)

****         Now output the command file to apriori file.
             call sumcmd(unv,cmdfile)
             write(unv,260)
 260         format('*  Long         Lat        Evel    Nvel    dEv  ',
     .           '   dNv    E +-    N +-    Rne      Hvel     dHv    ',
     .           'H +-  Site',/,
     .           '*  deg          deg       mm/yr   mm/yr   mm/yr   ',
     .           'mm/yr   mm/yr   mm/yr            mm/yr   mm/yr',
     .           '   mm/yr')
          else
             unv = 0
          endif
      endif

      if( trimlen(outapr_file).gt.0 ) then
          una = 203
          open(una,file=outapr_file,iostat=ierr, status='unknown')
          call report_error('IOSTAT',ierr,'open',outapr_file,0,
     .        'tsfit_head')
          if( ierr.eq.0 ) then
             write(una,320) cmdfile(1:trimlen(cmdfile))
 320         format('* TSFIT Apriori coordinate file from command file ',
     .               a)
****         Now output the command file to apriori file.
             call sumcmd(una,cmdfile)
          else
             una = 0
          endif
      endif

      if( report_renames ) then
          unt = 204
          open(unt,file=rename_file,iostat=ierr, status='unknown')
          call report_error('IOSTAT',ierr,'open',rename_file,0,
     .        'tsfit_head')
          if( ierr.eq.0 ) then
             write(unt,410) cmdfile(1:trimlen(cmdfile))
 410         format('* TSFIT rename edit file from command file ',
     .               a,/,'* Only type 1 (sigma too large) and 4 ',
     .                   ' (outlier) are reported' )
          else
             unt = 0
          endif
      endif

      if( trimlen(out_posf).gt.0 ) then
          unf = 207
          open(unf,file=out_posf,iostat=ierr, status='unknown')
          call report_error('IOSTAT',ierr,'open',out_posf,0,
     .                      'tsfit_head')
          if( ierr.eq.0 ) then
             write(unf,520) cmdfile(1:trimlen(cmdfile))
 520         format('* TSFIT Position Adjustments from command file ',a)

****         Now output the command file to apriori file.
             call sumcmd(unf,cmdfile)
             write(unf,560)
 560         format('*  Long         Lat        EdP     NdP     dEp  ',
     .           '   dNp    E +-    N +-    Rne      HdP      dHp    ',
     .           'H +-  Site',/,
     .           '*  deg          deg          mm      mm      mm   ',
     .           '   mm      mm      mm               mm      mm',
     .           '      mm')
          else
             unf = 0
          endif
      endif

      return
      end


CTITLE EDIT_TS

      subroutine edit_ts


      implicit none

*     Routine to apply edits to the time series based on
*     rename commands and max_sigma


      include 'tsfit.h'

      integer*4 i,j,k  ! Loop counters

      logical rnm      ! Logical true is solution string from 
                       ! P067.pbo.final
                       !      ^^^       matches hfile rename

*     Check editing
      do i = 1, num_ts
      
*        Add predefined offsets for this site
         do k = 1, num_off
C            print*,'off_codes ts_code ts_mjd off_time-',
C    .         off_codes(k),':',ts_code,ts_mjd(i),off_times(k)
	     if ( off_codes(k) .eq. ts_code .and. 
     .            ts_mjd(i) .ge. off_times(k) ) then
	        do j = 1,3 
                  ts_neu(j,i) = ts_neu(j,i) + off_dpos(j,k)
	        end do
	     endif
         end do
	 
*        Check max_sigma
         do j = 1, 3
            ts_edt(j,i) = 0
            if(ts_neu_std(j,i).gt.max_sigma(j) ) then
               ts_edt(j,i) = 1
            endif
         end do

*        See if outside time range
         if( ts_mjd(i).lt.time_range(1) .or. 
     .       ts_mjd(i).gt.time_range(2) ) then
            do j = 1,3
               call sbit(ts_edt(j,i),4,1)
            end do
         endif

*        Now see if a rename will edit this point
         do k = 1, num_rn
            if( ts_mjd(i).ge.rn_times(1,k) .and.
     .          ts_mjd(i).le.rn_times(2,k) .and.
     .          ts_code.eq.rn_codes(1,k)(1:4)  ) then
 
*               Add any rename offsets found in eq file     
	        do j = 1,3 
                   ts_neu(j,i) = ts_neu(j,i) + rn_dpos(j,k)
	        end do     

*               See if we edit based on solution type
                rnm = .true.
                if( lrn_hf(k).gt.0 ) then
                    if( solnstr(6:5+lrn_hf(k)).ne.
     .                  rn_hfiles(k)(1:lrn_hf(k)) ) then
                        rnm = .false.
                    end if
                end if
                    
                if( (rn_codes(2,k)(6:8).eq.'XCL' .or.
     .              rn_codes(2,k)(6:8).eq.'XPS') .and. rnm ) then
                    do j = 1,3 
                       call sbit(ts_edt(j,i),2,1)
                       call sbit(ts_edt(j,i),1,0)   ! Turn off maxsigma
                    end do
                endif
            endif
         end do
      end do

****  Thats all 
      return
      end
    

CTITLE FIT_TSPOS

      subroutine fit_tspos( edits, iter)

      implicit none

*     Routine to fit parameters to time series

      include 'tsfit.h' 

      logical edits  ! Set true while there edits based on n-sigma
      integer*4 iter

      if( iter.eq.1 ) then
         call count_par
      end if

      call fit_accum( edits ) 

      return
      end

CTITLE COUNT_PAR 

      subroutine count_par

      implicit none

*     Routine to count the number of parameters to be estimated

      include 'tsfit.h' 

      integer*4 np ! Parameter number count
     .,         i, j, k, l   ! Loop counter
     .,         num_rin   ! Number of break-renames for this site
     .,         rin(max_rn)  ! Time sorted list of renames needed
     .,         ind(max_rn)  ! Index slots where the renames are 
                             ! for this site

      integer*4 ltm_hf(max_rn)  ! Temp length of hfile name
      integer*4 date_start(5), date_end(5)
      real*8 sectag
      real*8 tm_times(2,max_rn) ! Times
     .,      tm_dpos(3,max_rn)   ! Positions
     .,      dVlog  ! Estimate velocity from log term at start
                    ! of data span based on distance dependent term.

      character*8 tm_codes(2,max_rn)  ! Codes
     .,           tm_types(max_rn)  ! type
      character*16 tm_hfiles(max_rn)  ! H-file mask

      logical sorted  ! Stays true if renames are sorted
     .,   changed     ! Set true is renames are changed
     .,   nodata      ! Set true when no data on a site for a
                      ! earthquake but we want to keep parameter
                      ! name (not really sure we want to do this).
     .,   log_OK      ! Set false if log term likely to be less
                      ! than 1 mm/yr at start of data

      real*8 dist  ! Cord distance of site from EQ.

****  Initialize the apriori constraints
      do i = 1,max_par
         do j = 1,3
             apr_con(i,j) = 0.d0  ! only non-zero values are applied
         end do
      end do

****  Before we start, make sure that the renames for breaks at this
*     site are in time order
      j = 0
      changed = .false.
      do i = 1, num_rn
         if( rn_codes(1,i)(1:4).eq. ts_code .and.
     .      (rn_codes(2,i)(6:8).ne.'XPS' .and.
     .       rn_codes(2,i)(6:8).ne.'XCL') )  then
*           This rename will be used, add to list
*           and see if sorted.
*           See if time order is OK
            sorted = .true.
            do k = 1,j   ! j is current number of entries
*               See if sorted and whether we need to insert
*               eariler in list
                if( rn_times(1,i).lt.rn_times(1,rin(k)) .and.
     .              sorted ) then
*                   Time is out of order so insert this one
*                   rather than add at ended
                    sorted = .false.
                    changed = .true.
                    do l = j,k,-1
                       rin(l+1) = rin(l)
                    end do
*                   Now add this value
                    rin(k) = i
                end if
            enddo
*           If list is sorted so far then just add
            j = j + 1
            if( sorted ) rin(j) = i
            ind(j) = i   ! This is staight list of slots used.
         end if
      end do 
*     Save number of used renames
      num_rin = j                                   
    
****  Now check that renames are OK and output if we change them
*     Copy site values to tmp storage
      do i = 1,num_rin
         do j = 1,3
            tm_dpos(j,i) = rn_dpos(j,rin(i))
         end do
         do j = 1,2
            tm_times(j,i) = rn_times(j,rin(i))
            tm_codes(j,i) = rn_codes(j,rin(i))
         end do
         ltm_hf(i) = lrn_hf(rin(i))
         tm_types(i) = rn_types(rin(i))
         tm_hfiles(i) = rn_hfiles(rin(i))
      enddo

****  Everying is not sorted; See if any changes are needed
      do i = 2,num_rin
         if( tm_times(1,i).ne.tm_times(2,i-1)) then
*            End of of i-1 entry not correct; FIX
             tm_times(2,i-1) = tm_times(1,i)
             changed = .true.
         end if
      end do

****  If we haved changes any values then output new entries
      if( changed ) then
         write(*,100) tm_codes(1,1)(1:4)
 100     format('* Updated rename entries for site ',a)
         do i = 1, num_rin
             call jd_to_ymdhms(tm_times(1,i), date_start, sectag)
             call jd_to_ymdhms(tm_times(2,i), date_end, sectag)
             if( (tm_dpos(1,i)**2+tm_dpos(2,i)**2+
     .            tm_dpos(3,i)**2).gt. 0.d0 ) then
                write(*,120) tm_codes(1,i), tm_codes(2,i),
     .               tm_hfiles(i),
     .               date_start, date_end, (tm_dpos(j,i), j=1,3),
     .               tm_types(i)
             else
                write(*,120) tm_codes(1,i), tm_codes(2,i),
     .               tm_hfiles(i),
     .               date_start, date_end
             endif
 120         format(' RENAME ',1x,a8,1x,a8,1x,a16,i4,1x,4(i2,1x),3x,
     .           i4,1x,4(i2,1x),2x,
     .           3(F10.4,1x),a3)
*            Now put entry back in list so the all is OK
             do j = 1,3
                rn_dpos(j,ind(i)) = tm_dpos(j,i) 
             end do
             do j = 1,2
                rn_times(j,ind(i)) = tm_times(j,i)
                rn_codes(j,ind(i)) = tm_codes(j,i)
             end do
             lrn_hf(ind(i)) = ltm_hf(i)
             rn_types(ind(i)) = tm_types(i) 
             rn_hfiles(ind(i)) = tm_hfiles(i)
         end do
      end if 

****  OK: Start counting
      np = 1   ! Offset term
      pcode(np) = 1
      pindx(np) = 0
      pflag(np) = 1 
      if( .not. mean_only ) then 
         np = np + 1
         pcode(np) = 2
         pindx(np) = 0
         pflag(np) = 1  
      endif 

*     Loop over periodic terms
      do i = 1,num_per
         np = np + 1
         pcode(np) = 3
         pindx(np) = i
         pflag(np) = 1  
         np = np + 1
         pcode(np) = 4
         pindx(np) = i
         pflag(np) = 1  
      end do

****  Loop over the breaks
      do i = 1, num_rn
         if( rn_codes(1,i)(1:4).eq. ts_code .and.
     .       (rn_codes(2,i)(6:8).ne.'XPS' .and.
     .        rn_codes(2,i)(6:8).ne.'XCL') )  then
*            Make sure time falls in our interval
C            if( rn_times(1,i).gt.ts_mjd(1) .and. 
C    .           rn_times(1,i).lt.ts_mjd(num_ts)  ) then
             np = np + 1
             pcode(np) = 5
             pindx(np) = i
*            If the data overlap this rename indicate
*            this
             if( rn_times(1,i).gt.ts_mjd(1) .and. 
     .           rn_times(1,i).lt.ts_mjd(num_ts)  ) then
                pflag(np) = 1 
             else
                pflag(np) = 0
*               Set apriori variance if no data
                do k = 1,3
                   apr_con(np,k) = 1.d-12
                end do
             endif  
             
*            Check to see if break is before the start
*            time we are using
             if( rn_times(1,i).lt.time_range(1) .or.
     .           rn_times(1,i).gt.time_range(2) ) then
*                We can't determine and so contrain
                 pflag(np) = 0  
                 do k = 1,3
                    apr_con(np,k) = 1.d-12
                 end do
             endif
         endif
      enddo

****  Now test earthquakes
      do i = 1, num_eq
          call eval_dist( ref_xyz, eq_pos(1,i), dist )

* MOD TAH 130312: Compute likely velocity for site form log to 
*         see if we should estimate.
          log_OK = .true.
          if( max(ts_mjd(1),time_range(1)).gt.eq_epoch(i) .and.
     .        log_eq(i) .and. dist.le.eq_rad(i)*eq_rad_scale ) then  ! Data starts after EQ and log, so
*             see if dV contribution is too small. (mm/yr)
              dVlog =  (eq_apr_coseismic(1+3,i)*(eq_depth(i)/dist)**2)/
     .           ((max(ts_mjd(1),time_range(1))-eq_epoch(i))/365.)*1000

              if( dVlog < 1.0 ) then
                  write(*,220) i, eq_codes(i), dVlog, 
     .                (max(ts_mjd(1),time_range(1))-eq_epoch(i))/365.
  220             format('Log too small EQ # ',i3,' Code ',a,
     .                ' dVlog ',F5.3,' mm/yr, dT ',F6.2,' yr')
                  log_OK = .false.
              end if
   
          end if
          if( dist.le.eq_rad(i)*eq_rad_scale .and. log_OK ) then
*             Site is within distance of earthquake; see how times
*             match
C             if( eq_epoch(i).ge. ts_mjd(1) .and. 
C    .            eq_epoch(i).lt. ts_mjd(num_ts) ) then
*             Set the offset parameter
              np = np + 1
              pcode(np) = 6
              pindx(np) = i
              pflag(np) = 1  
              nodata = .true.
              if( eq_epoch(i).ge.max(ts_mjd(1),time_range(1)) .and. 
     .            eq_epoch(i).lt.min(ts_mjd(num_ts),time_range(2))) then
                 nodata = .false.
                 if( use_constraints ) then
* MOD TAH Compute the distance dependent sigma (same as globk)
                     do j = 1,3
                        apr_con(np,j) = eq_apr_coseismic(j,i)**2 + 
     .                      (eq_apr_coseismic(j+3,i)*
     .                                 (eq_depth(i)/dist)**2)**2
                     enddo
                 endif
              else
*                if the data ends after the earthquake we can still estimate
*                log and exp terms so set no_data back to false
                 if( eq_epoch(i).lt.min(ts_mjd(num_ts),time_range(2)) )
     .                 nodata = .false.
                 pflag(np) = 0  ! Still no data for offset
                 do j = 1,3
                    apr_con(np,j) = 1.d-12
                 end do
              endif
* MOD TAH 190723: Scan times series and edit out data point on day of
*             earthquake (needed for 2019/07/04 and 2019/07/06 events in
*             california.
              do j = 1, num_ts
                  if( int(ts_mjd(j))-int(eq_epoch(i)).eq.0 ) then
                     call sbit(ts_edt(1,j),5,1) 
                     call sbit(ts_edt(2,j),5,1) 
                     call sbit(ts_edt(3,j),5,1) 
                     exit
                  endif
              end do

*             If we have data after this earthquake add the 
*             log or exp term
              if( log_eq(i) ) then
                 np = np + 1
                 pcode(np) = 7
                 pindx(np) = i
                 pflag(np) = 1
                 if( use_constraints ) then
* MOD TAH Compute the distance dependent sigma (same as globk)
                     do j = 1,3
                        apr_con(np,j) = eq_log_sig(j,i)**2 + 
     .                      (eq_log_sig(j+3,i)*(eq_depth(i)/dist)**2)**2
                     enddo
                 endif
                 if( nodata ) then
                    pflag(np) = 0  
                    do j = 1,3
                        apr_con(np,j) = 1.d-12
                    end do
                 endif
              endif
              if( exp_eq(i) ) then
                 np = np + 1
                 pcode(np) = 8
                 pindx(np) = i
                 pflag(np) = 1  
                 if( use_constraints ) then
                     do j = 1,3
                        apr_con(np,j) = eq_exp_sig(j,i)**2 + 
     .                      (eq_exp_sig(j+3,i)*(eq_depth(i)/dist)**2)**2
                     enddo
                 endif
                 if( nodata ) then
                    pflag(np) = 0  
                    do j = 1,3
                        apr_con(np,j) = 1.d-12
                    end do
                 endif
              endif
*
*             See if rate offset
              if( dtr_eq(i) ) then
                 np = np + 1
                 pcode(np) = 9
                 pindx(np) = i
                 pflag(np) = 1  
                 if( nodata ) then
                    pflag(np) = 0  
                    do j = 1,3
                        apr_con(np,j) = 1.d-12
                    end do
                 endif
              endif
          endif
      enddo
*
*     Thats all 
      num_par = np

      return
      end


CTITLE FIT_INIT

      subroutine fit_init( cmp )

      implicit none

*     Routine to initialize the fit norm equations

      include 'tsfit.h' 

      integer*4 cmp   ! Component being initialized (123 for NEU)

      integer*4 i,j  ! Loop counters

      do i = 1, num_par
         bvec(i) = 0.d0
         do j = 1, num_par
            norm_eq(i,j) = 0.d0
         end do
*        Now add contraints to parameters as needed
         norm_eq(i,i) = 1.d-2
         if( apr_con(i,cmp).ne.0 ) then
            norm_eq(i,i) = 1/apr_con(i,cmp)
         endif
      end do

      do i = 1,4
         stats(i) = 0
      end do

      return
      end


CTITLE FIT_PART

      subroutine fit_part( mjd_no )

      implicit none

*     Routine to set the partials for a specific observation

      include 'tsfit.h' 
      include '../includes/const_param.h'

      real*8 mjd_no  ! MJD of the observation being processed.

* npf  -- First parameter number to counted as offset and now
*     optional rate.

      integer*4 i, npf   ! Loop counter


***** Loop over the parameters and get the partials
      cen_mjd = (ts_mjd(1)+ts_mjd(num_ts))/2
      do i = 1, num_par
         apart(i) = 0.d0
      end do
  
      apart(1) = 1.d0
      if( .not.mean_only ) then
          npf = 3
          apart(2) = (mjd_no - cen_mjd)/365.25d0
      else
          npf = 2
      endif

****  Now see which parameters are on
      do i = npf, num_par

*        Check periodic
         if( pcode(i).eq.3 ) then
             apart(i) = cos(2*pi*(mjd_no-51544.0d0)/
     .                      per_per(pindx(i)))
         endif
         if( pcode(i).eq.4 ) then
             apart(i) = sin(2*pi*(mjd_no-51544.0d0)/
     .                      per_per(pindx(i)))
         endif

*        Check break
         if( pcode(i).eq.5 ) then
             if( mjd_no.ge.rn_times(1,pindx(i)) ) then
                apart(i) = 1.d0
             end if
         end if


*        Check for earthquake break
         if( pcode(i).eq.6 ) then
             if( mjd_no.ge.eq_epoch(pindx(i)) ) then
                apart(i) = 1.d0
             end if
         end if

*        Check for log
         if( pcode(i).eq.7 ) then
             if( mjd_no.ge.eq_epoch(pindx(i)) ) then
                 apart(i) = log(1+(mjd_no-eq_epoch(pindx(i)))/
     .                             eq_log_tau(pindx(i)))
             end if
         end if

*        Check for exponent
         if( pcode(i).eq.8 ) then
             if( mjd_no.ge.eq_epoch(pindx(i)) ) then
                 apart(i) = 1-exp(-(mjd_no-eq_epoch(pindx(i)))/
     .                             eq_exp_tau(pindx(i)))
             end if
         end if

*        Check for rate change
         if( pcode(i).eq.9 ) then
             if( mjd_no.ge.eq_epoch(pindx(i)) ) then
                 apart(i) = (mjd_no - eq_epoch(pindx(i)))/365.25d0
             end if
         end if
      end do

      return
      end


CTITLE FIT_ACCUM

      subroutine fit_accum( edits )

      implicit none

*     Routine to accumulate the normal eqations and statistics
*     N = 1; E = 2 ; U = 3

      include 'tsfit.h' 

      integer*4 i, j, l  ! Loop counters
     .,         k     ! Component (NEU)
     .,         no    ! Number of observation.
     .,         npa   ! Actual number of parameters
     .,         ipivot(max_par) ! Pivot elements
     .,         unts    ! Unit number for ts_<site> file output
     .,         trimlen ! Length of string.

      real*8 wgh      ! Weight of obs
     .,      res      ! Pre and postfit residual
     .,      ures     ! Pre and postfit residual without user specified estimates applied
     .,      scale(max_par)  ! Scale vector
     .,      fact     ! Factor for sigma comparison

* MOD TAH 190503: Continued refining of algorithm used to remove parameters that are
*     not well determined.
      real*8 max_chi  ! Maximum chi value for the parameter estimates
                      ! as each component is evaluted.  This is used in test 
                      ! of max_chi in stopping a "high" sigma parameter from 
                      ! being removed.

      logical edits   ! Set true while there are edits 
      logical rempar  ! Set true when parameters are removed
      logical OK(max_par)  ! Set false if parameter to be removed
                      ! due to being too large
      logical apply_est ! set true if user specified this parameter to be included in residual estimates	      
   
****  Loop over components
      edits  = .false.
      rempar = .false.
      stats(4) = 0.d0
      max_chi = 0.d0

      do k = 1, 3   ! NEU

*        Initalize
         call fit_init( k )
          
****     Loop over observations
         do no = 1, num_ts

****        If data point is OK, increment equations and
*           falls in time range
            if( ts_edt(k,no).eq.0  ) then 

****            Form the apartials
                call fit_part(ts_mjd(no))

****            Increment b-vec
                wgh = 1.d0/ts_neu_std(k,no)**2
                stats(4) = stats(4) + 1
                do i = 1, num_par
                   bvec(i) = bvec(i) + 
     .                       apart(i)*(ts_neu(k,no)-ref_neu(k))*wgh
                   do j = 1, num_par
                      norm_eq(i,j) = norm_eq(i,j) +
     .                       apart(i)*apart(j)*wgh
                   enddo
               end do
            end if
          end do

****      Compute mean epoch before inversion
!         mean_mjd(k) = (ts_mjd(1)+ts_mjd(num_ts))/2 + 
!    .                  (norm_eq(2,1)/norm_eq(1,1))*365.25d0
****      Finished accumulating solution; now solve
          call invert_vis(norm_eq, bvec, scale, ipivot, num_par, 
     .                    max_par, 1) 

****      Compute mean epoch after inversion
          mean_mjd(k) = (ts_mjd(1)+ts_mjd(num_ts))/2 + 
     .                  (norm_eq(2,1)/norm_eq(1,1))*365.25d0

****      Now compute the statistics
          stats(4) = 0.d0
          do no = 1, num_ts
             call fit_part(ts_mjd(no))
             wgh = 1.d0/ts_neu_std(k,no)**2
             res = ts_neu(k,no)-ref_neu(k)
             ures = ts_neu(k,no)-ref_neu(k)
             do i = 1,num_par
	     
**** Increment residuals by parameter estimates
	        
                res = res - apart(i)*bvec(i)
		
**** Increment residuals witout user specified parameters estimates (rcode) applied
                apply_est = .true.
                do j = 1, max_ptyp
		  if ( pcode(i) .eq. rcode(j) ) then
		     apply_est = .false.
		  endif
		enddo
		if ( apply_est ) then
                  ures = ures - apart(i)*bvec(i)   
		endif   
             enddo
	     
             ts_neu_res(no,k) = res
             ts_neu_ures(no,k) = ures
             ts_neu_sig(no,k) = ts_neu_std(k,no)
             if( ts_edt(k,no).eq.0 ) then 
                ts_neu_sig(no,k) = ts_neu_std(k,no)
                stats(1) = stats(1) + res*wgh
                stats(2) = stats(2) + wgh
                stats(3) = stats(3) + res**2*wgh
                stats(4) = stats(4) + 1
             else
                ts_neu_sig(no,k) = -ts_neu_std(k,no)
             endif
          end do

****      See if realistic-sigma requested
*         Compute actual number of parameters
          npa = 0
          do i = 1,num_par
             if( pflag(i).gt.0 ) npa = npa + 1
          end do
          if( real_sigma .and. stats(4).gt.min_rsig ) then
              call real_stats(ts_mjd, ts_neu_res(1,k),ts_neu_sig(1,k),
     .             num_ts, sig_scale(k), taufin(k) )
          else
              if( stats(4)-npa.gt. 0 ) then
                  sig_scale(k) = sqrt(stats(3)/(stats(4)-npa))
              else
                  sig_scale(k) = 1.d0
              endif
          endif
*****     If there is not enough data, set scale to 1.d0 (Ver 1.0.3)
          if( min_sigscale.gt.stats(4)-npa ) sig_scale(k) = 1.d0

          if( stats(4)-npa.gt.0 ) then
             wrms(k) = sqrt(stats(4)/stats(2)*
     .                  (stats(3)/(stats(4)-npa)))*1.e3
          else
             wrms(k) = 100.0
          end if

          ts_used(k) = nint(stats(4))

*         Complete the solution
          do i = 1,num_par
            soln(i,k) = bvec(i)*1000.d0
            smet(i,k) = bvec(i)
            solsig(i,k) = sqrt(norm_eq(i,i))*sig_scale(k)*1000.d0
            smtsig(i,k) = sqrt(norm_eq(i,i))*sig_scale(k)
            do j = 1,num_par
               solcov(i,j,k) = norm_eq(i,j)*sig_scale(k)**2
            enddo
          enddo
          if( (stats(4)-npa).gt. 0 ) then
             wn_nrms(k) = sqrt(stats(3)/(stats(4)-npa))
          else
             wn_nrms(k) = 1.d0
          endif

****      See if we need to edit data
          if( nsigma.gt.0 ) then
              call edit_ns( k, edits)
          else
              edits = .false.
          end if

          if( 1.eq.2 ) then 
             write(*,120) num_par, k, num_ts, nint(stats(4))
 120         format('Parms ',i4,' CMP ',i1,' TS All ',i5,' TS Used ',i5)
             write(*,*) 'Sig_scale ',sig_scale(1:3)
             write(*,140) (i,pcode(i),soln(i,k),solsig(i,k),
     .                       i=1,num_par)
 140         format(50('P ',i2,' PCode ',I3,' Est ',F6.2,' +- ',f6.2))
             write(*,160) stats(1)/stats(2)*1000, 
     .           sqrt(stats(3)/(stats(4)-num_par)),
     .           sqrt(stats(4)/stats(2)*(stats(3)/(stats(4)-num_par)))*1.e3
 160         format('Stats ',3F10.2)
          end if

****      Now see if any parameters are not well enough determined
*         based on max_persig and max_eqsig (Ver 1.0.3)
          fact = 1.d0
          if( k.eq.3 ) fact = htscale
          do i = 1,num_par
             OK(i) = .true.
             if( pcode(i).eq. 3 .or. pcode(i).eq.4 ) then
*               Check sigma
                if( solsig(i,k).gt.max_persig*fact .and.
     .              abs(soln(i,k)/solsig(i,k)).lt.max_rmchi ) then
*                   Sigma is too large, mark to remove from estimates
                    OK(i) = .false.
                    if( pcode(i).eq.4 ) OK(i-1) = .false.  ! Turn off cos
                endif
             endif

             if( pcode(i).eq. 7 .or. pcode(i).eq.8 ) then
*               Check sigma
                max_chi = 0 
                do j = 1,k
                   max_chi = max(abs(soln(i,j)/solsig(i,j)),max_chi)
                end do 
                if( solsig(i,k).gt.max_eqsig*fact .and.
     .              max_chi.lt.max_rmchi ) then
*                   Sigma is too large, mark to remove from estimates
                    OK(i) = .false.
                 endif
                 if ( .not. OK(i) ) then
                    write(*,205) i, k, pcode(i), solsig(i,k), 
     .                  max_eqsig*fact, abs(soln(i,k)/solsig(i,k)),
     .                  max_rmchi, max_chi
 205                format('RM PAR ',I3,' NEU ',I1, ' p-CODE ',i2, 
     .                     ' Sig/Tol ',2F12.4,' Chi/Tol/Max '32F12.4)
                 endif
             endif
          enddo
*         Now eliminate parameters
          if( 1.eq.2 ) then
             write(*,210) num_par,(i,OK(i),i=1,num_par)
 210         format('OK Array ',I3,50(' P'I2.2,1x,L1)) 
          endif 
          i = 2
          do while ( i.lt. num_par )
              i = i + 1
              if ( .not.OK(i) ) then
                 rempar = .true.
                 do j = i+1, num_par
                     pcode(j-1) = pcode(j)
                     pindx(j-1) = pindx(j)
                     pflag(j-1) = pflag(j)
                     apr_con(j-1,:) = apr_con(j,:)
                     OK(j-1) = OK(j)
*                    MOD TAH 190503: Move solution as well for testing
*                    max_chi.
                     do l = 1,k-1
                        soln(j-1,l) = soln(j,l)
                        solsig(j-1,l) = solsig(j,l)
                     end do
                 end do
                 num_par = num_par - 1
              end if
          end do 

          if( 1.eq.2 ) then 
             write(*,220) num_par, k, num_ts, nint(stats(4))
 220         format('Parms ',i4,' CMP ',i1,' TS All ',i5,' TS Used ',i5)
             write(*,240) (i,pcode(i),bvec(i)*1000,
     .                       sqrt(norm_eq(i,i))*1000,i=1,num_par)
 240         format(50('P ',i2,' PCode ',I3, ' Est ',F6.2,' +- ',f6.2))
             write(*,260) stats(1)/stats(2)*1000, 
     .           sqrt(stats(3)/(stats(4)-num_par)),
     .           sqrt(stats(4)/stats(2)*(stats(3)/(stats(4)-num_par)))*1.e3
 260         format('Stats ',F10.2)
          endif
     
      end do

****  If no edits, output results
      if( .not.edits .and. rempar ) edits = rempar
      if( .not.edits ) then
*         Generate the names of sites that should be used
*         Here we take into account the jumps and eq_renames
*         that have occurred.
          call gen_pnm
          unts = 0

*         Now reset the sigmas for the parameters than are not
*         determined (use large values)
          do i = 1, num_par
             do k = 1,3
                if( pflag(i).eq.0 ) then
                    solsig(i,k) = 1.0d3
                    smtsig(i,k) = 1.0d0
                    solcov(i,i,k) = 1.0d6
                end if
             end do
          end do

* MOD TAH 141206: Only output the solution if we are not going to
*         the Kalman filter code
          if( trimlen(kfopt).eq.0 ) then
             do k = 1,3
                 call fit_out(k,unts)
             end do
          endif
      endif
        


      return
      end

CTITLE GEN_PNM

      subroutine gen_pnm

      implicit none

*     Routine to generate site names that are needed for output
*     of parameters 

      include 'tsfit.h' 

      integer*4 i, j, k, n
      real*8 dist
      real*8 save_start   ! Save of the start of the segment we
                          ! are dividing.  Needed during move
                          ! of segments.

      character*1 cc      ! Character to replace ?
      integer*4   icc     ! Integer ASCII of cc
      logical     dup     ! Set true is duplicate

****  Start with base name
      n = 1
      pnm_site(n) = ts_code // '_GPS'
      pnm_indx(n) = 1   ! Point to first parameter
c     time_pnm(1,n) = min(min(ts_first,rn_times(1,1)),eq_epoch(1))
      time_pnm(1,n) = 15020.0d0  ! Make 1900/1/1 as min so that all events
                          ! are added (needed for correct name generation).
      time_pnm(2,n) = 88069.0d0  ! Make 2100/1/1 as end (was ts_last) 

****  OK; Start with renames. (Code assumes that renames are
*     time sorted).
      do i = 1, num_rn
*        See if site effected
         if( rn_codes(1,i)(1:4).eq.ts_code .and.
     .      (rn_codes(2,i)(5:8).ne.'_XPS' .and.
     .       rn_codes(2,i)(5:8).ne.'_XCL') ) then
*           Name matches and not a delete of data;
*           end the current site name and make new
*           name
            time_pnm(2,n) = rn_times(1,i)
*           Add next site name
            n = n + 1
            pnm_site(n) = rn_codes(2,i)
            pnm_type(n) = rn_codes(2,i)
            pnm_rn(n) = i
            time_pnm(1,n) = rn_times(1,i)
            time_pnm(2,n) = 88069.0d0  ! Make 2100/1/1 as end (was ts_last)
*           Find the parameter for this offsett
            do j = 2, num_par
               if( pcode(j).eq.5 .and. pindx(j).eq.i ) then
                  pnm_indx(n) = j  ! Points to the parameter for this break
               endif
            end do
        end if
      end do

****  Now update any ? names with an appropropriate name
      cc = '0' 
      do i = 1, n
*        See if we have ? char
         if ( pnm_site(i)(6:6).eq.'?' ) then
*            We need a new name here.  Scan to see what we
*            can use
             dup = .true.
             do while (  dup )
                icc = ichar(cc)+1
                if( icc.eq.58 ) icc = 65
                if( icc.gt.90 ) then
                    call report_stat('FATAL','TSFIT',
     .                   'GEN_PNM',rn_codes(1,i),
     .                   'Too many BREAKS at site',0)
                endif
                cc = char(icc)
                dup = .false.
                do j = 1, n
                   if( pnm_site(j)(6:6).eq.cc ) dup = .true.
                end do
             end do
*            Assign the name
             pnm_site(i)(6:6) = cc(1:1)
             rn_codes(2,pnm_rn(i))(6:6) = cc(1:1)
         end if
      end do  
         
****  After renames we have now added n-1 new
*     names.  Now see what happens as earthquakes happen
*

* MOD TAH 140227: Set the time of the first measurement to the
*     actual time.  At 1900 date had been set to make sure all
*     renames were accounted for even though there may have been
*     no data in the processing.
      time_pnm(1,1) = ts_first

      do i = 1,num_eq
         call eval_dist( ref_xyz, eq_pos(1,i), dist )
         if( dist.le.eq_rad(i)*eq_rad_scale ) then
*            Site is within distance of earthquake.  Insert
*            earthquake effected name into the list of site
*            names
C            do j = 1, n  ! Do loop but n can change during loop
             j = 0
             do while ( j.lt.n )
                j = j + 1
                if( eq_epoch(i).lt.time_pnm(1,j) ) then
*                   Add eq_code to name
                    pnm_site(j)(7:8) = eq_codes(i)(1:2)
*                   See if this EQ may also break the next segment
*                   as well. (Strictly should go to next earthquake
*                   but this should be OK -- problem if two EQ-codes
*                   at same time.  The second one will be missed.
                elseif( eq_epoch(i).lt.time_pnm(2,j) .and.
     .                  eq_epoch(i).gt.time_pnm(1,j) ) then
*                   Earthquake fall before the end of this
*                   rename, so we will need to add a new
*                   name.  Add eq code to all future names
*                   and move up values
*                   Here j is the section to be split with
*                   j becoming the start and j+1 the end of
*                   segment.
*                   Save start time of current segment because
*                   we will over write for move
                    save_start = time_pnm(1,j)
*                   Now j will become j+1 segment, so set 
*                   correct start time (we saved in save_start above).
                    time_pnm(1,j) = eq_epoch(i) 
                    do k = n,j,-1  ! Need to move backward
*                                  ! so not over written
                       pnm_site(k+1) = pnm_site(k)(1:6) //
     .                                 eq_codes(i)(1:2)
                       pnm_type(k+1) = pnm_type(k)
                       time_pnm(1,k+1) = time_pnm(1,k) 
                       time_pnm(2,k+1) = time_pnm(2,k)
                       pnm_indx(k+1) = pnm_indx(k)
                    end do
*                   Fix up times in current j segment. Replace start
*                   time and add eq_time as end.
                    time_pnm(1,j) = save_start
                    time_pnm(2,j) = eq_epoch(i)
                    pnm_type(j+1) = 'EQ ' // eq_codes(i)(1:2)
*                   Find the parameter number for this offset
                    do k = 2, num_par
                         if( pcode(k).eq.6 .and. pindx(k).eq.i ) then
                             pnm_indx(j+1) = k  ! Points to the parameter
                                              ! for this break
                         endif
                    end do

                    n = n + 1
                endif   ! EQ after rename so no change
             enddo
         endif
      enddo
*     Save number of site names
      num_pnm = n

C     do j = 1, n
C        write(*,320) j,   pnm_site(j),  pnm_type(j), pnm_rn(j), 
C    .                time_pnm(1,j),  time_pnm(2,j) 
C320     format('TIMES ',i3,1x,a,1x,a,1x,i6,' Epoch ',2F12.2)
C     end do

*     Now sum up all the NEU offsets due to breaks and earthquakes.
*     These will be used later for extended entries in apriori coordinates
      do j = 1,3
          soff_pnm(j,1) = 0.d0
      end do
      do i = 2, num_pnm
          do j = 1,3
             soff_pnm(j,i) = soff_pnm(j,i-1)+smet(pnm_indx(i),j)
          enddo
      end do

* MOD TAH 161220: Now scan through data to find actual start and stop
*     times
      n = 1
      time_act(1,1) = time_pnm(1,1)
      time_act(2,1) = time_pnm(2,1)
      if( num_pnm.gt.1 ) then 
         do j = 1, num_ts
          
             if ( ts_mjd(j).ge.time_pnm(1,n+1) .and. 
     .            n+1.le. num_pnm ) then
*               We past the time of next data span so save
*               epoch and increment n
                n = n + 1
                if( n.gt.num_pnm ) then
                    write(*,420) num_pnm, j
 420                format('*** ERROR ** Too many values > ',i4,
     .                     ' getting time_act at count ',i6)
                    n = num_pnm
                endif
                time_act(1,n) = ts_mjd(j)
* MOD TAH 200701: Save an end time as well because in some cases
*               we may have only one day (Ridge-crest sequence).
                time_act(2,n) = ts_mjd(j)
             else
                time_act(2,n) = ts_mjd(j)
             endif
         enddo
         time_act(2,num_pnm) = ts_mjd(num_ts)
      endif
*      

      return
      end

CTITLE EDIT_NS

      subroutine edit_ns(k, edits)

      implicit none

      include 'tsfit.h' 

      integer*4 k  ! NEU (1,2,3 component) 

      logical edits   ! Set true when edits are applied

      integer*4 ne,nr     ! Number of edits
     .,         i     ! Loop counter
      real*8 err      ! Abs(residual/sigma scaled by white NRMS)
      logical kbit    ! tests bit



***** Loop over residuals to see how large and remove large ones.
      ne = 0
      nr = 0
      do i = 1, num_ts
         err = abs(ts_neu_res(i,k)/(ts_neu_sig(i,k)*wn_nrms(k)))
         if( ts_edt(k,i).eq.0 .and. err.gt. nsigma ) then
             call sbit(ts_edt(k,i),3,1)
             ne = ne + 1
         elseif( kbit(ts_edt(k,i),3) .and. err.le. nsigma*0.9d0 ) then
             call sbit(ts_edt(k,i),3,0)
             ne = ne + 1
             nr = nr + 1    ! Number of restored values
         endif
      enddo 
      if( ne.gt .0 ) then
         edits = .true.
      end if

****  Thats all
      return
      end

 

CTITLE FIT_OUT

      subroutine fit_out( cmp, unts )

      implicit none

*     Routine to generate output ts_file

      include 'tsfit.h' 
      include '../includes/const_param.h'

      integer*4 cmp  ! NEU (1,2,3 component) 

      integer*4 date(5), dats(5), datf(5)
     .,         i,j, k, l
     .,         ierr, jerr
     .,         trimlen
     .,         ne, un  ! Earthquake number and file unit number
     .,         cnt     ! count for site name 
     .,         lr, lm, ls ! reference frame index
     .,         ps      ! Start parameter number for periodic terms
                        ! (normally 3 but with mean_only starts at 2)

      logical kbit, out

      real*8 per
     .,      sec
     .,      tau
     .,      decyrs
     .,      loc_coord(3)  ! Geod colat, long and height
* MOD TAH 00701: Added upto to lat/long.ht so they match XYZ values.
     .,      loc_update(3) ! Geod colat, long and height based on updated
                           ! XYZ coordinates.
     .,      xyz_coord(3),xyz_coord_var(3,3),xyz_vel(3),xyz_vel_var(3,3)  ! XYZ coordinates and sigmas
     .,      neu_coord(3),neu_coord_var(3,3),neu_vel(3),neu_vel_var(3,3)  ! Final NEU and sigmas
     .,      rot_mat(3,3),temp_covar(36)
     .,      dec_yrs
     .,      llu_coord(3)  ! Lat/long/ht
     .,      dgdt(3)       ! Rates of change in deg/yr and m/yr
     .,      sum_soln(3)  ! Summed solution due to offsets
     .,      osc   ! Scale factor on log to give values at
                   ! outlog_days after event
     .,      chi   ! Chi of jumps
     .,      tend  ! End epoch for PBO Velocity format (either data or time_range(2)

      integer*4 unts  ! Unit number of site dependent output
                      ! If detroot is NONE not output

      logical wrt_pbov    ! Set true if we are to output value (might be skipped
                          ! if no data)
 
      integer*4 zone      ! Zone for UTM coordinates
      character*1 hemi    ! Hemisphere for UTM coordinates

      character*150 fmt   ! Format for vel line outputs
     .,             fmtsml   ! Assigned to fmt when velocities are small
     .,             fmtlrg   ! Assigned to fmt when velocities are large.


      character*5 comp(3)
      character*256 of    ! Output file name
     .,             eqf   ! Earthquake file name
      character*8 name    ! Name of EQ or rename site in OFFSET lines.
      character*64 lab    ! Label with edit types in renames

      data comp / 'North','East ','Up   ' /

* MOD TAH 190629: Format for .pos lines for small and large velocities.
*     Used when pbovel == .true.

      fmtsml = '(1x,a4,1x,a16,1x,i4,4i2.2,"00",1x,F10.4,1x,3F15.5,'//
     .           '2F16.10,1x,F10.5,2x,3F9.5,3F8.5,3F7.3,2x,3F9.5,' //
     .           '3F8.5,3F7.3,1x,i4,4i2.2,"00",1x,i4,4i2.2,"00")'
      fmtlrg = '(1x,a4,1x,a16,1x,i4,4i2.2,"00",1x,F10.4,1x,3F15.5,'//
     .           '2F16.10,1x,F10.5,2x,3F9.3,3F8.4,3F7.3,2x,3F9.3,' //
     .           '3F8.3,3F7.3,1x,i4,4i2.2,"00",1x,i4,4i2.2,"00")'

***** OK: If this is component 1 then open the output file
      k = cmp
      ls = trimlen(solnstr)
      if( cmp.eq.1 .and. detroot(1:trimlen(detroot)).ne.'NONE' ) then
          unts = 201
          call systime( date, sec )
          of = detroot(1:trimlen(detroot)) // solnstr(1:ls) // '.det'
          open(unts,file=of,iostat=ierr,status='unknown')
          call report_error('IOSTAT',ierr,'open',of,1,'FIT_OUT')

***       Ouput header
          write(unts,120) date, nint(sec)
 120      format('------- ',I4,'/',i2.2,'/',i2.2,1x,i2.2,':',i2.2,1x,
     .           I2.2,' -------')
      endif

      if( unts.ne.0 ) then
         write(unts,140) solnstr(1:ls),comp(k)
 140     format('Detrend of ',a,1x,a)

*        Output statistics
         write(unts,160) wrms(k), sig_scale(k), ts_used(k), num_ts
 160     format('WRMS: ',F10.2,' mm NRMS: ',F10.2,' #: ',i5,
     .          ' used of ',i5,' data')
         if( real_sigma ) then
            write(unts,180) wn_nrms(k), sig_scale(k), taufin(k) 
 180        format('Real_Sigma: White Noise NRMS ',F6.2,
     .             ' RS NRMS ',F6.2,' Tau ',F6.1,' days')
         endif

         if( trimlen(kfopt).gt.0 ) then 
            write(unts,190) kfopt, RWvar(cmp)*1.e6
 190        format('Kalman Filter Option: ',a,' RW Var ',F10.2,
     .             ' mm^2/yr')
         endif
   
*        Now loop over the parameters
         do i = 1, num_par
            if( solsig(i,k).gt.9999.99 ) then
                solsig(i,k) = 9999.99d0
                solcov(i,i,k) = 1.d8
                smtsig(i,k) = solsig(i,k)/1000.d0
            end if
            if( pcode(i).eq. 1 ) then
                write(unts,210) soln(i,k), solsig(i,k)
 210            format('Offset',21x,F15.2,' +- ',F7.2,' mm')
            endif
            if( pcode(i).eq. 2 ) then
                write(unts,220) soln(i,k), solsig(i,k)
 220            format('Rate  ',21x,F15.2,' +- ',F7.2,' mm/yr')
            endif
            if( pcode(i).eq. 3 ) then
*               Cos periodic
                per = per_per(pindx(i))
                write(unts,230) per, soln(i,k), solsig(i,k)
 230            format('Cos ',F8.2,' d ',20x,F7.2,' +- ',F7.2,' mm')
            endif

            if( pcode(i).eq. 4 ) then
*               Cos periodic
                per = per_per(pindx(i))
                write(unts,240) per, soln(i,k), solsig(i,k)
 240            format('Sin ',F8.2,' d ',20x,F7.2,' +- ',F7.2,' mm')
            endif

            if( pcode(i).eq.5 ) then
*               Break
                call mjd_to_ymdhms(rn_times(1,pindx(i)),date,sec)
                write(unts,250) date,soln(i,k), solsig(i,k)
 250            format('Break ',I4,4(1x,I2.2),10x,F10.2,
     .                                     ' +- ',F7.2,' mm')
            endif
   
            if( pcode(i).eq.6 ) then
*               Break
                call mjd_to_ymdhms(eq_epoch(pindx(i)),date,sec)
                write(unts,260) date,soln(i,k), solsig(i,k)
 260            format('OffEq ',I4,4(1x,I2.2),10x,F10.2,
     .                                     ' +- ',F7.2,' mm')
            endif

            if( pcode(i).eq.7 ) then
*               Break
                call mjd_to_ymdhms(eq_epoch(pindx(i)),date,sec)
                tau = eq_log_tau(pindx(i))
                write(unts,270) date,tau, soln(i,k), solsig(i,k)
 270            format('Log   ',I4,4(1x,I2.2),' Tau ',F7.1,1x,
     .                  F7.2,' +- ',F7.2,' mm')
            endif

            if( pcode(i).eq.8 ) then
*               Break
                call mjd_to_ymdhms(eq_epoch(pindx(i)),date,sec)
                tau = eq_exp_tau(pindx(i))
                write(unts,280) date,tau, soln(i,k), solsig(i,k)
 280            format('Exp   ',I4,4(1x,I2.2),' Tau ',F7.1,1x,
     .                  F7.2,' +- ',F7.2,' mm')
            endif

         enddo
         write(unts,'(1x)')
      endif

****  See if last component
      if( cmp.eq.3 ) then

***       See if we need to output edits
          if( report_edits  ) then
             if( unts.gt.0 )
     .       write(unts,320)
 320         format('R DELETED DATA',/,
     .              'R  YYYY MM DD HH MN    DecYr     MJD           ',
     .              'dN       +-   F        dE       +-   F',
     .              '        dU       +-   F')
             do i = 1, num_ts
                if( ts_edt(1,i).ne.0 .or. ts_edt(2,i).ne.0 .or.
     .              ts_edt(3,i).ne.0  ) then 
                   call mjd_to_ymdhms(ts_mjd(i),date,sec)
                   call jd_to_decyrs(ts_mjd(i),decyrs) 
                   if( unts.gt.0 )
     .             write(unts,340) date, decyrs,ts_mjd(i), 
     .                  (ts_neu_res(i,j)*1000,  
     .                   abs(ts_neu_sig(i,j)*1000),
     .                   ts_edt(j,i), j = 1,3)
 340               format('R ',1x,I4,4(1x,I2.2),1x,F10.4,1x,F8.2,
     .                    3(1x,F9.1,1x,F8.1,1x,I3))
                   if( unt.gt.0  ) then
*                     See if we need to output (bit 1 or 3 set)
                      out = .false.
                      lab = ''
                      do j = 1,3
                         if( kbit(ts_edt(j,i),1)) then
                            out = .true.
                            if( index(lab,'MAXSIG').eq.0 )
     .                      lab = trim(lab) // '+' // 'MAXSIG' 
                         endif
                         if( kbit(ts_edt(j,i),3)) then
                            out = .true.
                            if( index(lab,'OUTLIER').eq.0 )
     .                      lab = trim(lab) // '+' // 'OUTLIER' 
                         endif
* MOD TAH 190724: Output for day of earthquake too.
                         if( kbit(ts_edt(j,i),5)) then
                            out = .true.
                            if( index(lab,'EQDAY').eq.0 )
     .                      lab = trim(lab) // '+' // 'EQDAY'
                         endif
                      enddo
                      if( out ) 
     .                write(unt, 360) pnm_site(1)(1:4),
     .                   pnm_site(1)(1:4), (date(j),j=1,3),
     .                   (date(j),j=1,3),(ts_edt(j,i),j=1,3),
     .                   lab
 360                  format(' RENAME ',a4,6x,a4,'_XPS',2x,
     .                   I4,1x,I2,1x,I2,' 00 00 ',2x,
     .                   I4,1x,I2,1x,I2,' 24 00 ',1x,
     .                   ' ! Edit flags ',3o4,' ',a)
                   endif

                endif
             end do
*            Now see if we should output breaks with sigmas and chi
             do i = 2, num_par
                if( pcode(i).eq.5 ) then  ! Regular break
                    chi = 0.d0
                    do j = 1, 3
                       chi = chi + (soln(i,j)/solsig(i,j))**2
                    end do
                    call mjd_to_ymdhms(rn_times(1,pindx(i)),dats,sec)
                    call mjd_to_ymdhms(rn_times(2,pindx(i)),date,sec)
                    if( unt.gt.0 )
     .              write(unt, 365) pnm_site(1)(1:4),
     .                   rn_codes(2,pindx(i)), (dats(j),j=1,5),
     .                   (date(j),j=1,5), sqrt(chi), 
     .                   (soln(i,j), solsig(i,j),j=1,3)
 365                format(' RENAME ',a4,6x,a8,2x,
     .                   I4,1x,I2,1x,I2,1x,I2,1x,I2,3x,
     .                   I4,1x,I2,1x,I2,1x,I2,1x,I2,2x,' ! Chi ',
     .                   F12.2,' Break (NEU) ',3(F8.1,' +- ',F7.1,1x),
     .                   ' (mm)' )
                endif
             end do
             
          end if
   
          if( unts.gt.0 ) close(unts)

***       See if we need to output residuals
          unts = 0
          if( resroot(1:trimlen(resroot)).ne.'NONE' ) then
             unts = 205
             call systime( date, sec )

             of=resroot(1:trimlen(resroot))//'.'//solnstr(1:ls)//'.res'
             open(unts,file=of,iostat=ierr,status='unknown')
             call report_error('IOSTAT',ierr,'open',of,1,'RESROOT')

***          Output files header
             write(unts,370) date, nint(sec)
 370         format('------- ',I4,'/',i2.2,'/',i2.2,1x,i2.2,':',i2.2,1x,
     .           I2.2,' -------')
          
             write(unts,375,iostat=ierr) ts_code
 375         format('4-character ID: ',a4)
             write(unts,380,iostat=ierr) ts_full
 380         format('Station name  : ',a16)
             call jd_to_ymdhms(ts_first+1d-6, date,sec)
             write(unts,385,iostat=ierr) date, nint(sec)
 385         format('First Epoch   : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
             call jd_to_ymdhms(ts_last+1d-6, date,sec)
             write(unts,390,iostat=ierr) date, nint(sec)
 390         format('Last Epoch    : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
             call systime( date, sec )
             write(unts,395,iostat=ierr) date, nint(sec)
 395         format('Release Date  : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
             lr = trimlen(reference_frame) 
             if( lr.gt.0 ) then
                lm = 1
                do while ( lm.lt.lr .and. 
     .              reference_frame(lm:lm).ne.' ' )
                    lm = lm + 1
                end do
                if( reference_frame(lm:lm).eq.' ' ) lm = lm - 1
                write(unts,400,iostat=ierr)ref_xyz,reference_frame(1:lm)
 400            format('XYZ Reference position : ',3F15.5,' (',a,')')
                write(unts,405,iostat=ierr) (pi/2-ref_llu(1))*180/pi,
     .          ref_llu(2)*180/pi,ref_llu(3),
     .          reference_frame(1:lm)
 405            format('NEU Reference position : ',2F16.10,1x,F10.5,
     .          ' (',a,'/WGS84)')
             else
                write(unts,410,iostat=ierr) ref_xyz
 410            format('XYZ Reference position : ',3F15.5,' (',a,')')
                write(unts,415,iostat=ierr) (pi/2-ref_llu(1))*180/pi,
     .                            ref_llu(2)*180/pi,ref_llu(3)
 415            format('NEU Reference position : ',2F16.10,1x,F10.5)
             endif

*** Write fit summary info to header of the .res file
             do k = 1, 3
               if( .not.mean_only ) then
                   write(unts,420) comp(k),soln(2,k),solsig(2,k),
     .                 wrms(k), sig_scale(k),
     .                 ref_neu(k)+soln(1,k)/1000,solsig(1,k)/1000,
     .                 ts_used(k),(ts_mjd(num_ts)-ts_mjd(1))/365.25d0,
     .                 decyrs
               else
                   write(unts,420) comp(k),0.d0,0.d0,
     .                 wrms(k), sig_scale(k),
     .                 ref_neu(k)+soln(1,k)/1000,solsig(1,k)/1000,
     .                 ts_used(k),(ts_mjd(num_ts)-ts_mjd(1))/365.25d0,
     .                 decyrs
               endif 
 420           format(a5,' stats: vel= '
     .         f7.2,' +- ',f5.2,' mm/yr wrms= ',
     .         f6.1,' mm nrms= ',f6.2,' len= '
     .         f15.4,' +- ',f6.4,' m ',' #= ',i4,' dur= ',
     .         f5.2,' yrs mean= ',f7.2,' yr')
             end do

*** Write out parameter fit estimates to the residual file
             write(unts,422) 
 422         format('Parameter Estimates : ',11x,'N          Sig      E',
     .              '          Sig      U          Sig')
             do i = 1, num_par
               if( pcode(i).eq. 1 ) then
                  write(unts,425) (soln(i,k), solsig(i,k), k=1,3)
 425              format('Offsets : ',17x,3(1x,F8.2,' +- ',F7.2),' mm')
               endif
               if( pcode(i).eq. 2 ) then
                 write(unts,430) (soln(i,k), solsig(i,k), k=1,3)
 430             format('Rates :  ',18x,3(1x,F8.2,' +- ',F7.2),' mm/yr')
               endif
               if( pcode(i).eq. 3 ) then
*                Cos periodic
                 per = per_per(pindx(i))
                 write(unts,435) per, (soln(i,k), solsig(i,k), k=1,3)
 435             format('Periodic : Cos ',F8.2,' d ',1x,3(1x,F8.2,
     .           ' +- ',F7.2),' mm')
               endif

               if( pcode(i).eq. 4 ) then
*                Sin periodic
                 per = per_per(pindx(i))
                 write(unts,440) per, (soln(i,k), solsig(i,k), k=1,3)
 440             format('Periodic : Sin ',F8.2,' d ',1x,3(1x,F8.2,
     .           ' +- ',F7.2),' mm')
               endif

               if( pcode(i).eq.5 ) then
*                Break
                 chi = 0.d0
                 do k = 1, 3
                    chi = chi + (soln(i,k)/solsig(i,k))**2
                 end do
                 call mjd_to_ymdhms(rn_times(1,pindx(i)),dats,sec)
c                 call mjd_to_ymdhms(rn_times(2,pindx(i)),date,sec)
                 write(unts, 445) dats,nint(sec),
     .                   (soln(i,k), solsig(i,k),k=1,3),sqrt(chi)
 445             format('Break : ',1x,I4,2(I2.2),1x,3(I2.2),3x,
     .                  3(1x,F8.2,' +- ',F7.2),' mm',
     .                  ' ! Chi ',F12.2)
               endif
     
               if( pcode(i).eq.6 ) then
*                Break
                 chi = 0.d0
                 do k = 1, 3
                    chi = chi + (soln(i,k)/solsig(i,k))**2
                 end do
                 call mjd_to_ymdhms(eq_epoch(pindx(i)),date,sec)
                 write(unts,455) date,nint(sec),(soln(i,k), 
     .           solsig(i,k), k=1,3),sqrt(chi)
 455             format('OffEq : ',1x,I4,2(I2.2),1x,3(I2.2),3x,
     .                  3(1x,F8.2,' +- ',F7.2),' mm',
     .                  ' ! Chi ',F12.2)                  
               endif
 
               if( pcode(i).eq.7 ) then
*                 Break
                  call mjd_to_ymdhms(eq_epoch(pindx(i)),date,sec)
                  tau = eq_log_tau(pindx(i))
                  write(unts,460) date,nint(sec),
     .            (soln(i,k),solsig(i,k),k=1,3),tau
 460              format('EqLog : ',1x,I4,2(I2.2),1x,3(I2.2),
     .            3x,3(1x,F8.2,' +- ',F7.2),' mm ! Tau ',F7.1,1x,'d')
               endif

               if( pcode(i).eq.8 ) then
*                 Break
                  call mjd_to_ymdhms(eq_epoch(pindx(i)),date,sec)
                  tau = eq_exp_tau(pindx(i))
                  write(unts,465) date,nint(sec),
     .            (soln(i,k),solsig(i,k),k=1,3),tau
 465              format('EqExp : ',1x,I4,2(I2.2),1x,3(I2.2),
     .            3x,3(1x,F8.2,' +- ',F7.2),' mm ! Tau ',F7.1,1x,'d')
               endif

             end do
             
             write(unts,470)
 470         format('# YYYYMMDD HHMNSC    DecYr     MJD              ',
     .              'N         E         U        ',
     .              'dN       +-    F       dE       +-    F        ',
     .              'dU       +-   F')

***          Output residuals
             do i = 1, num_ts
                call mjd_to_ymdhms(ts_mjd(i),date,sec)
                call jd_to_decyrs(ts_mjd(i),decyrs) 
                write(unts,495) date, nint(sec), decyrs, ts_mjd(i), 
     .                  ((ts_neu(k,i)-ref_neu(k))*1000.d0, k=1,3),
     .                   (ts_neu_ures(i,j)*1000.d0,  
     .                   abs(ts_neu_sig(i,j)*1000.d0),
     .                   ts_edt(j,i), j = 1,3)
 495            format(' ',1x,I4,2(I2.2),1x,3(I2.2),1x,F10.4,1x,F11.4,
     .                    3(1x,F9.1),3(1x,F9.1,1x,F8.1,1x,I3))
             end do
          endif

          if( unts.gt.0 ) close(unts)

****      See if we are output velocity file
          if( unv.gt.0 .and. .not.mean_only ) then
              if( solsig(2,1)+solsig(2,2)+solsig(2,3).lt.100 )
     .        write(unv,500) ref_llu(2)*180/pi, 
     .           (pi/2-ref_llu(1))*180/pi,
     .           soln(2,2), soln(2,1),soln(2,2), soln(2,1), 
     .           solsig(2,2), solsig(2,1), 0.00, 
     .           soln(2,3), soln(2,3), solsig(2,3),pnm_site(1)
 500          format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .               3(1x,f7.2), 1x,a8,' ')
          endif

****      See if apriori coordimate file
          if( una.gt.0 ) then
              if( .not.rf_out ) then
*                Reference frame is only known after first
*                time series is read, so wait til to write line
                 write(una,510) reference_frame
 510             format('+REFERENCE_FRAME',1x,a,/,'#')
                 rf_out = .true.
              end if

*             Compute and store final XYZ coordinates and velocities
              do i = 1,3
                 neu_coord(i) = ref_neu(i) + smet(1,i)
                 if( .not. mean_only ) then
                     neu_vel(i) = smet(2,i)
                 else
                     neu_vel(i) = 0.0d0
                 endif

		 do j = 1,3
*                  Extract the NEU variances for rotation to XYZ		 
		   if ( i .eq. j ) then
		     neu_coord_var(i,i) = smtsig(1,i)**2
                     if( .not. mean_only ) then
		         neu_vel_var(i,i) = smtsig(2,i)**2
                     else
                         neu_vel_var(i,i) = 0.d0
                     endif 
		   else
		     neu_coord_var(i,j) = 0
		     neu_vel_var(i,j) = 0
		   endif
		 end do		   
              end do
	      
*             Convert NEU coordinates to geodetic and XYZ
* MOD TAH 191217: Added routine to account for possible reference latitude
*             change affecting the estimated neu_coord.
*             call geod_to_loc(neu_coord,loc_coord,xyz_coord)
              call convert_gtol( ref_llu, neu_coord,loc_coord,xyz_coord)

*             Convert NEU velocities to XYZ
              call rotate_geod(neu_vel,xyz_vel,'NEU','XYZ',xyz_coord,
     .            loc_coord,rot_mat)

*             Convert NEU coordinate uncertainties to XYZ  uncertainties
              call var_comp(rot_mat, neu_coord_var, xyz_coord_var,
     .                      temp_covar, 3,3,1)
     
*             Convert NEU velocity uncertainties to XYZ  uncertainties
              call var_comp(rot_mat, neu_vel_var, xyz_vel_var,
     .                      temp_covar, 3,3,1)

*             OK Output apr line for this site
              call jd_to_decyrs(cen_mjd, dec_yrs)
              write(una,520) pnm_site(1), xyz_coord, xyz_vel, dec_yrs,
     .  	    sqrt(abs(xyz_coord_var(1,1))),
     .  	    sqrt(abs(xyz_coord_var(2,2))),
     .  	    sqrt(abs(xyz_coord_var(3,3))),
     .  	    sqrt(abs(xyz_vel_var(1,1))),
     .  	    sqrt(abs(xyz_vel_var(2,2))),
     .  	    sqrt(abs(xyz_vel_var(3,3))),
     .  	    xyz_vel_var(1,2),
     .  	    xyz_vel_var(1,3),
     .  	    xyz_vel_var(2,3)

 520          format(1x,a8,1x,3(f14.5,1x),3(f10.5,1x),f9.4,
     .                1x,  ' ! XYZ +- ',6(f8.5,1x),3(f14.11,1x))

* MOD TAH 091008: Added geodetic coordinate output for database population
              call geod_to_geod(xyz_coord, llu_coord, 'XYZ', 'GEOD',
     .              'WGS84','WGS84',zone,hemi)
*             Convert the rates to angular rates
              dgdt(1) = neu_vel(1)/Earth_rad*180/pi
              dgdt(2) = neu_vel(2)/Earth_rad*180/pi/
     .                  sin(llu_coord(1))
              dgdt(3) = neu_vel(3)
              write(una,530) pnm_site(1), 
     .             (pi/2-llu_coord(1))*180/pi,llu_coord(2)*180/pi,
     .             llu_coord(3), dgdt, dec_yrs, 
     .            (smtsig(1,i),i=1,3),(smtsig(2,i),i=1,3)
 530          format('GEOD',2x,a8,2F15.9,F10.4,1x,2F15.9,F10.5,1x,F9.4,
     .               1x,' ! NEU +- ', 3(f8.5,1x), 1x, 3(f8.5,1x) )


****          See if need extended lines for jumps and EQs.
*  MOD TAH 090925: Removed OFFSET extended.  Not really needed and
*             has errors when EQ's intersperced with renames offsets
C             do i = 3,num_par
****             See if offsets from rename or earthquake
C                if( pcode(i).eq.5 .or. pcode(i).eq. 6) then
*                    Now sum up the previous values so that this is 
*                    total offset
C                    if( pcode(i).eq.5 ) then
C                       call accum_offsets(i, sum_soln, cnt, 'BREAK') 
C                       call jd_to_ymdhms(rn_times(1,pindx(i)),date,sec)
C                       name = rn_codes(2,pindx(i))
C                    else
C                       call accum_offsets(i, sum_soln, cnt, 'EQOFF')
C                       call jd_to_ymdhms(eq_epoch(pindx(i)),date,sec)
C                       name = eq_codes(pindx(i))
C                    endif
* MOD TAH 090928: Replaced with prior summation of terms.  Differences
*                    at each epoch are output here
* MOD TAH 140225: Commented out OFFSET term in output since not needed
*             when renames are used and it reduces number of extended terms.
              do i = 2, num_pnm
                 do j = 1,3
                     sum_soln(j) = soff_pnm(j,i)-soff_pnm(j,i-1)
                 end do
                 call jd_to_ymdhms(time_pnm(1,i),date,sec)
                 write(una,540) pnm_site(1), 'OFFSET', date, 0.0, 
     .              sum_soln(1), 0.0, sum_soln(2), 0.0, 
     .              sum_soln(3), 0.0, (solsig(pnm_indx(i),j),j=1,3),
     .                                 pnm_type(i)
 540             format('!EXTENDED ',a8,1x,a8,1x,i5,4i3,F8.3,
     .                6F9.5,1x,' !  m +- ',3F8.2,' mm, Type ',a)
              end do

*             Convert NEU coordinates to geodetic and XYZ coords
C             call geod_to_loc(neu_coord,loc_coord,xyz_coord)
* MOD TAH 090928: Compute and output the position by name
              do i = 2, num_pnm
*                Compute final XYZ coordinates
                 do j = 1,3
                    neu_coord(j) = ref_neu(j) + smet(1,j) + 
     .                                          soff_pnm(j,i)
                 end do

*                Convert NEU coordinates to geodetic and XYZ coords
                 call geod_to_loc(neu_coord,loc_coord,xyz_coord)
* MOD TAH 191217: Added routine to account for possible reference latitude
*                change affecting the estimated neu_coord.
*                call geod_to_loc(neu_coord,loc_coord,xyz_coord)
                 call convert_gtol( ref_llu, neu_coord,loc_coord,
     .                                                 xyz_coord)
 

*                OK Output line
                 write(una,580) pnm_site(i), xyz_coord, xyz_vel,
     .               dec_yrs,(smtsig(pnm_indx(i),j),j=1,3),
     .                       (smtsig(2,j),j=1,3) 
 580             format(1x,a8,1x,3(f14.5,1x),3(f10.5,1x),f9.4,
     .               1x,' ! NEU +- ',3(f7.4,1x), 1x, 3(f7.4,1x) )
               
              end do


****          See if periodic: Get start parameter number for mean_only case
              if( .not. mean_only ) then
                 ps = 3
              else
                 ps = 2
              endif 
              do i = ps,num_par
                 if( pcode(i).eq.3 ) then
*                    Output lines for each renamed site so that all
*                    have the same periodic signals.  Output periodics
*                    for all site names
                     do j = 1,num_pnm
                         write(una,620) pnm_site(j), 'PERIODIC',
     .                       2000,1,1,0,0, per_per(pindx(i)),
     .                       smet(i,1), smet(i+1,1), smet(i,2),
     .                       smet(i+1,2), smet(i,3), smet(i+1,3),
     .                      (solsig(i,l),solsig(i+1,l),l=1,3)
 620                     format(' EXTENDED ',a8,1x,a8,1x,i5,4i3,
     .                       F8.3, 6F9.5,' ! m +- ',
     .                       6F6.2,' mm')
                     end do
                 endif
              end do


****          Now do log/exp terms: These are EXTENDED lines
              do i = ps,num_par 
                 if( pcode(i).eq.7 .or. pcode(i).eq.8 ) then 
*                   Convert time and date
                    call jd_to_ymdhms(eq_epoch(pindx(i)),date,sec)

*                   Log term.  Get the site name.  Find first occurance
*                   of this code and use all future names
                    cnt = 0
                    do j = 1, num_pnm
                       if( pnm_site(j)(7:8).eq.
     .                     eq_codes(pindx(i))(1:2) .and. cnt.eq.0 ) then
                           cnt = j
                       endif
                    end do

                    do j = cnt,num_pnm
                       if( pcode(i).eq.7 )
     .                 write(una,640) pnm_site(j),'LOG', date,
     .                      eq_log_tau(pindx(i)), (smet(i,l),l=1,3),
     .                     (solsig(i,l), l=1,3)
                       if( pcode(i).eq.8 )
     .                 write(una,640) pnm_site(j),'EXP', date,
     .                      eq_exp_tau(pindx(i)), (smet(i,l),l=1,3),
     .                     (solsig(i,l), l=1,3)
640                   format(' EXTENDED ',a8,1x,a8,1x,i5,4i3,1x,F8.3,
     .                    3F9.5,1x,'  !  m; +- ',3F8.2,' mm')
                    enddo

                 end if

              end do
          end if

*         See if Log output have been requested
          if( eq_out ) then
*             Loop over earthquakes and see if we are writing files
*             yet
              do i = 3, num_par
                 if( pcode(i).eq.6 .and. 
     .               apr_con(i,1).ne.1.d-12 ) then  ! Earthquake offset
*                    See which earthquake and if file is open
                     ne = pindx(i)
                     if( une(1,ne).eq.0 ) then
*                        Need to open file
                         une(1,ne) = 300 + (ne-1)*4
                         un = une(1,ne)
                         eqf = eqo_root(1:trimlen(eqo_root)) //
     .                         eq_codes(ne)(1:2) // '_Off.off'
                         open(un,file=eqf,iostat=jerr, 
     .                        status='unknown')
                         call report_error('IOSTAT',jerr,'open',
     .                        eqf,0,'tsfit/fit_out/eq_off')
                         write(un,820)eq_codes(ne)(1:2), 
     .                         cmdfile(1:trimlen(cmdfile)), 
     .                         eq_rad(ne)/1000, eq_rad_scale

 820                     format('* TSFIT Earthquake offset for EQ ',a2,
     .                          ' from command file ',a,/,
     .                          '* EQ Radius ',F8.2,' km, Scale ',F6.2)

                         write(un,840)
 840                     format('*  Long         Lat        Eoff    ',
     .                        'Noff    dEo     dNo    E +-    N +-   ',
     .                        ' Rne      Hoff     dHo    H +-  Site',/,
     .                        '*  deg          deg        mm      ',
     .                        'mm      mm      mm      mm      mm   ',
     .                        '            mm      mm      mm  ')
                     endif
*                    Output line
                     un = une(1,ne)
                     if( solsig(i,1)+solsig(i,2)+solsig(i,3).lt.
     .                  100 .and.
     .                  solsig(i,1)+solsig(i,2)+solsig(i,3).gt.
     .                  0.d0 )
     .               write(un,860) ref_llu(2)*180/pi, 
     .                    (pi/2-ref_llu(1))*180/pi,
     .                     soln(i,2), soln(i,1),soln(i,2), soln(i,1), 
     .                     solsig(i,2), solsig(i,1), 0.00, 
     .                     soln(i,3), soln(i,3), solsig(i,3),
     .                     pnm_site(1)(1:6) // eq_codes(ne)(1:2)
 860                 format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .                      3(1x,f7.2), 1x,a8,' ')
                 endif
              end do

****          Now see if log output
              do i = 3, num_par
                 if( pcode(i).eq.7  .and. 
     .               apr_con(i,1).ne.1.d-12 ) then  ! Earthquake log term
*                    See which earthquake and if file is open
                     ne = pindx(i)
                     if( une(2,ne).eq.0 ) then
*                        Need to open file
                         une(2,ne) = 301 + (ne-1)*4
                         un = une(2,ne)
*                        Compute scale for output at time after
*                        earthquake. 
C                        if( outlog_days.gt.0 ) then
C                          osc = log(1+outlog_days/eq_log_tau(ne))
C                        else
C                          osc = 1.d0
C                        endif
                         osc = 1
                         eqf = eqo_root(1:trimlen(eqo_root)) //
     .                         eq_codes(ne)(1:2) // '_Log.ln'
                         open(un,file=eqf,iostat=jerr, 
     .                        status='unknown')
                         call report_error('IOSTAT',jerr,'open',
     .                        eqf,0,'tsfit/fit_out/eq_off')
                         write(un,920)eq_codes(ne)(1:2), outlog_days,
     .                         cmdfile(1:trimlen(cmdfile)),
     .                         eq_log_tau(ne), eq_rad(ne)/1000, 
     .                         eq_rad_scale
 920                     format('* TSFIT Earthquake Log for EQ ',a2,
     .                          1x,F8.2,' days after ',
     .                          ' from command file ',a,/,
     .                          '* Log tau ',F8.2,' days, Radius ',
     .                          F8.2,' km, Scale ',F6.2)
                         write(un,940)
 940                     format('*  Long         Lat        Elog    ',
     .                        'Nlog    dEl     dNl    E +-    N +-    ',
     .                        'Rne      Hlog     dHl    H +-  Site',/,
     .                        '*  deg          deg        mm      mm',
     .                        '      mm      mm      mm      mm     ',
     .                        '          mm      mm      mm  ')
                     endif
*                    Output line
                     un = une(2,ne)
                     if( solsig(i,1)+solsig(i,2)+solsig(i,3).lt.
     .                  100/osc .and.
     .                  solsig(i,1)+solsig(i,2)+solsig(i,3).gt.0 )
     .               write(un,960) ref_llu(2)*180/pi, 
     .                    (pi/2-ref_llu(1))*180/pi,
     .                     soln(i,2)*osc, soln(i,1)*osc,soln(i,2)*osc, 
     .                     soln(i,1)*osc, solsig(i,2)*osc,  
     .                     solsig(i,1)*osc, 0.00,soln(i,3)*osc, 
     .                     soln(i,3)*osc, solsig(i,3)*osc, 
     .                     pnm_site(1)(1:6) // eq_codes(ne)(1:2)
 960                 format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .                      3(1x,f7.2), 1x,a8,' ')
                 endif
              end do

****          Now see if Exp output
              do i = 3, num_par
                 if( pcode(i).eq.8  .and. 
     .               apr_con(i,1).ne.1.d-12 ) then  ! Earthquake log term
*                    See which earthquake and if file is open
                     ne = pindx(i)
                     if( une(3,ne).eq.0 ) then
*                        Need to open file
                         une(3,ne) = 302 + (ne-1)*4
                         un = une(3,ne)
*                        Compute scale for output at time after
*                        earthquake. 
C                        if( outlog_days.gt.0 ) then
C                          osc = 1-exp(-outlog/eq_exp_tau(ne))
C                        else
C                          osc = 1.d0
C                        endif
                         osc = 1
                         eqf = eqo_root(1:trimlen(eqo_root)) //
     .                         eq_codes(ne)(1:2) // '_Exp.exp'
                         open(un,file=eqf,iostat=jerr, 
     .                        status='unknown')
                         call report_error('IOSTAT',jerr,'open',
     .                        eqf,0,'tsfit/fit_out/eq_off')
                         write(un,1020)eq_codes(ne)(1:2), outlog_days,
     .                         cmdfile(1:trimlen(cmdfile)),
     .                         eq_exp_tau(ne), eq_rad(ne)/1000,
     .                         eq_rad_scale
1020                     format('* TSFIT Earthquake Log for EQ ',a2,
     .                          1x,F8.2,' days after ',
     .                          ' from command file ',a,/,
     .                          '* Exp tau ',F8.2,' days, Radius ',
     .                          F8.2,' km, Scale ',F6.2)
                         write(un,1040)
1040                     format('*  Long         Lat        Elog    ',
     .                        'Nlog    dEl     dNl    E +-    N +-    ',
     .                        'Rne      Hlog     dHl    H +-  Site',/,
     .                        '*  deg          deg        mm      mm',
     .                        '      mm      mm      mm      mm     ',
     .                        '          mm      mm      mm  ')
                     endif
*                    Output line
                     un = une(3,ne)
                     if( solsig(i,1)+solsig(i,2)+solsig(i,3).lt.
     .                  100/osc )
     .               write(un,1060) ref_llu(2)*180/pi, 
     .                    (pi/2-ref_llu(1))*180/pi,
     .                     soln(i,2)*osc, soln(i,1)*osc,soln(i,2)*osc, 
     .                     soln(i,1)*osc, solsig(i,2)*osc,  
     .                     solsig(i,1)*osc, 0.00,soln(i,3)*osc, 
     .                     soln(i,3)*osc, solsig(i,3)*osc, 
     .                     pnm_site(1)(1:6) // eq_codes(ne)(1:2)
1060                 format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .                      3(1x,f7.2), 1x,a8,' ')
                 endif
              end do

*             See if total accumulated motion
              if( outlog_days.gt.0 ) then 
                  call total_accum
              end if
          end if

****      See if we are output velocity file
          if( unf.gt.0 ) then
* MOD TAH 160407: Increased sig-limit to 1000 (should still fit format)
              if( solsig(1,1)+solsig(1,2)+solsig(1,3).lt.1000 )
     .        write(unf,1110) ref_llu(2)*180/pi, 
     .           (pi/2-ref_llu(1))*180/pi,
     .           soln(1,2), soln(1,1),soln(1,2), soln(1,1), 
     .           solsig(1,2), solsig(1,1), 0.00, 
     .           soln(1,3), soln(1,3), solsig(1,3),pnm_site(1)
1110          format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .               3(1x,f7.2), 1x,a8,' ')
          endif


******    Write PBO-format velocity file. Added by MAF (2014-10-23, MIT)

          if( trimlen(out_pbovel).gt.0 ) then
            unp = 206
            open(unp,file=out_pbovel,iostat=ierr, status='unknown')
            call report_error('IOSTAT',ierr,'open',out_pbovel,0,
     .          'fit_out')
            if( ierr.eq.0 ) then

              if( .not.pbov_header ) then
*               Reference frame is only known after first
*               time series is read, so wait til to write line
                write(unp,511) sumfile(1:trimlen(sumfile)),
     .                         reference_frame
* Changed label from PBO to GAGE
 511            format('GAGE Velocity file from ',a,
     .                 ' Reference Frame : ',a)

*               Write standard PBO-format velocity file header
                call systime( date, sec )
                write(unp,417) date, int(sec)
 417            format('Format Version: 1.1.0',/,'Release Date  : ',i4,
     .                 i2.2,i2.2,i2.2,i2.2,i2.2,/,
     .               'Start Field Description',/,'Dot#        4-'
     .               'character identifier for a given station',/,'Name'
     .               '        16-character station name',/,'Ref_epoch  '
     .               ' Date and time at which the station position is '
     .               'as given in ref_XYZ and ref_NEU. Format is '
     .               'YYYYMMDDhhmmss.',/,'Ref_jday    Reference epoch, '
     .               'represented as Modified Julian day',/,'Ref_X     '
     .               '  Reference X coordinate at Ref_epoch, meters',/,
     .               'Ref_Y       Reference Y coordinate at Ref_epoch, '
     .               'meters',/,'Ref_Z       Reference Z coordinate at '
     .               'Ref_epoch, meters',/,'Ref_Nlat    Reference North'
     .               ' latitude WGS-84 ellipsoid, decimal degrees',/,
     .               'Ref_Elong   Reference East Longitude WGS-84 '
     .               'ellipsoid, decimal degrees',/,'Ref_Up      '
     .               'Reference Height  WGS-84 ellipsoid, meters ',/,
     .               'dX/dt       X component of station '
     .               'velocity, meters/yr',/,'dY/dt       Y component '
     .               'of station velocity, meters/yr',/,'dZ/dt       Z '
     .               'component of station velocity, meters/yr',/,'SXd '
     .               '        Standard deviation of X velocity, '
     .               'meters/yr',/,'SYd         Standard deviation of Y'
     .               ' velocity,  '
     .               'meters/yr',/,'SZd         Standard deviation of Z'
     .               ' velocity,  meters/yr',/,'Rxy         Correlation'
     .               ' of X and Y velocity',/,'Rxz         '
     .               'Correlation of X and Z velocity',/,
     .               'Ryz         Correlation of Y and Z velocity, ',/,
     .               'dN/dt       North component of '
     .               'station velocity, meters/yr',/,'dE/dt       East '
     .               'component of station velocity, meters/yr',/,
     .               'dU/dt       Vertical component of station '
     .               'velocity, meters/yr',/,'SNd         Standard '
     .               'deviation of North velocity, meters/yr',/,'SEd   '
     .               '      Standard deviation of East velocity, '
     .               'meters/yr',/,'SUd         Standard deviation of '
     .               'vertical velocity, meters/yr',/,'Rne         '
     .               'Correlation of North and East velocity',/,
     .               'Rnu         Correlation of North '
     .               'and vertical velocity',/,'Reu         '
     .               'Correlation of East and vertical velocity',/,
     .               'first_epoch Epoch of first data '
     .               'used to derive the station velocity, in the same '
     .               'format as ref_epoch.',/,'last_epoch  Epoch of '
     .               'last data used to derive the station velocity, in'
     .               ' the same format as ref_epoch.',/,'End Field '
     .               'Description ',/,'*Dot#     Name           '
     .               'Ref_epoch      Ref_jday      Ref_X          Ref_Y'
     .               '           Ref_Z         Ref_Nlat        '
     .               'Ref_Elong       Ref_Up     dX/dt    dY/dt   '
     .               'dZ/dt    SXd     SYd     SZd    Rxy     Rxz    '
     .               'Rzy      dN/dt     dE/dt    dU/dt   SNd     SEd  '
     .               '   SUd     Rne    Rnu    Reu   first_epoch    '
     .               'last_epoch ')
                pbov_header = .true.
              end if

*             Compute and store final XYZ coordinates and velocities
              do i = 1,3
                neu_coord(i) = ref_neu(i) + smet(1,i)
                neu_vel(i) = smet(2,i)
                do j = 1,3
*                 Extract the NEU variances for rotation to XYZ                 
                  if ( i .eq. j ) then
                    neu_coord_var(i,j) = smtsig(1,i)**2
                    neu_vel_var(i,j) = smtsig(2,i)**2
                  else
                    neu_coord_var(i,j) = 0
                    neu_vel_var(i,j) = 0
                  endif
                end do         
              end do
       
*             Convert NEU coordinates to geodetic and XYZ
* MOD TAH 191217: Added routine to account for possible reference latitude
*             change affecting the estimated neu_coord.
*             call geod_to_loc(neu_coord,loc_coord,xyz_coord)
              call convert_gtol( ref_llu, neu_coord,loc_coord,xyz_coord)

*             Convert NEU velocities to XYZ
c              call rotate_geod(neu_vel,xyz_vel,'NEU','XYZ',ref_xyz,
c     .                         loc_coord,rot_mat)
              call rotate_geod(neu_vel,xyz_vel,'NEU','XYZ',xyz_coord,
     .                         loc_coord,rot_mat)

*             Convert NEU velocity uncertainties to XYZ  uncertainties
              call var_comp(rot_mat, neu_vel_var, xyz_vel_var,
     .                      temp_covar, 3,3,1)

*             Write velocity record
              call mjd_to_ymdhms(cen_mjd,date,sec)
              call mjd_to_ymdhms(ts_first,dats,sec)
              if ( num_pnm.eq.1 ) then  ! No renames
                call mjd_to_ymdhms(ts_last,datf,sec)
              else
C               call mjd_to_ymdhms(time_pnm(2,1),datf,sec)
                call mjd_to_ymdhms(time_act(2,1),datf,sec)
              end if

* MOD TAH 141119: Make sure end time matchs data processed (replaced
*             ts_last with tend in code below).
              tend = ts_last
              if( ts_last.gt.time_range(2) ) tend = time_range(2)

* MOD TAH 141119: Before output line, make sure the start time is 
*             before the end time (break before data start causes this
*             problem).
              wrt_pbov = .true.
              if( time_pnm(2,1).lt.ts_first ) then
                  wrt_pbov = .false.
                  if( num_pnm.gt.1 ) time_pnm(1,2) = ts_first
              endif

* MOD TAH 190629: Change velocity format if values are large (>999 or <-99)
              fmt = fmtsml 
              do i = 1,3
                 if( xyz_vel(i).gt.99. or. xyz_vel(i).lt.-9.99 ) 
     .              fmt = fmtlrg
                 if( neu_vel(i).gt.99. or. neu_vel(i).lt.-9.99 ) 
     .              fmt = fmtlrg
              end do 
              
              if( wrt_pbov )
     .        write(unp,fmt) ts_code,ts_full(1:16),date,cen_mjd,
     .             xyz_coord,(pi/2-loc_coord(1))*180/pi,
     .             loc_coord(2)*180/pi,loc_coord(3),xyz_vel,
     .             sqrt(abs(xyz_vel_var(1,1))),
     .             sqrt(abs(xyz_vel_var(2,2))),
     .             sqrt(abs(xyz_vel_var(3,3))),
     .             xyz_vel_var(1,2)/sqrt(abs(xyz_vel_var(1,1)
     .             )*abs(xyz_vel_var(2,2))),
     .             xyz_vel_var(1,3)/sqrt(abs(xyz_vel_var(1,1)
     .             )*abs(xyz_vel_var(3,3))),
     .             xyz_vel_var(2,3)/sqrt(abs(xyz_vel_var(2,2)
     .             )*abs(xyz_vel_var(3,3))),neu_vel,
     .             sqrt(abs(neu_vel_var(1,1))),
     .             sqrt(abs(neu_vel_var(2,2))),
     .             sqrt(abs(neu_vel_var(3,3))),
     .             neu_vel_var(1,2)/(sqrt(abs(neu_vel_var(1,1)))*
     .             sqrt(abs(neu_vel_var(2,2)))),
     .             neu_vel_var(1,3)/(sqrt(abs(neu_vel_var(1,1)))*
     .             sqrt(abs(neu_vel_var(3,3)))),
     .             neu_vel_var(2,3)/(sqrt(abs(neu_vel_var(2,2)))*
     .             sqrt(abs(neu_vel_var(3,3)))), dats, datf
C MOD TAH 190629: Format now variable depending on velocity size.
C 418         format(1x,a4,1x,a16,1x,i4,4i2.2,'00',1x,F10.4,1x,3F15.5,
C    .               2F16.10,1x,F10.5,2x,3F9.5,3F8.5,3F7.3,2x,3F9.5,
C    .               3F8.5,3F7.3,1x,i4,4i2.2,'00',1x,i4,4i2.2,'00')

*             Loop over site renames
              do i = 2,num_pnm

*               Compute final XYZ coordinates
                do j = 1,3
                  neu_coord(j) = ref_neu(j) + smet(1,j) + soff_pnm(j,i)
                end do

*               Convert NEU coordinates to geodetic and XYZ coords
* MOD TAH 191217: Added routine to account for possible reference latitude
*               change affecting the estimated neu_coord.
*               call geod_to_loc(neu_coord,loc_coord,xyz_coord)
                call convert_gtol( ref_llu, neu_coord,loc_coord,
     .                                                xyz_coord)
* MOD TAH 200701: Compute updated geodetic coordinates based on XYZ.
                call XYZ_to_GEOD(rot_mat, xyz_coord, loc_update )

*               Write velocity record
                call mjd_to_ymdhms(cen_mjd,date,sec)
C               call mjd_to_ymdhms(time_pnm(1,i),dats,sec)
                call mjd_to_ymdhms(time_act(1,i),dats,sec)
                if ( i.eq.num_pnm ) then  ! Last site rename
                  call mjd_to_ymdhms(tend,datf,sec)
                else
C                 call mjd_to_ymdhms(time_pnm(2,i),datf,sec)
                  call mjd_to_ymdhms(time_act(2,i),datf,sec)
                end if
                wrt_pbov = .true.
                if( time_pnm(2,i).lt.ts_first ) then
                    wrt_pbov = .false.
                    if( num_pnm.gt.i ) time_pnm(1,i+1) = ts_first
                endif
                if( time_pnm(2,i).gt.tend ) then
*                   Break after end of data so set end time back to actual end
                    call mjd_to_ymdhms(tend,datf,sec)
                endif
*               Make sure we have data 
                if( time_pnm(1,i).gt.tend ) then
*                   Break after end of data, so don't output
                    wrt_pbov = .false.
                endif
 
* MOD TAH Used updated latitude and longitude based on updated XYZ coordinates.
                if( wrt_pbov )
     .          write(unp,fmt) ts_code,ts_full(1:16),date,cen_mjd,
     .                xyz_coord,(pi/2-loc_update(1))*180/pi,
     .                loc_update(2)*180/pi,loc_update(3),xyz_vel,
     .                sqrt(abs(xyz_vel_var(1,1))),
     .                sqrt(abs(xyz_vel_var(2,2))),
     .                sqrt(abs(xyz_vel_var(3,3))),
     .                xyz_vel_var(1,2)/(sqrt(abs(xyz_vel_var(1,1)))*
     .                sqrt(abs(xyz_vel_var(2,2)))),
     .                xyz_vel_var(1,3)/(sqrt(abs(xyz_vel_var(1,1)))*
     .                sqrt(abs(xyz_vel_var(3,3)))),
     .                xyz_vel_var(2,3)/(sqrt(abs(xyz_vel_var(2,2)))*
     .                sqrt(abs(xyz_vel_var(3,3)))),
     .                neu_vel,sqrt(abs(neu_vel_var(1,1))),
     .                sqrt(abs(neu_vel_var(2,2))),
     .                sqrt(abs(neu_vel_var(3,3))),
     .                neu_vel_var(1,2)/(sqrt(abs(neu_vel_var(1,1)))*
     .                sqrt(abs(neu_vel_var(2,2)))),
     .                neu_vel_var(1,3)/(sqrt(abs(neu_vel_var(1,1)))*
     .                sqrt(abs(neu_vel_var(3,3)))),
     .                neu_vel_var(2,3)/(sqrt(abs(neu_vel_var(2,2)))*
     .                sqrt(abs(neu_vel_var(3,3)))),
     .                dats,datf
C419            format(1x,a4,1x,a16,1x,i4,4i2.2,'00',1x,F10.4,1x,3F15.5,
C    .                 2F16.10,1x,F10.5,2x,3F9.5,3F8.5,3F7.3,2x,3F9.5,
C    .                 3F8.5,3F7.3,1x,i4,4i2.2,'00',1x,i4,4i2.2,'00')

C 418           format(1x,a4,1x,a16,1x,i4,4i2.2,'00',1x,F10.4,1x,3F15.5,
C    .               2F16.10,1x,F10.5,2x,3F9.5,3F8.5,3F7.3,2x,3F9.5,
C    .               3F8.5,3F7.3,1x,i4,4i2.2,'00',1x,i4,4i2.2,'00')
              end do

            else
               unp = 0
            endif

          end if
****** END: Added by MAF (2014-10-23, MIT)
                
      endif

****  Now write the summary information in the ENSUM format
* MOD TAH 101207: Set k before call to decyrs
      k = cmp
      call jd_to_decyrs(mean_mjd(k),decyrs)
      if( .not.mean_only ) then 
          if( ts_used(k).gt.0 )
     .    write(uns,100) pnm_site(1),comp(k)(1:1), ts_used(k),
     .         ref_neu(k)+soln(1,k)/1000,solsig(1,k)/1000, 
     .         wrms(k), wn_nrms(k), soln(2,k),solsig(2,k),
     .         wrms(k), sig_scale(k), 
     .        (ts_mjd(num_ts)-ts_mjd(1))/365.25d0, decyrs
      else
          if( ts_used(k).gt.0 .and. .not. resid_sum )
     .    write(uns,100) pnm_site(1),comp(k)(1:1), ts_used(k),
     .         ref_neu(k)+soln(1,k)/1000,solsig(1,k)/1000, 
     .         wrms(k), wn_nrms(k),0.d0 ,0.d0,
     .         wrms(k), sig_scale(k), 
     .        (ts_mjd(num_ts)-ts_mjd(1))/365.25d0, decyrs
          if( ts_used(k).gt.0 .and. resid_sum )
     .    write(uns,110) pnm_site(1),comp(k)(1:1), ts_used(k),
     .         soln(1,k),solsig(1,k), 
     .         wrms(k), wn_nrms(k),0.d0 ,0.d0,
     .         wrms(k), sig_scale(k), 
     .        (ts_mjd(num_ts)-ts_mjd(1))/365.25d0, decyrs
      end if

 100  format(a4'_GPS-',7x,a1,1x,' 1',1x,i4,1x,
     .      f16.5,1x,f7.5,1x,f6.2,1x,f6.2,3x,
     .      f8.2,1x,f8.2,1x,f6.2,1x,f6.2,1x,
     .      f5.2,1x,f7.2,1x,i7,"-",i7)
 110  format(a4'_GPS-',7x,a1,1x,' 1',1x,i4,1x,
     .      f14.2,1x,f7.2,1x,f6.2,1x,f6.2,3x,
     .      f8.2,1x,f8.2,1x,f6.2,1x,f6.2,1x,
     .      f5.2,1x,f7.2,1x,i7,"-",i7)

      return
      end

CTITLE TOTAL_ACCUM

      subroutine total_accum

      implicit none

*     Routine to output total accumulation of log
*     and exp parts of earthquake

      include 'tsfit.h'
      include '../includes/const_param.h'

* LOCAL
      integer*4 i,j,k, c  ! counters
     .,    pardP(3)       ! Parameter numbers for log, exp, drate
     .,    ne             ! EQ number
     .,    nc             ! number of components contributing
     .,    un             ! Unit number for output
     .,    jerr           ! IOSTAT error
     .,    trimlen        ! String length

      real*8 osc    ! Scale factor from coefficien to offset including
                    ! conversion to mm
     .,    comdP(3) ! Partials for conversion
     .,    dP(3)    ! position offsets for NEU
     .,    sig(3)   ! Sigmas of offsets
     .,    cov(3,3) ! Covariance for log, exp, drate
     .,    var      ! Variance and temp value for var_comp

      character*256 eqf  ! File name


****  Loop over earthquakes to see which one applies here
      do i = 1, num_eq

         nc = 0
*        See if we have estimate
         do j = 3, num_par
            if( pcode(j).eq.7 .and. pindx(j).eq.i ) then
*               Log coefficient for this earthquake
                nc = nc + 1
                osc = log(1+outlog_days/eq_log_tau(i))*1000 ! mm
                comdP(nc) = osc
                pardP(nc) = j
            elseif( pcode(j).eq.8 .and. pindx(j).eq.i ) then
*               Exp coefficient for this earthquake
                nc = nc + 1
                osc = (1-exp(-outlog_days/eq_exp_tau(i)))*1000 ! mm
                comdP(nc) = osc
                pardP(nc) = j
           endif
         end do

****     See we have any output
         if( nc.gt.0 ) then
            ne = i
            if( une(4,ne).eq.0 ) then
*               Need to open file
                une(4,ne) = 303 + (ne-1)*4
                un = une(4,ne)
                eqf = eqo_root(1:trimlen(eqo_root)) //
     .                eq_codes(ne)(1:2) // '_Eq.tot'
                open(un,file=eqf,iostat=jerr, 
     .               status='unknown')
                call report_error('IOSTAT',jerr,'open',
     .               eqf,0,'tsfit/fit_out/eq_off')
                write(un,220)eq_codes(ne)(1:2), outlog_days,
     .                cmdfile(1:trimlen(cmdfile)),
     .                eq_log_tau(ne), eq_exp_tau(ne),
     .                eq_rad(ne)/1000, eq_rad_scale 

 220            format('* TSFIT Earthquake motion for EQ ',a2,
     .                 1x,F8.2,' days after ',
     .                 ' from command file ',a,/,
     .                 '* Log tau ',F8.2,' days, Exp tau ',F8.2,
     .                 ' days, Radius ', F8.2,' km, Scale ',F6.2)
                write(un,240)
 240            format('*  Long         Lat        ETot    NTot',
     .               '    dEt     dNt    E +-    N +-    Rne',
     .               '      HTot     dHt    H +-  Site',/,
     .               '*  deg          deg        mm      mm   ',
     .               '   mm      mm      mm      mm        ',
     .               '       mm      mm      mm  ')
            endif
            un = une(4,ne)

****        Now compute the total motion and sigma
            do c = 1,3    ! Loop over NEU
               dP(c) = 0.0
               do j = 1,nc
                  dP(c) = dP(c) + comdP(j)*smet(pardP(j),c)
                  do k = 1,nc
                     cov(j,k) = solcov(pardP(j),pardP(k),c)
                  end do
               enddo
*              Compute varaiance
               var = 0.d0
               do j = 1, nc
                  do k = 1,nc
                     var = var + comdP(j)*cov(j,k)*comdP(k)
                  end do
               end do
c               call var_comp(comdP,cov, var,tcov,1,3,1)
               sig(c) = sqrt(var)
            end do
*           Now output values
            if( sig(1)+sig(2)+sig(3).lt.100 )
     .      write(un,260) ref_llu(2)*180/pi, 
     .           (pi/2-ref_llu(1))*180/pi,
     .            dP(2), dP(1), dP(2), dP(1), sig(2),sig(1),0.0,
     .            dP(3), dP(3), sig(3), 
     ,            pnm_site(1)(1:6) // eq_codes(ne)(1:2)

 260        format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .             3(1x,f7.2), 1x,a8,' ')
         end if
      end do

****  Thats all
      return
      end

CTITLE sort_tsfeq
 
      subroutine sort_tsfeq

      implicit none
  
*     This routine will sort the earthquake list into ascending
*     time order.
 
      include 'tsfit.h'
 
*   i,j,k       - Loop counters
*   temp_i4     - I*4 variaable
 
      integer*4 i,j,k, temp_i4
 
*   temp_8      - Temporary Real*8 for switching
 
      real*8 temp_8
 
*   temp_ch     - Temporary character for switching.
 
      character*8 temp_ch
 
***** Now do sort
      do i = 1, num_eq-1
          do j = 1, num_eq-i
*                                                ! Time sort 
              if( eq_epoch(j).gt.eq_epoch(j+1) ) then
                  call switch_ch(eq_codes(j), eq_codes(j+1),
     .                           temp_ch )
                  call switch_8(eq_epoch(j), eq_epoch(j+1), temp_8)
                  call switch_8(eq_rad(j), eq_rad(j+1), temp_8)
                  call switch_8(eq_depth(j), eq_depth(j+1), temp_8)
                  call switch_8(eq_log_tau(j), eq_log_tau(j+1), temp_8)
                  call switch_8(eq_exp_tau(j), eq_exp_tau(j+1), temp_8)

                  do k = 1,2
                      call switch_8(eq_dur(k,j),
     .                              eq_dur(k,j+1),temp_8)
                  end do
                  do k = 1,3
                      call switch_8(eq_pos(k,j),
     .                              eq_pos(k,j+1),temp_8)

                  end do
                  do k = 1,6
                      call switch_8(eq_apr_coseismic(k,j),
     .                              eq_apr_coseismic(k,j+1),temp_8)
                      call switch_8(eq_mar_pre(k,j),
     .                              eq_mar_pre(k,j+1),temp_8)
                      call switch_8(eq_mar_post(k,j),
     .                              eq_mar_post(k,j+1),temp_8)
                      call switch_8(eq_log_sig(k,j),
     .                              eq_log_sig(k,j+1),temp_8)
                       call switch_8(eq_exp_sig(k,j),
     .                              eq_exp_sig(k,j+1),temp_8)
                 end do
                 call switch_I4(eq_rename(j), eq_rename(j+1), temp_i4)
                 call switch_I4(log_eq(j), log_eq(j+1), temp_I4)
                 call switch_I4(exp_eq(j), exp_eq(j+1), temp_I4)
                 call switch_I4(dtr_eq(j), dtr_eq(j+1), temp_I4)
 
              end if
          end do
      end do
 
***** Thats all
      return
      end
 
CTITLE ACCUM_OFFSETS

      subroutine accum_offsets(np, sum_soln, cnt, type) 

      implicit none

*     Routine to sum up all the offsets to get final 
*     positions.  Type indicates if BREAK or EQOFFSET.

      include 'tsfit.h'

      integer*4 np  ! Parameter number being evaluted
      integer*4 cnt ! Counter for number of breaks and earthqaukes
                    ! to generate name.

      real*8 sum_soln(3)  ! NEU sum of offsets

      character*(*) type

* LOCAL
      integer*4 j,k


*     Now sum up the previous values so that this is 
*     total offset
      do k = 1,3
         sum_soln(k) = 0.0
      end do

***   See if simple break
      cnt = 1
      if( type(1:1).eq.'B' ) then

****      Add up earilier offsets and earthquakes
         do j = 3, num_par
            if( pcode(j).eq.5 .and.
     .           rn_times(1,pindx(j)).le.
     .           rn_times(1,pindx(np)) ) then
                 do k = 1,3
                    sum_soln(k) = sum_soln(k) + smet(j,k)
                 end do
*                Add count to name
                 cnt = cnt + 1
             endif
*            See if EQ offsets need to be included
*            as well
             if( pcode(j).eq.6 .and.
     .           eq_epoch(pindx(j)).le.
     .           rn_times(1,pindx(np)) ) then
                 do k = 1,3
                    sum_soln(k) = sum_soln(k) + smet(j,k)
                 end do
*                Add count for name
                 cnt = cnt + 1
             endif
          end do
      elseif( type(1:1).eq.'E' ) then

****      Add up earilier earthquakes and offsets
          do j = 3, num_par
            if( pcode(j).eq.5 .and.
     .           rn_times(1,pindx(j)).le.
     .           eq_epoch(pindx(np)) ) then
                 do k = 1,3
                    sum_soln(k) = sum_soln(k) + smet(j,k)
                 end do
                 cnt = cnt + 1
             endif
*            See if EQ offsets need to be included
*            as well
             if( pcode(j).eq.6 .and.
     .           eq_epoch(pindx(j)).le.
     .           eq_epoch(pindx(np)) ) then
                 do k = 1,3
                    sum_soln(k) = sum_soln(k) + smet(j,k)
                 end do
                 cnt = cnt + 1
             endif
          end do
      endif

****  Thats all
      return
      end

CTITLE SUMCMD

      subroutine sumcmd(un,cmdfile)

      implicit none

*     Summarize command file to unit un.

      integer*4 un   ! Unit for output
      character*(*) cmdfile  ! Name of command file

      integer*4 ierr, j, trimlen
      character*256 line

****  re-Open command file and output
      open(300,file=cmdfile,iostat=ierr,status='old')
      do while ( ierr.eq.0 ) 
          read(300,'(a)',iostat=ierr) line
          if( ierr.eq.0 .and. trimlen(line).gt.0 ) then
             write(un,120) line(1:trimlen(line))
 120         format('# ',a)
          endif
      end do
****  Thats all
      close(300)
      return
      end

CTITLE KF_TSPOS

      subroutine kf_tspos

      implicit none

*     Routine to run the Kalman filter on the timeseries and
*     parameters set up by the least squares solution.
*     N = 1; E = 2 ; U = 3

      include 'tsfit.h' 

      integer*4 i, j, n  ! Loop counters
     .,         k     ! Component (NEU)
     .,         no    ! Number of observation.
     .,         np    ! Actual number of parameters
     .,         unts    ! Unit number for ts_<site> file output

      real*8 Kgain(max_par)   ! Kalman gain
     .,     dnorm_eq(max_par,max_par)   ! Decrement to covariance matrix
     .,     dt        ! Time step (years)
     .,     S         ! State transition for FOGM parameter
     .,     modvar    ! Variance of model projected + data_var
     .,     xe        ! Estimate of data from current state

* MOD 161121: Needed to compute residuals    
       real*8 wgh      ! Weight of obs
     .,      res      ! Pre and postfit residual
     .,      ures     ! Pre and postfit residual without user specified estimates applied

      logical apply_est ! set true if user specified this parameter to be included in residual estimates
      integer*4 npa     ! Number of parameters applied for residuals	      

****  Loop over NEU   
      np = num_par   
      do k = 1, 3   ! NEU

          ts_used(k) = 0 

*         Initalize
          call kf_init( k )
           
****      Loop over observations
          do no = 1, num_ts

****        If data point is OK, increment the Kalman filter
            if( ts_edt(k,no).eq.0  ) then 

                ts_used(k) = ts_used(k) + 1

****            Form the apartials
                call fit_part(ts_mjd(no))
*               Add a 1 at the end for the RW/FOGM term
                apart(num_par+1) = 1.d0
                dt = 0.d0
                if ( no.gt.1 ) dt = (ts_mjd(no)-ts_mjd(no-1))/365.25d0
                if( kfopt(1:2).eq.'FO' ) then
                   S = exp(-dt*365.25d0/taufin(k))
                else 
                   S = 1.d0
                endif

*               Propagate state forward.  Only the np+1 FOGM/RW
*               state entry has any change (Use norm_eq as covariance
*               matrix
                bvec(np+1) = bvec(np+1)*S
                norm_eq(np+1,:) =  norm_eq(np+1,:)*S
                norm_eq(:,np+1) =  norm_eq(:,np+1)*S
                norm_eq(np+1,np+1) = norm_eq(np+1,np+1)+RWvar(k)*dt

*               NOTE: This multiplies the diagonal by S^2
*               Now compute apart*norm_eq*apart'
                modvar = 0.0d0
                do i = 1, np+1
                   do j = 1, np+1
                      modvar = modvar + apart(i)*norm_eq(i,j)*apart(j)
                   end do
                end do
                modvar = modvar + (ts_neu_std(k,no))**2
                do i = 1, np+1
                   Kgain(i) = 0.0d0
                   do j = 1, np+1
                      Kgain(i) = Kgain(i) + norm_eq(i,j)*apart(j)/modvar
                   end do
                enddo
                xe = 0
                do i = 1, np+1
                   xe = xe + apart(i)*bvec(i)
                end do
                prechi(k) = prechi(k) + 
     .                      (ts_neu(k,no)-ref_neu(k)-xe)**2/modvar
                do i = 1, np+1
                    bvec(i) = bvec(i) + Kgain(i)*
     .                            (ts_neu(k,no)-ref_neu(k)-xe)
                    do j = 1,np+1
                       dnorm_eq(i,j) = 0.d0
                       do n = 1, np+1
                          dnorm_eq(i,j) = dnorm_eq(i,j) + 
     .                             Kgain(i)*apart(n)*norm_eq(n,j)
                       end do 
                    enddo
                end do
                norm_eq(1:np+1,1:np+1) = norm_eq(1:np+1,1:np+1) -
     .                                  dnorm_eq(1:np+1,1:np+1)
            endif   ! Data OK
          end do     ! Looping over the data points.

                       
****      Compute mean epoch before inversion
          mean_mjd(k) = (ts_mjd(1)+ts_mjd(num_ts))/2 + 
     .                    (norm_eq(2,1)/norm_eq(1,1))*365.25d0
          prechi(k) = prechi(k)/ts_used(k)
          ! print *,'MEAN PREF ',k,mean_mjd(k), prechi(k), ts_used(k)

*         Complete the solution
          sig_scale(k) = 1.0
          do i = 1,num_par
            soln(i,k) = bvec(i)*1000.d0
            smet(i,k) = bvec(i)
            solsig(i,k) = sqrt(norm_eq(i,i))*sig_scale(k)*1000.d0
            smtsig(i,k) = sqrt(norm_eq(i,i))*sig_scale(k)
            do j = 1,num_par
               solcov(i,j,k) = norm_eq(i,j)*sig_scale(k)**2
            enddo
          enddo

* MOD TAH 161121: Compute residuals (WRMS will be greater than WLS
*         solutuion.
****      Now compute the statistics
          stats = 0.d0
          npa = np  ! Applied number is the same.
          do no = 1, num_ts
             call fit_part(ts_mjd(no))
             wgh = 1.d0/ts_neu_std(k,no)**2
             res = ts_neu(k,no)-ref_neu(k)
             ures = ts_neu(k,no)-ref_neu(k)
             do i = 1,num_par
	     
**** Increment residuals by parameter estimates
	        
                res = res - apart(i)*bvec(i)
		
**** Increment residuals witout user specified parameters estimates (rcode) applied
                apply_est = .true.
                do j = 1, max_ptyp
		  if ( pcode(i) .eq. rcode(j) ) then
		     apply_est = .false.
		  endif
		enddo
		if ( apply_est ) then
                  ures = ures - apart(i)*bvec(i)   
		endif   
             enddo
	     
             ts_neu_res(no,k) = res
             ts_neu_ures(no,k) = ures
             ts_neu_sig(no,k) = ts_neu_std(k,no)
             if( ts_edt(k,no).eq.0 ) then 
                ts_neu_sig(no,k) = ts_neu_std(k,no)
                stats(1) = stats(1) + res*wgh
                stats(2) = stats(2) + wgh
                stats(3) = stats(3) + res**2*wgh
                stats(4) = stats(4) + 1
             else
                ts_neu_sig(no,k) = -ts_neu_std(k,no)
             endif
          end do
* MOD TAH 161121: Finish up the statistcs.
          if( stats(4)-npa.gt.0 ) then
             wrms(k) = sqrt(stats(4)/stats(2)*
     .                  (stats(3)/(stats(4)-npa)))*1.e3
          else
             wrms(k) = 100.0
          end if

          ts_used(k) = nint(stats(4))
    
      end do

****  If no edits, output results
      do k = 1,3
          call fit_out(k,unts)
      end do        

      return
      end

CTITLE KF_INIT 

      subroutine kf_init( k )

      implicit none

      include 'tsfit.h' 

*     Initialize KF run.
      integer*4 k   ! Component NEU

* LOCAL
      integer*4 i,j

      RWVar(k) = (solsig(2,k)/1000)**2*
     .           (ts_mjd(num_ts)-ts_mjd(1))/365.d0 + 
     .            min_rwvar*1.d-6  

      if( kfopt(1:2).eq.'WN' ) then
          RWVar(k) = 0.0d0
      endif 
      do i = 1, num_par+1
         bvec(i) = 0.d0
         do j = 1, num_par+1
            norm_eq(i,j) = 0.d0
         end do
*        Now add contraints to parameters as needed
         norm_eq(i,i) = 1.d0
         if( i.lt.num_par .and. apr_con(i,k).ne.0 ) then
            norm_eq(i,i) = apr_con(i,k)
         endif
      end do

      norm_eq( num_par+1, num_par+1) =  RWVar(k)/(2*taufin(k)/365.25d0); 

      return
      end




