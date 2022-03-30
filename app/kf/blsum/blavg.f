      program blavg

      implicit none 
c
c     This program will read the output file from a SOLVK solution
c     and produce a baseline summary.  The summary will be output
c     to one file and the values to another file.
c     The program can be run in batch so that it can be scheduled
c     after a batch solvk run.  In batch mode the runstring is
c     :RU,blavg,crt,printer,summary_file,values_file,solution1,solution2,
c          ...... up to maximum number of solutions
c
c Include files
c -------------
*                                 ! the blavg control common block
      include 'blavg.h'
c
c Variables
c ---------
c  ema storage of baselines lengths from solution file
c data -- storage to data from file.  When this array is used it is
c     passed to subroutines where it is split up as
c        data(1) = site_1
c        data(2) = site_2
c        data(3) = solution number
c        data(5:6)   = epoch (julian date)
c        data(7:8)  = length (m)
c        data(9:10) = sigma (m)
c
      integer*4 data(num_field,max_ent)
 
c
      common /values/ data
 
c
c file names
c ----------
c more_files    -- logical to indicate that we are still reading files
c new_name      -- dummy name for solution file
c status_old  -- variable which gives the status old to file opening
c status_unk  -- varaible which gives the status unknown to file opening
c unit  -- the unit number to be opened
c
      logical more_files
 
c
      character status_old*4, status_unk*8
 
c
c
c Runstring parameters
c --------------------
c irm -- rmpar parameters
c run_string -- the runstring as a character array
c irun_string -- integer array equivalenced to run_string (so that rhpar
c     will work.)
c count -- used to count our way through the runstring parameters
* i     -- Loop counter
* len_run -- Length of runstring
* ierr    -- IOSTAT error flag
* rcpar   -- Reads the rsunstring
c
      integer*4 count, i, len_run, ierr, rcpar
c
      character*128 run_string(max_runstring)
 
      character*10 run_num
 
      integer*4 start_date(3)
 
*       fjldy_8     - Julian date function
      real*8 fjldy_8
 
c
c Data statements for status
c --------------------------
      data status_old / 'old '/, status_unk / 'unknown' /
c
c.... Get the runstring and decide if we are in batch mode
      icrt = 6      
*     Get the minumnum of values needed for average
      len_run = rcpar(1, run_num)
      read(run_num,*, iostat=ierr) min_av_num
      if( len_run.eq.0 .or. ierr.ne.0 ) then
          call proper_runstring('blavg.hlp', 'blavg/min', 1)
      end if
*                              ! read rest of runstring
      do i = 1,3             
         len_run = rcpar(i+1, run_string(i))
c
         if( len_run.eq.0 ) then 
             call proper_runstring('blavg.hlp', 'blavg/names', 1)
         end if
c
      end do
c
***** Get the rest of the information from the runstring
      do i = 1,3
         len_run = rcpar(4+i,run_num)
         if( len_run.gt.0 ) then
             read(run_num,*) start_date(i)
         else
             call proper_runstring('blavg.hlp', 'blavg/date', 1)
         end if
      end do
      start_ep = fjldy_8(start_date(2), start_date(3), start_date(1))
 
      len_run = rcpar(8,run_num)
      if( len_run.gt.0 ) then
          read(run_num,*) stepep
      else
          stepep = 90.d0
      end if 

*     FORCE batch mode run 
      batch = .true.
c
c.... write headr message
      write(icrt,50)
  50  format(//" Program to aveage baseline lengths from BLSUM",
     .   " values file",/,
     .         " ----------------------------------------------",
     .   "------------",/)
 
      write(icrt,60) start_date, start_ep, stepep
  60  format(' Start Date ',i4,'/',i2,'/',i2,' (',F10.1,'). Step ',
     .       f5.1,' days')
c
c.... if we are in batch mode assign the file names from the runstring
      if( batch ) then
         summary_file     = run_string(1)
         values_file      = run_string(2)
c
c....    Now get the solution files from the runstring
         num_files = 0
         count = 3
         do while ( run_string(count)(1:1).ne.' ' .and.
     .      count.le.max_runstring)
            call check_soln_file(run_string(count),ierr,100,
     .         status_old)
            close(100)
c
c....       see if file found
*                                  ! file found ok
            if( ierr.eq.0 ) then
               num_files = num_files + 1
               solution_file(num_files) = run_string(count)
            end if
c
c....       goto next string
            count = count + 1
         end do
c
c....    open the summary file and the values files
         call check_soln_file(summary_file, ierr,200 ,status_unk)
         if ( ierr.ne.0 ) stop
     .      ' Error opening summary file in batch mode'
 
         call check_soln_file(values_file, ierr,201,status_unk)
         if ( ierr.ne.0 ) stop
     .       ' Error opening values file in batch mode'
c
*                           ! get files from user
      else
c
         num_files = 1
         more_files = .true.
*                                   ! get more file names
         do while ( more_files )
                                    
            write(icrt,'(a,i2,$)') 
     .           ' Enter the name of the solution file # ',num_files
c           the above replaces below
c            write(icrt,'(" Enter the name of solution file #",i2,
c     .         $)') num_files
c
            call get_soln_file(solution_file(num_files))
c
*                                                             ! see if file
            if( solution_file(num_files)(1:2).ne.'::' ) then
c                                                              exists
               call check_soln_file(solution_file(num_files),ierr,
     .            100,status_old)
               close(100)
c
c....          increment the file counter if file exists
               if( ierr.eq.0 ) num_files = num_files + 1
c
c....          check maximum of solutions
*                                                 ! this end of list
               if( num_files.gt.max_files ) then
                  write(icrt,'(" Maximum number of files reached")')
                  more_files = .false.
                  stop 'Too many entries: Update max_ent in blavg.h'
                  num_files  = num_files - 1
               end if
c
*                             ! this is the end of the input
            else
               more_files = .false.
               num_files = num_files - 1
 
            end if
*                           ! looping over solution files
         end do
c
c....    see if we have any files
         if( num_files.le.0 ) then
            stop ' No solution files given'
         end if
c
c....    Now get the name of the summary file
         ierr = 1
c
         do while ( ierr.ne.0 )
            write(icrt,'(/" Name of output summary file ",$)')
            call get_soln_file(summary_file)
c
            if( summary_file(1:2).eq.'::' ) then
               stop ' :: in response to summary file name'
            end if
c
            call check_soln_file(summary_file,ierr,200,
     .         status_unk)
         end do
c
c....    get the values file name
         ierr = 1
c
         do while ( ierr.ne.0 )
            write(icrt,'(" Name of output values file ",$)')
            call get_soln_file(values_file)
c
            if( values_file(1:2).eq.'::' ) then
               stop ' :: in response to values file name'
            end if
c
            call check_soln_file(values_file,ierr,201,
     .         status_unk)
         end do
c
      end if
c
c
c.... Now read in the data from the solution files
      num_ent = 0
      do i = 1, num_files
         write(icrt,'(" Reading baselines from ",a)') solution_file(i)
c
c....    use the check subroutine to open the solution file
         call check_soln_file(solution_file(i), ierr, 100,'old ')
         call read_soln(i,data)
         close(100)
c
      end do
c
*     Now sort the data into time order
C     call blsort(data)
 
c.... Now produce the summaries of the baseline lengths
      call sum_solution(data)
c
c.... Thats all
      close(200)
      close(201)
      end
c
c.....

      subroutine clear_av

      implicit none  
 
*     Clears information about averages
      include 'blavg.h'
 
      av_num = 0
      av_norm(1) = 0
      av_norm(2) = 0
      av_norm(3) = 0
      av_b(1) = 0
      av_b(2) = 0

      av_rms  = 0
      av_bm   = 0 
 
*     Thats all
      return
      end
 
CTITLE ACC_AV
 
      subroutine acc_av(epoch, length, sigma)

      implicit none 
 
 
*     Increments the averaging values
      include 'blavg.h'
 
*   epoch     - Epoch
*   wgh      - weight
 
      real*8 dt, epoch, length, sigma, wgh
 
 
      if( av_norm(1).eq.0 ) then
          ref_av_ep = epoch
          ref_av_bl = length
      end if
 
***** Incrememt normal quations
      wgh = 1.d0/sigma**2
      dt = (epoch - ref_av_ep)/365.25d0
 
      av_num = av_num + 1
 
      av_norm(1) = av_norm(1) + wgh
      av_norm(2) = av_norm(2) + wgh*dt
      av_norm(3) = av_norm(3) + wgh*dt**2
 
      av_b(1) = av_b(1) + (length-ref_av_bl)*wgh
      av_b(2) = av_b(2) + (length-ref_av_bl)*dt*wgh

      av_bm   = av_bm + (length-ref_av_bl-intercept-
     .                   slope*(epoch-ref_epoch)/365.25d0)*wgh

      av_rms = av_rms + (length-ref_av_bl-intercept-
     .                   slope*(epoch-ref_epoch)/365.25d0)**2*wgh
 
***** Thats all
      return
      end
 
CTITLE FINISH_AV
 
      subroutine finish_av
 
      implicit none 

 
*     Finished computation of avergaes
      include 'blavg.h'
 
*       det   - determinate
*       inverse(3)   - Inverse of matrix
 
 
      real*8 det, inverse(3)
 
*                                ! No data
      if( av_num.eq.0 ) return
 
      av_length = av_b(1)/av_norm(1) + ref_av_bl
      av_lsig   = sqrt(1.d0/av_norm(1))
 
***** Now get slope
      det = av_norm(1)*av_norm(3) - av_norm(2)**2
 
      av_epoch = ref_av_ep + (av_norm(2)/av_norm(1))*365.25d0
 
***** See if we have slope
*                             ! No slope
      if( det.lt.tol ) then
          av_slope = 0
          av_ssig = 0
      else
          inverse(1) = av_norm(3)/det
          inverse(2) = -av_norm(2)/det
          inverse(3) = av_norm(1)/det
 
*****     Get slope
          av_slope = inverse(2)*av_b(1) + inverse(3)*av_b(2)
          av_ssig = sqrt(inverse(3))
      end if

      if( av_num.gt.1 ) then
          av_nrms = sqrt(abs(av_rms - av_bm*(av_bm/av_norm(1)))/
     .                   (av_num-1))
          av_wrms = av_nrms*sqrt(av_num/av_norm(1))
      else
          av_nrms = 1.d0
          av_wrms = av_nrms*sqrt(av_num/av_norm(1))
      end if
 
***** Thats all
      return
      end
 
CTITLE WRITE_AV
 
      subroutine write_av(unit)

      implicit none 
 
c     Routine to write the data record to the values file
c
      include 'blavg.h'
c
c variables
c ---------
c unit -- the output unit number
c
      integer*4 unit

* date(5) -- Date of this value

      integer*4 date(5)

* sectag  -- Seconds tag

      real*8 sectag, dtime, residual
 
c
c.... set up for output, convert epoch to yr,mon,day,hr,min
      if( av_num.eq.0 ) return
 
      call jd_to_ymdhms(av_epoch, date, sectag)
      dtime = (av_epoch - ref_epoch)/365.25d0
      residual = av_length - (ref_length+intercept) - dtime*slope
c
c.... write the record
      write(unit,100) date,av_length, av_lsig, residual, av_lsig,
     .                av_wrms, av_nrms, 
     .                av_slope, av_ssig, av_num
 
  100 format(i5,4i3,f16.5,1x,f10.5,1x,6(f10.5,1x),1x,i4)
c
      return
      end
c
 
CTITLE    ......................................................
c
      subroutine read_soln(soln_num, data)

      implicit none 
c
c     routine to read in the solution file and save the data in ema. As
c     the solution file is read, the list of sites involved in these
c     experiments is incremented.
c
c     The solution file will already be open to unit 100 when this
c     routine is called.
c
c Include files
c -------------
c
*                                    ! the control common block
      include 'blavg.h'
c
c
c passed variables
c ----------------
c data -- ema storage of the baselines from the solution file
c
      integer*4 data(num_field,max_ent)
 
c
c
c soln_num -- the number of the solution file
c
      integer*4 soln_num
 
c
c Local variables
c ---------------
c buffer -- the buffer used to read the solutiom file
c new_buffer -- a second copy of the buffer
c
      character*150 buffer, new_buffer
 
c
c new_soln -- a new solution indentifier from SOLVK
c new_time -- a new epoch indentifier
c bas_to   -- baseline indentifier
c
      character new_soln*10, new_time*18, bas_to*4
 
c
c clk_date -- epoch of the experiment
c values -- used to get baseline length and sigma
c
      integer*4 clk_date(5)
 
c
      real*8 values(3)
 
c
*   indx        - Pointer for use in multiread
*   err         - Error flag returned from multiread
*   iel         - index position of string
*   iclk        - Index position for clock time
*   ibas        - index position for baseline length
*   in          - length of string
 
      integer*4 indx, err, iel, iclk, ibas, in
 
*   cdummy      - Dummy character string for multiread
 
      character*1 cdummy
 

*   new_gtime*24    - solution commenced with line
 
      character new_gtime*24
 
c.... Define the identifier strings
c
      data new_soln / ' SOLUTION ' /,
     .     new_time / 'EXPERIMENT DATE :'    /,
     .     bas_to   / ' TO '     /
c
      data new_gtime / 'SOLUTION REFERS TO     :' /
 
c.... Firstly read the header from solution
      read(100,'(a)',iostat=err) header(soln_num)
c
* MOD TAH 881031 Check to see if this is VALUES_FILE
      if( header(soln_num)(1:11).eq.'VALUES_FILE' ) then
          call read_values_file( soln_num, data )
*                      ! Do not need to do the rest
          return
      end if

*     See if name says that it is a values files
      if( solution_file(soln_num)(1:2).eq.'va' ) then
          rewind(100)
          call read_values_file( soln_num, data )
          return
      end if
 
 
      call trailing(in,header(soln_num))
      write(*,100) header(soln_num)(:in)
 100  format(" Solution: ",a)
c
c.... Now start reading the solution file extracting the baseline
c     lengths
      err = 0
c
      num_soln = 1
      max_soln = 1
 
*                              ! read until error
      do while ( err.eq.0 )
c
         read(100,'(a)',iostat=err) buffer
         call casefold(buffer)
c
c....    continue processing if there was no error
         if( err.eq.0 ) then
c
c....       see if start of a new solution
* MOD TAH 861008 Checked number after character 59 to aviod finding
*           clock break warning message.
 
            iel = index(buffer,new_soln(1:10))
 
*                                 ! Solution found
            if( iel.eq.47 ) then
               indx = iel + 10
               call multiread(buffer,indx,'I4',err,num_soln,cdummy,1)
 
               if( num_soln.eq.0 .or. err.ne.0 ) then 
                  write(*,'(a,a,/,a,4i8)') ' Solution number error '
     .                 ,buffer,' iel, indx, num_soln, err = '
     .                 , iel, indx, num_soln, err
c                 the above replaces the below to avoid Hollerith continuation
c                  write(*,'(" Solution number error in ",a,/,
c     .                    " iel, indx, num_soln, err = ",4i8)')
c     .                buffer, iel, indx, num_soln, err
                  num_soln = 1
               end if
 
*                                                              ! keep track
               if( num_soln.gt.max_soln ) max_soln = num_soln
c                                of the number of solutions we need to
c                                monitor the perferomance of
c
            end if
c
c....       see if a new epoch has been encountered
            iclk = index(buffer, new_time)
*                                 ! new clock epoch -- get value
            if( iclk.gt.0 ) then
*                                            ! 18  is the length of new_time
               read( buffer(iclk+18:),140, iostat=err) clk_date
 140           format(i4,1x,4(i2,1x))
c
c....          convert to julian date
               call ymdhms_to_jd( clk_date, 0.d0, sol_epoch )
            end if
 
*           Check clock time for a global solution
            iclk = index(buffer, new_gtime)
*                                 ! New clock epoch
            if( iclk.gt.0 ) then
               read( buffer(iclk+27:),150, iostat=err)
     .            clk_date
 150           format(5(i2,1x))
               call ymdhms_to_jd( clk_date, 0.d0, sol_epoch )
            end if
c
c....       Now see if we have a baseline length
            ibas = index(buffer, bas_to)
*                                   ! we found a baseline length
!            if( ibas.eq.14 ) then
            if( ibas.eq.14 .and. buffer(1:1).eq.' ' ) then
c
c....          get the first site in the baseline. If there are no
c              sites yet then skip the call to strip_name.
*                                         ! just set iel=-1
               if( num_site.eq.0 ) then
                  iel = -1
*                                         ! see if in list
               else
                  call strip_name(buffer,names,num_site,iel,0)
               end if
c
c....          see if we already know this site
*                                    ! site name not found
               if( iel.eq.-1 ) then
                  call add_name(buffer)
                  call strip_name(buffer,names,num_site,iel,0)
               end if
c
               site(1) = iel
c
c....          Now strip out the ' to ' from the buffer.  The first
c              name has already been stripped out by site_name.
               new_buffer = buffer(4:)
c
c....          Get second site in baseline
               call strip_name(new_buffer,names,num_site,iel,0)
c
c....          see if we know this site
*                                    ! site name not found
               if( iel.eq.-1 ) then
                  call add_name(new_buffer)
                  call strip_name(new_buffer,names,num_site,iel,0)
               end if
c
               site(2) = iel
c
c....          Now get the baseline data
               indx = 1
               call multiread(new_buffer,indx,'R8',err,values,
     .                         cdummy,3)
C              call get_real(new_buffer, values, 3)
               base_data(1) = values(1)
               base_data(2) = values(3)
               if( base_data(2).eq.0 ) base_data(2) = 0.0001
c
c....          Make sure that the site numbers are in ascending order
               if( site(2).lt.site(1) ) then
                  iel = site(1)
                  site(1) = site(2)
                  site(2) = iel
               end if
c
c.....         Save this data in ema
               num_ent = num_ent + 1
c
c.....         see if maximun number of entries exceeded
               if( num_ent.ge.max_ent ) then
                  write(*,200) num_ent
  200             format(/" Maximum number of baselines ",i8,
     .               " has been exceeded")
*                              ! Force exit on next read
                  err = -99
               end if
c
               call save_ema(data(1,num_ent), data(5,num_ent),
     .           data(7,num_ent),data(9,num_ent))
c
c....       Thats all
*                         ! baseline found
            end if
*                         ! no error on file read
         end if
c
*                         ! looping over the file to be read
      end do
c
c.... Write out number of entries
      write(*,300) num_ent, num_site
  300 format(" There are ",i8," baseline entries between ",i3," sites")
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine check_soln_file(file_name, ierr, unit, status)

      implicit none 
c
c     Routine to open a file 'file_name' on unit 'unit' with status
c     'status'.
c
c Variables
c ---------
c file_name -- file to be opened
c ierr   -- iostat error opening file
c unit   -- the unit number of the device
c status -- the status with which the file should be opened
c in     -- index for character counts in trailing
c
      integer*4  ierr, unit, in
 
c
      character*(*) file_name, status
 
c
c
c.... open the file
c
      open(unit, file=file_name, iostat=ierr, status=status)
c
c.... check error
      if( ierr.ne.0 ) then
         call trailing(in,file_name)
         write(*,100) ierr, file_name(:in), status
  100    format(" IOSTAT error ",i4," occurred opening ",a,
     .      " status ",a)
      end if
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine get_soln_file(file_name)

      implicit none 
c
c     Routine to ask for a file name.
c
c Variables
c ---------
c file_name -- the name of the file selected
c
c
      character*(*) file_name
 
c
c.... get the file name, the message asking for file name has already
c     been printed
      read(*,'(a)') file_name
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine sum_solution(data)
c
c     Routine to read through the ema data and accumualte statistics
c     for each baseline.
c
c
c Include files
c -------------
c
*                                     ! the control common block
      include 'blavg.h'
 
* i,j,k,l,m  - Loop counters
* ierr       - IOSTAT error
* inf, inh   - lengths of file and header
 
      integer*4 i,j,k,m, inf, inh, ierr

* printed   - Indicates that some values have been pronted

      logical printed

* epoch    - JD of current data
* ema8     - Function to return real value from integer array

      real*8 epoch, ema8

* MOD TAH 980326: Check for ENSUM type solution.  Make output
*     more readable

* istart, iend -- Start and stop of i - loop
* jstart, jend -- Start and stop of j - loop
* ensum_input -- Ensum input

      integer*4 istart, iend, jstart, jend
      logical ensum_input
c
c EMA variables
c -------------
c data -- the stored baseline data
c
      integer*4 data(num_field,max_ent)
c
c
c.... start looping over all baselines.  Write out the headers to the
c     summary file
c
      write(200,100)
  100 format(/" SOLVK baseline summary: includes the following",
     .   " solutions")
c
      do i = 1, num_files
         call trailing(inf, solution_file(i))
         call trailing(inh,header(i))
         write(200,120) solution_file(i)(1:inf), header(i)(:inh)
  120    format(" File ",a,": ",a)
      end do
c
c.... Write out descriptive header
      write(200,140)
 140  format("    Baseline       S#   #   mean length   sig   wrms",
     .  "  nrms    slope   sig    wrms nrms   dur   mean    line #",
     . /,    "                                (m)       (m)   (mm)",
     .  "        (mm/yr) (mm/yr)  (mm)       (yrs)  (yrs)")
 
c.... start the line counter for tplot file extraction
      line_st = 1
      line_en = 0
      if( max_soln.eq.0 ) max_soln = 1

*     See if ensum type input
      if( names(2).eq.'N SOLUTI' ) then
         ensum_input = .true.
         istart = 1
         iend = num_site
         names(2) = 'NORTH'
         names(3) = 'EAST'
         names(4) = 'HEIGHT'
      else
         ensum_input = .false.
         istart = 1
         iend = num_site - 1
      endif
c
c.... loop over all baseline combinations
      do i = istart, iend
c
         if( ensum_input ) then
            jstart = 2
            jend = 4
            if( i.eq.1 .or. i.gt.4 ) then
               write( *,9000) names(i)
            end if
         else
            write( *,9000) names(i)
 9000       format(" Processing data from site ",a)
            jstart = i+1
            jend = num_site
         endif

         do j = jstart, jend 
c
c
c....       loop over the solution
            do k = 1, max_soln
c
               call clear_stat
c
c....          Now get all the entries
               nrms_mean = 1
               nrms_slope = 1
               curr_ep = start_ep
               call clear_av

               printed = .false.

c
c....          Now get all the entries
               do m = 1,num_ent
c
c....             see if baseline and solution matches
                  if( ((data(1,m).eq.i .and. data(2,m).eq.j).or.
     .                 (data(1,m).eq.j .and. data(2,m).eq.i)) .and.
     .                data(3,m).eq.k ) then
c
c....                accumulate statistics
                     epoch = ema8(data(5,m))
C                                                          ! Finished
                     if( epoch.gt.curr_ep + stepep ) then  

                          call finish_av
                          curr_ep = start_ep +
     .                        int((epoch-start_ep)/stepep)*stepep
                          if( av_num.gt.0 ) then
                             call slope_acc(av_epoch,av_length,
     .                                      av_lsig)
                          end if
                          call clear_av
                     end if
                     call acc_av(data(5,m), data(7,m), data(9,m))
                  end if
               end do
c
*****          See of any data left over
               if( av_num.ne.0 ) then
                   call finish_av
                   call slope_acc(av_epoch,av_length, av_lsig)
               end if

c....          complete the computation of the statistics
               call compute_slope(ierr)

*              Now loop over the data again and output values
*              with their residuals and rms
               curr_ep = start_ep
               call clear_av

c....          Now get all the entries
               do m = 1,num_ent
c
c....             see if baseline and solution matches
                  if( ((data(1,m).eq.i .and. data(2,m).eq.j).or.
     .                 (data(1,m).eq.j .and. data(2,m).eq.i)) .and.
     .                data(3,m).eq.k ) then
c
c....                accumulate statistics
                     epoch = ema8(data(5,m))
C                                                          ! Finished
                     if( epoch.gt.curr_ep + stepep ) then  

*                         Since we foumd something
                          if( .not.printed ) then
                              call write_header(201,i,j,k)
                              printed = .true.
                          end if

                          call finish_av
                          call write_av(201)
                          call clear_av
                          curr_ep = start_ep +
     .                        int((epoch-start_ep)/stepep)*stepep
                     end if
                     call acc_av(data(5,m), data(7,m), data(9,m))
                  end if
               end do
c
*****          See of any data left over
               if( av_num.ge.min_av_num ) then
                   call finish_av
                   call write_av(201)
               end if

               if( printed ) then
                  call write_summary(200,i,j,k,'short ')
                  call write_summary(201,i,j,k,'long')
c
c....             save line numbers
                  line_st = line_en + 5
                  line_en = line_st - 1
               end if

c
*                     ! looping over solutions   -- k
            end do
*                     ! looping over second site -- j
         end do
*                     ! looping over first site  -- i
      end do
c
c.... Thats all
      return
      end
c
CTITLE    ......................................................
c
      subroutine write_header(unit,i,j,k)

      implicit none 
c
c     Routine to write the header records to the values file
c
      include 'blavg.h'
c
c Variables
c ---------
c unit -- output unit number
* in   -- Length of string
c
      integer*4 unit, i, j, k, in
 
c
c.... The header record to be used by tplot is the solution
c     description, baseline + solution #, blank record
c
c.... delete trailing blanks from header
      call trailing(in,header(1))
      write(unit,'(a)') header(1)(:in)
c
      write(unit,'(a," to ",a," Solution ",i2)') names(i),names(j),k
c
      write(unit,'(1x)')
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine compute_slope(ierr)

      implicit none 
c
c     Routine to compute the slope of best fit straight line
c
      include 'blavg.h'
c
c Local variables
c ---------------
c det -- the determinate of the normal equations
c inverse -- the inverse of the normal equations
c ierr -- an error number returned to main routine
c      = 0 if all OK
c      = -1 if no data
c      = +1 if slope inversion is singular
* sug_inter - sigma of intersect
c
      real*8 det, inverse(3), sig_inter
 
* ierr  - IOSTAT error
 
      integer*4 ierr
 
c
c.... see if we have enough data to do a mean value
      ierr = 0
      if( norm_eq(1).lt. tol ) then
         wgh_mean = 0.d0
         sig_mean = 0.d0
c
         intercept = 0.d0
         sig_inter = 0.d0
         slope     = 0.d0
         sig_slope = 0.d0
c
*                        ! probably no data
         ierr = -1
c
*                 ! compute mean
      else
c
c....    Now get weighted mean
         wgh_mean  = b_vec(1)/norm_eq(1)
         sig_mean  = sqrt(1.d0/norm_eq(1))
c
c....    Invert the normal equations
         det = norm_eq(1)*norm_eq(3) - norm_eq(2)**2
c
c....    compute the mean epoch of the data  ! julian date
         mean_epoch = ref_epoch + (norm_eq(2)/norm_eq(1))*365.25d0

*        Finish up the statistics.
         if( num_tot.gt.1 ) then
             nrms_mean = sqrt((rms - b_vec(1)*wgh_mean)/(num_tot-1))
         else
             nrms_mean = 1.0
         end if
         wrms_mean = nrms_mean*sqrt(num_tot/norm_eq(1))

c
c....    see if singular
*                                ! not enough data to do slope
         if( det.lt. tol ) then
            intercept = wgh_mean
            sig_inter = sig_mean
c
            slope     = 0.d0
            sig_slope = 0.d0
c
*                        ! matrix singular -- only one data
            ierr = +1
c
*                   ! compute the slope
         else
c
            inverse(1) =  norm_eq(3)/det
            inverse(2) = -norm_eq(2)/det
            inverse(3) =  norm_eq(1)/det
c
c....       compute intercept and slope
*                                                                  ! meters
            intercept = inverse(1)*b_vec(1) + inverse(2)*b_vec(2)
*                                                                  ! m/year
            slope     = inverse(2)*b_vec(1) + inverse(3)*b_vec(2)
c
            sig_inter = sqrt(inverse(1))
            sig_slope = sqrt(inverse(3))

            if( num_tot.gt.2 ) then
                nrms_slope = sqrt((rms - b_vec(1)*intercept 
     .                                 - b_vec(2)*slope )/(num_tot-2))
            else
                nrms_slope = 1.0
            end if
            wrms_slope = nrms_slope*sqrt(num_tot/norm_eq(1))
c
         end if
      end if
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine write_summary(unit,site_1,site_2, soln, format)

      implicit none 
c
c     Routine to write the summary information.  This can be written
c     in two forms.. the short where all results are put on a single
c     line or the long form which takes multiple lines
c
      include 'blavg.h'
c
c variables
c ---------
c unit -- the output unit number
c site_1 -- first site in baseline
c site_2 -- second site in baseline
c soln -- solution number
c format -- determines if long or short form
c
      integer*4 unit, site_1, site_2, soln
 
c
      character*(*) format
 
c
c Local variables
c ---------------
c sig_mean_scaled -- scaled sigma of mean (m)
c sig_slope_scaled -- scaled sigma of slope (mm/yr)
c slope_out -- output slope (mm/yr)
c mean_out  -- the output weighted mean (m)
c wmrs_mean_out -- wrms_mean output (mm)
c wmrs_slope_out -- wrms_slope ouput (mm)
c duration -- output duration of data (yrs)
c mean_year -- output mean epoch (yrs)
* date      -- Date for JD
* sec_tag   -- sec_tag of time
c
      real*8 sig_mean_scaled, sig_slope_scaled, slope_out, mean_out,
     .    wrms_mean_out, wrms_slope_out, duration, mean_year, sec_tag
 
* num_len -- number of values in mean length estimate (used so we can output
*            with an I4 format
 
      integer*4 date(5), num_len
 
c
c.... Convert the units form meters to millimeters where necessary
c     and scale the sigmas of the various quanities
c
*                                  ! scale wgh_mean sigma
      if( nrms_mean.gt.1 ) then
         sig_mean_scaled = sig_mean*nrms_mean
      else
         sig_mean_scaled = sig_mean
      end if
c
c.... Now do the slope quantities
*                                  ! scale slope sigma
      if( nrms_slope.gt.1 ) then
*                                                         ! mm/yr
         sig_slope_scaled = sig_slope*nrms_slope*1000.d0
      else
*                                                         ! mm/yr
         sig_slope_scaled = sig_slope*1000.d0
      end if
c
c.... Convert units on slope
      slope_out = slope*1000.d0
      mean_out  = ref_length + wgh_mean
c
      wrms_mean_out = wrms_mean*1000.d0
      wrms_slope_out = wrms_slope*1000.d0
c
c.... get the duration of the data
      duration = (max_epoch - min_epoch)/365.25
      call jd_to_ymdhms( mean_epoch, date, sec_tag)
      mean_year = date(1) + date(2)/12.d0 + date(3)/365.25
c
c.... Now do the output
      num_len = stat(4,1)
      num_len = num_tot
      if( format.eq.'short' ) then
c
         write(unit,100) names(site_1),names(site_2), soln, num_len,
     .      mean_out,sig_mean_scaled, wrms_mean_out, nrms_mean,
     .      slope_out, sig_slope_scaled, wrms_slope_out, nrms_slope,
     .      duration, mean_year, line_st, line_en
 100     format(a8,"-",a8,1x,i2,1x,i4,1x,
     .      f13.4,1x,f6.4,1x,f6.2,1x,f7.2,3x,
     .      f8.2,1x,f8.2,1x,f5.1,1x,f5.2,1x,
     .      f6.2,1x,f7.2,1x,i5,"-",i5)
c
*                  ! do the long output
      else
c
         write(unit,200) mean_out, sig_mean_scaled, num_len,
     .      wrms_mean_out, nrms_mean
 200     format(/"Wmean ",f13.4," m +- ",f6.4," from ",i4," data.",
     .      " WRMS ",f5.1," mm, NRMS ",f5.2)
         write(unit,220) slope_out, sig_slope_scaled, wrms_slope_out,
     .      nrms_slope, duration, mean_year
 220     format("Slope ",f8.2," +- ",f8.2," mm/yr, WRMS ",f5.1,
     .      " mm, NRMS ",f5.2,", dur ",f6.2," <> ",f7.2," yr"/)
c
      end if
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine slope_acc(epoch,length,sigma)

      implicit none 
c
c     Routine to accumulate the data needed to compute the slope
c
      include 'blavg.h'
c
c Variables
c ---------
c epoch -- the epoch of this measurement
c length -- the length of the baseline
c sigma  -- the sigma of the baseline
c
      real*8 epoch, length, sigma
 
c
c
c Local variables
c ---------------
c wgh  -- the length measurement weight
c dtime -- time difference between epoch and ref_epoch in years
c
      real*8 wgh, dtime
 
c
c.... If this is the first value save the epoch and length
c
      if ( norm_eq(1).eq.0.d0 ) then
         ref_epoch = epoch
         ref_length = length
      end if
c
c.... see if this this value is before or after all others
      if( epoch.lt. min_epoch ) min_epoch = epoch
      if( epoch.gt. max_epoch ) max_epoch = epoch
c
c.... Now accumulate normal equations
      wgh = 1.d0/sigma**2
      dtime = (epoch-ref_epoch)/365.25d0
      norm_eq(1) = norm_eq(1) + wgh
      norm_eq(2) = norm_eq(2) + dtime*wgh
      norm_eq(3) = norm_eq(3) + dtime**2*wgh
c
      b_vec(1)   = b_vec(1) + (length-ref_length)*wgh
      b_vec(2)   = b_vec(2) + (length-ref_length)*dtime*wgh

      rms = rms + (length-ref_length)**2*wgh
      num_tot = num_tot + 1
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine clear_stat

      implicit none 
c
c     This routine will clear the statistics accumulation and
c     slope computation  arrays
c
      include 'blavg.h'
 
* i,j  -- loop counters
      integer*4 i,j
c
c.... clear all of the arrays
      do i = 1,4
         do j = 1,2
*                              ! summation variables for mean and rate
            stat(i,j) = 0.d0
c                                statistics.
         end do
      end do
c
      do i = 1,3
*                              ! normal equations for slope computation
         norm_eq(i) = 0.d0
      end do
c
      do i = 1,2
*                              ! b_vector slope computation
         b_vec(i)   = 0.d0
      end do

      rms = 0
      num_tot = 0
c
c.... initialise the min and max times
      min_epoch = 1.d30
      max_epoch = -1.d30
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine accum_stat(epoch, length, sigma)

      implicit none 
c
c     Routine to accumulate the statistics from a set of baseline
c     results.  The statistics are accumulated for both the residuals
c     about the mean and about the best fit straight line
c
      include 'blavg.h'
c
c EMA variables
c -------------
c epoch -- epoch of measuarement
c length -- baseline length
c sigma  -- sigma of length
c
      real*8 epoch, length, sigma
 
c
c
c Local variables
c ---------------
c wgh -- the weight of the observation
c residual -- the residual
c dtime -- epoch difference in years
c
      real*8 wgh, residual, dtime
 
c
c.... Compute the residuals and sum into the statistics
      wgh = 1.d0/sigma**2
      dtime = (epoch-ref_epoch)/365.25d0
c
c.... do statistics about the mean
      residual = length - (ref_length+wgh_mean)
c
      stat(1,1) = stat(1,1) + residual*wgh
      stat(2,1) = stat(2,1) + 1.d0*wgh
      stat(3,1) = stat(3,1) + residual**2*wgh
      stat(4,1) = stat(4,1) + 1.d0
c
c.... statistics about the slope
      residual = length - (ref_length+intercept) - dtime*slope
c
      stat(1,2) = stat(1,2) + residual*wgh
      stat(2,2) = stat(2,2) + 1.d0*wgh
      stat(3,2) = stat(3,2) + residual**2*wgh
      stat(4,2) = stat(4,2) + 1.d0
c
c.... Thats all
      return
      end
c
CTITLE READ_VALUES_FILE
 
      subroutine read_values_file( soln_num, data )

      implicit none 
 
 
*     Routine to read a values file to extract information about the
*     baseline lengths.
 
      include 'blavg.h'
 
*         soln_num      - Number of current file
*         data(num_field, max_ent)   - Place where data is actually
*                       - saved.
 
      integer*4 soln_num, data(num_field, max_ent)
 
c
c Local variables
c ---------------
c buffer -- the buffer used to read the solutiom file
c new_buffer -- a second copy of the buffer
c
      character*150 buffer, new_buffer
 
c
c clk_date -- epoch of the experiment
c values -- used to get baseline length and sigma
c
      integer*4 clk_date(5)
 
c
      real*8 values(3)
 
c
*   indx        - Pointer for use in multiread
*   err         - Error flag returned from multiread
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   i           - Loop counter
*   iel, in     - Index positon temporary storgae
 
      integer*4 indx, ierr, trimlen, i, iel, in
 
*        still_reading   - Indicates that we are still reading entries
 
      logical still_reading
 
*   cdummy      - Dummy character string for multiread
 
      character*1 cdummy
 
*     Read the actual header for this file
 
      read(100,'(a)', iostat=ierr) header(soln_num)
      call trailing(in,header(soln_num))
      write(*,100) header(soln_num)(:in)
 100  format(" Solution: ",a)
c
c.... Now start reading the solution file extracting the baseline
c     lengths
      ierr = 0
c
*                              ! read until error
      do while ( ierr.eq.0 )
c
         read(100,'(a)',iostat=ierr) buffer
         call casefold(buffer)
c
c....    continue processing if there was no error
         if( ierr.eq.0 ) then
 
*        Get the baseline from buffer
         if( num_site.eq.0 ) then
             iel = -1
         else
             call strip_name(buffer,names, num_site, iel,0)
         end if
 
         if( iel.eq.-1) then
             call add_name(buffer)
             call strip_name(buffer,names, num_site, iel,0)
         end if
 
         site(1) = iel
 
*****    Now get rid of 'to '
         new_buffer = buffer( 4: )
 
*        Get second name
         call strip_name(new_buffer,names,num_site, iel,0)
         if( iel.eq.-1) then
             call add_name(new_buffer)
             call strip_name(new_buffer,names, num_site, iel,0)
         end if
 
         site(2) = iel
 
*        Get sites in correct order
         if( site(2).lt.site(1) ) then
             iel = site(1)
             site(1) = site(2)
             site(2) = iel
         end if
 
*****    Get out the solution number
         indx = 10
         call read_line(new_buffer,indx,'I4',ierr,num_soln,cdummy)
         if( ierr.ne.0 ) then
*            check for LNEU
             if( index(new_buffer(10:),'L').gt.0 ) num_soln = 1
             if( index(new_buffer(10:),'N').gt.0 ) num_soln = 2
             if( index(new_buffer(10:),'E').gt.0 ) num_soln = 3
             if( index(new_buffer(10:),'U').gt.0 ) num_soln = 4
         end if
         if( num_soln.eq.0 ) num_soln = 1

         if( num_soln.gt.max_soln ) max_soln = num_soln
 
****     Skip line and read in all the baselines
         read(100,'(a)') buffer
         still_reading = .true.
 
*                                         ! Loop getting baseline values
         do while ( still_reading )
             read(100,'(a)', iostat=ierr ) buffer
*                                    ! Then procesws
             if( ierr.eq.0 ) then
*                                                   ! End
                 if ( trimlen(buffer).eq.0 ) then
                      still_reading = .false.
                 end if
             else
                 still_reading = .false.
             end if
 
****         See if we should decode
             if( still_reading ) then
                 indx = 1
                 call multiread(buffer, indx, 'I4', ierr, clk_date,
     .                           cdummy,5)
                 call ymdhms_to_jd(clk_date, 0.d0, sol_epoch )
 
*                Now read the values, length and sigma
                 call multiread(buffer, indx,'R8',ierr,values,
     .                           cdummy,2)
 
                 base_data(1) = values(1)
                 base_data(2) = values(2)
                 if( base_data(2).eq.0 ) base_data(2) = 0.0001
 
c.....           Save this data in ema
                 num_ent = num_ent + 1
c
c.....           see if maximun number of entries exceeded
                 if( num_ent.ge.max_ent ) then
                    write(*,200) num_ent
  200               format(/" Maximum number of baselines ",i8,
     .                 " has been exceeded")
*                                ! Force exit on next read
                    ierr = -99
                 end if
c
                 call save_ema(data(1,num_ent), data(5,num_ent),
     .                        data(7,num_ent), data(9,num_ent) )
 
****          Thats all
*                        ! Still reading data
              end if
*                        ! Looping over the data for this baseline
          end do
 
****      Now skip the lines with statistics and header
          do i = 1, 4
             read(100,'(a)', iostat=ierr) buffer
          end do
 
*                        ! No file read error
          end if
*                        ! Looping until we run out of baselines
      end do
 
***** Write out number of entries
      write(*,300) num_ent, num_site
  300 format(' There are ',i8,' baselines entries between ',
     .                     i3,' sites' )
*
      return
      end
 
CTITLE    ......................................................
c
      subroutine save_ema(site_dat, epoch, length, sigma)

      implicit none 
c
c     routine to save the baseline values in ema
c
c Include files
c -------------
c
*                                   ! the control common block
      include 'blavg.h'
c
c Variables
c ---------
c site_dat -- integer values for sites and solution number
c epoch     -- real*8 epoch of measurement
c length    -- real*8 length value
c sigma     -- real*8 sigma
c
      integer*4 site_dat(3)
 
c
      real*8 epoch, length, sigma
 
c
c
c.... Save the values
      site_dat(1) = site(1)
      site_dat(2) = site(2)
      site_dat(3) = num_soln
c
      epoch = sol_epoch
      length = base_data(1)
      sigma  = base_data(2)
c
c
c
      return
      end
c
CTITLE    ......................................................
c
      subroutine add_name(buffer)

      implicit none 
c
c     Routine to add the first name in buffer to the list in names
c
c Include file
c ------------
c
*                                       ! the BLSUM parameter file
      include 'blavg.h'
c
c Variables
c ---------
c buffer -- the buffer read from file
* i      -- loop counter
c
      integer*4 i
      character*(*) buffer
 
c
c.... Firstly kill any leading blanks
      i = 1
      do while ( buffer(i:i).eq.' ' )
         i = i + 1
      end do
c
c.... copy the next 8 characters to names
      num_site = num_site + 1
      names(num_site) = buffer(i:i+8)
c
      return
      end
c
CTITLE EMA8    ......................................................
c
      real*8 function ema8(data)

      implicit none 
 
      real*8 data
 
****  Just set ema8 equal to data
 
      ema8 = data
 
      return
      end
 
CTITLE    ......................................................
c
      subroutine trailing(j,title)

      implicit none 
c
c     routine to kill trailing blanks -- actually return the index of
c     the last non-blank character
c
c Variables
c ---------
c j -- counter
c title -- general character string
c
      integer*4 j
 
c
      character*(*) title
 
c
c.... loop over string until last blank character found
      j = len(title)
      do while (title(j:j).eq.' ' )
*                                     ! get out if null string
         if( j.eq.1 ) return
         j = j - 1
      end do
c
      return
      end
 
