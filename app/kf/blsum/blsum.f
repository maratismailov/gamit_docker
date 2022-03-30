      program BLSUM

      implicit none 
c
c     This program will read the output file from a SOLVK solution
c     and produce a baseline summary.  The summary will be output
c     to one file and the values to another file.
c     The program can be run in batch so that it can be scheduled
c     after a batch solvk run.  In batch mode the runstring is
c     :RU,BLSUM,crt,printer,summary_file,values_file,solution1,solution2,
c          ...... up to maximum number of solutions
c
c Include files
c -------------
*                                  ! the BLSUM control common block
      include 'blsum.h'
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
* MOD TAH 990101: Changed data storage to use real*8 and integer*4 arrays
c
      integer*4 data(3,max_ent)
      real*8 vdata(3,max_ent)
 
c
      common /values/ vdata, data
 
c
c file names
c ----------
c more_files    -- logical to indicate that we are still reading files
c status_old  -- variable which gives the status old to file opening
c status_unk  -- varaible which gives the status unknown to file opening
c
      logical more_files
 
c
      character status_old*4, status_unk*8
 
c
*  i    - loop counter
*  ierr - iostat error on read/write
	
      integer*4 i, ierr
 
c
c Runstring parameters
c --------------------
c run_string -- the runstring as a character array
c count -- used to count our way through the runstring parameters
c len_run(max_runstring) - length of runstrings
* sort_string  -- Indicates we should sort data
* len_sort     -- length of sort
* rcpar        -- Gets the rustring
c
      integer*4 count, len_run(max_runstring), len_sort, rcpar
 
c
      character*128 run_string(max_runstring)

      character*16 sort_string
 
c
c
c Data statements for status
c --------------------------
      data status_old / 'old '/, status_unk / 'unknown' /
c
c.... Get the runstring and decide if we are in batch mode
*     Check for sort option
      len_sort = rcpar(1, sort_string)
      call casefold( sort_string )
      if( len_sort.le.0 ) then
          call proper_runstring('blsum.hlp', 'blsum', 1)
      else
          read(sort_string,*, iostat=ierr) min_out
          if( ierr.eq.0 ) then

*             if negative then tell program to sort and save
*             absolute values
              if( min_out.lt.0 ) then
                  sort = .true.
              else
                  sort = .false.
              end if
              min_out = abs(min_out)
          else     
*             check old version of sort
              if( sort_string(1:2).eq.'SO' ) then
                  sort = .true.
              else
                  sort = .false.
              end if
          end if
      end if

*                              ! read rest of runstring
      do i = 1,max_runstring
         len_run(i) = rcpar(i+1, run_string(i))

c...     see if zero string
         if( len_run(i).eq.0 ) run_string(i) = '  '
c
      end do
c
c.... see if runstring was passed
*                              ! we are in batch mode
      if( len_run(1).ne.0 ) then 
         batch = .true.
*                              ! interactive run
      else
         batch = .false.
      end if
c
c.... write headr message
      write(*,50)
  50  format(//" BLSUM: Program to extract baseline lengths from GLOBK",
     .   " solution",/,
     .         " -----------------------------------------------------",
     .   "---------",/)
c
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
                  
            write(*,'(a,i2,a,$)') 
     .         'Enter the name of the solution file #',num_files
     .        ,' (:: to end) '   
c           above replaces below to avoid Hollerith split for DEC (rwk 970920)
c            write(*,'(" Enter the name of solution file #",i2,
c     .         " (:: to end) ",$)') num_files
c
            call get_soln_file(solution_file(num_files))
c
*                                                             ! see if file
            if( solution_file(num_files)(1:2).ne.'::' ) then
c                                                              exists
               call check_soln_file( solution_file(num_files),ierr,
     .            100,status_old)
               close(100)
c
c....          increment the file counter if file exists
               if( ierr.eq.0 ) num_files = num_files + 1
c
c....          check maximum of solutions
*                                                 ! this end of list
               if( num_files.gt.max_files ) then
                  write(*,'(" Maximum number of files reached")')
                  more_files = .false.
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
            write(*,'(/" Name of output summary file ",$)')
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
            write(*,'(" Name of output values file  ",$)')
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
         write(*,'(" Reading baselines from ",a)') solution_file(i)
c
c....    use the check subroutine to open the solution file
         call check_soln_file(solution_file(i), ierr, 100,'old ')
         call read_soln(i,vdata, data)
         close(100)
c
      end do
c
*     Now sort the data into time order
      if( sort ) then
          call blsort(vdata, data)
      end if
 
c.... Now produce the summaries of the baseline lengths
      call sum_solution(vdata, data)
c
c.... Thats all
      close(200)
      close(201)
      end
c
c........................................................................
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
c........................................................................
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
c........................................................................
c
      subroutine read_soln(soln_num, vdata, data)

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
      include 'blsum.h'
c
c
c passed variables
c ----------------
c data -- ema storage of the baselines from the solution file
c
      integer*4 data(3,max_ent)
      real*8 vdata(3,max_ent)
 
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
 
      character  new_gtime*24
 
c.... Define the identifier strings
c
      data new_soln / ' SOLUTION ' /, 
     .     new_time / 'EXPERIMENT DATE :'    /,
     .     bas_to   / ' TO '     /
c
      data  new_gtime / 'SOLUTION REFERS TO     :' /
 
c.... Firstly read the header from solution
      read(100,'(a)',iostat=err) header(soln_num)
c
* MOD TAH 881031 Check to see if this is VALUES_FILE
      if( header(soln_num)(1:11).eq.'VALUES_FILE' ) then
          call read_values_file( soln_num, vdata, data )
*                      ! Do not need to do the rest
          return
      end if

*     See if name says that it is a values files
      if( solution_file(soln_num)(1:2).eq.'va' .or.
     .    index(solution_file(soln_num),'/va').gt.0 ) then
          rewind(100)
          call read_values_file( soln_num, vdata, data )
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
               call multiread(buffer,indx,'i4',err,num_soln,cdummy,1)
 
               if( num_soln.eq.0 .or. err.ne.0 ) then
C                 write(*,'(" Solution number error in ",a,/,
C    .                    " iel, indx, num_soln, err = ",4i8)')
C    .                buffer, iel, indx, num_soln, err
                  write(*,'(a,a,/,a,4i8)' ) 
     .                    ' Solution number error in ', buffer,
     .                    ' iel, indx, num_soln, err = ',
     .                    iel, indx, num_soln, err
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
* MOD TAH 150930: Allowed range for ' TO ' location due to changing
*           formats.
            if( ibas.ge.15 .and. ibas.le.16 
     .          .and. buffer(1:1).eq.' ' ) then
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
               call save_ema(data(1,num_ent), vdata(1,num_ent))
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
c........................................................................
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
      include 'blsum.h'
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
c........................................................................
c
      subroutine save_ema(site_dat, real_dat)

      implicit none 
c
c     routine to save the baseline values in ema
c
c Include files
c -------------
c
*                                   ! the control common block
      include 'blsum.h'
c
c Variables
c ---------
c site_dat -- integer values for sites and solution number
c real_dat -- Real*8 values to save

      integer*4 site_dat(3)
      real*8 real_dat(3)
 
c
c.... Save the values
      site_dat(1) = site(1)
      site_dat(2) = site(2)
      site_dat(3) = num_soln
c
      real_dat(1) = sol_epoch
      real_dat(2) = base_data(1)
      real_dat(3) = base_data(2)
c
      return
      end
c
c........................................................................
c
      subroutine sum_solution(vdata, data)

      implicit none 
c
c     Routine to read through the ema data and accumualte statistics
c     for each baseline.
c
c
c Include files
c -------------
c
*                                     ! the control common block
      include 'blsum.h'

* i,j,k,l,m  - Loop counters
* ierr       - IOSTAT error
* inf, inh   - lengths of file and header

      integer*4 i,j,k,m, inf, inh, ierr
c
c EMA variables
c -------------
c data -- the stored baseline data
c
      integer*4 data(3,max_ent)
      real*8 vdata(3,max_ent)
c
c
c.... start looping over all baselines.  Write out the headers to the
c     summary file
c
      write(200,100)
  100 format(/" BLSUM: GLOBK baseline summary: includes the following",
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
     .  "  nrms     slope    sig    wrms nrms   dur   mean    line #",
     . /,    "                                (m)       (m)   (mm)",
     .  "         (mm/yr)  (mm/yr)  (mm)       (yrs)  (yrs)")
 
c.... start the line counter for tplot file extraction
      line_st = 1
      line_en = 0
c
c.... loop over all baseline combinations
      do i = 1, num_site-1
         write( *,9000) names(i)
 9000    format(" Processing data from site ",a)
c
         do j = i+1, num_site
c
c
c....       loop over the solution
            do k = 1, max_soln
c
               call clear_norm
c
c....          Now get all the entries
               do m = 1,num_ent
c
c....             see if baseline and solution matches
                  if( data(1,m).eq.i .and. data(2,m).eq.j .and.
*                                             ! we found a value
     .                data(3,m).eq.k ) then
c
c....                accumulate statistics
                     call slope_acc(vdata(1,m),vdata(2,m),vdata(3,m))
c
                  end if
               end do
c
c....          complete the computation of the statistics
               call compute_slope(ierr)
c
c....          see if any data
*                                     ! yes there is data, so continue
               if( ierr.ge.0 .and. stat(4,1).ge.min_out ) then

*                 Clear the statistics
                  call clear_stat
c
c....             write the header    ! now that we know there is data
                  call write_header(201,i,j,k)
                  line_en = line_en + 3
                  stat(4,1) = 0.d0
c
c....             Now accumlate residual statistics
                  do m = 1,num_ent
c
c....                see if baseline and solution matches
                     if( data(1,m).eq.i .and. data(2,m).eq.j .and.
*                                               ! we found a value
     .                  data(3,m).eq.k ) then
c
c....                   write record to file
                        call write_rec(201,vdata(1,m),vdata(2,m),
     .                      vdata(3,m))
                        line_en = line_en + 1
c
c....                   accumulate statistics
                        call accum_stat(vdata(1,m),vdata(2,m),
     .                                  vdata(3,m))
c
                     end if
                  end do
c
                  call complete_stat
c
c....             Now write out the summary
                  if( stat(4,1).ge.min_out ) then
                       call write_summary(200,i,j,k,'short ')
                  end if 
                  call write_summary(201,i,j,k,'long')
c
c....             save line numbers
                  line_st = line_en + 5
                  line_en = line_st - 1
c
*                     ! producing summary if there was data
               end if
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
c........................................................................
c
      subroutine clear_stat

      implicit none 
c
c     This routine will clear the statistics accumulation and
c     slope computation  arrays
c
      include 'blsum.h'

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

      return
      end
c
c........................................................................
c
      subroutine clear_norm 

      implicit none 
c
c     This routine will clear the statistics accumulation and
c     slope computation  arrays
c
      include 'blsum.h'

* i,j  -- loop counters
      integer*4 i
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
c
c.... initialise the min and max times
      min_epoch = 1.d30
      max_epoch = -1.d30
c
      return
      end
c
c........................................................................
c
      subroutine write_header(unit,i,j,k)

      implicit none 
c
c     Routine to write the header records to the values file
c
      include 'blsum.h'
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
c........................................................................
c
      subroutine slope_acc(epoch,length,sigma)

      implicit none 
c
c     Routine to accumulate the data needed to compute the slope
c
      include 'blsum.h'
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
         stat(4,1) = 0
      end if
c
c.... see if this this value is before or after all others
      if( epoch.lt. min_epoch ) min_epoch = epoch
      if( epoch.gt. max_epoch ) max_epoch = epoch
c
c.... Now accumulate normal equations
      stat(4,1) = stat(4,1) + 1
      wgh = 1.d0/sigma**2
      dtime = (epoch-ref_epoch)/365.25d0
      norm_eq(1) = norm_eq(1) + wgh
      norm_eq(2) = norm_eq(2) + dtime*wgh
      norm_eq(3) = norm_eq(3) + dtime**2*wgh
c
      b_vec(1)   = b_vec(1) + (length-ref_length)*wgh
      b_vec(2)   = b_vec(2) + (length-ref_length)*dtime*wgh
c
      return
      end
c
c........................................................................
c
      subroutine compute_slope(ierr)

      implicit none 
c
c     Routine to compute the slope of best fit straight line
c
      include 'blsum.h'
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
c
         end if
      end if
c
      return
      end
c
c........................................................................
c
      subroutine accum_stat(epoch, length, sigma)

      implicit none 
c
c     Routine to accumulate the statistics from a set of baseline
c     results.  The statistics are accumulated for both the residuals
c     about the mean and about the best fit straight line
c
      include 'blsum.h'
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
c........................................................................
c
      subroutine complete_stat

      implicit none 
c
c     routine to complete the computation of the statistics
c
      include 'blsum.h'
c
c
c.... Firstly do the statistics about the weighted mean, if there were
c     any data at all.
c
      if( stat(4,1).gt.1 ) then
c
         nrms_mean = sqrt( stat(3,1)/(stat(4,1)-1) )
         wrms_mean = sqrt( stat(4,1)/stat(2,1) )*nrms_mean
c
*              ! only one data so we can not compute nrms scatter
      else
c
         nrms_mean = 1.d0
         wrms_mean = sqrt( stat(4,1)/stat(2,1) )*nrms_mean
c
      end if
c
c.... Now do slope statistics
      if( stat(4,2).gt.2 ) then
c
         nrms_slope = sqrt( stat(3,2)/(stat(4,2)-2) )
         wrms_slope = sqrt( stat(4,2)/stat(2,2) )*nrms_slope
c
*              ! only one data so we can not compute nrms scatter
      else
c
         nrms_slope = 1.d0
         wrms_slope = sqrt( stat(4,2)/stat(2,2) )*nrms_slope
c
      end if
c
c.... Thats all
      return
      end
c
c........................................................................
c
      subroutine write_rec(unit, epoch, length, sigma)

      implicit none 
c
c     Routine to write the data record to the values file
c
      include 'blsum.h'
c
c EMA values
c ----------
c epoch -- the epoch of the measurement
c length -- the measured baseline
c sigma  -- the sigma of the length
c
      real*8 epoch, length, sigma
 
c
c
c variables
c ---------
c unit -- the output unit number
* date -- Date of this measurement
c
      integer*4 unit, date(5)
 
c
c rel_epoch  -- main memory value of epoch
c sec_tag    -- seconds tag
c
      real*8 rel_epoch, sec_tag, dtime, residual
 
c
c.... set up for output, convert epoch to yr,mon,day,hr,min
      rel_epoch = epoch
      call jd_to_ymdhms(rel_epoch, date, sec_tag)
      dtime = (epoch-ref_epoch)/365.25d0

c
c.... do statistics about the mean
      residual = length - (ref_length+intercept) - dtime*slope
c
c
c.... write the record
      write(unit,100) date,length,sigma, residual, sigma
  100 format(i5,4i3,f15.4,1x,f10.4,1x,2f10.4)
c
      return
      end
c
c........................................................................
c
      subroutine write_summary(unit,site_1,site_2, soln, format)

      implicit none 
c
c     Routine to write the summary information.  This can be written
c     in two forms.. the short where all results are put on a single
c     line or the long form which takes multiple lines
c
      include 'blsum.h'
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
      if( format.eq.'short' ) then
c
         write(unit,100) names(site_1),names(site_2), soln, num_len,
     .      mean_out,sig_mean_scaled, wrms_mean_out, nrms_mean,
     .      slope_out, sig_slope_scaled, wrms_slope_out, nrms_slope,
     .      duration, mean_year, line_st, line_en
 100     format(a8,"-",a8,1x,i2,1x,i4,1x,
     .      f13.4,1x,f6.4,1x,f5.1,1x,f5.2,3x,
     .      f8.2,1x,f8.2,1x,f5.1,1x,f5.2,1x,
     .      f6.2,1x,f7.2,1x,i7,"-",i7)
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
c........................................................................
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
 
CTITLE READ_VALUES_FILE
 
      subroutine read_values_file( soln_num, vdata, data )

      implicit none 
 
 
*     Routine to read a values file to extract information about the
*     baseline lengths.
 
      include 'blsum.h'
 
*         soln_num      - Number of current file
*         data(num_field, max_ent)   - Place where data is actually
*                       - saved.
 
      integer*4 soln_num, data(3, max_ent)
      real*8 vdata(3,max_ent)
 
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
         if( new_buffer(1:3).eq.' N ' ) call sub_char(new_buffer,' N ',
     .                                               ' NORTH   ')
         if( new_buffer(1:3).eq.' E ' ) call sub_char(new_buffer,' E ',
     .                                               ' EAST    ')
         if( new_buffer(1:3).eq.' U ' ) call sub_char(new_buffer,' U ',
     .                                               ' HEIGHT  ')
 
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
         call read_line(new_buffer,indx,'i4',ierr,num_soln,cdummy)
         if( ierr.ne.0 ) then
*            check for LNEU
             if( index(new_buffer(10:),'L').gt.0 ) num_soln = 1
             if( index(new_buffer(10:),'N').gt.0 ) num_soln = 2
             if( index(new_buffer(10:),'E').gt.0 ) num_soln = 3
             if( index(new_buffer(10:),'U').gt.0 ) num_soln = 4
         end if
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
                 call multiread(buffer, indx, 'i4', ierr, clk_date,
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
                 call save_ema(data(1,num_ent), vdata(1,num_ent))
 
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
      write(*,300) num_ent, num_site, max_soln
  300 format(' There are ',i8,' baselines entries between ',
     .                     i3,' sites, and ',I3,' solns' )
      if( max_soln.eq.0 ) max_soln = 1
*
      return
      end
 
CTITLE BLSORT
 
      subroutine blsort(vdata, data)

      implicit none 
 
 
*     Routine to sort the baseline data into time order.
 
      include 'blsum.h'
 
*         data(num_field,max_ent)   - Data to be time sorted. Time
*                                   - is in elemnet 4-7.
*   i,j,k           - Loop counters
*   temp(num_field) - Temporary storage for switch
 
      integer*4 data(3,max_ent), i,j,k, temp(3)
      real*8 vdata(3,max_ent), vtemp(3)
      
 
*   bepochj, bepochnext  - Epochs of the two test values
*   ema8             - Returns a real*8 value from EMA
 
      real*8 bepochj, bepochnext
 
 
****  Loop doing a bubble sort
      write(*,100) num_ent
  100 format(' Starting bubble sort with ',i5,' entries')
 
      do i = 1, num_ent - 1
          do j = 1, num_ent - 1
 
****          Get the two epochsc
              bepochj = vdata(1,j)
              bepochnext = vdata(1,j+1)
 
*                                                ! Switch values
              if( bepochj.gt.bepochnext ) then
 
                  do k = 1, 3
                      temp(k) = data(k,j)
                      vtemp(k) = vdata(k,j)
                      data(k,j) = data(k,j+1)
                      vdata(k,j) = vdata(k,j+1)
                      data(k,j+1) = temp(k)
                      vdata(k,j+1) = vtemp(k)
                  end do
              end if
          end do
      end do
 
***** Thats all
      return
      end
 
 
