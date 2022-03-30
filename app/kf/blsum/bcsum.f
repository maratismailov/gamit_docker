      program BCSUM

      implicit none 
c
*     This program is an extension of BLSUM.  Here baseline components
*	are also read and a new file (the sumcomp) file is produced.  The
*	values file also has its definition modified to include the
*	components as well as the lengths.

c     This program will read the output file from a SOLVK solution
c     and produce a baseline summary.  The summary will be output
c     to one file and the values to another file.
c     The program can be run in batch so that it can be scheduled
c     after a batch solvk run.  In batch mode the runstring is
c     % BCSUM, SO summary_file values_file sum_comp solution1,solution2,
c          ...... up to maximum number of solutions
c
c Include files
c -------------
*                                  ! the BCSUM control common block
      include 'bcsum.h'
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
c        data(7:8 )  = length (m)
c        data(9:10) = sigma (m)
c        data(11:12) = North (m)
c        data(13:14) = Sigma (m)
c        data(15:16) = East (m)
c        data(17:18) = sigma (m)
c        data(19:20) = Up (m)
c        data(21:22) = sigma (m)
c
      integer*4 data(num_field,max_ent)
 
c
      common /values/ data
 
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

      character*4 sort_string
 
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
          call proper_runstring('bcsum.hlp', 'bcsum', 1)
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
  50  format(//' BCSUM: Program to extract baseline lengths and ',
     .         ' components from GLOBK solution',/,
     .         ' -----------------------------------------------',
     .         '-------------------------------',/)
c
c
c.... if we are in batch mode assign the file names from the runstring
      if( batch ) then
         summary_file     = run_string(1)
         values_file      = run_string(2)
         sumcomp_file     = run_string(3)
c
c....    Now get the solution files from the runstring
         num_files = 0
         count = 4
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
         call check_soln_file(sumcomp_file, ierr,202,status_unk)
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
                    
            write(*,'(a,i2,$)') 
     .           ' Enter the name of the solution file # ',num_files
c           the above replaces below
c            write(*,'(" Enter the name of solution file #",i2,
c     .         " ",$)') num_files
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

c....    get the sumcomp file name
         ierr = 1
c
         do while ( ierr.ne.0 )
            write(*,'(" Name of output component summary file  ",$)')
            call get_soln_file(sumcomp_file)
c
            if( sumcomp_file(1:2).eq.'::' ) then
               stop ' :: in response to values file name'
            end if
c
            call check_soln_file(sumcomp_file,ierr,202,
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
         call read_soln(i,data)
         close(100)
c
      end do
c
*     Now sort the data into time order
      if( sort ) then
          call blsort(data)
      end if
 
c.... Now produce the summaries of the baseline lengths
      call sum_solution(data)
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
      include 'bcsum.h'
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
      real*8 values(10)
 
c
*   indx        - Pointer for use in multiread
*   err         - Error flag returned from multiread
*   iel         - index position of string
*   iclk        - Index position for clock time
*   ibas        - index position for baseline length
*   in          - length of string
*   prev_num_ent - Number of entries at the time that a new solution
*                  is found.  (used so that as components are found
*                  they can be assigned to the correct baseline)
*   i           - Loop counter
*   posm        - Position of '-' in line
*   lenb        - Length of string read from file

      integer*4 indx, err, iel, iclk, ibas, in, prev_num_ent, i,
     .          posm, lenb, jel, trimlen
 
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

      prev_num_ent = 0
 
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

*              save the current number of entries
               prev_num_ent = num_ent
            end if
 
*           Check clock time for a global solution
            iclk = index(buffer, new_gtime)
*                                 ! New clock epoch
            if( iclk.gt.0 ) then
               read( buffer(iclk+27:),150, iostat=err)
     .            clk_date
 150           format(5(i2,1x))
               call ymdhms_to_jd( clk_date, 0.d0, sol_epoch )
               prev_num_ent = num_ent
            end if
c
c....       Now see if we have a baseline length
            ibas = index(buffer, bas_to)
*                                   ! we found a baseline length
!            if( ibas.eq.15 ) then
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
               base_data(1) = values(1)
               base_data(2) = values(3)
               if( base_data(2).eq.0 ) base_data(2) = 0.1d0    
c
c....          Make sure that the site numbers are in ascending order
               if( site(2).lt.site(1) ) then
                  iel = site(1)
                  site(1) = site(2)
                  site(2) = iel
               end if
c
c.....         Save this data in ema
* MOD TAH 930219.  If sigma is too big just ignore this line.
               if( base_data(2).lt.0.5d0 ) then
                  num_ent = num_ent + 1
c
c.....            see if maximun number of entries exceeded
                  if( num_ent.ge.max_ent ) then
                     write(*,200) num_ent
  200                format(/," Maximum number of baselines ",i10,
     .                  " has been exceeded")
*                                 ! Force exit on next read
                     err = -99
                  end if
c
                  call save_ema(data(1,num_ent), data(5,num_ent),
     .              data(7,num_ent),data( 9,num_ent))
               end if
c
c....       Thats all
*                         ! baseline found
            end if

*           See if we can find the baseline components
            if( index(buffer,'BASELINE COMPONENTS').gt.0 .and.
     .          err.eq.0 ) then

*               we have found the line for baseline components
*               Skip the next three lines
                read(100,'(a)', iostat=err) buffer
                read(100,'(a)', iostat=err) buffer
                read(100,'(a)', iostat=err) buffer

*               Loop until we find a blank line
                lenb = 1
                do while ( lenb.gt.0 )
                    read(100,'(a)', iostat=err) buffer
                    lenb = trimlen(buffer)
                    if ( err.ne.0 ) lenb = 0
*                   Clear the - from the baseline name
                    posm = index(buffer(8:),'-')
                    if( posm.gt.0 ) buffer(posm+7:posm+7) = ' '
                    call strip_name(buffer, names, num_site, jel,0)
                    call strip_name(buffer, names, num_site, iel,0)

*                   now read out the values
                    indx = 1
                    if( lenb.gt.0 ) then
                       call multiread(buffer,indx,'R8',err,values,
     .                                cdummy,10)

*                      Now see if we can the baseline
                       i = prev_num_ent
                       do while (i.lt. num_ent ) 
                          i = i + 1
                          if( iel.eq.data(1,i).and.
     .                        jel.eq.data(2,i) ) then
                              call save_comp(data(11,i), values, 1)
*                                            ! Force loop exit
                             i = num_ent
                          end if
                          if( jel.eq.data(1,i).and.
     .                        iel.eq.data(2,i) ) then
                              call save_comp(data(11,i), values,-1)
*                                            ! Force loop exit
                              i = num_ent
                         end if
                       end do
*                         ! Buffer is non-zero length
                    end if
*                         ! Loop over reading file for components
                 end do
*                         ! components found

             end if
                
*                         ! no error on file read
         end if
c
*                         ! looping over the file to be read
      end do
c
c.... Write out number of entries
      write(*,300) num_ent, num_site
  300 format(" There are ",i10," baseline entries between ",i6," sites")
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
*                                       ! the BCSUM parameter file
      include 'bcsum.h'
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
      if( num_site.gt. max_site ) then
          write(*,100) max_site
 100      format(/,'*** DISASTER *** Maximum number of sites ',
     .             I8,' has bee exceeded',/,
     .             'Increase max_site in bcsum.h and remake')
          stop 'Too Many Sites'
      endif

      names(num_site) = buffer(i:i+8)
c
      return
      end
c
c........................................................................
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
      include 'bcsum.h'
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
c........................................................................
c
      subroutine sum_solution(data)

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
      include 'bcsum.h'

* i,j,k,l,m  - Loop counters
* ierr       - IOSTAT error
* inf, inh   - lengths of file and header

      integer*4 i,j,k,m, inf, inh, ierr
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
      write(202,100)
  100 format(/" BCSUM: GLOBK baseline component summary: includes the",
     .   " following solutions")
c
      do i = 1, num_files
         call trailing(inf, solution_file(i))
         call trailing(inh,header(i))
         write(200,120) solution_file(i)(1:inf), header(i)(:inh)
         write(202,120) solution_file(i)(1:inf), header(i)(:inh)
  120    format(" File ",a,": ",a)
      end do
c
c.... Write out descriptive header
      write(200,140)
 140  format("    Baseline       S#   #    mean length   sig   wrms",
     .  "  nrms     slope    sig    wrms nrms   dur   mean    line #",
     . /,    "                                (m)       (m)   (mm)",
     .  "         (mm/yr)  (mm/yr)  (mm)       (yrs)  (yrs)")
      write(202,150)
 150  format('    Baseline      Type  #      Length       Type Mean ',
     .  '   +-    wrms  nrms        slope    +-    wrms  nrms   dur',
     .  ' mean',/,
     .  '                                (km)         (mm)      (mm)',
     .  '    (mm)          (mm/yr)  (mm/yr)   (mm)          (yrs)')
 
c.... start the line counter for tplot file extraction
      line_st = 1
      line_en = 0
c
c.... loop over all baseline combinations
      do i = 1, num_site-1
         write( *,200 ) names(i)
 200     format(" Processing data from site ",a)
c
         do j = i+1, num_site
c
c
            call clear_norm

c....       Now get all the entries
            do m = 1,num_ent
c
c....          see if baseline and solution matches
               if( data(1,m).eq.i .and. data(2,m).eq.j ) then 
c
c....             accumulate statistics
                  call slope_acc(data(5,m),data(7 ,m),data( 9,m),1)
                  call slope_acc(data(5,m),data(11,m),data(13,m),2)
                  call slope_acc(data(5,m),data(15,m),data(17,m),3)
                  call slope_acc(data(5,m),data(19,m),data(21,m),4)
c
               end if
            end do
c
c....       complete the computation of the statistics
            call compute_slope(ierr,1)
            call compute_slope(ierr,2)
            call compute_slope(ierr,3)
            call compute_slope(ierr,4)
c
c....       see if any data
*                                     ! yes there is data, so continue
            if( ierr.ge.0 .and. stat(4,1,1).ge.min_out ) then

*               Reset the number of statistics
                call clear_stat    

*               Loop over all four data types
                do k = 1, 4
c
c....             write the header    ! now that we know there is data
                  call write_header(201,i,j,k)
                  line_en = line_en + 3
c
c....             Now accumlate residual statistics
                  do m = 1,num_ent
c
c....                see if baseline and solution matches
                     if( data(1,m).eq.i .and. data(2,m).eq.j ) then
c
c....                   write record to file
                        call write_rec(201,data(5,m),
     .                                 data(7+(k-1)*4,m),k)
                        line_en = line_en + 1
c
c....                   accumulate statistics
                        call accum_stat(data(5,m),
     .                                  data(7+(k-1)*4,m),
     .                                  data(9+(k-1)*4,m),k)
c
                     end if
                  end do
c
                  call complete_stat(k)
c
c....             Now write out the summary
                  if( stat(4,1,1).ge.min_out ) then
                      call write_summary(200,i,j,k,'short ',k)
                      if( k.eq.4 ) call write_sumcomp(202,i,j,k)
                  end if
                  call write_summary(201,i,j,k,'long',k)
c
c....             save line numbers
                  line_st = line_en + 5
                  line_en = line_st - 1
c
*                     ! producing summary if there was data
*                     ! looping over type   -- k
               end do 
*                     ! ierr was zero on compute slope
            end if 
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
      include 'bcsum.h'

* i,j  -- loop counters
      integer*4 i,j,k
c
c.... clear all of the arrays
      do k = 1,4
         do i = 1,4
            do j = 1,2
*                              ! summation variables for mean and rate
               stat(i,j,k) = 0.d0
c                                statistics.
            end do
         end do
      end do
c
      return
      end
c
c........................................................................
c
      subroutine clear_norm

      implicit none 
c
c     This routine will clear the normal equations and
c     slope computation  arrays
c
      include 'bcsum.h'

* i,j  -- loop counters
      integer*4 i,k
c
c.... clear all of the arrays
      do k = 1,4
         do i = 1,3
*                              ! normal equations for slope computation
            norm_eq(i,k) = 0.d0
         end do
c
         do i = 1,2
*                              ! b_vector slope computation
            b_vec(i,k)   = 0.d0
         end do
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
      include 'bcsum.h'
c
c Variables
c ---------
c unit -- output unit number
* in   -- Length of string
c
      integer*4 unit, i, j, k, in

* type -- Type of quanity output (L,N,E,U -- Based on k value)

      character*4 type(4)

      data type / 'L', 'N', 'E', 'U' /
 
c
c.... The header record to be used by tplot is the solution
c     description, baseline + solution #, blank record
c
c.... delete trailing blanks from header
      call trailing(in,header(1))
      write(unit,'(a)') header(1)(:in)
c
      write(unit,'(a," to ",a," Solution ",a1)') names(i),names(j),
     .          type(k)
c
      write(unit,'(1x)')
c
      return
      end
c
c........................................................................
c
      subroutine slope_acc(epoch,length,sigma, it)

      implicit none 
c
c     Routine to accumulate the data needed to compute the slope
c
      include 'bcsum.h'
c
c Variables
c ---------
c epoch -- the epoch of this measurement
c length -- the length of the baseline
c sigma  -- the sigma of the baseline
* it     -- Index in stats array to use 1 -baseline length
*                                       2 -North
*                                       3 -East
*                                       4 -Up
c
      integer*4 it
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
      if ( norm_eq(1,it).eq.0.d0 ) then
         ref_epoch = epoch
         ref_length(it) = length
         stat(4,1,it) = 0
      end if

*     Check to see if this value is discordant by about 10.0 m
*     if it is then we we have slipped a 10 m interger in the
*     compoent.  NOTE: Units are mm at the moment.  
      if( length - ref_length(it).lt.-9.0d3) length = length + 10.d3
      if( length - ref_length(it).gt. 9.0d3) length = length - 10.d3
c
c.... see if this this value is before or after all others
      stat(4,1,it) = stat(4,1,it) + 1
      if( epoch.lt. min_epoch ) min_epoch = epoch
      if( epoch.gt. max_epoch ) max_epoch = epoch
c
c.... Now accumulate normal equations if the sigma is greater than
*     zero (In case we have some missing data)
      if( sigma.gt.0 ) then
         wgh = 1.d0/sigma**2
         dtime = (epoch-ref_epoch)/365.25d0
         norm_eq(1,it) = norm_eq(1,it) + wgh
         norm_eq(2,it) = norm_eq(2,it) + dtime*wgh
         norm_eq(3,it) = norm_eq(3,it) + dtime**2*wgh
c
         b_vec(1,it)   = b_vec(1,it) + (length-ref_length(it))*wgh
         b_vec(2,it)   = b_vec(2,it) + (length-ref_length(it))*dtime*wgh
      end if
c
      return
      end
c
c........................................................................
c
      subroutine compute_slope(ierr,it)

      implicit none 
c
c     Routine to compute the slope of best fit straight line
c
      include 'bcsum.h'
c
c Local variables
c ---------------
c det -- the determinate of the normal equations
c inverse -- the inverse of the normal equations
c ierr -- an error number returned to main routine
c      = 0 if all OK
c      = -1 if no data
c      = +1 if slope inversion is singular
* it   - Type of quantity to compute slope with
* sig_inter - sigma of intersect
c
      real*8 det, inverse(3), sig_inter

* ierr  - IOSTAT error

      integer*4 ierr, it
 
c
c.... see if we have enough data to do a mean value
      ierr = 0
      if( norm_eq(1,it).lt. tol ) then
         wgh_mean(it) = 0.d0
         sig_mean(it) = 0.d0
c
         intercept(it) = 0.d0
         sig_inter     = 0.d0
         slope(it)     = 0.d0
         sig_slope(it) = 0.d0
c
*                        ! probably no data
         ierr = -1
c
*                 ! compute mean
      else
c
c....    Now get weighted mean
         wgh_mean(it)  = b_vec(1,it)/norm_eq(1,it)
         sig_mean(it)  = sqrt(1.d0/norm_eq(1,it))
c
c....    Invert the normal equations
         det = norm_eq(1,it)*norm_eq(3,it) - norm_eq(2,it)**2
c
c....    compute the mean epoch of the data  ! julian date
         mean_epoch = ref_epoch + 
     .                (norm_eq(2,it)/norm_eq(1,it))*365.25d0
c
c....    see if singular
*                                ! not enough data to do slope
         if( det.lt. tol ) then
            intercept(it) = wgh_mean(it)
            sig_inter     = sig_mean(it)
c
            slope(it)     = 0.d0
            sig_slope(it) = 0.d0
c
*                        ! matrix singular -- only one data
            ierr = +1
c
*                   ! compute the slope
         else
c
            inverse(1) =  norm_eq(3,it)/det
            inverse(2) = -norm_eq(2,it)/det
            inverse(3) =  norm_eq(1,it)/det
c
c....       compute intercept and slope
*                                                                  ! meters
            intercept(it) = inverse(1)*b_vec(1,it) 
     .                    + inverse(2)*b_vec(2,it)
*                                                                  ! m/year
            slope(it)     = inverse(2)*b_vec(1,it) 
     .                    + inverse(3)*b_vec(2,it)
c
            sig_inter     = sqrt(inverse(1))
            sig_slope(it) = sqrt(inverse(3))
c
         end if
      end if
c
      return
      end
c
c........................................................................
c
      subroutine accum_stat(epoch, length, sigma, it)

      implicit none 
c
c     Routine to accumulate the statistics from a set of baseline
c     results.  The statistics are accumulated for both the residuals
c     about the mean and about the best fit straight line
c
      include 'bcsum.h'
c
c EMA variables
c -------------
c epoch -- epoch of measuarement
c length -- baseline length
c sigma  -- sigma of length
* it     -- Type of quantity to accumuate

      integer*4 it
 
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

* use_point -- Logical which indicates whether or not we
*              use this point.
C     logical use_point

*     See if we should use this point
C     use_point = .false.
C     if( it.eq.1 ) then
C         if( sigma.gt.0 .and. sigma.lt.0.5 ) use_point = .true.
C     else
C         if( sigma.gt.0 .and. sigma.lt.500.) use_point = .true.
C     end if 
c
c.... Compute the residuals and sum into the statistics
      if( sigma.gt.0 ) then
C     if( use_point ) then
         wgh = 1.d0/sigma**2
         dtime = (epoch-ref_epoch)/365.25d0
c
c....    do statistics about the mean
         residual = length - (ref_length(it)+wgh_mean(it))
c
         stat(1,1,it) = stat(1,1,it) + residual*wgh
         stat(2,1,it) = stat(2,1,it) + 1.d0*wgh
         stat(3,1,it) = stat(3,1,it) + residual**2*wgh
         stat(4,1,it) = stat(4,1,it) + 1.d0
c
c....    statistics about the slope
         residual = length - (ref_length(it)+intercept(it)) 
     .                     - dtime*slope(it)
c
         stat(1,2,it) = stat(1,2,it) + residual*wgh
         stat(2,2,it) = stat(2,2,it) + 1.d0*wgh
         stat(3,2,it) = stat(3,2,it) + residual**2*wgh
         stat(4,2,it) = stat(4,2,it) + 1.d0
      end if
c
c.... Thats all
      return
      end
c
c........................................................................
c
      subroutine complete_stat(it)

      implicit none 
c
c     routine to complete the computation of the statistics
c
      include 'bcsum.h'

* it - Type of quantity
      integer*4 it
c
c
c.... Firstly do the statistics about the weighted mean, if there were
c     any data at all.
c
      if( stat(4,1,it).gt.1 ) then
c
         nrms_mean(it) = sqrt( stat(3,1,it)/(stat(4,1,it)-1) )
         wrms_mean(it) = sqrt( stat(4,1,it)/stat(2,1,it) )
     .                           *nrms_mean(it)
c
*              ! only one data so we can not compute nrms scatter
      else
c
         nrms_mean(it) = 1.d0
         wrms_mean(it) = sqrt( stat(4,1,it)/stat(2,1,it) )
     .                         *nrms_mean(it)
c
      end if
c
c.... Now do slope statistics
      if( stat(4,2,it).gt.2 ) then
c
         nrms_slope(it) = sqrt( stat(3,2,it)/(stat(4,2,it)-2) )
         wrms_slope(it) = sqrt( stat(4,2,it)/stat(2,2,it) )
     .                         *nrms_slope(it)
c
*              ! only one data so we can not compute nrms scatter
      else
c
         nrms_slope(it) = 1.d0
         wrms_slope(it) = sqrt( stat(4,2,it)/stat(2,2,it) )
     .                         *nrms_slope(it)
c
      end if
c
c.... Thats all
      return
      end
c
c........................................................................
c
      subroutine write_rec(unit, epoch, data, type)

      implicit none 
c
c     Routine to write the data record to the values file
c
      include 'bcsum.h'
c
c EMA values
c ----------
c epoch -- the epoch of the measurement
c data(8) -- Length +- and components NEU with +-
c
      real*8 epoch, data(2)
 
c
c
c variables
c ---------
c unit -- the output unit number
* date -- Date of this measurement
* type -- Set type of data (k=1 is len in m and k2,3,4 are mm).  mm
*         are converted to meters
c
      integer*4 unit, date(5), type
 
c
c rel_epoch  -- main memory value of epoch
c sec_tag    -- seconds tag
* data_out(2) -- Output data (converted from mm to m)
c
      real*8 rel_epoch, sec_tag, data_out(2), dtime, residual 
 
c
c.... set up for output, convert epoch to yr,mon,day,hr,min
      rel_epoch = epoch
      call jd_to_ymdhms(rel_epoch, date, sec_tag)
      dtime = (epoch-ref_epoch)/365.25d0
c
c.... do statistics about the mean
      residual = data(1) - (ref_length(type)+intercept(type)) 
     .                   - dtime*slope(type)
c
c.... write the record
      if( type.eq.1 ) then
          data_out(1) = data(1) 
          data_out(2) = data(2) 
      else
          data_out(1) = data(1) / 1000.d0
          data_out(2) = data(2) / 1000.d0
          residual    = residual /1000.d0
      end if

      write(unit,100) date,data_out, residual, data_out(2)
  100 format(i5,4i3,f15.4,1x,f10.4,1x,2f10.4)
c 100 format(i5,4i3,f15.4,1x,f10.4,1x,3(F8.1,1x,f8.1,2x))
c
      return
      end
c
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
 
      subroutine read_values_file( soln_num, data )

      implicit none 
 
 
*     Routine to read a values file to extract information about the
*     baseline lengths.
 
      include 'bcsum.h'
 
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
      real*8 values(10)
 
c
*   indx        - Pointer for use in multiread
*   err         - Error flag returned from multiread
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   i           - Loop counter
*   iel, in     - Index positon temporary storgae
*   prev_num_ent - Previous number entries.   USed so that when
*                  components are found they can be put in the
*                  correct places
 
      integer*4 indx, ierr, trimlen, i, iel, in, prev_num_ent
 
*        still_reading   - Indicates that we are still reading entries
 
      logical still_reading
 
*   cdummy      - Dummy character string for multiread
 
      character*1 cdummy

      prev_num_ent = 0

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
*        If the length has been found, then this is a new
*        baseline so save the prev_num_ent
         if( num_soln.eq.1 ) then
             prev_num_ent = num_ent
         end if
         num_ent = prev_num_ent

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
 
c.....           Save this data in ema. 
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
*                If more than the length then convert to 
*                millimeters now.
                 if( num_soln.gt.1 ) then
                     base_data(1) = base_data(1)*1000.d0
                     base_data(2) = base_data(2)*1000.d0
                 end if

                 iel = 7 + (num_soln-1)*4
                 call save_ema(data(1,num_ent), data(5,num_ent),
     .                        data(iel,  num_ent), 
     .                        data(iel+2,num_ent) )
 
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
  300 format(' There are ',i10,' baselines entries between ',
     .                     i6,' sites' )
*
      return
      end
 
CTITLE BLSORT
 
      subroutine blsort(data)
 
      implicit none 

 
*     Routine to sort the baseline data into time order.
 
      include 'bcsum.h'
 
*         data(num_field,max_ent)   - Data to be time sorted. Time
*                                   - is in elemnet 4-7.
*   i,j,k           - Loop counters
*   temp(num_field) - Temporary storage for switch
 
      integer*4 data(num_field,max_ent), i,j,k, temp(num_field)
 
*   bepochj, bepochnext  - Epochs of the two test values
*   ema8             - Returns a real*8 value from EMA
 
      real*8 bepochj, bepochnext, ema8
 
 
****  Loop doing a bubble sort
      write(*,100) num_ent
  100 format(' Starting bubble sort with ',i5,' entries')
 
      do i = 1, num_ent - 1
          do j = 1, num_ent - 1
 
****          Get the two epochsc
              bepochj = ema8(data(5,j))
              bepochnext = ema8(data(5,j+1))
 
*                                                ! Switch values
              if( bepochj.gt.bepochnext ) then
 
                  do k = 1, num_field
                      temp(k) = data(k,j)
                      data(k,j) = data(k,j+1)
                      data(k,j+1) = temp(k)
                  end do
              end if
          end do
      end do
 
***** Thats all
      return
      end
 
CTITLE EMA8
 
 
      real*8 function ema8(data)

      implicit none 
  
      real*8 data
 
 
****  Just set ema8 equal to data

      ema8 = data

      return
      end

CTITLE SAVE_COMP

      subroutine save_comp( data, values, sign )

      implicit none 

*     Save the component values in the data array
 
* sign - sign to be applied (+-1) (Depends on direction of baseline)
* data(6) - The row of the data array for this baseline
* values(10) - 10 values read from component line

      integer*4 sign

      real*8 data(6), values(10)

*     Asign the values (North +- East +- Up +-)
      data(1) = values(1)*sign
      data(2) = values(3)
      if( data(2).eq.0 ) data(2) = 0.1d0
      data(3) = values(4)*sign
      data(4) = values(6)
      if( data(4).eq.0 ) data(4) = 0.1d0
      data(5) = values(8)*sign
      data(6) = values(10) 
      if( data(6).eq.0 ) data(6) = 0.1d0

****  Thats all
      return
      end

      subroutine write_sumcomp(unit,site_1,site_2, soln)

      implicit none 
c
c     Routine to write the summary information.  This can be written
c     in two forms.. the short where all results are put on a single
c     line or the long form which takes multiple lines
c
      include 'bcsum.h'
c
c variables
c ---------
c unit -- the output unit number
c site_1 -- first site in baseline
c site_2 -- second site in baseline
c soln -- solution number
* it  -- Output type
c
      integer*4 unit, site_1, site_2, soln, it
 
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
* mean_len  -- Mean length of the baseline
c
      real*8 sig_mean_scaled, sig_slope_scaled, slope_out, mean_out,
     .    wrms_mean_out, wrms_slope_out, duration, mean_year, sec_tag,
     .    mean_len

* num_len -- number of values in mean length estimate (used so we can output
*            with an I4 format
 
      integer*4 date(5), num_len 

* ctype(4) -- Type of output as a string (LNEU)

      character*1 ctype(4)

      data ctype / 'L','N','E','U' /   

      mean_len = 0
c
c.... Convert the units form meters to millimeters where necessary
c     and scale the sigmas of the various quanities
c
*                                  ! scale wgh_mean sigma
***** Loop over the four component types Length, N, E, U
      do it = 1, 4
         if( nrms_mean(it).gt.1 ) then
            sig_mean_scaled = sig_mean(it)*nrms_mean(it)
         else
            sig_mean_scaled = sig_mean(it)
         end if
c
c....    Now do the slope quantities
*                                  ! scale slope sigma
         if( nrms_slope(it).gt.1 ) then
*                                                         ! mm/yr
            sig_slope_scaled = sig_slope(it)*nrms_slope(it)
         else
*                                                         ! mm/yr
            sig_slope_scaled = sig_slope(it)
         end if
c
c....    Convert units on slope
         slope_out = slope(it)
         mean_out  = ref_length(it) + wgh_mean(it)
         if( it.eq.1 ) mean_len = mean_out
c
         wrms_mean_out = wrms_mean(it)
         wrms_slope_out = wrms_slope(it)
c
c....    get the duration of the data
         duration = (max_epoch - min_epoch)/365.25
         call jd_to_ymdhms( mean_epoch, date, sec_tag)
         mean_year = date(1) + date(2)/12.d0 + date(3)/365.25

*        Special code for length to convert from m to mm
         if( it.eq.1 ) then
            mean_out         = mean_out         * 1000.d0
            sig_mean_scaled  = sig_mean_scaled  * 1000.d0
            wrms_mean_out    = wrms_mean_out    * 1000.d0
            slope_out        = slope_out        * 1000.d0
            sig_slope_scaled = sig_slope_scaled * 1000.d0
            wrms_slope_out   = wrms_slope_out   * 1000.d0
         end if
c
c....    Now do the output
         num_len = stat(4,1,it)
c
         write(unit,100) names(site_1),names(site_2), it,ctype(it), 
     .      num_len, mean_len/1.d3, mean_out,sig_mean_scaled, 
     .      wrms_mean_out, nrms_mean(it),
     .      slope_out, sig_slope_scaled, wrms_slope_out, nrms_slope(it),
     .      duration, mean_year
 100     format(a8,"-",a8,1x,i2,1x,a1,1x,i4,1x,F10.3,
     .      f15.1,1x,f6.1,1x,f6.1,1x,f5.2,3x,
     .      f8.2,1x,f8.2,1x,f6.1,1x,f5.2,1x,
     .      f6.2,1x,f7.2)
c
      end do
c
      return
      end

CTITLE WRITE_SUMMARY

      subroutine write_summary(unit,site_1,site_2, soln, format,it)

      implicit none 
c
c     Routine to write the summary information.  This can be written
c     in two forms.. the short where all results are put on a single
c     line or the long form which takes multiple lines
c
      include 'bcsum.h'

c variables
c ---------
c unit -- the output unit number
c site_1 -- first site in baseline
c site_2 -- second site in baseline
c soln -- solution number
c format -- determines if long or short form
c it     -- Quantity type
c
      integer*4 unit, site_1, site_2, soln, it
 
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
* scale     -- Scale to be used for conversions from m to mm or to keep
*              as mm.
c
      real*8 sig_mean_scaled, sig_slope_scaled, slope_out, mean_out,
     .    wrms_mean_out, wrms_slope_out, duration, mean_year, sec_tag, 
     .    scale 

* num_len -- number of values in mean length estimate (used so we can output
*            with an I4 format
 
      integer*4 date(5), num_len 

*     Depending on type set the scale converiosn
      if( it.eq.1 ) then
          scale = 1000.0d0
      else
          scale = 1.d0
      end if
 
c
c.... Convert the units form meters to millimeters where necessary
c     and scale the sigmas of the various quanities
c
*                                  ! scale wgh_mean sigma
      if( nrms_mean(it).gt.1 ) then
         sig_mean_scaled = sig_mean(it)*nrms_mean(it)*scale/1000.d0
      else
         sig_mean_scaled = sig_mean(it)*scale/1000.d0
      end if
c
c.... Now do the slope quantities
*                                  ! scale slope sigma
      if( nrms_slope(it).gt.1 ) then
*                                                         ! mm/yr
         sig_slope_scaled = sig_slope(it)*nrms_slope(it)*scale
      else
*                                                         ! mm/yr
         sig_slope_scaled = sig_slope(it)*scale
      end if
c
c.... Convert units on slope
      slope_out = slope(it)*scale
      mean_out  = (ref_length(it) + wgh_mean(it))*scale/1000.d0
c
      wrms_mean_out = wrms_mean(it)*scale
      wrms_slope_out = wrms_slope(it)*scale
c
c.... get the duration of the data
      duration = (max_epoch - min_epoch)/365.25
      call jd_to_ymdhms( mean_epoch, date, sec_tag)
      mean_year = date(1) + date(2)/12.d0 + date(3)/365.25
c
c.... Now do the output
      num_len = stat(4,1,it)
      if( format.eq.'short' ) then
c
         write(unit,100) names(site_1),names(site_2), soln, num_len,
     .      mean_out,sig_mean_scaled, wrms_mean_out, nrms_mean(it),
     .      slope_out, sig_slope_scaled, wrms_slope_out, nrms_slope(it),
     .      duration, mean_year, line_st, line_en
 100     format(a8,"-",a8,1x,i2,1x,i4,1x,
     .      f14.4,1x,f6.4,1x,f6.1,1x,f5.2,3x,
     .      f8.2,1x,f8.2,1x,f6.1,1x,f5.2,1x,
     .      f6.2,1x,f7.2,1x,i7,"-",i7)
c
*                  ! do the long output
      else
c
*        Reconvert to meters for the long output
         if( it.eq.1 ) then
             scale = 1.0d0
         else
             scale = 1.0d-0
         end if
         write(unit,200) mean_out*scale, sig_mean_scaled*scale, num_len,
     .      wrms_mean_out, nrms_mean(it)
 200     format(/"Wmean ",f14.4," m +- ",f6.4," from ",i4," data.",
     .      " WRMS ",f6.1," mm, NRMS ",f5.2)
         write(unit,220) slope_out, sig_slope_scaled, wrms_slope_out,
     .      nrms_slope(it), duration, mean_year
 220     format("Slope ",f8.2," +- ",f7.2," mm/yr, WRMS ",f6.1,
     .      " mm, NRMS ",f5.2,", dur ",f5.2," <> ",f7.2," yr"/)
c
      end if
c
      return
      end
