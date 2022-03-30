      program ENFIT

      implicit none 
c
c     This program will read the output file from a SOLVK solution
c     and produce a ENU      summary.  This program is derived from
c     ensum and differs in that extra parameters can be estimated.
c     These parameters
c
c Include files
c -------------
*                                  ! the BLSUM control common block
      include 'enfit.h'
c
c Variables
c ---------
c  ema storage of baselines lengths from solution file
c data -- storage to data from file.  When this array is used it is
c     passed to subroutines where it is split up as
c        data(1) = site_1
c        data(2) = enu type (1=>N,2=>E,3=>U)
c        data(3) = solution number
c        data(5:6)   = epoch (julian date)
c        data(7:8)   = length (m)
c        data(9:10)  = sigma (m)
c
      integer*4 data(4,max_ent)
      real*8 vdata(4,max_ent)
c
      common /values/  vdata, data
 
c
c file names
c ----------
c more_files    -- logical to indicate that we are still reading files
c status_old  -- variable which gives the status old to file opening
c status_unk  -- varaible which gives the status unknown to file opening
c
c
      character status_old*4, status_unk*8
 
c
*  i    - loop counter
*  ierr - iostat error on read/write
	
      integer*4 i, j, ierr
 
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
c Data statements for status
c --------------------------
      data status_old / 'old '/, status_unk / 'unknown' /
c
c.... Get the runstring and decide if we are in batch mode
*     Check for sort option
      len_sort = rcpar(1, sort_string)
      call casefold( sort_string )
      if( len_sort.le.0 ) then
          call proper_runstring('enfit.hlp', 'enfit', 1)
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

* MOD TAH 9812129: See if separation option passed
      len_run(1) = rcpar(2, run_string(1))
      call check_num(run_string(1), ierr )
      if( ierr.eq.0 ) then
          read(run_string(1),*) match_sep
          j = 2
      else
          match_sep = -1.d0
          j = 1
      endif

* MOD TAH 030316: See if option for realistic sigmas passed
      len_run(1) = rcpar(j+1,run_string(1))
      if( run_string(1)(1:3).eq.'-RS' .or.
     .    run_string(1)(1:3).eq.'-rs' ) then
          j = j + 1
          realsigma = .true.
      else
          realsigma = .false.
      endif 

* MOD TAH 990504: See if fit command file passed
      fit_cmd_file = ' '
      len_run(1) = rcpar(j+1, run_string(1))
      if( run_string(1)(1:2).eq.'-f' .or.
     .    run_string(1)(1:2).eq.'-F'  ) then
          len_run(1) = rcpar(j+2,fit_cmd_file)
          j = j + 2
      end if


*                              ! read rest of runstring
      do i = 1,max_runstring
         len_run(i) = rcpar(i+j, run_string(i))

c...     see if zero string
         if( len_run(i).eq.0 ) run_string(i) = '  '
c
      end do
c
c.... see if runstring was passed
      batch = .true.

c.... write headr message
      write(*,50)
  50  format(//" ENFIT: Program to fit ENU components from GLOBK",
     .   " solution",/,
     .         " -----------------------------------------------",
     .   "---------",/)
c 
*     Save names of files

      summary_file     = run_string(1)
      values_file      = run_string(2)

***   Open and read the fit file
      call read_fit_file
c
c.... if we are in batch mode assign the file names from the runstring
      if( batch ) then
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
      end if

* MOD TAH 981229: Clear all site coordinates
      do i = 1, max_site
         do j = 1, 3
            coords(j,i) = 0.d0
         end do
      end do
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

* MOD TAH 981229: Now see if we should map sites
      if( match_sep.gt.0 ) then
          call match_sites( vdata, data )
      endif 
 
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
      include 'enfit.h'
c
c
c passed variables
c ----------------
c data -- ema storage of the baselines from the solution file
c
      integer*4 data(4,max_ent)
      real*8 vdata(4,max_ent)
 
c
c
c soln_num -- the number of the solution file
c
      integer*4 soln_num
 
c
c Local variables
c ---------------
c buffer -- the buffer used to read the solutiom file
      character*150 buffer
 
c
c new_soln -- a new solution indentifier from SOLVK
c new_time -- a new epoch indentifier
c bas_to   -- baseline indentifier -- not used
c l1,l2,l3 -- identifier of enu coordinates  - not used
c meter    -- Meter indentifier 
c
      character new_soln*10, new_time*18, meter*3
c     character l1*1,l2*1,l3*1, bas_to*4  -- not used
c
c clk_date -- epoch of the experiment
c
      integer*4 clk_date(5)
 
c
*   indx        - Pointer for use in multiread
*   err         - Error flag returned from multiread
*   iel         - index position of string
*   iclk        - Index position for clock time
*   ibas        - index position for baseline length
*   in          - length of string
*   ferr        = File reading error

      integer*4 indx, err, iel, iclk, ibas, in, ferr

*   lenu        - Indicates that we found NEU

      logical lenu
 
*   cdummy      - Dummy character string for multiread
 
      character*1 cdummy
 
*   new_gsoln*5     - Name GLOBK --not used
*   new_gtime*24    - solution commenced with line 
 
      character new_gtime*24, message*256
c     character new_gsoln*5 -- not used

c.... Define the identifier strings
c
      data new_soln / ' SOLUTION ' /, 
     .     new_time / 'EXPERIMENT DATE :'    /,
     .     meter / '(M)' /
c     data l1 / 'E' /, l2 / 'N' / , l3 / 'U' /,  -- not used
c    .     bas_to   / ' TO '     /  -- not used

      data  new_gtime / 'SOLUTION REFERS TO     :' /
c     data new_gsoln / 'GLOBK' /   -- not used

c.... Firstly read the header from solution
      read(100,'(a)',iostat=ferr) header(soln_num)
c
* MOD TAH 881031 Check to see if this is VALUES_FILE
      if( header(soln_num)(1:11).eq.'VALUES_FILE' ) then
          call read_values_file( soln_num, vdata, data )
*                      ! Do not need to do the rest
          return
      end if
*     See if name says that it is a values files
      if( solution_file(soln_num)(1:2).eq.'va' .or. 
     .    solution_file(soln_num)(1:4).eq.'VAL.' ) then
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
      ferr = 0
c
      num_soln = 1
      max_soln = 1
 
*                              ! read until error
      do while ( ferr.eq.0 )
c
         read(100,'(a)',iostat=ferr) buffer
         call casefold(buffer)
c
c....    continue processing if there was no error
         if( ferr.eq.0 ) then
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
                  write(message,'(a,a,4i8)') 'Solution number error in '
     .                 , buffer,' iel indx num_soln err '
     .                 , iel,indx,num_soln,err 
                  call report_stat('WARNING','BLSUM','blsum', ' '
     .                            ,message,0)
c                 above replaces below -- rwk 970920
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
c....       Now see if we have a local coordinate 
            lenu=.false.
            ibas=index(buffer,meter)

            if( ibas.gt.0 ) then

               call find_enu( buffer, ibas, lenu )
               
               if( lenu ) call save_ema(data(1,num_ent), 
     .                                  vdata(1,num_ent))

* MOD TAH 981229: See if we have coordinate already
               if ( coords(site(2),site(1)).eq.0.d0 .and.
     .              vdata(3,num_ent).lt.100.d0 ) then
                    coords(site(2),site(1)) = vdata(2,num_ent)
               end if
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
  300 format(" There are ",i10," baseline entries between ",i5," sites")
c
      return
      end
c
CTITLE FIND_ENU

      subroutine find_enu(buffer,indx,ok)

      implicit none 
c
c     subroutine to get the ENU site coordinates from the buffer
c
c     input:
c     buffer - char*150, data buffer read from the SOLVK result file
c     output to common:
c     site(2)     - int*2, site code and type of the coordinate
c                   (1->N, 2->E, 3->U)
c     base_dat(2) - real*8, coordinate and its rms error
c     num_ent     - number of entries
c
c     Programmer:  TGQ, 1989-03-09, CfA
c
c  include files
c  -------------
c
*                                   !the control common block
      include 'enfit.h'
c  variables defined
c  -----------------
      real*8 values(3)

      integer*4 i1,i2,i3,indx, ll, iel,  ierr 

      character buffer*(*),new_buffer*150

      character l1*1,l2*1,l3*1, cdummy*1

      logical ok,t,f

      data l1,l2,l3/'N','E','U'/,t,f/.true.,.false./
c  initiate the pointers
      i1=0
      i2=0
      i3=0
      ll=0
c to check if it is a enu site coordinate line
c and determine which coordinate is
      ll=indx-15
      
c     write(1,101)buffer(ll:ll)
c101  format('coordinate (',a1,')')   !debug
      if (buffer(ll:ll).eq.l1 .or.
     .    buffer(ll+1:ll+1).eq.l1 ) then
*                                     !coordinate type
         if(  buffer(ll+1:ll+1).eq.l1 ) ll = ll+1
         site(2)=1
         ok=t
      else if(buffer(ll:ll).eq.l2 .or.
     .    buffer(ll+1:ll+1).eq.l2 ) then
          if(  buffer(ll+1:ll+1).eq.l2 ) ll = ll+1
         site(2)=2
         ok=t
      else if(buffer(ll:ll).eq.l3 .or.
     .    buffer(ll+1:ll+1).eq.l3 ) then
          if(  buffer(ll+1:ll+1).eq.l2 ) ll = ll+1
         site(2)=3
         ok=t
      else
         ok=f
         return
      endif
c
c     an enu data record is caught
c.... If it is the first we have ever found
c
*                               ! just set iel=-1
      if (num_site.eq.0) then
             iel = -1
c     If we have already known some, see if it's a new comer
*                                ! see if in list
      else
          new_buffer=buffer(ll-9:)
c         write(1,'(A)')new_buffer
          call strip_name(new_buffer,names,num_site,iel,0)
      end if
c
c.... if it is the first site, file it
*                                  ! site name unknown
      if( iel.eq.-1 ) then
          new_buffer=buffer(ll-9:)
          call add_name(new_buffer)
          call strip_name(new_buffer,names,num_site,iel,0)
      end if
c
      site(1) = iel

* MOD TAH 000123: See if too many sites
      if( num_site.gt.max_site ) then
          call report_stat('FATAL','ENFIT','read_soln',
     .         ' ','Too many sites: Max allowed',max_site)
      endif

c
c....          Now get the site coordinate data
      indx = indx-9
c     write(1,'(A)')new_buffer(indx:)
      call multiread(new_buffer,indx,'R8',ierr,values,
     .                         cdummy,3)
      base_data(1) = values(1)
      base_data(2) = values(3)
      if( base_data(2).eq.0 ) base_data(2) = 0.0001
c
c.....   number of data in ema
      num_ent = num_ent + 1
c
c.....         see if maximun number of entries exceeded
      if( num_ent.ge.max_ent ) then
         write( *  ,200) num_ent
  200    format(/" Maximum number of baselines ",i8,
     .               " has been exceeded")
*                     ! Force exit on next read
         ierr = -99
      end if
c
c....       Thats all
      return
      end

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
      include 'enfit.h'
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
      subroutine save_ema(site_dat, real_dat )

      implicit none 
c
c     routine to save the baseline values in ema
c
c Include files
c -------------
c
*                                   ! the control common block
      include 'enfit.h'
c
c Variables
c ---------
c site_dat -- integer values for sites and solution number
c epoch     -- real*8 epoch of measurement
c length    -- real*8 length value
c sigma     -- real*8 sigma
c
      integer*4 site_dat(4)
      real*8 real_dat(4)
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
c
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
      include 'enfit.h'

* i,j,k,l,m  - Loop counters
* ierr       - IOSTAT error
* inf, inh   - lengths of file and header

      integer*4 i,j,k,m, inf, inh, ierr, trimlen, iter

* converged  -- Site true when time constant estimates have converged or reached
*     max iterations
* selected -- Set true is site matches selections
      logical converged, selected
c
c EMA variables
c -------------
c data -- the stored baseline data
c
      integer*4 data(4,max_ent)
      real*8 vdata(4,max_ent)
      real*8 residual ! Residual value
      real*8 sig_scale  ! Replacement for Nrms_slope when realistic 
                        ! sigams are computed
      real*8 taufin     ! Time constant for time series
c
c
c.... start looping over all baselines.  Write out the headers to the
c     summary file
c
      write(200,100)
  100 format(/"*ENFIT: GLOBK NEU component summary: includes the",
     .   " following solutions")
c
      do i = 1, num_files
         call trailing(inf, solution_file(i))
         call trailing(inh,header(i))
         write(200,120) solution_file(i)(1:inf), header(i)(:inh)
  120    format("*File ",a,": ",a)
      end do
      write(200,130) fit_cmd_file(1:max(1,trimlen(fit_cmd_file)))
 130  format('* Command File: ',a)
c
c.... Write out descriptive header
      call write_sum_head(200)

*     Open and write headers to fit_val_file
      if( trimlen(fit_val_file).gt.0 ) then
          call write_out_head(202)
      endif
 
c.... start the line counter for tplot file extraction
      line_st = 1
      line_en = 0
c
c.... loop over all baseline combinations
      if( max_soln.eq.0  ) max_soln = 1
      do i = 1, num_site
         call check_ss( names(i), selected)
         if ( selected ) then
            write( *,9000) names(i)
 9000       format(" Processing data from site ",a)
c
            do j = 1, 3        
c
c....          loop over the solution
               do k = 1, max_soln
c
*                 Set the initial estimates of tau
                  do m = 1, num_exp
                     tau_est(m) = tau_exp(m)
                     exp_est(m) = 0.d0
                  end do
                  do m = 1, num_log
                     dtl_est(m) = dtl_log(m)
                     log_est(m) = 0.d0
                  end do
c
c....             Now get all the entries
                  converged = .false. 
                  iter = 0
                  do while ( .not.converged )
                     call clear_norm

                     do m = 1,num_ent
c
c....                   see if baseline and solution matches
                        if( data(1,m).eq.i .and. data(2,m).eq.j .and.
*                                                   ! we found a value
     .                      data(3,m).eq.k ) then
c
c....                      accumulate statistics
                           call fit_acc(vdata(1,m),vdata(2,m),
     .                                  vdata(3,m))
c
                        end if
                     end do

c....                complete the computation of the statistics
                     call compute_fit(iter, converged)
                  end do
c
c....             see if any data ! yes there is data, so continue
                  if( ierr.ge.0 .and. stat(4,1).ge.min_out ) then

*                    Clear the statistics
                     call clear_stat
c
c....                write the header    ! now that we know there is data
                     call write_header(201,i,j,k)
                     line_en = line_en + 3
c
c....                Now accumlate residual statistics
                     num_ts = 0
                     do m = 1,num_ent
c
c....                   see if baseline and solution matches
                        if( data(1,m).eq.i .and. data(2,m).eq.j .and.
*                                                  ! we found a value
     .                     data(3,m).eq.k ) then
c
c....                      write record to file
                           call write_rec(201,vdata(1,m),vdata(2,m),
     .                                        vdata(3,m))
                           line_en = line_en + 1
c
c....                      accumulate statistics
                           call accum_stat(vdata(1,m),vdata(2,m),
     .                                     vdata(3,m), residual)
*                          Save the time series
                           if( vdata(3,m).gt.0 ) then
                              num_ts = num_ts + 1
                              times(num_ts) = vdata(1,m)
                              res(num_ts) = residual
                              Rerr(num_ts) = vdata(3,m)
                           end if
c
                       end if
                     end do
c
                     call complete_stat

* MOD TAH 030316: Compute realistic sigma
                     if( realsigma ) then 
                         call real_stats(times,res,Rerr,num_ts,
     .                                   sig_scale, taufin)
                         nrms_mean = sig_scale 
                     end if
c
c....                Now write out the summary
                     if( stat(4,1).ge.min_out ) then
                         call write_summary(200,i,j,k,'short ')
                     end if
                     call write_summary(201,i,j,k,'long')
c
c....                save line numbers
                     line_st = line_en + 5
                     line_en = line_st - 1

*                    See if we are generating output values
                     if( trimlen(fit_val_file).gt.0 ) then
                        call gen_out(202, i,j,k)
                     end if
*                        ! producing summary if there was data
                  end if
*                        ! looping over solutions   -- k
               end do
*                        ! looping over second site -- j
            end do
         end if       ! Site selected
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
      include 'enfit.h'

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
      return
      end
c
c........................................................................
      subroutine clear_norm

      implicit none 
c
c     This routine will clear the statistics accumulation and
c     slope computation  arrays
c
      include 'enfit.h'

* i,j  -- loop counters
      integer*4 i,j, np
c
      do i = 1,num_parm
         do j = 1, num_parm
            norm_eq(i,j) = 0.d0
         end do
      end do
*
*     Assign apriori weights
      norm_eq(1,1) = 1.d-6
      norm_eq(2,2) = 1.d-6
      np = 2
      do i = 1, num_exp
         np = np + 1
         norm_eq(np,np) = 1.d0/con_exp(i)**2
         if( est_tau(i) ) then
            np = np + 1
            norm_eq(np,np) = 1.d0/tau_sig(i)**2
         end if
      end do
      do i = 1, num_log
         np = np + 1
         norm_eq(np,np) = 1.d0/con_log(i)**2
         if( est_dtl(i) ) then
            np = np + 1
            norm_eq(np,np) = 1.d0/dtl_sig(i)**2
         end if
      end do
      do i = 1, num_per
         do j = 1,2
            np = np+1
            norm_eq(np,np) = 1.d0/con_per(i)**2
         end do
      end do

c
      do i = 1,num_parm
         bvec(i)   = 0.d0
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
      include 'enfit.h'
c
c Variables
c ---------
c unit -- output unit number
* in   -- Length of string
c
      integer*4 unit, i, j, k, in

* enu(3) - Denote NEU

      character*1 enu(3)

      data enu / 'N', 'E', 'U' /
 
c
c.... The header record to be used by tplot is the solution
c     description, baseline + solution #, blank record
c
c.... delete trailing blanks from header
      call trailing(in,header(1))
      write(unit,'(a)') header(1)(:in)
c
      write(unit,'(a," to ",a," Solution ",i2)') names(i),enu(j),k
c
      write(unit,'(1x)')
c
      return
      end
c
c........................................................................
c
      subroutine fit_acc(epoch,length,sigma)
c
c     Routine to accumulate the data needed to compute the slope
c
      include 'enfit.h'
c
c Variables
c ---------
c epoch -- the epoch of this measurement
c length -- the length of the baseline
c sigma  -- the sigma of the baseline
c
      real*8 epoch, length, sigma
      integer*4 i, j
c
c
c Local variables
c ---------------
c wgh  -- the length measurement weight
c
      real*8 wgh, resnl
 
c
c.... If this is the first value save the epoch and length
c
      if ( norm_eq(1,1).eq.1.d-6 ) then
         ref_epoch = ep_exp(1)
         if( num_log.gt.0 ) ref_epoch = ep_log(1)
         ref_length = length
         stat(4,1) = 0
      end if
c
c.... see if this this value is before or after all others
      stat(4,1) = stat(4,1) + 1
      if( epoch.lt. min_epoch ) min_epoch = epoch
      if( epoch.gt. max_epoch ) max_epoch = epoch
c
c.... Now accumulate normal equations
      wgh = 1.d0/sigma**2
      call form_apart(epoch)
      
      resnl = (length-ref_length)
      call nonlin(resnl, epoch)
   
      do i = 1, num_parm
         bvec(i) = bvec(i) +    resnl*wgh*apart(i)

         do j = 1, num_parm
            norm_eq(i,j) = norm_eq(i,j) + apart(i)*wgh*apart(j)
         end do
      end do

c
      return
      end
c
c........................................................................
c
      subroutine compute_fit(iter, converged)

      implicit none 
c
c     Finish solving the least-squares problem

      include 'enfit.h'
c
c Local variables
c ---------------
c scale(max_parm) -- Scaling vector

      real*8 scale(max_parm)

* iter -- Iteratiion count
* ipivot(max_parm)  -- Pivot elements

      integer*4 iter, ipivot(max_parm), np, i

      logical converged, loc_conv
 
c
c.... see if we have enough data to do a mean value
 
      call invert_vis(norm_eq, bvec, scale, ipivot, num_parm,
     .                max_parm, 1 )

*     Now check the estimated tau paramters
      converged = .true. 
      iter = iter + 1
      np = 2
*     Update and check exponent time constant
      do i = 1, num_exp
         np = np + 1
         loc_conv = .true.
         if( est_tau(i) ) then
             np = np + 1

*            Now check convergence and that the time offset can not 
*            go negative.
*            Check the size of the adjustment for convergence
             loc_conv = .false.
             if( (abs(bvec(np)).lt.1d-2 .and.
     .           abs(bvec(np-1)).lt.1.d-6) .or. iter.gt.20 ) then
                 loc_conv = .true.
             end if
*            Check for negative time offsets first
             if ((tau_est(i) + bvec(np)).lt.0 ) then
*                Set the offset correction back to zero
                 bvec(np) = 0.0d0
                 bvec(np-1) = 0.0d0
             endif

*            Now update the values
             tau_est(i) = tau_est(i) + bvec(np)
             exp_est(i) = exp_est(i) + bvec(np-1)
             bvec(np) = 0.d0
             bvec(np-1) = 0.d0
         end if
         if( loc_conv .and. est_tau(i)) then
             write(*,220) np, iter, 
     .                   tau_est(i), sqrt(norm_eq(np,np))
220          format(' Nonlinear Exp time # ',i3,' Converged ',i3,
     .              ' Iterations; Estimate ',F6.2,' +- ',F6.2,' days')
         endif
         converged = converged .and. loc_conv
      end do
*     Update and check convergence of log time offsets
      do i = 1, num_log
         np = np + 1
         loc_conv = .true.
         if( est_dtl(i) ) then
             np = np + 1

*            Now check convergence and that the time offset can not 
*            go negative.
*            Check the size of the adjustment for convergence
             loc_conv = .false.
             if( (abs(bvec(np)).lt.1d-2 .and.
     .           abs(bvec(np-1)).lt.1.d-6) .or. iter.ge.20 ) then
                 loc_conv = .true.
             end if
*            Check for negative time offsets first
             if ((dtl_est(i) + bvec(np)).lt.0 ) then
*                Set the offset correction back to zero
                 bvec(np) = 0.0d0
                 bvec(np-1) = 0.0d0
             endif

*            Now update the values
             dtl_est(i) = dtl_est(i) + bvec(np)
             log_est(i) = log_est(i) + bvec(np-1)
             bvec(np) = 0.d0
             bvec(np-1) = 0.d0
         end if
         if( loc_conv .and. est_dtl(i) ) then
             write(*,320) np, iter, 
     .                   dtl_est(i), sqrt(norm_eq(np,np))
320          format(' Nonlinear log time # ',i3,' Converged ',i3,
     .              ' Iterations; Estimate ',F6.2,' +- ',F6.2,' days')
         endif
         converged = converged .and. loc_conv
      end do
c
      return
      end
c
c........................................................................
c
      subroutine accum_stat(epoch, length, sigma, residual)

      implicit none 
c
c     Routine to accumulate the statistics from a set of baseline
c     results.  The statistics are accumulated for both the residuals
c     about the mean and about the best fit straight line
c
      include 'enfit.h'
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
c dlen     -- Model calculation based on fit to data
c
      real*8 wgh, residual, dlen

* i   -- Loop counter
      integer*4 i
 
c
c.... Compute the residuals and sum into the statistics
      wgh = 1.d0/sigma**2
c
c.... do statistics about the mean
      call form_apart(epoch)
      dlen = 0.d0
      do i = 1, num_parm
         dlen = dlen + apart(i)*bvec(i)
      end do

      residual = (length-ref_length)
      call nonlin(residual, epoch)

      residual = residual - dlen
c
      stat(1,1) = stat(1,1) + residual*wgh
      stat(2,1) = stat(2,1) + 1.d0*wgh
      stat(3,1) = stat(3,1) + residual**2*wgh
      stat(4,1) = stat(4,1) + 1.d0
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
      include 'enfit.h'
c
c
c.... Firstly do the statistics about the weighted mean, if there were
c     any data at all.
c
      if( stat(4,1).gt.num_parm ) then
c
         nrms_mean = sqrt( stat(3,1)/(stat(4,1)-num_parm) )
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
      if( stat(4,2).gt.num_parm ) then
c
         nrms_slope = sqrt( stat(3,2)/(stat(4,2)-num_parm) )
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
      include 'enfit.h'
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
* dlen       -- Correction to length 
* cov        -- Variance of length estimate
* est        -- Estimate of length

      real*8 rel_epoch, sec_tag,  residual, dlen, cov, est
      integer*4 i, j
c
c.... set up for output, convert epoch to yr,mon,day,hr,min
      rel_epoch = epoch
      call jd_to_ymdhms(rel_epoch, date, sec_tag)
      call form_apart( epoch )
      dlen = 0
      call nonlin(dlen, epoch)
      dlen = -dlen
      cov  = 0
      do i = 1, num_parm
          dlen = dlen + bvec(i)*apart(i)
          do j = 1, num_parm
             cov = cov + apart(i)*norm_eq(i,j)*apart(j)
          end do
      end do
      est = ref_length + dlen

c
c.... do statistics about the mean
      residual = length - est
c
c.... write the record
      write(unit,100) date,length,sigma, residual, sigma, 
     .                est, sqrt(cov)
  100 format(i5,4i3,f15.5,1x,f10.5,1x,2f10.5,1x,f15.5,1x,f10.5)

      return
      end
c
c........................................................................
c
      subroutine write_summary(unit,site_1,enu_type, soln, format)

      implicit none 
c
c     Routine to write the summary information.  This can be written
c     in two forms.. the short where all results are put on a single
c     line or the long form which takes multiple lines
c
      include 'enfit.h'
c
c variables
c ---------
c unit -- the output unit number
c site_1 -- first site in baseline
c enu_type -- type of coordinate being used.
c soln -- solution number
c format -- determines if long or short form
c
      integer*4 unit, site_1, enu_type, soln
 
c
      character*(*) format
 
c
c Local variables
c ---------------
c mean_out  -- the output weighted mean (m)
c
      real*8 mean_out

* num_len -- number of values in mean length estimate (used so we can output
*            with an I4 format
 
      integer*4 num_len, i, np 
 
* enu(3) -- Denote NEU

      character*1 enu(3)

* out_format  -- String containing output format
      character*256 out_format

      data enu / 'N', 'E', 'U' /
c
c.... Convert the units form meters to millimeters where necessary
c     and scale the sigmas of the various quanities
c
*                                  ! scale wgh_mean sigma
      if( nrms_mean.gt.1 ) then
         sig_mul = nrms_mean
      else
         sig_mul = 1.d0
      end if

*     Compute and scale parameter sigmas
      do i = 1, num_parm
         sig_est(i) = sqrt(norm_eq(i,i))*sig_mul
      end do

      do i = 1, num_parm
         finsol(i) = bvec(i)
      end do
      
      np = 2
      do i = 1, num_exp
         np = np + 1 
         finsol(np) = finsol(np) + exp_est(i)
         if( est_tau(i) ) then
             np = np + 1
             finsol(np) = finsol(np) + tau_est(i)
         end if
      end do
      do i = 1, num_log
         np = np + 1 
         finsol(np) = finsol(np) + log_est(i)
         if( est_dtl(i) ) then
             np = np + 1
             finsol(np) = finsol(np) + dtl_est(i)
         end if
      end do

c
      mean_out  = ref_length + bvec(1)
c
c.... Now do the output
      num_len = stat(4,1)
      if( format.eq.'short' ) then
c
         write(out_format,100) num_parm-1
 100     format('(a8,1x,a1,1x,i2,1x,i4,1x, f13.4,1x,f6.4,',
     .          i2.2,'(1x,F7.2,1x,F6.2),F6.2,F6.2,1x,i7,1x,i7)')

         write(unit,out_format) names(site_1),enu(enu_type), 
     .      soln, num_len,  mean_out, sig_est(1), 
     .      (finsol(i), sig_est(i), i = 2, num_parm), 
     .      wrms_mean*1000.d0, nrms_mean,  line_st, line_en
c
*                  ! do the long output
      else
c
         write(unit,200) mean_out, sig_est(1), num_len,
     .      wrms_mean*1000, nrms_mean
 200     format(/,"Offset ",f13.4," m +- ",f6.4," from ",i4," data.",
     .      " WRMS ",f5.1," mm, NRMS ",f5.2)
         write(unit,220) (finsol(i), sig_est(i), i = 2, num_parm)
 220     format('Other Parameters ',F7.2,' +- ',F6.2,' mm/yr',
     .          20(1x,F7.2,' +- ',f6.2,' mm; '))
         write(unit,'(1x)')
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
 
      include 'enfit.h'
 
*         soln_num      - Number of current file
*         data(4, max_ent)   - Place where data is actually
*                       - saved.
 
      integer*4 soln_num, data(4, max_ent)
      real*8 vdata(4,max_ent) 
 
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
      num_soln = 1
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
 
         if( new_buffer(2:2).eq.'N' ) site(2) = 1
         if( new_buffer(2:2).eq.'E' ) site(2) = 2
         if( new_buffer(2:2).eq.'U' ) site(2) = 3
 
*****    Get out the solution number
         indx = 10
         call read_line(new_buffer,indx,'I4',ierr,num_soln,cdummy)
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
                 call save_ema(data(1,num_ent), vdata(1,num_ent))
* MOD TAH 981229: See if we have coordinate already
                 if ( coords(site(2),site(1)).eq.0.d0 .and.
     .                vdata(3,num_ent).lt.0.5d0 ) then
                      coords(site(2),site(1)) = vdata(2,num_ent)
                 end if
 
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
  300 format(' There are ',i10,' baselines entries between ',
     .                     i5,' sites, Solns ',i3 )
*
      return
      end
 
CTITLE BLSORT
 
      subroutine blsort(vdata, data)
 
      implicit none 
 
*     Routine to sort the baseline data into time order.
 
      include 'enfit.h'
 
*         data(4,max_ent)   - Data to be time sorted. Time
*                                   - is in elemnet 4-7.
*   i,j,k           - Loop counters
*   temp(num_field) - Temporary storage for switch
 
      integer*4 data(4,max_ent), i,j,k, temp(num_field)
      real*8 vdata(4,max_ent), vtemp(4) 
 
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
 
                  do k = 1, 4  
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
 
CTITLE MATCH_SITES

      subroutine match_sites( vdata, data )

      implicit none 

*     Routine to map from site to another based on separation being
*     less than specified max.

c Include files
c -------------
c
*                                    ! the control common block
      include 'enfit.h'
c
c
c passed variables
c ----------------
c data -- ema storage of the baselines from the solution file
c
      integer*4 data(4,max_ent)
      real*8 vdata(4,max_ent)

* Local variables 

* i,j,k, l, m -- Loop counters
* num_ent_org -- Original number of entries.  We use this stop
*     the new mapped values from being used to map other sites.

      integer*4 i,j,k,l,m, li, num_ent_org

* stats(2)   -- Accumulation statistics
* sep        -- Site separation
* wgh        -- Weight used to get mean

      real*8 stats(2), sep, wgh

***** Now loop over sites see which are close enough
      do i = 1, num_site
        write(*,320) names(i), (coords(j,i),j=1,3)
 320    format(1x,a8,3F20.4)
      end do 

****  Save the original number of entries
      num_ent_org = num_ent

      write(*,120) num_site, match_sep
 120  format('Checking ',i4,' sites for separations < ',F8.2,' m')
      do i = 1, num_site-1
         do j = i+1, num_site

*           Get the separation:
            sep = 0
            do k = 1,3
               sep = sep + (coords(k,i)-coords(k,j))**2
            end do

            if( sep.le.match_sep**2 ) then
                write(*,150) names(j), names(i)
 150            format('Mapping site ',a8,' to ',a8)

*               Now see if we more than one repeat and get
*               mean difference of the positions.
                ni = 0
*               First find all the site i entries
                do m = 1, num_ent_org
                   if ( data(1,m).eq.i .and. ni.le.max_match ) then
                      k = data(2,m)
                      li = 0
                      if( k.eq.1 ) then
                         ni = ni + 1
                         li = ni
                      else
                         do l = 1, ni
                            if( vdata(1,m).eq.si_ep(l)) li = l
                         end do
                      endif
                      if( li.gt.0 ) then
                         si_ep(li) = vdata(1,m)
                         si_cd(k,li) = vdata(2,m)
                         si_sg(k,li) = vdata(3,m) 
                      end if
                   end if
                end do

*               Now scan again to find the matching site j entries
                do m = 1,3
                   nj(m) = 0
                   na(m) = 0
                enddo
                do m = 1, num_ent_org
                   if( data(1,m).eq.j  ) then

*                     Scan the site i epochs to see if match
                      k = data(2,m)
                      do l = 1, ni
                         if( vdata(1,m).eq.si_ep(l) ) then
*                           Epoch match.  Save the difference and
*                           sigma.  Second part of 'if' is to handle
*                           case of duplicate epochs where we would
*                           not have stored the second and greater
*                           E and U values. 
                            if( nj(k).lt.max_match .and.
     .                          abs(si_cd(k,l)-vdata(2,m)).lt.
     .                              match_sep) then 
                                nj(k) = nj(k) + 1
                                sj_cd(k,nj(k)) = si_cd(k,l)-vdata(2,m)
                                sj_sg(k,nj(k)) = sqrt(si_sg(k,l)**2+
     .                                           vdata(3,m)**2)
                            end if
                         end if
                      end do
*                     Now save all the entries for site j
                      if( na(k).lt.max_match ) then
                          na(k) = na(k) + 1
                          sa_ep(k,na(k)) = vdata(1,m)
                          sa_cd(k,na(k)) = vdata(2,m)
                          sa_sg(k,na(k)) = vdata(3,m)
                      end if
                   end if   
                end do 

                write(*,180) ni, nj, na
 180            format('Found ',i4,' and ',3i4,' match, ',3i4,' total')

*               Now get the mean difference
                do m = 1, 3
                   if( nj(m).ge.1 .and. na(m).gt.1 ) then
                      do l = 1, 2
                         stats(l) = 0.d0
                      end do
                      do l = 1, nj(m)
                         wgh = 1.d0/sj_sg(m,l)**2

                         stats(1) = stats(1) + sj_cd(m,l)*wgh
                         stats(2) = stats(2) + wgh
                      end do

                      if( stats(2).gt.0 ) then
                         stats(1) = stats(1)/stats(2)
                      else
                         stats(1) = 0.d0
                      end if

*                     Now create entries for mapped sites
                      if( stats(1).ne.0 ) then
                          write(*,220) na(m),m, stats(1)
 220                      format('Adding ',i4,' Entries for coord ',
     .                           i2,' Adjusting ',F10.4, 'm')
  
                          do l = 1, na(m)
                             num_ent = num_ent + 1
                             data(1,num_ent) = i
                             data(2,num_ent) = m
                             data(3,num_ent) = 1
                             vdata(1,num_ent) = sa_ep(m,l)
                             vdata(2,num_ent) = sa_cd(m,l) + stats(1)
                             vdata(3,num_ent) = sa_sg(m,l)
                          end do
                      end if
                   end if
                end do

*            End of checks and loops for match_sep, and looping
*            over stations.
             end if
          end do
      end do

*     Now tell user number of entries
      write(*,300) num_ent, num_ent_org
 300  format('After site mapping there are ',i8,' entries (',
     .        I8,' Originally)')

      return
      end

CTITLE READ_FIT_FILE

      subroutine read_fit_file

      implicit none 

*     This routine will read the enfit command file to find out
*     what needs to be estimated and output.

      include 'enfit.h'

* LOCAL VARIABLES
* ierr, jerr -- IOSTAT errors
* indx, jndx -- Pointers to positions in string
* ival(2)    -- Interger values read from line
* date(5)    -- Calender date
* i          -- Loop counter

      integer*4 ierr, jerr, indx, jndx, ival(2), trimlen, date(5),
     .          i

* sectag     -- Seconds tag for conversion to JD
* rval(2)    -- R*8 values read from line

      real*8 sectag, rval(2)
* cdum       -- Dummy for readline
* word       -- Word extracted from command line

      character*64 cdum, word 
* line       -- Command line read from file.

      character*256 line


*     Initialize before we start
      num_exp = 0
      num_log = 0
      num_per = 0
      num_out = 0
      num_parm = 2
      num_ss = 0
      fit_val_file = ' '
      ep_exp(1) =  2451179.50d0
      ep_log(1) =  2451179.50d0

*     See if command file used
      if( trimlen(fit_cmd_file).eq.0 ) RETURN

****  Open the command file
      open(103, file=fit_cmd_file, iostat=ierr,status='old')
      call report_error('IOSTAT',ierr,'open',fit_cmd_file,1,'enfit')


*     Now read the file
      do while ( ierr.eq.0 )
         read(103,'(a)',iostat=ierr) line
         if( ierr.eq.0 .and. trimlen(line).gt.0 .and. 
     .       line(1:1).eq.' ' ) then

*           Get the first word
            indx = 1
            call getword(line, word, indx )
            call casefold(word)

*           Check out commands:
* EXP  -- Exponential Fit command 
            if( word(1:3).eq.'EXP' ) then
*              read the date and decay time
               do i = 1, 5
                  call read_line(line,indx,'I4',jerr,date(i),cdum)
               end do
               call read_line(line,indx,'R8',jerr,rval(1),cdum)
               call read_line(line,indx,'R8',jerr,rval(2),cdum)
               call report_error('IOSTAT',ierr,'decod',line,0,
     .              'EXP command')
               if( jerr.eq.0 .and. num_exp.lt.max_exp ) then
                   num_exp = num_exp + 1
                   sectag = 0.d0
                   call ymdhms_to_jd(date, sectag, ep_exp(num_exp))
                   tau_exp(num_exp) = rval(1)
                   con_exp(num_exp) = rval(2)

*                  See if tau is to be estimated
                   call read_line(line,indx,'R8',jerr,rval(1),cdum)
                   if( jerr.eq.0 .and. rval(1).gt.0.d0 ) then
                       tau_sig(num_exp) = rval(1)
                       est_tau(num_exp) = .true.
                   else
                       est_tau(num_exp) = .false.
                   end if

               else
                   if( num_exp.eq.max_exp ) then
                      write(*,120) max_exp
 120                  format('**WARNING** Max number of exponentials ',
     .                       i4,' exceeded.  New entries ignored')
                   end if
               end if
* LOG  -- logonential Fit command 
            elseif( word(1:3).eq.'LOG' ) then
*              read the date and time offset
               do i = 1, 5
                  call read_line(line,indx,'I4',jerr,date(i),cdum)
               end do
               call read_line(line,indx,'R8',jerr,rval(1),cdum)
               call read_line(line,indx,'R8',jerr,rval(2),cdum)
               call report_error('IOSTAT',ierr,'decod',line,0,
     .              'log command')
               if( jerr.eq.0 .and. num_log.lt.max_log ) then
                   num_log = num_log + 1
                   sectag = 0.d0
                   call ymdhms_to_jd(date, sectag, ep_log(num_log))
                   dtl_log(num_log) = rval(1)
                   con_log(num_log) = rval(2)

*                  See if dtl is to be estimated
                   call read_line(line,indx,'R8',jerr,rval(1),cdum)
                   if( jerr.eq.0 .and. rval(1).gt.0.d0 ) then
                       dtl_sig(num_log) = rval(1)
                       est_dtl(num_log) = .true.
                   else
                       est_dtl(num_log) = .false.
                   end if

               else
                   if( num_log.eq.max_log ) then
                      write(*,130) max_log
 130                  format('**WARNING** Max number of logarithms ',
     .                       i4,' exceeded.  New entries ignored')
                   end if
               end if
* PER  -- Periodic function fit
            else if( word(1:3).eq.'PER' ) then

*              get the Period (days)
               call read_line(line,indx,'R8',jerr,rval(1),cdum)
               call read_line(line,indx,'R8',jerr,rval(2),cdum)
               call report_error('IOSTAT',ierr,'decod',line,0,
     .              'PER command')
               if( jerr.eq.0 .and. num_per.lt.max_per ) then
                   num_per = num_per + 1
                   per_per(num_per) = rval(1)
                   con_per(num_per) = rval(2)
               else
                   if( num_per.eq.max_per ) then
                      write(*,140) max_per
 140                  format('**WARNING** Max number of periodic',
     .                       ' terms', i4,' exceeded. ',
     .                       ' New entries ignored')
                   end if
               end if

* OUT  -- Output days for function fit.
            else if( word(1:3).eq.'OUT' ) then

*              Get the file name
               call getword(line,fit_val_file, indx)
               call wild_card(fit_val_file,summary_file)

*              Read accross line getting the days to output values
               jerr = 0
               do while ( jerr.eq.0 )
                  call read_line(line,indx,'CH',jerr,rval,word)
                  if( index(word,'x').gt.0 .and. jerr.eq.0 ) then
*                     Then form is 10x100 which is 10 values spaced
*                     by 100 days
                      call sub_char(word,'x',' ')
                      jndx = 1
                      call read_line(word,jndx,'I4',jerr,ival(1),cdum)
                      call read_line(word,jndx,'R8',jerr,rval(1),cdum)
                      do i = 1, min(ival(1),max_out-num_out)
                         out_exp(num_out+i) = i*rval(1)
                      end do
                      num_out = num_out + min(ival(1),max_out-num_out)

                  else if ( jerr.eq.0 ) then
                      if( num_out.lt.max_out ) then
                         num_out = num_out + 1
                         jndx = 0
                         call read_line(word,jndx,'R8',jerr,rval(1),
     .                                  cdum)
                         out_exp(num_out) = rval(1)
                      else
                         write(*,160) max_out
 160                     format('**WARNING** Max number of output',
     .                          ' values', i4,' exceeded. ',
     .                          ' New entries ignored')
                      endif
                  end if
               end do 
            elseif( word(1:3).eq.'SEL' ) then
               call extract_ss( line,indx, sel_site, num_ss, 
     .                          max_site_sel)

* NoSuch Command -- No command found
            else
               call report_error('NoSuchCommand',-1,'decod',line,0,
     .              'enfit/no such command')
            end if
         end if
      end do

*     Summarize what we found out
      num_parm = 2 + num_exp +num_log + 2*num_per
*     See how many additional parameters for tau estimation
      do i = 1, num_exp
         if( est_tau(i) ) num_parm = num_parm + 1
      end do
      do i = 1, num_log
         if( est_dtl(i) ) num_parm = num_parm + 1
      end do

      write(*,210) num_exp, num_log, num_per, num_out
 210  format('Parameters: ',i4,' Expontenitals; ',
     .       i4, ' Logarithm;', I4,' Periodic; ',
     .       i4,' Output values')
      if( trimlen(fit_val_file).gt.0 ) then
          write(*,220) fit_val_file(1:trimlen(fit_val_file))
 220      format('Values output to ',a)
      endif
      if( num_ss.gt.0 ) then

         write(*,240) num_ss, (sel_site(i),i=1,num_ss)
 240     format('Site Selections: ',i4,' wild cards',/,
     .       100(10(a8,:,1x),/))
      endif

****  Thats all
      return
      end

CTITLE WRITE_SUM_HEAD

      subroutine write_sum_head(unit)

      implicit none 

*     Routine to write the the header lines for the output

      include 'enfit.h'

* PASSED VARIABLES
* unit  -- Unit number
      integer*4 unit

* LOCAL VARIABLES

* i   -- Loop counter
* trimlen  -- Length of string
* pos      -- Position in string
* date(5)  -- date for output

      integer*4 i, trimlen, pos, date(5)

* sectag   -- Seconds tag
      real*8 sectag

* line     -- Line to be written
      character*256 line 


***** Make the first line
      line = 'Site    ENU #         Offset at'
      write(line(41:),'(a)') '      Rate'
      pos = 58 
      do i = 1, num_exp

         write(line(pos:),120) i
 120     format('Exponential',i2)
         if( est_tau(i) ) then
            pos = pos + 15
            write(line(pos:),125) i
 125        format('Tau Est',i2)
         end if
         pos = pos + 15
      end do 

      do i = 1, num_log

         write(line(pos:),130) i
 130     format('Logarithm',i2)
         if( est_dtl(i) ) then
            pos = pos + 15
            write(line(pos:),135) i
 135        format('dtl Est',i2)
         end if
         pos = pos + 15
      end do 

      do i = 1, num_per

         write(line(pos:),140) i, per_per(i)
 140     format('Periodic ',i2,1x,F7.1,' d')
         pos = pos + 30
      end do

      write(line(pos:),160) 
 160  format('wrms  nrms  line #')
      write(unit,'(a)') line(1:trimlen(line))

*     Make the second line
      if( num_exp.gt.0 ) then
          call jd_to_ymdhms(ep_exp(1), date, sectag)
      elseif ( num_log.gt.0 ) then
          call jd_to_ymdhms(ep_log(1), date, sectag)
      else
          call jd_to_ymdhms(ref_epoch, date, sectag)
      endif
 
      line = ' '
      write(line(21:),210) date
 210  format(i4,4(1x,i2.2),' +-')
      write(line(45:), 220)
 220  format('mm/yr    +-')
      pos = 58 
      do i = 1, num_exp
         write(line(pos:),240) tau_exp(i)
 240     format('Tau ',F6.1,' days')
         if( est_tau(i) ) then
             pos = pos + 15
             write(line(pos:),245)
 245         format(' Days ')
         end if
         pos = pos + 15
      end do
      do i = 1, num_log
         write(line(pos:),250) dtl_log(i)
 250     format('dTl ',F6.1,' days')
         if( est_dtl(i) ) then
             pos = pos + 15
             write(line(pos:),255)
 255         format(' Days ')
         end if
         pos = pos + 15
      end do

      do i = 1, num_per
         write(line(pos:),260) 
 260     format(' Cos     +-      Sin     +-')
         pos = pos + 30
      end do

      write(unit,'(a)') line(1:trimlen(line))

*     Make the third and final line
      line = ' '
      pos = 58 
      do i = 1, num_exp
         call jd_to_ymdhms(ep_exp(i),date, sectag)
         write(line(pos:),320) date
 320     format(i4,4(1x,i2.2))
         if( est_tau(i) ) pos = pos + 15
         pos = pos + 15
      end do

      do i = 1, num_log 
         call jd_to_ymdhms(ep_log(i),date, sectag)
         write(line(pos:),330) date
 330     format(i4,4(1x,i2.2))
         if( est_dtl(i) ) pos = pos + 15
         pos = pos + 15
      end do

      do i = 1, num_per

         write(line(pos:),340)
 340     format('  mm     mm       mm     mm')
         pos = pos + 30
      end do

      write(line(pos:),360) 
 360  format(' mm')

      if( trimlen(line).gt.0 ) then
         write(unit,'(a)') line(1:trimlen(line))
      end if

***** Thats all
      return
      end

CTITLE FORM_APART

      subroutine form_apart( epoch )

      implicit none 

*     Routine to form partial derivatives

      include 'enfit.h'
      include '../includes/const_param.h'

* PASSED VARIABLE
* epoch  -- JD of this measurement

      real*8 epoch

* LOCAL VARIABLES
* i, np    -- Loop counter and parameter number
      integer*4 i, np

* dt       -- time difference (days)
* jd_2000  -- JD of 2000.0

      real*8 dt, jd_2000

      data jd_2000 / 2451544.5d0 /

*     Start forming the partials
      dt = epoch - ref_epoch
      apart(1) = 1.d0
      apart(2) = dt / (365.25d0*1000.d0 )
      np = 2

*     Now do exponentials
      do i = 1, num_exp
         np = np + 1
         if( epoch.ge.ep_exp(i) ) then
* MOD TAH 030612: Check exp to 1-exp.  Time constant partial changed as well
             apart(np) = (1-exp(-(epoch-ep_exp(i))/tau_est(i)))/1000.d0
         else
             apart(np) = 0.d0
         end if
         if( est_tau(i) ) then
             np = np + 1
             if( epoch.ge.ep_exp(i) ) then
                 apart(np) = -exp_est(i)*
     .                       exp(-(epoch-ep_exp(i))/tau_est(i))**2*
     .                       (epoch-ep_exp(i))/tau_est(i)**2/1000.d0
             else
                apart(np) = 0.d0
             end if
          end if
      end do

*     Now do logiarithms
      do i = 1, num_log
         np = np + 1
         if( epoch.ge.ep_log(i) ) then
* MOD TAH 030612: Changed log(dt-tau) to log(1+dt/tau).  
C            apart(np) = log((epoch-ep_log(i)+dtl_est(i))/365.25d0)/1.d3
             apart(np) = log(1+(epoch-ep_log(i))/dtl_est(i))/1.d3
         else
             apart(np) = 0.d0
         end if
         if( est_dtl(i) ) then
             np = np + 1
             if( epoch.ge.ep_log(i) ) then
C               apart(np) = log_est(i)/
C    .             ((epoch-ep_log(i)+dtl_est(i))/365.25d0)/365.25d3
                apart(np) = -log_est(i)/
     .                       (1+(epoch-ep_log(i))/dtl_est(i))*
     .                       (epoch-ep_log(i))/dtl_est(i)**2/1.d3
             else
                apart(np) = 0.d0
             end if
          end if
      end do

*     Now do periodic
      do i = 1, num_per
         np = np + 1
         apart(np) = cos((epoch-jd_2000)/per_per(i)*2*pi)/1000.d0
         np = np + 1
         apart(np) = sin((epoch-jd_2000)/per_per(i)*2*pi)/1000.d0
      end do

****  OK All done
      return 
      end

CTITLE GEN_OUT

      subroutine gen_out(unit, is,jc,ks)

      implicit none 

*     Routine to generate results at specified times

      include 'enfit.h'

* PASSED VARIABLES
* is -- Station number
* jc -- Component number
* ks -- Solution number
* unit -- Unit number

      integer*4 is,jc, ks, unit

* LOCAL VARIABLES
* i,j,k  -- Loop counters
* date(5) -- Date

      integer*4 i,j,k, date(5) 

* dlen, var  -- DLen from expotential epoch and variance
* adiff(max_parm) -- Partials for difference
* jd_out     -- JD of output values
* sectag     -- Seconds tag

      real*8 dlen, var, adiff(max_parm), jd_out, sectag


* enu(3) - Denote NEU

      character*1 enu(3)

      data enu / 'N', 'E', 'U' /
 

****  Loop over the times needed
      do i = 1, num_out

*        Compute the output JD
         jd_out = ep_exp(1) + out_exp(i)
         if( num_log.gt.0 ) jd_out = ep_log(1) + out_exp(i)
*        Compute partials and save
         call form_apart(jd_out)
         do j = 1, num_parm
            adiff(j) = apart(j)
         end do

*        Now form negative partials at reference epoch
         if( num_exp.gt.0 ) then
            call form_apart(ep_exp(1))
         else if ( num_log.gt.0 ) then
            call form_apart(ep_log(1))
         else
            call form_apart(ref_epoch)
         end if

         do j = 1, num_parm
             adiff(j) = adiff(j) - apart(j)
         end do

*        Now get dlen and variance
         dlen = 0.d0

*        Get nonlinear term (reverse sign since this is removed from data)
         call nonlin(dlen, jd_out)
         dlen = -dlen

         var = 0.d0
         do j = 1, num_parm
             dlen = dlen + adiff(j)*bvec(j)
             do k = 1, num_parm
                var = var + adiff(j)*norm_eq(j,k)*adiff(k)
             end do
         end do

*        Get the date
         call jd_to_ymdhms(jd_out, date, sectag)

*        Now write the line:
         write(unit,120) names(is), enu(jc), ks, date, out_exp(i),
     .                   dlen*1000.d0, sqrt(var)*1000.d0*sig_mul
 120     format(a8,1x,a1,1x,i2,1x,i4,4(1x,i2),F9.1,f8.1,1x,f9.1)
      end do

****  Thats all 
      return
      end

CTITLE WRITE_OUT_HEAD

      subroutine write_out_head( unit ) 

      implicit none 

*     Routine to write out output value file header

      include 'enfit.h'

* PASSED VARIABLES
* unit  -- Unit number

      integer*4 unit

* LOCAL VARIABLES
* ierr  -- IOSTAT errors
* date(5) -- Calender date
* trimlen -- Length of string

      integer*4 ierr, date(5), trimlen, i 
* sectag  -- Seconds tag

      real*8 sectag

****  Open file:  If error, set name to blank
      open(unit,file=fit_val_file, iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'open',fit_val_file,0,'enfit')
      if( ierr.eq.0 ) then
         write(unit,120) values_file(1:trimlen(values_file))
 120     format('ENFIT: Spaced value estimates from ',a)
         call jd_to_ymdhms(ep_exp(1), date, sectag)
         write(unit,140) date
 140     format('Values given as difference from ',i4,4(1x,i2))
         write(unit,160)
 160     format(' SITE   NEU #        Date',8x,'Days',4x,
     .          'Change (mm)  +- ') 
         write(*,180) num_out, (out_exp(i),i=1,num_out)
 180     format('Output of ',i4,' values at ',50f10.1)
      else
         fit_val_file = ' '
      end if

****  Thats all
      return
      end 

CTITLE NONLIN

      subroutine nonlin( resnl, epoch )

      implicit none 

*     Routine to apply the non-linear part of the model to the residual

      include 'enfit.h'

* PASSED VARIABLES
* res  -- residual to be corrected for non-linear exponential term
* epoch -- JD of measurement

      real*8 resnl, epoch

* LOCAL VARIABLES
* i -- Loop counter
      integer*4 i

* dt -- Delta time over tau
      real*8 dt

      do i = 1, num_exp
         dt = (epoch-ep_exp(i))/tau_est(i)
         if( dt.ge.0 ) then
             resnl = resnl - exp_est(i)*(1-exp(-dt))/1000.d0
         endif
      end do
      do i = 1, num_log
         dt = (epoch-ep_log(i))/dtl_est(i)
         if( dt.ge.0 ) then
             resnl = resnl - log_est(i)*log(1+dt)/1000.d0
         endif
      end do

****  Thats all
      return
      end

CTILTE CHECK_SS

      subroutine check_ss( name, selected)

      implicit none 

*     Routine to check if site names has been selected
      include 'enfit.h'

* PASSED VARIABLES
      logical selected  ! Set true if site selected
      character*(*) name

* LOCAL VARIABALES
      integer*4 i
      character*8 newname


      selected = .true.
      if ( num_ss.eq.0 ) RETURN

*     OK Loop over selection criteria and see if matched
      selected = .false.
      do i = 1, num_ss
         newname = sel_site(i)
         call wild_card(newname, name)
         if ( name.eq.newname ) selected = .true.
      end do

      return
      end

CTITLE EXTRCAT_SS

      subroutine extract_ss( line,indx, sel_site, num_ss,
     .                       max_site_sel)

      implicit none 
 
*     Routine to extract pattern selections

* PASSED VARIABLES
      integer*4 indx,     ! Current position in buffer
     .          num_ss,    ! Current number of patterns
     .          max_site_sel  ! Maximum number allowed

      character*(*) sel_site(max_site_sel),  ! Patterns
     .          line      ! Line read from command file    


* LOCAL VARIABLES
      integer*4 trimlen, jndx
      character*8 next_word
      logical done

      done = .false.
      do while ( .not.done )
         jndx = indx
         call GetWord( line, next_word, indx)
         if( trimlen(next_word).gt.0 ) then
             if( next_word(1:1).eq.'!' .or.
     .           next_word(1:1).eq.'#'     ) then
                 done = .true.
             else
                 num_ss = num_ss + 1
                 if( num_ss.gt.max_site_sel ) then
                     write(*,120) max_site_sel, line(jndx:trimlen(line))
 120                 format('Too many sel_site: Max ',i4,': ',
     .                       a,' ignored')
                     done = .true.
                 else
                     sel_site(num_ss) = next_word
                 end if
             end if
          else
             done = .true.
          end if
       end do

       return
       end

             

      



