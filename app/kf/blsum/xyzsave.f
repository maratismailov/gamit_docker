      program XYZSAVE

      implicit none
c
*     This prpgram will read a series of org files and save the 
*     XYZ coordinates in separate files for each site.  
* MOD TAH 090205: New version of this program to output XYZ files in the
*     Reason/SIO format.  The valules file is replaced with a directory
*     for the files to sent.
c
c Include files
c -------------
*                                  ! the BLSUM control common block
      include 'ensum.h'
c
c Variables
c ---------
c  ema storage of baselines lengths from solution file
c data -- storage to data from file.  When this array is used it is
c     passed to subroutines where it is split up as
c        data(1) = site_1
c        data(2) = enu type (1=>N,2=>E,3=>U)
c        data(3) = solution number
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

      character*256 xyzdir   ! Directory for XYZ files
 
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
          call proper_runstring('xyzsave.hlp', 'xyzsave', 1)
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
  50  format(//" XYZSAVE: Program to extract XYZ components from GLOBK",
     .   " solution",/,
     .         " ---------------------------------------------------",
     .   "---------",/)
c
c
c.... if we are in batch mode assign the file names from the runstring
      if( batch ) then
         xyzdir     = run_string(1)
c
c....    Now get the solution files from the runstring
         num_files = 0
         count = 2
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
      end if

c
c.... Now read in the data from the solution files
      num_ent = 0
      do i = 1, num_files
         write(*,'(" Reading stations from ",a)') solution_file(i)
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
      write(*,*) 'Generating solution summary'
      call out_xyz(vdata, data, xyzdir)
c
c.... Thats all
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
      include 'ensum.h'
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
      character new_soln*10, new_time*18, meter*16
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
 
      integer*4 indx, err, ierr, iel, iclk, ibas, in

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
     .     meter / 'COORDINATE  (M)' /

      data  new_gtime / 'SOLUTION REFERS TO     :' /
c     data new_gsoln / 'GLOBK' /   -- not used

c.... Firstly read the header from solution
      read(100,'(a)',iostat=ierr) header(soln_num)
c 
      call trailing(in,header(soln_num))
      write(*,100) header(soln_num)(:in)
 100  format(" Solution: ",a)
c
c.... Now start reading the solution file extracting the baseline
c     lengths
      ierr = 0
c
      num_soln = 1
      max_soln = 1
 
*                              ! read until error
      do while ( ierr.eq.0 )
c
         read(100,'(a)',iostat=ierr) buffer
         call casefold(buffer)
c
c....    continue processing if there was no error
         if( ierr.eq.0 ) then
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
                  call report_stat('WARNING','XYSUM','xysum', ' '
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
            if( index(buffer,'TRANSLTN').gt.0 .or.
     .          index(buffer,'TRANSLATION').gt.0 ) ibas = 0

            if( ibas.gt.0 ) then

               call find_xyz( buffer, ibas, lenu )
               
               if( lenu ) call save_ema(data(1,num_ent), 
     .                                  vdata(1,num_ent))

*                         ! baseline found
            end if
*                         ! no error on file read
         end if
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
CTITLE FIND_XYZ

      subroutine find_xyz(buffer,indx,ok)

      implicit none
c
c     subroutine to get the XYZ site coordinates from the buffer
c
c     input:
c     buffer - char*150, data buffer read from the SOLVK result file
c     output to common:
c     site(2)     - int*2, site code and type of the coordinate
c                   (1->X, 2->Y, 3->Z)
c     base_dat(2) - real*8, coordinate and its rms error
c     num_ent     - number of entries
c
c     Programmer:  TGQ, 1989-03-09, CfA
c
c  include files
c  -------------
c
*                                   !the control common block
      include 'ensum.h'
c  variables defined
c  -----------------
      real*8 values(3)

      integer*4 i1,i2,i3,indx, ll, iel,  ierr 

      character buffer*(*),new_buffer*150

      character l1*1,l2*1,l3*1, cdummy*1

      logical ok,t,f

      data l1,l2,l3/'X','Y','Z'/,t,f/.true.,.false./
c  initiate the pointers
      i1=0
      i2=0
      i3=0
      ll=0
c to check if it is a enu site coordinate line
c and determine which coordinate is
      indx = index(buffer,'(M)')
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
          if(  buffer(ll+1:ll+1).eq.l3 ) ll = ll+1
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
*         Make all names _GPS
          new_buffer(5:8) = '_GPS'
c         write(1,'(A)')new_buffer
          call strip_name(new_buffer,names,num_site,iel,0)
      end if
c
c.... if it is the first site, file it
*                                  ! site name unknown
      if( iel.eq.-1 ) then
          new_buffer=buffer(ll-9:)

*         Make all names _GPS
          new_buffer(5:8) = '_GPS'

          call add_name(new_buffer)

* MOD TAH 000123: See if too many sites
          if( num_site.ge.max_site ) then
              call report_stat('FATAL','ENSUM','read_soln',
     .             ' ','Too many sites: Max allowed',max_site)
          endif
          call strip_name(new_buffer,names,num_site,iel,0)
      end if
c
      site(1) = iel
c
c....          Now get the site coordinate data
      indx = indx-9
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
      include 'ensum.h'
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
      if( num_site.gt.max_site ) then 
          call report_stat('FATAL','XYZSAVE','add_name','',
     .         'Too many sites',max_site)
      endif
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
      include 'ensum.h'
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
      subroutine out_xyz(vdata, data, xyzdir)

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
      include 'ensum.h'

* i,j,k,l,m  - Loop counters
* ierr       - IOSTAT error
* inf, inh   - lengths of file and header

      integer*4 i,m, ierr
c
c EMA variables
c -------------
c data -- the stored baseline data
c
      integer*4 data(4,max_ent)
      real*8 vdata(4,max_ent)

      integer*4 trimlen, lendir, date(5), doy
      real*8 decyrs, sectag

      character*(*) xyzdir   ! Directory to output results

      character*256 outfile  ! File name
      character*4 lcname     ! Lower case site name
c
c
c.... start looping over all baselines.  Write out the headers to the
c     summary file
c
*     Make sure directory has / at end of name
      lendir = trimlen(xyzdir)
      if( xyzdir(lendir:lendir).ne.'/' ) then
          lendir = lendir + 1
          xyzdir(lendir:lendir) = '/'
      endif
      
*     Loop over sites
      do i = 1, num_site
         write( *,120) i, names(i)
 120     format(" Processing data from site ",I5,2x,a)

****     Create output file
         lcname = names(i)(1:4)
         call caseunfold(lcname)
         outfile = xyzdir(1:lendir) // lcname // 'Raw.xyz'
         open(201,file=outfile,iostat=ierr, status='unknown')
*        If error make fatal
         call report_error('IOSTAT',ierr,'open',outfile,1,'xyzsave')

c....    Now get all the entries
         do m = 1,num_ent
c
c....       see if baseline and solution matches
            if( data(1,m).eq.i .and. data(2,m).eq.1 ) then
*               Found an X coordinates, next entries will be Y and Z
*               Get the dates we need
                call jd_to_decyrs( vdata(1,m), decyrs)
                call jd_to_ymdhms( vdata(1,m), date, sectag)
                call ymd_to_doy( date, doy)
*               Output line
                write(201,220) decyrs, date(1), doy,
     .                vdata(2,m), vdata(2,m+1), vdata(2,m+2),              
     .                vdata(3,m), vdata(3,m+1), vdata(3,m+2)
 220            format(f9.4,1x,I4,1x,I3,1x,3(F14.4,1x),3(F8.4,1x))

C1995.2726 1995 100 -2341332.8851 -3539049.5113 4745791.3847  0.0015  0.0020  0.0021
c
            end if
         end do
*        Close the output file
         close(201)
      end do
c
c.... Thats all
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
 
 
CTITLE BLSORT
 
      subroutine blsort(vdata, data)

      implicit none
 
 
*     Routine to sort the baseline data into time order.
 
      include 'ensum.h'
 
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
 

 


                  


 


