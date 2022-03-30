      program merge_igs_clk

      implicit none

*     Program to merge multiple igs clock files into a single
*     clock file 

      include 'merge_igs_clk.h'

* LOCAL VARIABLES

      integer*4 nt   ! Counter for netwerk number
     .,         i    ! Loop over iterations
     .,         ep   ! Epoch counter
     .,         num_it   ! Number of iterations (only 1 if
                     ! if 1 input


***** Start: Decode the runstring to get the input networks
* MOD TAH 200928: Option to pass sampling rate
      space_ep = 0.d0  ! Maybe set from runsring 
      call mic_runstring

***** Read the data files
      num_site = 0
      num_sat = 0
      num_ref = 0
      num_ep = 0
      start_ep = 0.d0
      call init_flags
      clk_gnss = ""  ! Default (none and set when data read)
*     Start reading clock files.
      do nt = 1, num_net
         call read_clk(nt)
      end do

****  Iterate the allign of the clocks
* MOD TAH 201006: Set just 1 iteration if num_net == 1
      num_it = 5
      if( num_net.eq.1 ) num_it = 1

      do i = 1, num_it
         write(*,120) i
 120     format(/,'ALIGNMENT PASS ',i2)


****     Now fit all the clocks that we have read.  This determines the
*        linear trends for each 
         call fit_clk(.true.)


*        For those clocks with good fits to the linear trends, use these 
*         to allign the different clock files.
         call allign_clk

*        Now check continuity of clocks and fix any offsets
         call offset_clk( i ) 

****     Now that we have the allignment, run through and average clocks
*        from over lapping networks
C        call gen_mesh
         call fit_clk(.false.)


*****    Based on list of reference sites, compute average residual
*        at each epoch and remove a constant (by epoch) from all
*        clocks
* MOD TAH 201006: Only reset_refs if num_net > 1
         if( num_net.gt.1 ) then
            call reset_refs
         endif

*****    Once we have iterated a few times, check the quality of the 
*        individual values to the fits and remove any outlier values
         if( i.gt.2 ) then
             call edit_clk(i)
         end if

*        TEST cleaning befow final runs
         if( i.gt.3 ) then
            write(*,'(a,1x,I2)') 'AVERAGING NETS Iter ',i
            do ep = 1,num_ep
               call average_net(ep)
            end do
         end if 
     

      end do

****  Average the common clocks and form the final results
      debug = .false.
* MOD TAH 201006: If only one-net copy clocks across.
      if( num_net.gt.1 ) then
         call gen_fin
      else
         call copy_fin
      endif

****  Now write out the clock file
      call wr_clk


****  Thats all
      end

CTITLE MIC_RUNSTRING

      subroutine mic_runstring

      implicit none

*     Routine to decode the Merge_igs_clock runstring

      include 'merge_igs_clk.h'

* LOCAL VARIABLES

      integer*4 rcpar, len_run, n,m
     .,         trimlen ! Length of string
      integer*4 jerr    ! ISTAT error

      logical done

      character*256 runstring
 
****  Loop over the runstring
      len_run = 1
      n = 0
      out_file = ' '
      debug = .false.
* MOD TAH 200723: Changed default to abs model.
      abs_mod = .true.
      num_usr_ref = 0

      do while (len_run.gt.0 )
         n = n + 1
         len_run = rcpar(n,runstring)

*        See what option has been passed
         if( runstring(1:2).eq.'-o' ) then
*            Output file
             n = n + 1
             len_run = rcpar(n,out_file)
         else if ( runstring(1:2).eq.'-i' ) then
*            input files to follow.  We do this as
*            a loop
             m = 0
             done = .false.
             do while ( .not.done )
                n = n + 1
                len_run = rcpar(n,runstring)
                if( len_run.gt.0 .and. runstring(1:1).ne.'-' ) then
                    m = m + 1
                    clk_files(m) = runstring
                    num_net = m
                    if( num_net.gt.max_net) then
                        write(*,'(a,1x,I3)') 
     .                     '** FATAL ** Too many nets Max ',max_net
                       stop 'MERGE_IGS_CLK: Too many Networks'
                    endif
                else
                    done = .true.
                    n = n - 1
                endif
             end do
         else if( runstring(1:2).eq.'-a' ) then
             abs_mod = .true.
         else if( runstring(1:2).eq.'-r' ) then
* MOD TAH 200623: Added supplying reference site list in runstring.
*        rather than using list in .clk files
*        Loop over the list
             m = 0
             done = .false.
             do while ( .not.done )
                n = n + 1
                len_run = rcpar(n,runstring)
                if( len_run.gt.0 .and. runstring(1:1).ne.'-' ) then
                    m = m + 1
                    ref_codes(m) = runstring
                    num_usr_ref = m
                else
                    done = .true.
                    n = n - 1
                endif
             end do

* MOD TAH 200928: Added -s <sampling> option to output clock rate
*       
         else if( runstring(1:2).eq.'-s' ) then
             n = n + 1
             len_run = rcpar(n,runstring)
             if( len_run.gt.0 .and. runstring(1:1).ne.'-' ) then
                read(runstring,*,iostat=jerr) space_ep
                if( jerr.ne.0 ) then
                   print *,'IOSTAT error ',jerr,' decoding space_ep'
                   space_ep = 0.0d0
                else
                   space_ep = space_ep/86400.d0  ! Convert to days
                   write(*,'(a,1x,F6.2," sec")')'Clock Sample Rate ',
     .                space_ep*86400
                endif
             end if
               
         end if

      end do

****  Tell user what is happening
      write(*,120)
 120  format('MERGE_IGS_CLK: Program parameters')
      if( trimlen(out_file).gt.0 ) then
          write(*,140) out_file(1:trimlen(out_file))
 140      format('Merged file Name: ',a)
      else
          call proper_runstring('merge_igs_clk.hlp','merge_igs_clk',1)
      end if

      do m = 1, num_net
         write(*,160) m, clk_files(m)(1:trimlen(clk_files(m)))
 160     format('Input file ',i3,': ',a)
      end do

***** Thats all
      return 
      end

CTITLE READ_CLK

      subroutine read_clk(nt)

      implicit none

*     Routine to read the clock files

      include 'merge_igs_clk.h'

* PASSED VARIABLES

      integer*4 nt   ! File/network number to read

* LOCAL VARIABLES

      integer*4 ierr, jerr  ! IOSTAT error
     .,         indx  ! Position in string
     .,         num   ! Number of estimates of clock (always 1)
     .,         ns    ! Site or satellite nu,ber
     .,         ep    ! Epoch number
     .,         date(5) ! YMDHM of clock estimate
     .,         prn    ! PRN number from file
     .,         tn(2)    ! Total number of clock values read
     .,         i,j      ! Loop counters
c      logical   OK     ! Set true is ref_site is in standard list

* Functions used
      integer*4 get_clk_ep   ! Returns epoch number
     .,         ptol    ! PRN number to entry in list

* MOD TAH 200627: GNSS offset function
      integer*4 conoff   ! Returns numeric offset based on GNSS type

      logical head_done ! Set true when header finished

      real*8 clk_vers  ! RINEX Version of clock file
     .,      sectag    ! Seconds tag of measurement
     .,      clk, sig  ! Clock value and sigma from file

      character*128 line ! Line read from file
      character*4 site   ! site code
      character*1 gnss   ! GNSS type


***** Start by opening file 
      open(100, file=clk_files(nt), status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open', clk_files(nt),
     .                   1,'READ_CLK')

*     OK: Loop over the file reading
      read(100,'(a)',iostat=ierr) line
*     Verify clock file and read OK
      indx = index(line,'C')
      if( indx.ne.21 .or. ierr.ne.0 ) then
         call report_error('IOSTAT/TYPE',ierr,'read',line,1,
     .                     'READ_CLK/First line')
      end if
      read(line,'(f9.2)',iostat=ierr) clk_vers

*     Read next line to get progr information
      read(100,'(a)',iostat=ierr) line
      autcln_ver = line(1:20)

*     Now start reading the file.  Base action on type of line
      head_done = .false.
      do while ( .not.head_done )
         read(100,'(a)',iostat=ierr) line

*        Check type of line
         indx = index(line,'# OF CLK REF')
         if( indx.eq.61 ) then
             call read_ref(100,nt,line)
         endif
         indx = index(line,'# OF SOLN STA')
         if( indx.eq.61 ) then
             call read_sta(100, nt,line)
         end if
         indx = index(line,'# OF SOLN SATS')
         if( indx.eq.61 ) then
             call read_sats(100,nt, line)
         end if
         indx = index(line,'SYS / PCVS APPLIED')
         if( indx.eq.61 ) then
             call read_pcvs(100,nt, line)
         end if

         indx = index(line,'END OF HEADER')
         if( indx.eq.61 ) then
             head_done = .true.
         end if
         if ( ierr.ne.0 ) then
            call report_error('IOSTAT',ierr,'read',clk_files(nt),
     .                        1,'READ_CLK/Unexpected end of header')
         end if
      end do         ! Finished reading header

*     Reference sites used in autcln:
*     algo albh nrc2 wtzr amc2 brus hob2 drao wsrt usno nya1 pie1 mate brew nlib
*     Remove any sites from reference list that are not in the std_ref list
*     Write final list of reference sites
* MOD TAH 200623: See if we are using user supplied list
      if( num_usr_ref.gt.0 ) then
         j = 0
         do i = 1, num_usr_ref
            indx = 0
            call get_cmd(ref_codes(i),site_codes, num_site, ns, indx)
            if( ns.lt.0 ) then
               write(*,110) i, ref_codes(i)
 110           format('User Reference site ',i3,1x,a,' not found')
            else
*              Add site into index list
               j = j + 1
               ref_sites(j) = ns
            endif
         enddo
         num_ref = j
         write(*,115) num_ref
 115     format('USING ',i2,' user reference clocks')
      endif

      write(*,120) num_ref, (site_codes(ref_sites(i)),i=1,num_ref)
 120  format('REFERENCE SITE LIST ',i3,20(1x,a4))


****  OK Now start reading the clock values themselves
      tn(1) = 0
      tn(2) = 0
      do while ( ierr.eq.0 )
         read(100,'(a)',iostat=ierr) line

*        See if site or satellite
         if( line(1:2).eq.'AR' .and. ierr.eq.0 ) then
*           Receiver clock
            read(line, 220,iostat=jerr) site, date, sectag, num, 
     .                                  clk, sig
 220        format(3x,a4,1x,i4,4i3.2,1x,f9.6,1x,i2,
     .             2x,2(1x,E19.12)) 
            if( jerr.eq. -1 .or. sig.le.0.d0) then
                sig = 10.d-9
                jerr = 0
            endif
* MOD TAH 200624: Make the sigmas more uniform to avoid weighting
            if( sig.le.1.d-9 ) sig = 1.d-9
*
            indx = 1
            call get_cmd(site,site_codes, num_site, ns, indx)
c Decimation code
c           if( abs(nint((date(5)*60.d0+sectag)/300.d0)*300.d0-
c     .          (date(5)*60.d0+sectag) ).lt.1.d0 ) then   
c              ep = get_clk_ep(date,sectag)
c           else
c              ep = 0
c           endif

            ep = get_clk_ep(date,sectag)

*           If indices are OK, save the value
            if( ns.gt.0 .and. ep.gt.0 .and. jerr.eq.0 ) then
                rcv_clk(ns,nt,ep) = clk
                rcv_sgc(ns,nt,ep) = sig
                rcv_flgc(ns,nt,ep) = 1
                tn(1) = tn(1) + 1
                num_ep = ep
            elseif( jerr.ne.0 ) then
* MOD TAH 200928: Report if only an error in read. Data can be skipped
*               when new -s option set.
                write(*,230) ns, ep, jerr, line(1:80)
 230            format('**WARNING** Site/Epoch/jerr',3i4,
     .                 ' read error ',a)
            end if
         elseif( line(1:2).eq.'AS' .and. ierr.eq.0 ) then
*           Satellite clock
* MOD TAH 200628: Add GNSS read and PRN update
            read(line, 320, iostat=jerr) gnss, prn, date, sectag, num, 
     .                                   clk, sig
 320        format(3x,a1,i2,1x,1x,i4,4i3.2,1x,f9.6,1x,i2,
     .             2x,2(1x,E19.12))
            prn = prn + conoff(gnss)  ! Change PRN for GNSS type
            if( jerr.eq. -1 .or. sig.le.0.d0) then
                sig = 10.d-9
                jerr = 0
            endif
* MOD TAH 200624: Make the sigmas more uniform to avoid weighting
            if( sig.le.1.d-9 ) sig = 1.d-9

            ns = ptol(prn) 
c Decimation code
c           if( abs(nint((date(5)*60.d0+sectag)/300.d0)*300.d0-
c     .          (date(5)*60.d0+sectag) ).lt.1.d0 ) then   
c              ep = get_clk_ep(date,sectag)
c           else
c              ep = 0
c           endif

            ep = get_clk_ep(date,sectag)

*           If indices are OK, save the value
            if( ns.gt.0 .and. ep.gt.0 .and. jerr.eq.0 ) then
                sat_clk(ns,nt,ep) = clk
                sat_sgc(ns,nt,ep) = sig
                sat_flgc(ns,nt,ep) = 1
                tn(2) = tn(2) + 1
            elseif( jerr.ne.0 ) then
* MOD TAH 209928: Only report due to error (-s option added)
                write(*,330) ns, ep, jerr, line(1:80)
 330            format('**WARNING** Sat/Epoch/jerr ',3i4,
     .                 ' Read Error ',a)
            end if
         end if
      end do
      

*     OK, all done reading file
      close(100)

      write(*,420) nt, tn
 420  format('For Network ',i3,' Total site values ',i6,' Sats ',i6)
C     write(*,440) start_ep, space_ep*1440.d0, ep
C440  format('Epoch Start ',F12.5,' Step ',F6.1,' mins, Last Ep ',i5)

*     Save the number of epochs found
      if ( ep.gt. num_ep ) num_ep = ep
      if( num_ep.gt.max_ep) then
         write(*,'(a,1x,I6)') 
     .       '** FATAL ** Too epochs Max ',max_ep
         stop 'MERGE_IGS_CLK: Too many epochs'
      endif


*     Tell user what we read
c      call report_clk_read(nt)

*     Thats all
      return
      end

CTITLE INIT_FLAGS

      subroutine init_flags

      implicit none

*     Routine to set all clock flags to 0, when values read they
*     are initially set to 1.

      include 'merge_igs_clk.h'

* LOCAL VARIABLES
      integer*4 i,j,k

*     Loop and set everything to zero
      do i = 1, max_ep
         do j = 1, max_site
            rcv_flag(j,i) = 0
            do k = 1, max_net
               rcv_flgc(j,k,i) = 0
            end do
         end do
      end do

      do i = 1, max_ep
         do j = 1, max_sat
            sat_flag(j,i) = 0
            do k = 1, max_net
               sat_flgc(j,k,i) = 0
            end do
         end do
      end do

      do i = 1, max_site
         do j = 1, max_net
            do k = 1,4
               rcv_fit(k,i,j) = 0.d0
            end do
         end do
      end do

      do i = 1, max_sat
         do j = 1, max_net
            do k = 1,4
               sat_fit(k,i,j) = 0.d0
            end do
         end do
      end do
            

****  Thats all
      return 
      end

CTITLE READ_REF

      subroutine read_ref(unit, nt, line)

      implicit none

*     Routine to read the reference sites from the clock file
*     and to add to an existing list


      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 nt  ! Network number
     .,         unit ! Unit being read from 
      character*(*) line  ! Line read from file

* LOCAL VARIABLES
      integer*4 i,j,k   ! Loop counters
     .,         nr      ! Number of reference sites from file
     .,         indx    ! Position in string
     .,         ierr, jerr    ! IOSTAT error
     .,         ns      ! Site number from list of sites
     .,         trimlen ! Length of string

      character*4 site  ! Site code
      character*20 name ! Long name of site
      character*6 lnum



****  Start: get the number of reference sites
      jerr = 0
      lnum = line
      read(lnum,'(i6)',iostat=jerr) nr
      call report_error('IOSTAT',jerr,'read',lnum,1,'READ_REF')

*     Now loop reading the reference site information
      i = 0
      do while ( i.lt. nr )
         read(unit,'(a)',iostat=ierr) line
         call report_error('IOSTAT',ierr,'read','Reference clocks',
     .                     1, 'READ_REF')
         if( index(line,'ANALYSIS CLK REF').gt.0 ) then
            i = i + 1
            site = line(1:4)
            name = line(6:26)
*           see if we have this site name yet
            indx = 1
            call get_cmd(site,site_codes, num_site, ns, indx)
            if( ns.le.0 ) then
*               We do not know about this site, add to list
                num_site = num_site + 1
                if( num_site.gt.max_site) then
                    write(*,'(a,1x,I5)') 
     .                  '** FATAL ** Too many sites Max ',max_site
                    stop 'MERGE_IGS_CLK: Too many sites'
                endif

                call casefold(site)
                site_codes(num_site) = site
                ns = num_site
                site_names(num_site) = name
            endif

*           OK, we know site name; see if added to reference
*           list
            j = 0
            do k = 1, num_ref
               if( ref_sites(k).eq.ns) j = k
            enddo
            if( j.eq.0 ) then
                num_ref = num_ref + 1
                ref_sites(num_ref) = ns
            endif
         end if
      end do

*     Report to user
      write(*,420) clk_files(nt)(1:trimlen(clk_files(nt))),
     .             nr, num_ref, num_site
 420  format('REF SITES from ',a,' # ref new ',i4,' Total ',i4,
     .       ' # Sites ',i4)

****  Thats all
      return 
      end


CTITLE READ_STA

      subroutine read_sta(unit, nt, line)

      implicit none

*     Routine to read the site names and coords from the clock file
*     and to add to an existing list


      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 nt  ! Network number
     .,         unit ! Unit being read from 
      character*(*) line  ! Line read from file

* LOCAL VARIABLES
      integer*4 i,k   ! Loop counters
     .,         nr      ! Number of reference sites from file
     .,         indx    ! Position in string
     .,         ierr, jerr    ! IOSTAT error
     .,         ns      ! Site number from list of sites
     .,         trimlen ! Length of string

      real*8 xyz(3)     ! coordinates of site

      character*4 site  ! Site code
      character*20 name ! Long name of site
      character*6 lnum



****  Start: get the number of reference sites
      jerr = 0
      lnum = line
      read(line,'(i6,4x,a)',iostat=jerr) nr, clk_frame_ver
      call report_error('IOSTAT',jerr,'read',lnum,0,'READ_STA')

*     Now loop reading the reference site information
      i = 0
      do while ( i.lt. nr )
         read(unit,'(a)',iostat=ierr) line
         call report_error('IOSTAT',ierr,'read','Site names',
     .                     1, 'READ_STA')
         if( index(line,'SOLN STA NAME').gt.0 ) then
            i = i + 1
            read(line,120) site, name, xyz
 120        format(a4,1x,a20,f11.0,1x,f11.0,1x,f11.0)

*           see if we have this site name yet
            indx = 1
            call get_cmd(site,site_codes, num_site, ns, indx)
            if( ns.le.0 ) then
*               We do not know about this site, add to list
                num_site = num_site + 1
                ns = num_site
                call casefold(site)
                site_codes(num_site) = site
            endif
            site_names(ns) = name

*           Save the coordinates of the site
            do k = 1,3
               site_coords(k,ns) = xyz(k)
            end do

         end if
      end do

*     Report to user
      write(*,420) clk_files(nt)(1:trimlen(clk_files(nt))), 
     .             nr, num_site
 420  format('STA SITES from ',a,' # sta new ',i4,' Total ',i4)
 
****  Thats all
      return 
      end

CTITLE READ_SATS

      subroutine read_sats(unit, nt, line)

      implicit none

*     Routine to read the satellite numbers from the clock file
*     and to add to an existing list.
*     NOTE: This code only works for GPS satellites


      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 nt  ! Network number
     .,         unit ! Unit being read from 
      character*(*) line  ! Line read from file

* LOCAL VARIABLES
      integer*4 i,j,k   ! Loop counters
     .,         nr      ! Number of reference sites from file
     .,         ierr, jerr    ! IOSTAT error
     .,         trimlen ! Length of string
     .,         mxs     ! Number of entries on line
     .,         prn_loc(max_sat)  ! List of PRN's from this file
     .,         num_line ! Number of PRN lines to read

* MOD TAH 200628: Updated for GNSS and gnss prn
      integer*4 conoff   ! Function to return gnss offset
      character*1 gnss_loc(max_sat)  ! GNSS type

      character*6 lnum


****  Start: get the number of reference sites
      jerr = 0
      lnum = line
      read(lnum,'(i6)',iostat=jerr) nr
      call report_error('IOSTAT',jerr,'read',lnum,0,'READ_SATS')

*     Now loop reading the satellites (code assumes no comments in
*     lines)
*     We can get 15 entries on a line, see how many lines needed
      num_line = (nr-1)/15+1
      do i = 1, num_line
         mxs = 15
         if( i.eq.num_line ) mxs = nr-(i-1)*15

         read(unit,120,iostat=ierr) (gnss_loc(j+(i-1)*15), 
     .                                prn_loc(j+(i-1)*15),j=1,mxs) 
 120     format(15(a1,i2,1x))
         call report_error('IOSTAT',ierr,'read','PRN line',
     .                      1,'READ_SATS')

* MOD TAH 200723: See what system we have (change to M if changes)
         do j = 1,mxs
            if( clk_gnss.eq.' ' ) then
               clk_gnss = gnss_loc(j)
            elseif ( clk_gnss.ne.gnss_loc(j) ) then
               clk_gnss = 'M'
            endif
         enddo 
      
      end do

* MOD TAH 200628: Convert prn_loc to full gnss prn
      do i = 1, nr
         prn_loc(i) = prn_loc(i) + conoff(gnss_loc(i))
      end do

*     Now make sure we have these satellites
      do i = 1,nr
         k = 0
         do j = 1, num_sat
            if( prn_loc(i).eq.prn_list(j) ) k = j
         end do
         if( k.eq.0 ) then   ! We don't have this satellite
              num_sat = num_sat + 1
              if( num_sat.gt.max_sat) then
                  write(*,'(a,1x,I5)') 
     .                '** FATAL ** Too many sats Max ',max_sat
                  stop 'MERGE_IGS_CLK: Too many satellites'
              endif

              prn_list(num_sat) = prn_loc(i)
         end if
      end do

*     Report to user
      write(*,420) clk_files(nt)(1:trimlen(clk_files(nt))), 
     .             nr, num_sat
 420  format('PRNS      from ',a,' # prn new ',i4,' Total ',i4)

****  Thats all
      return 
      end

CTITLE READ_PCVS

      subroutine read_pcvs(unit, nt, line)

      implicit none

*     Routine to read the gamit version and PCVS model used


      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 nt  ! Network number
     .,         unit ! Unit being read from 
      character*(*) line  ! Line read from file


****  Start: get the number of reference sites
      read(line,140) clk_gamitver, clk_phsmod
 140  format(8x,a5,7x,a40)
!140  format('G GAMIT ',F5.2,'       ',a40,'SYS / PCVS APPLIED')

****  Thats all
      return 
      end

CTITLE GET_CLK_EP

      integer*4 function get_clk_ep(date,sectag)

      implicit none

*     Function to return the epoch number of this clock value.  Epoch
*     set up is based on first two epochs of first file

      include 'merge_igs_clk.h'

* PASSSED VARIABLES
      integer*4 date(5)    ! Date
      real*8 sectag        ! Seconds tag of date

* LOCAL VARIABLES

      integer*4 ep   ! Epoch number
      real*8 mjd     ! MJD of date
     .,      epf     ! Floating epoch to make sure fits



****  Convert the date to MJD
      call ymdhms_to_mjd( date,sectag, mjd)

****  See if we have the start epoch already
      if( start_ep.eq.0 ) then
          start_ep = mjd
      endif

****  See if we are start time
      if( abs(mjd-start_ep).lt.1.d-5 ) then
          ep = 1
      elseif ( space_ep.eq.0 ) then
          space_ep = mjd - start_ep
      endif

****  OK Compute epoch
      if( space_ep.gt.0 ) then
           epf = (mjd-start_ep)/space_ep + 1
           if( abs(epf - nint(epf)).lt.0.001 ) then
               ep = nint(epf)
           else
               ep = 0
           end if
      end if

****  See if epoch looks OK
      if( ep.lt.0 .or. ep.gt. max_ep ) then
         write(*,220) date, sectag, ep, mjd, start_ep, space_ep*1440
 220     format('**ERROR** Problem with epochs at ',i4,4i3,1x,f6.2,
     .          ' Computed EP ',i6,' MJD, Start ',2F12.5,' Space ',
     .           f10.1,' min')
         stop
      end if

      get_clk_ep = ep

***** Thats all
      return
      end

CTITLE PTOL

      integer*4 function ptol(prn)

      implicit none

*     Function to return list entry for a specific PRN

      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 prn

* LOCAL VARIABLES
      integer*4 i,j   ! Loop counters

***** Loop over list of PRNs and see which one this is
      j = 0
      do i = 1, num_sat
         if( prn.eq.prn_list(i) ) j = i
      end do
      if( j.eq.0 ) then
          write(*,220) prn, num_sat
 220      format('**ERROR** Unable to find PRN ',i4,' in list of ',
     .           i3,' statellites')
           stop
      end if

      ptol = j

****  Thats all
      return 
      end

CTITLE FIT_CLK

      subroutine fit_clk(use_comm)

      implicit none

*     Routine to fit linear clocks to the stations and satellites.  These clocks
*     are used to align the network clock fits.

      include 'merge_igs_clk.h'

* LOCAL VARIABLES
      integer*4 i,k, nt  ! Loop counters

* MOD TAH 200627: Added removal of reference frame sites
      integer*4 good_ rf(max_site) ! Set to 1 if good reference site, 0
                         ! if not (based on data being edited)

      logical use_comm ! Set true to use just common data
     .,       refs     ! Logical indicating this is a reference site

* MOD TAH 200628: Added GNSS capability
      character*1 offcon  ! Function to get GNSS code from full PRN Number
     .,       gnss        ! GNSS type
      integer*4 prn       ! Non-GNSS PRN number.

****  First get the linear trends and fits for the reference clocks
      do i = 1, num_site
         use_rcv(i) = .true.
      end do
      do i = 1, num_sat 
        use_sat(i) = .true.
      end do

      write(*,120) 'REFERENCE SITE FITS: USE_COMM',use_comm 
 120  format(/,a,1x,L1)
      good_rf = -1
      do i = 1, num_ref
         do nt = 1, num_net
            call fit_rcv('R',ref_sites(i),nt, .false., 1.d-9)
            if( nint(rcv_fit(1,ref_sites(i),nt)).gt.0  ) then
                if( nint(rcv_fit(1,ref_sites(i),nt)).eq.num_ep ) then
                   good_rf(i) = 1
                else
                   good_rf(i) = -nint(rcv_fit(1,ref_sites(i),nt))
                end if
            endif
         end do
      end do

* MOD TAH 200627: See if we should removed reference frame sites with not
*     enough data.
      k = 0
      do i = 1, num_ref
         if( good_rf(i).eq. 1 ) then
            k = k + 1
            ref_sites(k) = ref_sites(i)
         else
            write(*,140) i, site_codes(ref_sites(i)), -good_rf(i)
 140        format('REMOVING RF site ',i3,1x,a,' # data ',i5)
         endif
      end do
      num_ref = k

      if( use_comm ) then
          write(*,120) 'ALL COMMON SITES EXCEPT REFERENCE FITS' 
      else
          write(*,120) 'ALL SITES EXCEPT REFERENCE FITS' 
      end if 
      do i = 1, num_site

*        Check to see if reference site, do not run again if it
*        is.
         refs = .false.
         do k = 1,num_ref
             if( ref_sites(k).eq.i ) refs = .true.
         end do
         if ( .not.refs ) then
            do nt = 1, num_net
                call fit_rcv('R',i,nt, use_comm, 20.d-9)
             end do
         end if
      end do

      if( use_comm ) then 
          write(*,120) 'ALL COMMON SATELLITE FITS' 
      else
          write(*,120) 'ALL SATELLITE FITS' 
      endif 
      do i = 1, num_sat
         do nt = 1, num_net
            call fit_rcv('S',i,nt, use_comm, 1.d-9)
         end do
      end do

***** OK Now report on comparison
      if( use_comm ) then
         write(*,120) 'COMPARISON OF COMMON RCV FITS'
      else
         write(*,120) 'COMPARISON OF RCV FITS'
      end if 
      do i = 1, num_site
          fin_rcv(i) = .false.
          if( rcv_fit(1,i,1).gt.0.d0 .or. .not.use_comm ) then
             write(*,240) site_codes(i), 
     .            int(rcv_fit(1,i,1)),rcv_fit(2,i,1)*1.d9,
     .            rcv_fit(3,i,1)*1.d12,rcv_fit(4,i,1)*1.d12,
     .           (int(rcv_fit(1,i,k)-rcv_fit(1,i,1)),
     .            (rcv_fit(2,i,k)-rcv_fit(2,i,1))*1.d9,
     .            (rcv_fit(3,i,k)-rcv_fit(3,i,1))*1.d12,
     .             rcv_fit(4,i,k)*1.d12, k = 2,num_net)
 240         format(a4,1x,i5,2F13.3,1x,F8.0,10(1x,i5,2F13.3,1x,f8.0))
             if( int(rcv_fit(1,i,1)).gt.1 ) fin_rcv(i) = .true.

*            If number of values used exceeds one, remove from list
             do k = 2, num_net
                if( int(rcv_fit(1,i,k)).gt.1 ) fin_rcv(i) = .true.

                if( abs(rcv_fit(1,i,k)-rcv_fit(1,i,1)).gt.1 .or.
     .                  rcv_fit(1,i,k).lt.num_ep*0.9 )
     .                                      use_rcv(i) = .false.
             end do

          else
             use_rcv(i) = .false.
          end if
      end do
      write(*,120) 'COMPARISON OF SAT FITS'
      do i = 1, num_sat
          fin_sat(i) = .false.
          if( sat_fit(1,i,1).gt.0.d0 .or. .not.use_comm ) then
* MOD TAH Convert to GNSS:
             write(*,260) offcon(prn_list(i)), mod(prn_list(i),100),
     .            int(sat_fit(1,i,1)),sat_fit(2,i,1)*1.d9,
     .            sat_fit(3,i,1)*1.d12,sat_fit(4,i,1)*1.d12,
     .           (int(sat_fit(1,i,k)-sat_fit(1,i,1)),
     .            (sat_fit(2,i,k)-sat_fit(2,i,1))*1.d9,
     .            (sat_fit(3,i,k)-sat_fit(3,i,1))*1.d12,
     .             sat_fit(4,i,k)*1.d12, k = 2,num_net)
 260         format(a1,i2.2,2x,i5,2F13.3,1x,F8.0,
     .                                    10(1x,i5,2F13.3,1x,f8.0))
             if( int(sat_fit(1,i,1)).gt.1 ) fin_sat(i) = .true.
*            If number of values used exceeds one, remove from list
             do k = 2, num_net
                if( int(sat_fit(1,i,k)).gt.1 ) fin_sat(i) = .true.
                if( abs(sat_fit(1,i,k)-sat_fit(1,i,1)).gt.1 .or.
     .                 sat_fit(1,i,k).lt.num_ep*0.9 )
     .                                     use_sat(i) = .false.
 
             end do
          else
             use_sat(i) = .false.
          end if
      end do

      return
      end

CTITLE FIT_RCV

      subroutine fit_rcv(typ,nr,nt, use_comm, tol)

      implicit none

*     Routine to fit linear clocks ground receivers and Satellites


      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 nr    ! Receiver number
     .,         nt    ! Network number

      logical use_comm ! Set true to use just common data accross
                      ! all networks

      real*8 tol      ! Tolerance of fit of data 
      character*(*) typ  ! R for reciever/S for satellites

* LOCAL VARIABLES

      integer*4 num   ! Number of valid values
     .,         ep    ! Epoch number
     .,         ref_ep ! Epoch of reference values (usually 1)
     .,         nb    ! Number bad.  (If exceed 10% data, reset)


      real*8 res ! residual seconds
     .,      W   ! Weight
     .,      dt  ! Time difference in epochs
     .,      neq(2,2)  ! Normal eqautions
     .,      bv(2)     ! B-vector
     .,      sol(2)    ! Solution to linear fit
     .,      rms       ! RMS to fit (sec)
     .,      schi      ! Sum of weighted residual squared
     .,      ref_val   ! first value (to remove large offsets)
     .,      ref_dval  ! Change in value over 1-epoch
     .,      prev_res  ! Previous epoch residual
     .,      det       ! Determinate of neq
     .,      ddres     ! Difference between residual and expected value

      logical bad     ! Set true if clock does not fit tolerance
      logical comm_dat  ! Logical function to test if data is Common

* MOD TAH 200628: Added GNSS capability
      character*1 offcon  ! Function to get GNSS code from full PRN Number


****  Loop over epochs accumulating normal equations and stats
      ref_val = 0.d0
      ref_ep = 0
      ref_dval = 0.d0
      prev_res = 0.d0

*     Clear estimation values
      nb = 0

      call clear_neq( neq, bv, num, schi)

*     Now loop over times
      do ep = 1, num_ep
c         if( rcv_flgc(nr,nt,ep).ne. 0 ) then
          if( comm_dat(typ,nr, nt, ep, use_comm) ) then
             bad = .false.
*            OK; we have a value at this epoch
             if( ref_val.eq. 0.d0 ) then
                 if ( typ(1:1).eq.'R' ) then
                    ref_val = rcv_clk(nr,nt,ep)
                 else
                    ref_val = sat_clk(nr,nt,ep)
                 end if

                 ref_ep = ep
             end if

             if( typ(1:1).eq.'R' ) then
                res = rcv_clk(nr,nt,ep)-ref_val
             else
                res = sat_clk(nr,nt,ep)-ref_val
             end if

*            Check quality, initial points are a little difficult
*            because of rate not known yet
             if( prev_res.ne.0 ) then
                 if( ref_dval.ne.0 ) then
                     ddres = res - (prev_res+ref_dval*(ep-ref_ep))
                     if( abs(ddres).gt.tol ) then
                         bad = .true.
                         nb = nb + 1
                         if( debug )
     .                   write(*,220) ep, nr, nt, res*1e9,  
     .                        prev_res*1e9,ref_dval*1e9, ref_ep
 220                     format('BAD CLK Ep ',i4,' Site/net ',2i4,
     .                          ' Residuals ',3F10.3,' Ref ',i4)
                         if( nb.gt. 2*num ) then
*                            reset reference values
                             if( debug )
     .                           write(*,240) nb, num, ep
 240                         format('Resetting: Num Bad ',i4,
     .                              ' Num Good ',i4,' Ref Ep ',i5)
                             prev_res = 0.d0
                             ref_val  = 0.d0
                             ref_dval = 0.d0
                             ref_ep   = 0
                             call clear_neq( neq, bv, num, schi)
                             nb = 0
                         end if
                     else
                         ref_dval = res/(ep-ref_ep)
                     endif
                 else
                     ref_dval = res/(ep-ref_ep)
                 end if
             else
                 prev_res = res + 1.d-12
             end if

*            OK: If residual is not bad, add to normal equations
             if( .not.bad ) then
                if( typ(1:1).eq.'R' ) then
                    W = 1.d0/rcv_sgc(nr,nt,ep)**2
                else
                    W = 1.d0/sat_sgc(nr,nt,ep)**2
                end if
c               W = 1.d0
                dt = (ep-1440)*1.d0
                call inc_neq(res, W, dt, neq, bv, schi, num)
             end if
         end if
      end do

****  Finished accumlation, now solve solutions
      if( num.gt.num_ep*0.05 ) then
         det = (neq(1,1)*neq(2,2) - neq(1,2)**2)
         sol(1) = (bv(1)*neq(2,2)-bv(2)*neq(1,2))/det
         sol(2) = (bv(2)*neq(1,1)-bv(1)*neq(1,2))/det

         rms = sqrt(abs(schi-(sol(1)*bv(1)+sol(2)*bv(2)))/neq(1,1))
         if ( rms.lt.10d-12 ) then
             if( typ.eq.'R' ) then
                write(*,310) site_codes(nr),nt, num, rms*1.d12
 310            format('RMS too small Site ',a4,'  Net ',i3,' Num ',i5,
     .                 ' RMS ',F10.4, ' ps')
             else
                write(*,315) offcon(prn_list(nr)),mod(prn_list(nr),100),
     .                nt, num, rms*1.d12
 315            format('RMS too small Sat ',a1,I2.2,'  Net ',i3,' Num ',
     .                 i5,' RMS ',F10.4, ' ps')
             endif

             rms = 10d-12
         end if

*        If needed map refernce value back to epoch 1
         ref_val = ref_val + sol(2)*(1-ref_ep)

         if( typ.eq.'R' ) then
            write(*,320) nr, site_codes(nr), nt, num, 
     .                  (sol(1)+ref_val)*1.d9,
     .                   sol(2)*1.d12, rms*1.d9
 320        format('Ref: ',i3,1x,a4,' Net ',i3,' Num ',i5,' Off ',F12.3,
     .            ' Rate ',f12.3,' RMS ',f10.3,' ns')

*           Save the fit values
            rcv_fit(1,nr,nt) = num
            rcv_fit(2,nr,nt) = sol(1) + ref_val
            rcv_fit(3,nr,nt) = sol(2)
            rcv_fit(4,nr,nt) = rms

         else
            write(*,330) nr, offcon(prn_list(nr)), 
     .                   mod(prn_list(nr),100), nt, num, 
     .                  (sol(1)+ref_val)*1.d9,
     .                   sol(2)*1.d12, rms*1.d9
 330        format('Ref:  ',I2,1x,a1,i2.2,'  Net ',i3,' Num ',i5,' Off ',
     .              F12.3,' Rate ',f12.3,' RMS ',f10.3,' ns')

*           Save the fit values
            sat_fit(1,nr,nt) = num
            sat_fit(2,nr,nt) = sol(1) + ref_val
            sat_fit(3,nr,nt) = sol(2)
            sat_fit(4,nr,nt) = rms
         end if

       endif

****   Thats all
       return
       end

CTITLE INC_NEQ

      subroutine inc_neq(res, W, dt, neq, bv, schi, num)

      implicit none

*     Routine to increment the normal eqautions

      real*8 res ! residual seconds
     .,      W   ! Weight
     .,      dt  ! Time difference in epochs
     .,      neq(2,2)  ! Normal eqautions
     .,      bv(2)     ! B-vector
     .,      schi      ! Sum of weighted residual squared

      integer*4 num    ! Number of values in fit


****  Increment normal eqautions and statistics
      num = num + 1
      neq(1,1) = neq(1,1) + W
      neq(1,2) = neq(1,2) + W*dt
      neq(2,1) = neq(2,1) + W*dt
      neq(2,2) = neq(2,2) + W*dt**2
      bv(1) = bv(1) + res*W
      bv(2) = bv(2) + res*dt*W
      schi  = schi + res**2*W

****  Thats all
      return
      end

CTITLE CLEAR_NEQ

      subroutine clear_neq( neq, bv, num, schi)

      implicit none

*     Clear the normal eqautions
      real*8 neq(2,2), bv(2), schi
      integer*4 num, i, j

      num = 0
      do i = 1,2
         bv(i) = 0.d0
         do j = 1,2
            neq(i,j) = 0.d0
         end do
      end do
      schi = 0.d0

      return
      end

CTITLE COMM_DAT

      logical function comm_dat(typ, nr, nt, ep, use_comm)

      implicit none

*     Logical function to test if data is common to all receivers
*     or satellites  

      include 'merge_igs_clk.h'

* PASSED VARIABLES
      character*(*) typ   ! Set R for reciever, S for Satellite
      integer*4 nr     ! Receover or satellite number
     .,         nt     ! Network being processed
     .,         ep     ! Epoch being processed

      logical use_comm ! Set true to use common data

      integer*4 i

****  Start, see if receiver or satellites
      comm_dat = .true.
      if ( typ(1:1).eq.'R' ) then
          if( .not.use_comm ) then  ! Just check this site
              if( rcv_flgc(nr,nt,ep).ne. 1 .or. 
     .            rcv_sgc(nr,nt,ep).gt.100.d-9 ) comm_dat = .false. 
          else
              do i = 1, num_net
                  if( rcv_flgc(nr,i,ep).ne. 1 .or. 
     .                rcv_sgc(nr,i,ep).gt.100.d-9 ) comm_dat = .false.
              end do
          end if
      else
*         Check the satellite flags
          if( .not.use_comm ) then  ! Just check this site
              if( sat_flgc(nr,nt,ep).eq. 0 .or. 
     .            sat_sgc(nr,nt,ep).gt.100.d-9 ) comm_dat = .false. 
          else
              do i = 1, num_net
                  if( sat_flgc(nr,i,ep).eq. 0. or. 
     .                sat_sgc(nr,i,ep).gt.100.d-9  ) comm_dat = .false.
              end do
          end if
      end if

***** Thats all
      return
      end

CTITLE ALLIGN_CLK

      subroutine allign_clk

      implicit none

*     Routine to allign the offset and rates of each network by
*     computing average offset and rate for each network difference
*
 
      include 'merge_igs_clk.h'

* LOCAL VARIABLES
      integer*4 i,j,k   ! Loop counters
     .,         nt, nr  ! Network and site/sat
     .,         nums    ! Number of values in stats
     .,         ep      ! Epoch counter

      real*8 stats(3,2)  ! Stats for offset and rate
     .,      means(2), sig(2), rms(2)  ! Mean, +- and rms for offset/rate
     .,      doff, drat, W  ! Offset, rate and weight
     .,      summ        ! Summ of final Mean
     .,      offset      ! Offset to be applied as function of time


****  Loop over all the networks above 1 and get mean (weighted) offsets
*     and rates

      write(*,'(/,a)') 'FINAL ALLIGNMENT'
      do nt = 2, num_net
*        Clear the statistics and loop over sites
         nums = 0
         do i = 1, 3
            do j = 1, 2
               stats(i,j) = 0.d0
            end do
         end do

         do nr = 1, num_site 
*           Make sure number of estimates matchs
            if( abs(rcv_fit(1,nr,nt)-rcv_fit(1,nr,1)).le.1.d0 .and.
     .          rcv_fit(1,nr,1).gt.120 ) then
                doff = rcv_fit(2,nr,nt)-rcv_fit(2,nr,1)
                drat = rcv_fit(3,nr,nt)-rcv_fit(3,nr,1) 
                W = 1.d0/rcv_fit(4,nr,nt)**2

*               OK: Sum up statistics for offset and rate
                nums = nums + 1
                stats(1,1) = stats(1,1)+W
                stats(2,1) = stats(2,1)+doff*W
                stats(3,1) = stats(3,1)+doff**2*W
                stats(1,2) = stats(1,2)+W
                stats(2,2) = stats(2,2)+drat*W
                stats(3,2) = stats(3,2)+drat**2*W
            end if
         end do
*        Now add in satellites
         do nr = 1, num_sat 
*           Make sure number of estimates matchs
            if( abs(sat_fit(1,nr,nt)-sat_fit(1,nr,1)).le.1.d0 .and.
     .          sat_fit(1,nr,1).gt.num_ep*0.9 ) then
                doff = sat_fit(2,nr,nt)-sat_fit(2,nr,1)
                drat = sat_fit(3,nr,nt)-sat_fit(3,nr,1) 
                W = 1.d0/sat_fit(4,nr,nt)**2

*               OK: Sum up statistics for offset and rate
                nums = nums + 1
                stats(1,1) = stats(1,1)+W
                stats(2,1) = stats(2,1)+doff*W
                stats(3,1) = stats(3,1)+doff**2*W
                stats(1,2) = stats(1,2)+W
                stats(2,2) = stats(2,2)+drat*W
                stats(3,2) = stats(3,2)+drat**2*W
            end if
         end do

****     OK: Now finish up statistics
         if( nums.gt.1 ) then
             do i = 1,2
                 means(i) = stats(2,i)/stats(1,i)
                 rms(i) = sqrt((stats(3,i)-means(i)*stats(2,i))/
     .                          stats(1,i))
                 sig(i) = sqrt(rms(i)**2/(nums-1))
             end do
             write(*,220) nt, nums, 
     .                    means(1)*1.d9, sig(1)*1.d9, rms(1)*1.d9,
     .                    means(2)*1.d12,sig(2)*1.d12,rms(2)*1.d12

 220         format('For Net ',i2,' Num ',i3,
     .              ' Offset ',F6.3,' +- ',F6.3,' RMS ',f6.3,' ns',
     .              ' Rate   ',F6.3,' +- ',F6.3,' RMS ',f6.3,' ps/d')

*            Save the offsets
             net_mean(1,nt) = means(1)
             net_mean(2,nt) = means(2)
         end if
      end do

****  Finishup the offsets for each net
      do k = 1,2
         summ = 0.d0
         do nt = 2,num_net
            summ = summ + net_mean(k,nt)
         end do
         net_mean(k,1) = -summ/num_net
         do nt = 2, num_net
             net_mean(k,nt) = net_mean(k,nt)+net_mean(k,1)
         end do
      end do

      write(*,320) (net_mean(1,nt)*1.d9, nt=1,num_net)
      write(*,340) (net_mean(2,nt)*1.d12,nt=1,num_net)
 320  format(/,'BIASES BY NETWORK',/,
     .         'OFFSETS (ns)  ',100F10.3)
 340  format(  'RATES   (ps/ep) ',100F10.3)

****  Apply the allignment to all the clocks
      do nt = 1, num_net
         do ep = 1, num_ep
            offset = net_mean(1,nt) + net_mean(2,nt)*(ep-1440)
            do nr = 1, num_site
               rcv_clk(nr,nt,ep) = rcv_clk(nr,nt,ep) - offset
            end do
            do nr = 1, num_sat
               sat_clk(nr,nt,ep) = sat_clk(nr,nt,ep) - offset
            end do
         end do
      end do
        

****  Thats all
      return
      end
 
CITTLE GEN_MESH

      subroutine gen_mesh

      implicit none

*     Routine to generate meshed series that can simply be averaged         

      include 'merge_igs_clk.h'

* LOCAL VARIABLES
      integer*4 ep  ! Epoch counter
     .,         nt, nr  ! Loop counters
     .,         i,j 

      real*8 net_diffs(max_net)  ! Correction to be applied to
                                 ! each network
      logical test_com    ! Set true for common sites

      character*256 line

* MOD TAH 200628: Added GNSS capability
      character*1 offcon  ! Function to get GNSS code from full PRN Number


****  Here we loop over all the networks at each epoch and compute
*     and offset for each network.
      write(*,'(a)') 'ALIGNMENT PROCESS:'
      line = ' '
      j = 1
      do i = 1, num_site
          if( use_rcv(i) ) then
              line(j:) = site_codes(i)
              j = j + 5
           end if
      enddo

****  Check to see if we have any sites
      if ( j.eq.1 ) then
         write(*,120) 
 120     format('**WARNING** No common sites, Scanning')
*        There are no sites.  Scan for just common sites
         do nr = 1, num_site
            use_rcv(nr) = .false.
            do ep = 1,num_ep
               test_com = .true.
               do nt = 1, num_net
                   if( rcv_flgc(nr,nt,ep).ne.1 ) test_com = .false.
               end do
               if( test_com ) use_rcv(nr) = .true.
            end do
*           Now build list
            if( use_rcv(nr) ) then
                line(j:) = site_codes(nr)
                j = j + 5
            end if
         end do
      end if
            

      write(*,140) line(1:j)
 140  format('SITES: ',a)
      j = 1
      line = ' '
      do i = 1, num_sat
          if( use_sat(i) ) then
             write(line(j:),160) offcon(prn_list(i)), 
     .                              mod(prn_list(i),100)
 160         format(a1,i2.2)
             j = j + 5
          end if
      end do
      write(*,180) line(1:j)
 180  format('SATS  : ',a)

      do ep = 1, num_ep
*
*        For all common sites and satellites, get the differences
         call diff_net(ep, net_diffs)

*        For each network, apply the differences
         do nt = 1, num_net
            do nr = 1, num_site
                rcv_clk(nr,nt,ep) = rcv_clk(nr,nt,ep)-net_diffs(nt)
            end do
            do nr = 1, num_sat
                sat_clk(nr,nt,ep) = sat_clk(nr,nt,ep)-net_diffs(nt)
            end do
         end do

      end do

      end

CTITLE DIFF_NET

      subroutine diff_net(ep, net_diffs)

      implicit none

*     Routine to compute differences of networks at common sites/sats

      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 ep   ! Epoch being processed
      real*8 net_diffs(max_net)  ! Offsets to be applied to each network


* LOCAL VARIABLES
      integer*4 i,j, nt, nr      ! Loop counters
     .,         nums(max_net)    ! Number in stats

      logical comm_dat       ! Logical function to see if common data
     .,       use_comm       ! Set true to use common data

      real*8 stats(3,max_net)  ! Stats of differences between networds
     .,      doff, W           ! Offset and weight
     .,      mean(max_net), sig(max_net), rms(max_net)  ! Mean, sig and
                               ! rms of differences between networks
     .,      sumn              ! Sum of adjustments to network 


 
****  Find the common data
      use_comm = .true.

      do i = 1, 3
         do j = 1, num_net
             nums(j) = 0
             stats(i,j) = 0.d0
         end do
      end do
      do nr = 1, num_site     
          if( comm_dat('R',nr, 1, ep, use_comm) .and.
     .        use_rcv(nr) ) then
             do nt = 2,num_net
                doff = rcv_clk(nr,nt,ep)-rcv_clk(nr,1,ep)

                W = 1.d0/(rcv_sgc(nr,nt,ep)**2+
     .                    rcv_sgc(nr, 1,ep)**2)

*               OK: Sum up statistics for offset and rate
                if( abs(doff).lt.1.d-9 ) then 
                   nums(nt) = nums(nt) + 1
                   stats(1,nt) = stats(1,nt)+W
                   stats(2,nt) = stats(2,nt)+doff*W
                   stats(3,nt) = stats(3,nt)+doff**2*W
                end if
             end do
          end if
      end do 

      do nr = 1, num_sat     
          if( comm_dat('S',nr, 1, ep, use_comm) .and.
     .        use_sat(nr) ) then
             do nt = 2,num_net
                doff = sat_clk(nr,nt,ep)-sat_clk(nr,1,ep)
                W = 1.d0/(sat_sgc(nr,nt,ep)**2+
     .                    sat_sgc(nr, 1,ep)**2)

*               OK: Sum up statistics for offset and rate
                if( abs(doff).lt.1.d-9 ) then
                   nums(nt) = nums(nt) + 1
                   stats(1,nt) = stats(1,nt)+W
                   stats(2,nt) = stats(2,nt)+doff*W
                   stats(3,nt) = stats(3,nt)+doff**2*W
                end if
             end do
          end if
      end do 

****  Now finish up statistics
      do i = 2,num_net
          if( nums(i).gt.1 ) then
              mean(i) = stats(2,i)/stats(1,i)
              rms(i) = sqrt((stats(3,i)-mean(i)*stats(2,i))/
     .                   stats(1,i))
              sig(i) = sqrt(rms(i)**2/(nums(i)-1))
          else
              mean(i) = 0.d0
              rms(i) = 0.d0
              sig(i) = 0.d0
          end if
      end do


***** Now save the net diffs (making average zero)
      sumn = 0.d0
      do i = 2, num_net
         sumn = sumn + mean(i)
      end do
      net_diffs(1) = -sumn/num_net
      do i = 2,num_net
         net_diffs(i) = mean(i)+net_diffs(1)
      end do 

*     Write out result
      if( debug ) then
         write(*,220) ep, (nums(i),mean(i)*1.d12, sig(i)*1.d12, 
     .                rms(i)*1.d12, i=2,num_net)
 220     format('Ep ',i4,
     .        100(' Num ',i3,1x,F8.1,' +- ',f8.1,1x,f8.1:))
         write(*,240) ep, (net_diffs(i)*1.d12,i=1,num_net)
 240     format('EP ',i4,' Corrections ps ', 100(1x,F8.1))
      end if
 

****  Thats all
      Return
      end
 
CITTLE GEN_FIN

      subroutine gen_fin

      implicit none

*     Routine to generate the averaged final combinated clock estimates 

      include 'merge_igs_clk.h'

* LOCAL VARIABLES
      integer*4 ep  ! Epoch counter
     .,         i,j,k   ! Loop counters
     .,         num_fin  ! Final number of reference clocks
     .,         num      ! Number counter

* MOD TAH 200528: Add GNSS separation
      integer*4 prn        ! PRN with GNSS offset removed
      character*1 offcon   ! function gets GNSS type from PRN
     .,         gnss       ! GNSS type

      real*8 res, W     ! residual and weight
     .,      fin_fit(2,max_site)  ! Final linear fit
     .,      neq(2,2), bv(2), sol(2), det     ! Fit
     .,      dc, dcm           ! Offsets and mean
     .,      rms, schi, dt   ! Stats
     .,      stats(3)        ! Summation
     .,      rms_fin(max_site)
     .,      rms_min, rms_min_save
     .,      rms_best        ! RMS of best reference clock

      logical fin_ref(max_site)  ! Set true for final reference sites



****  Here we loop over all the networks at each epoch and compute
*     and offset for each network.
      do ep = 1, num_ep
*
*        Average the nets
         call average_net(ep)
         if( debug ) then
            do i = 1, num_site
               if( rcv_flag(i,ep).ne.0 ) then
                  write(*,120) site_codes(i), ep, 
     .                  rcv_fin(i,ep)*1.d9, rcv_sig(i,ep)*1.d9
 120              format('AR ',a4,1x,i5,1x, 2f16.3)
               end if
            end do
            do i = 1, num_sat
               if( sat_flag(i,ep).ne.0 ) then
* MOD TAH 200628: Convert PRN back to GNSS code and prn
                  gnss = offcon(prn_list(i))
                  prn = mod(prn_list(i),100)
                  write(*,140) gnss, prn, ep, 
     .                  sat_fin(i,ep)*1.d9, sat_sig(i,ep)*1.d9
 140              format('AS ',a1,i2.2,2x,i5,1x,2f16.3)
               end if
            end do
         end if

      end do

****  Now compute the statistics of the final clocks
      num_fin = 0
      rms_min = 1.d-6
      best_ref = 0
      rms_best = 1.d30 
      do i = 1, num_site
         fin_ref(i) = .false.
* MOD TAH 200627: Only out put receivers with RMS < 1.e-9
         fin_rcv(i) = .false.
         call clear_neq( neq, bv, num, schi) 

         do ep = 1, num_ep
            if( rcv_flag(i,ep).ne.0 .and.
     .          rcv_sig(i,ep).lt. 10.e-9 ) then
                W = 1.d0/rcv_sig(i,ep)**2
                dt = (ep-1440)*1.d0
                res = rcv_fin(i,ep)
                call inc_neq(res, W, dt, neq, bv, schi, num)
            end if
         end do
*        OK, Now do final solution 
         if( num.gt.num_ep*0.95 ) then
            det = (neq(1,1)*neq(2,2) - neq(1,2)**2)
            sol(1) = (bv(1)*neq(2,2)-bv(2)*neq(1,2))/det
            sol(2) = (bv(2)*neq(1,1)-bv(1)*neq(1,2))/det
            fin_fit(1,i) = sol(1)
            fin_fit(2,i) = sol(2)
         end if

*        Now re-do the residual calculation
         call clear_neq( neq, bv, num, schi) 
         do ep = 1, num_ep
            if( rcv_flag(i,ep).ne.0 .and.
     .          rcv_sig(i,ep).lt. 10.e-9 ) then
                W = 1.d0/rcv_sig(i,ep)**2
                dt = (ep-1440)*1.d0
                res = rcv_fin(i,ep) - (fin_fit(1,i)+dt*fin_fit(2,i))
                call inc_neq(res, W, dt, neq, bv, schi, num)
            end if
         end do
         if( num.gt.num_ep*0.95 ) then
            det = (neq(1,1)*neq(2,2) - neq(1,2)**2)
            sol(1) = (bv(1)*neq(2,2)-bv(2)*neq(1,2))/det
            sol(2) = (bv(2)*neq(1,1)-bv(1)*neq(1,2))/det
            fin_fit(1,i) = sol(1) + fin_fit(1,i)
            fin_fit(2,i) = sol(2) + fin_fit(2,i)     
* 
            rms = sqrt(abs(schi-(sol(1)*bv(1)+sol(2)*bv(2)))/neq(1,1))
            rms_fin(i) = rms
            if( rms.lt.rms_min ) rms_min = rms
            if( rms.lt.0.05d-9 ) then
               write(*,220) site_codes(i), rms*1.e9, num
 220           format('TEST CLOCK REF ',a4,' RMS ',f6.3,' ns. Num ',
     .                 I4)
               num_fin = num_fin + 1
               fin_ref(i) = .true.
               if( rms .lt. rms_best ) then
                   best_ref = i
                   rms_best = rms
               endif
            end if
*           Now output all "reasonable" clock values.  Every thing < 1 ns RMS
            if( rms.lt.1.0d-9 ) then
               write(*,240) site_codes(i), rms*1.e9, num
 240           format('OKCLK ',a4,' RMS ',f6.3,' ns. Num ', I4)
               fin_rcv(i) = .true.
            end if


         end if
      end do

****  OK See if we have any finals
      if( num_fin.lt.2 ) then
          write(*,280) rms_min*2.0e9
 280      format('NOT ENOUGH RMS VALUES LESS THAN 50 ps, using RMS < ',
     .            f8.3,' ns')
          do i = 1, num_site
             if( rms_fin(i).le.rms_min*2.0 .and. 
     .           rms_fin(i).gt.1.d-12) then
                write(*,290) site_codes(i), rms_fin(i)*1.e9
 290            format('TEST CLOCK MIN ',a4,' RMS ',f6.3,' ns')
                num_fin = num_fin + 1
                fin_ref(i) = .true.
                if( rms_fin(i) .lt. rms_best ) then
                    best_ref = i
                    rms_best = rms_fin(i)
                endif
             end if
          end do
      end if                  
 
*     OK Final adjustment
      write(*,320) num_fin
 320  format('FINAL ADJUSTMENT WITH ',i2,' CLOCKS')
c      debug = .true.
      do ep = 1, num_ep
         num = 0
         do k = 1,3
            stats(k) = 0.d0
         end do
         do j = 1, num_site
            if( fin_ref(j) .and. rcv_flag(j,ep).ne.0 .and.
     .          rcv_sig(j,ep).lt. 10.e-9 ) then
                dc = rcv_fin(j,ep) - (fin_fit(1,j)+fin_fit(2,j)*
     .                                             (ep-1440))
                W = 1.d0/rcv_sig(j,ep)**2
                stats(1) = stats(1)+W
                stats(2) = stats(2)+dc*W
                stats(3) = stats(3)+dc**2*W
                num = num + 1
            endif
         enddo
         if( num.gt.1 ) then
             dcm =  stats(2)/stats(1)
*            Now remove this value from all sites and satellites
             if( debug ) then
                write(*,330) ep, num, dcm*1.d9, 
     .                       sqrt(stats(3)/stats(2))*1.d9
 330            format('Fit EP ',i4,1x,i4,1x,2F10.2)
             end if
 
             do j = 1, num_site
                rcv_fin(j,ep) = rcv_fin(j,ep) - dcm
             end do
             do j = 1, num_sat
                sat_fin(j,ep) = sat_fin(j,ep) - dcm
             end do
         else
             write(*,340) ep
 340         format('**WARNING** No reference clocks at Epoch ',i4)
         endif
         if( debug ) then
           do i = 1, num_site
              if( rcv_flag(i,ep).ne.0 ) then
                 write(*,420) site_codes(i), ep, 
     .             rcv_fin(i,ep)*1.d9, rcv_sig(i,ep)*1.d9
 420             format('AR ',a4,1x,i5,1x, 2f16.3,' FIN')
              end if
           end do
           do i = 1, num_sat
              if( sat_flag(i,ep).ne.0 ) then
* MOD TAH 200628: Converted to get GNSS and PRN
                 gnss = offcon(prn_list(i))
                 prn = mod(prn_list(i),100)
                 write(*,440) gnss, prn, ep, 
     .              sat_fin(i,ep)*1.d9, sat_sig(i,ep)*1.d9
 440             format('AS ',a1,i2.2,2x,i5,1x,2f16.3,' FIN')
              end if
           end do
         end if

      end do

****  Now that we have done final adjustment-reestimate the RMS of 
*     the reference cites 
      num_fin = 0
      rms_min_save = max(rms_min,0.05d-9)  ! Old minimum value.
      rms_min = 1.d-6
      best_ref = 0
      rms_best = 1.d30 
      do i = 1, num_site
         if( fin_ref(i) ) then 
            call clear_neq( neq, bv, num, schi) 

            do ep = 1, num_ep
               if( rcv_flag(i,ep).ne.0 .and.
     .             rcv_sig(i,ep).lt. 10.e-9 ) then
                   W = 1.d0/rcv_sig(i,ep)**2
                   dt = (ep-1440)*1.d0
                   res = rcv_fin(i,ep) - (fin_fit(1,i)+dt*fin_fit(2,i))
                   call inc_neq(res, W, dt, neq, bv, schi, num)
               end if
            end do
            det = (neq(1,1)*neq(2,2) - neq(1,2)**2)
            sol(1) = (bv(1)*neq(2,2)-bv(2)*neq(1,2))/det
            sol(2) = (bv(2)*neq(1,1)-bv(1)*neq(1,2))/det
            fin_fit(1,i) = sol(1) + fin_fit(1,i)
            fin_fit(2,i) = sol(2) + fin_fit(2,i)     
* 
            rms = sqrt(abs(schi-(sol(1)*bv(1)+sol(2)*bv(2)))/neq(1,1))
            rms_fin(i) = rms
            if( rms.lt.rms_min ) rms_min = rms
            if( rms.lt.rms_min_save*2 ) then
               write(*,620) site_codes(i), rms*1.e9, num
 620           format('FINAL CLOCK REF ',a4,' RMS ',f6.3,' ns. Num ',
     .                 I4)
               num_fin = num_fin + 1
               if( rms .lt. rms_best ) then
                   best_ref = i
                   rms_best = rms
               endif
            end if
         end if
      end do

      end

CTITLE AVERAGE_NET

      subroutine average_net(ep)

      implicit none

*     Routine to average values from networks at common sites/sats

      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 ep   ! Epoch being processed


* LOCAL VARIABLES
      integer*4 i, nt, nr      ! Loop counters
     .,         nums    ! Number in stats
     .,         wnet    ! Network with worst residuals

      logical use_comm       ! Set true to use common data
     .,       convrgd        ! Set true when satellite iteration converged
     .,       edited         ! Set truw if data edited.

      real*8 stats(3)  ! Stats of differences between networds
     .,      doff, W           ! Offset and weight
     .,      mean, sig, rms    ! Mean, sig and
                               ! rms of differences between networks 
     .,      wres              ! Worst satellite residual
     .,      ref               ! reference value
      real*4 save_res(max_net) ! Saved residuals so they can be print
                               ! when data deleted.

* MOD TAH 200628: Added GNSS capability
      character*1 offcon  ! Function to get GNSS code from full PRN Number

 
****  Find the common data
      use_comm = .false.

      do nr = 1, num_site
          nums = 0
          do i = 1, 3
             stats(i) = 0.d0
          end do     
          do nt = 1,num_net

*            See if this value OK
             if( rcv_flgc(nr,nt,ep).eq.1 ) then
                doff = rcv_clk(nr,nt,ep)
                W = 1.d0/(rcv_sgc(nr,nt,ep)**2)

                nums = nums + 1
                stats(1) = stats(1)+W
                stats(2) = stats(2)+doff*W
                stats(3) = stats(3)+doff**2*W
             end if
          end do
*         Finish up and averaging that we need
          rcv_flag(nr,ep) = 1
          if ( nums.gt.1 ) then
             mean = stats(2)/stats(1)
             rms  = sqrt((stats(3)-mean*stats(2))/
     .                  stats(1))
             sig = max(sqrt(rms**2/(nums-1)),
     .              sqrt(1.d0/stats(1)))
          elseif ( nums.eq.1 ) then
             mean = stats(2)/stats(1)
             rms  = 0
             sig  = sqrt(1.d0/stats(1))
          else
*            No data at all.  Mark as bad
             mean = 0
             rms  = 0
             sig  = 0
             rcv_flag(nr,ep) = 0
          end if 
          rcv_fin(nr,ep) = mean
          rcv_sig(nr,ep) = max(sig*4,0.2d-9)  
      end do 

      do nr = 1, num_sat 
          convrgd = .false. 
          edited = .false.

          do while ( .not.convrgd )
             nums = 0
             do i = 1, 3
                stats(i) = 0.d0
             end do

             ref = sat_clk(nr,1,ep)
      
             do nt = 1,num_net

*               See if this value OK
                if( sat_flgc(nr,nt,ep).ne.0 ) then
                   doff = sat_clk(nr,nt,ep) - ref 
                   W = 1.d0/(sat_sgc(nr,nt,ep)**2)

                   nums = nums + 1
                   stats(1) = stats(1)+W
                   stats(2) = stats(2)+doff*W
                   stats(3) = stats(3)+doff**2*W
                end if
             end do

*            Finish up and averaging that we need
             sat_flag(nr,ep) = 1
             if ( nums.gt.1 ) then
                mean = stats(2)/stats(1)
                rms  = sqrt((stats(3)-mean*stats(2))/
     .                  stats(1))
                sig = sqrt(rms**2/(nums-1))
             elseif ( nums.eq.1 ) then
                mean = stats(2)/stats(1)
                rms  = 0
                sig  = sqrt(1.d0/stats(1))
             else
*               No data at all.  Mark as bad
                mean = 0
                rms  = 0
                sig  = 0
                sat_flag(nr,ep) = 0
             end if
*            Check quality of fit (RMS should be less than 0.25 ns)
             if( rms.gt.0.5d-9 .and. nums.gt.2 ) then
*                Remove the worst residual
                 wres = 0.d0
                 wnet = 0
                 save_res = 0.
                 do nt = 1, num_net
                    if( sat_flgc(nr,nt,ep).ne.0 ) then 
                       doff = sat_clk(nr,nt,ep)-mean-ref
                       save_res(nt) = doff
                       if( abs(doff).gt.wres ) then
                           wres = abs(doff)
                           wnet = nt
                       end if
                    end if
                 end do
                 if( wres.gt.0.25d-9 ) then
                     write(*,220) ep, offcon(prn_list(nr)), 
     .                      mod(prn_list(nr),100), wnet, wres*1.d9,
     .                      nums, rms*1.e9, save_res(1:num_net)*1.e9
 220                 format('DELETING EP ',i5,1x,a1,i2.2,' Net ',i1,
     .                      ' Res ',F8.2,' ns, NU ',I2,' RMS ', 
     .                      F8.2,' ns, Res (ns) ',50F8.2)
                     sat_flgc(nr,wnet,ep) = 0
                     sat_sgc(nr,wnet,ep) = 300.e-9
                     edited = .true.
                 end if
             else
                 convrgd = .true.
             end if                     
             if( edited ) write(*,250) ep,offcon(prn_list(nr)), 
     .                      mod(prn_list(nr),100),  
     .         (mean+ref)*1.e9,sig*1.e9
              
          end do

          sat_fin(nr,ep) = mean+ref
          sat_sig(nr,ep) = max(sig*4,0.2d-9)  
          if( edited ) write(*,250) ep, offcon(prn_list(nr)), 
     .        mod(prn_list(nr),100),
     .        sat_fin(nr,ep)*1.e9,sat_sig(nr,ep)*1.e9
 250      format('EDITED EP ',i5,1x,a1,I2.2,' Fin ',F12.2,1x,F8.2)
   
      end do 


****  Thats all
      Return
      end

CTITLE WR_CLK

      subroutine wr_clk

      implicit none

*     Routine to write the final clock file


      include 'merge_igs_clk.h'

* LOCAL VARIABLES
      integer*4 date(5)   ! Run date
     .,         ierr      ! IOSTAT error
     .,         i,j, ep   ! Loop counte
     .,         num_clk_out  ! Number of output clocks
     .,         ls        ! Line size
     .,         trimlen   ! Length of string
     .,         drs(5), dre(5)  ! Dates for start and stop of
                          ! reference clk.

* MOD TAH 200628: Modified for GNSS 
      integer*4 prn       ! PRN number from GNSS PRN
      character*1 offcon  ! Fuction to get gnss code from gnss prn
     .,      gnss         ! gnss code 

      real*8 sectag       ! Seconds tag
     .,      mjd          ! MJD of clock
     .,      srs, sre     ! Start and stop seconds
 
      integer*4 leap_sec  ! number of leap seconds

* prog, runby -- Character strings with program and runby
      character*20 prog, runby

* line  -- Line for outputing prn list
      character*80 line
* posl(3) -- Position written as floating point so that can be
*     translated to interger

      character*12 posl(3)
      character*40 phs_cen_mod

****  Open the output file
      open(200,file=out_file,iostat=ierr,status='unknown')
      call report_error('IOSTAT',ierr,'open',out_file,1,'WR_CLK')

****  Write the header information
      write(200,120) 3.00,'CLOCK DATA'
 120  format(F9.2,11x,a10,30x,'RINEX VERSION / TYPE')
      call systime( date, sectag)
      call getenv('INSTITUTE',runby)
      write(prog,130) autcln_ver(1:trimlen(autcln_ver))
 130  format(a,'+MIG')
      write(200,135) prog, runby, date
 135  format(a20,a20,i4,'-',i2.2,'-',i2.2,1x,i2.2,':',
     .       i2.2,4x,'PGM / RUN BY / DATE')
      if( .not.abs_mod ) then
         write(200,140)
 140     format('CLK ANT Z-OFFSET(M): II/IIA 1.0230; IIR 0.0000',
     .          '              COMMENT')
      else
c          write(200,145)
c 145     format('CLK ANT Z-OFFSET(M): Absolute model IGS05_1402',
c     .          '              COMMENT')
         phs_cen_mod = clk_phsmod
* MOD TAH 200723: Added clock GNSS
         write(200,145) clk_gnss, clk_gamitver, phs_cen_mod
  145    format(a1,' GAMIT ',a5,'       ',a40,'SYS / PCVS APPLIED')
      end if
      write(200,150)                                       
 150  format('     2    AS    AR                           ',
     .      '               # / TYPES OF DATA')
      write(200,160) runby
 160  format(a20,40x,'ANALYSIS CENTER')

*     See how many stations we will output
      num_clk_out = 0
      do i = 1, num_site
         if( fin_rcv(i) ) num_clk_out = num_clk_out + 1
      end do

****  Get the leap second information.
      call get_leapsec( start_ep+2400000.5d0, 1.d20, leap_sec)
      write(200,210) -(leap_sec+19)
 210  format(i6,54x,'LEAP SECONDS')

****  Get the domes numbers for the sites.
      call get_domes(101)
*     Write out the reference clock
*     Get the start and stop times for reference
      call mjd_to_ymdhms(start_ep, drs, srs)
      call mjd_to_ymdhms(start_ep+1-30/86400.d0, dre, sre)
      write(200,215) 1, drs,srs, dre, sre
 215  format(i6,1x,i4,4i3,F10.6,1x,i4,4i3,F10.6,'# OF CLK REF')

      write(200,220) site_codes(best_ref), site_names(best_ref)(1:20)
 220  format(A4,1x,A20,15x,19x   ,1x,'ANALYSIS CLK REF')
c220  format(A4,1x,A20,15x,e19.12,1x,'ANALYSIS CLK REF')
      write(200,240) num_clk_out,clk_frame_ver
 240  format(I6,4x,a5,'     ',40x,'# OF SOLN STA / TRF')
      do i = 1, num_site
         if( fin_rcv(i) ) then 
            do j = 1,3 
               write(posl(j),'(f12.0)') site_coords(j,i)              
            end do 
            write(200,250) site_codes(i), site_names(i)(1:20),
     .                     (posl(j),j=1,3)
         
 250        format(a4,1x,a20,A11,1x,A11,1x,A11, 'SOLN STA NAME / NUM')
         end if
      end do

****  Now output number of satellites
      write(200,260) num_sat 
 260  format(i6,54x,'# OF SOLN SATS')
      do i = 1, num_sat, 15
         if( i+15.le.num_sat ) then
            ls = 15
         else
            ls = mod(num_sat,15)
            if( ls .eq. 0 ) ls = 15
         end if
         write(line,270,iostat=ierr) (offcon(prn_list(i+j-1)),
     .                      mod(prn_list(i+j-1),100),j=1,ls)
         if( ierr.ne.0 ) print *,'Ierr ',ierr, line
 270     format(15(a1,i2.2,1x))
         line(61:) = 'PRN LIST'
         write(200,'(a)' ) line(1:69)
      end do

      write(200,280)
 280  format(60x,'END OF HEADER')

****  Now write out the final values
      do ep = 1, num_ep

*        Write the receiver clocks
         mjd = start_ep + (ep-1)*space_ep + 1.d-3/86400.d0
         call mjd_to_ymdhms(mjd, date, sectag)
         sectag = nint(sectag)
         if( sectag.eq.60.d0 ) then
             date(5) = date(5) + 1
             if( date(5).eq.60 ) then
                 date(4) = date(4) + 1
                 date(5) = 0
             end if
             sectag = 0
         endif

* MOD TAH 061104: See if should output sites: Only best_ref for 
*        5-minute clocks
         do i = 1, num_site
            if( fin_rcv(i) .and. rcv_flag(i,ep).ne.0 ) then
*               Added if to test epochs
                if( mod(nint((ep-1)*space_ep*86400.d0),300).eq.0 .or.
     .              i.eq.best_ref ) then
                   write(200,320) site_codes(i), date, sectag,
     .                       2,  rcv_fin(i,ep), rcv_sig(i,ep)
 320               format('AR',1x,a4,1x,i4,4i3.2,1x,f9.6,1x,i2,
     .                  2x,2(1x,E19.12))  
                end if 
            end if
         end do

*        Now do the satellite clocks
         do i = 1, num_sat
            if( sat_flag(i,ep).ne.0  ) then
               gnss = offcon(prn_list(i))
               prn = mod(prn_list(i),100)
               write(200,340) gnss, prn, date, sectag,
     .                        1, sat_fin(i,ep), sat_sig(i,ep)
 340           format('AS',1x,a1,i2.2,2x,i4,4i3.2,1x,f9.6,1x,i2,
     .                 2x,2(1x,E19.12))
            end if
         end do
      end do

****  We are done
      close(200)
      return
      end

CITTLE offset_clk

      subroutine offset_clk( it )

      implicit none

*     Routine to check breaks in clocks and re-allign them based on fitting
*     to non-broken clocks.        

      include 'merge_igs_clk.h'

* PASSED 
      integer*4 it   ! Iteration counter (satellite alligned bewteen 
                     ! networks on 2nd and above iteration.

* LOCAL VARIABLES
      integer*4 ep, eps, epe  ! Epoch counter
     .,         i,j,k         ! Loop counters
     .,         n, m          ! Number of values in estimates 
     .,         nd            ! Data counter

      real*8 stats(3)        ! statistic summation
     .,      dc, W, dcm, rms, sig ! Clock difference, weight, mean, rms,
                             ! and sigma of mean

      logical done, OK       ! Set when epochs finished and data OK
     .,       rcv_OK(max_site,max_net)
     .,       sat_OK(max_sat,max_net)  ! Set OK if no missing data
     .,       rcv_avail(max_site,max_net) ! Set true if data available in network


* MOD TAH 200628: Added GNSS capability
      character*1 offcon  ! Function to get GNSS code from full PRN Number


****  Here we loop over all the networks at each epoch and compute
*     and offset for each network.
      write(*,'(a,1x,I2)') 'OFFSET CHECKING PROCESS: Iter ', it

*     First check which sites and satellites seem to have complete
*     times series with no breaks.
      do i = 1, num_net
         do j = 1, num_site
            OK = .true.
            rcv_avail(j,i) = .false.
            nd = 0
            do k = 1,num_ep
*              See if data OK
               if ( rcv_flgc(j,i,k).ne.1 .or.
     .              rcv_sgc(j,i,k).gt. 10.e-9 ) then
                  OK = .false.
C                 write(*,110) 'AR',i,j,k, rcv_flgc(j,i,k), 
C    .                                     rcv_sgc(j,i,k)*1.e9
 110              format(a,' Net Site EP ',2I3,I6,' F ',i3,' Sig ',
     .                  F9.1,' ns')
               else
                  nd = nd + 1
               endif
               if( OK ) rcv_avail(j,i) = .true.
              
            end do
            rcv_OK(j,i) = OK        
C           write(*,115) i,site_codes(j), nd, OK
 115        format('Net ',i2,' Site ',a4,' Complete ',I5,1x,L1)
        end do
      end do
*
*     Now do the satellites
      do i = 1, num_net
         do j = 1, num_sat
            OK = .true.
            nd = 0
            do k = 1,num_ep
*              See if data OK
               if ( sat_flgc(j,i,k).eq.0 .or.
     .              sat_sgc(j,i,k).gt. 10.e-9 ) then
                  OK = .false.
C                 write(*,110) 'AS',i,j,k, sat_flgc(j,i,k), 
C    .                                     sat_sgc(j,i,k)*1.e9
               else
                  nd = nd + 1
               endif
               
            end do
            sat_OK(j,i) = OK
C           write(*,135) i,offcon(prn_list(j)),mod(prn_list(j),100),
C   .                    nd, OK
C135        format('Net ',i2,1x,a1,i2.2,1x,' Complete ',i5,1x,L1)
        end do
      end do
*
*     Now start check offsets for those that are not OK
      do i = 1, num_net
         do j = 1, num_site
*           See if OK
            if( .not.rcv_OK(j,i) .and.   rcv_avail(j,i) ) then
*              Compute and apply offsets for those parts that
*              are off.  For sites there will only be the small
*              number of over lapping ones
               done = .false.
               write(*,140) i,site_codes(j)
 140           format('Checking offsets of Net ',i2,' Site ',a4)
               eps = 0
               epe = 0
               ep = 0
               n = 0
               do k = 1, 3
                  stats(k) = 0.d0
               end do

               do while ( .not.done )
                   ep = ep + 1
                   if( rcv_flgc(j,i,ep).eq.1 .and.
     .                 rcv_sgc(j,i,ep).lt. 10.e-9 ) then
*                     This measurement is OK, so accumulate the differences
*                     from the good networks
                      if( eps.eq.0 ) eps = ep
                      do k = 1, num_net
* MOD TAH 200623: Added test k.ne.i so net compared with itself.
                         if( rcv_OK(j,k) .and. k.ne.i) then
                            dc = rcv_clk(j,i,ep) - rcv_clk(j,k,ep)
                            W = 1.d0/(rcv_sgc(j,i,ep)**2+
     .                                rcv_sgc(j,k,ep)**2)
                            n = n + 1
                            stats(1) = stats(1)+W
                            stats(2) = stats(2)+dc*W
                            stats(3) = stats(3)+dc**2*W
                            epe = ep
                         end if
                      end do
                   else
*                     Ok, we have reached a break in the clock.  
*                     Determine the offset to be applied if we have enough
*                     data
                      if ( n.gt.1 ) then 
                         dcm = stats(2)/stats(1)
                         rms  = sqrt((stats(3)-dcm*stats(2))/
     .                          stats(1))
                         sig = sqrt(rms**2/(n-1))
                         write(*,150) i,site_codes(j), eps, epe, 
     .                                dcm*1.d9, sig*1.d9, rms*1.d9, n
 150                     format('For Net ',i2,' ',A4,' Epochs ',i4,
     .                          ' to ',i4,' Mean offset ',F8.3,' +- ',
     .                          F8.3,' RMS ',F8.3,' ns, Num ',i5)
                         do k = eps, epe
                            rcv_clk(j,i,k) = rcv_clk(j,i,k)-dcm
                         end do
                         n = 0
                         do k = 1, 3
                             stats(k) = 0.d0
                         end do
                         eps = 0
                      end if
                   end if
*                  See if we are done yet
                   if( ep.eq.num_ep ) done = .true.
               end do
*              Finally see if we have some data left over
               if ( n.gt.1 ) then 
                  dcm = stats(2)/stats(1)
                  rms  = sqrt((stats(3)-dcm*stats(2))/
     .                    stats(1))
                  sig = sqrt(rms**2/(n-1))
                  write(*,160) i,site_codes(j), eps, epe, 
     .                         dcm*1.d9, sig*1.d9, rms*1.d9, n
 160              format('End Net ',i2,' ',A4,' Epochs ',i4,
     .                 ' to ',i4,' Mean offset ',F8.3,' +- ',
     .                  F8.3,' RMS ',F8.3,' ns, Num ',i5)
                  do k = eps, epe
                     rcv_clk(j,i,k) = rcv_clk(j,i,k)-dcm
                  end do
               end if
            end if

         end do
      end do

*     Now offset the satellites
      do i = 1, num_net
         do j = 1, num_sat
*           See if OK
            if( .not.sat_OK(j,i) ) then
*              Compute and apply offsets for those parts that
*              are off.
*              We do this as a do while loop
               done = .false.
               write(*,240) i,offcon(prn_list(j)),mod(prn_list(j),100)
 240           format('Checking offsets of Net ',i2,1x,a1,i2.2)
               eps = 0
               epe = 0
               ep = 0
               n = 0
               do k = 1, 3
                  stats(k) = 0.d0
               end do

               do while ( .not.done )
                   ep = ep + 1
                   if( sat_flgc(j,i,ep).ne.0 .and.
     .                 sat_sgc(j,i,ep).lt. 10.e-9 ) then
*                     This measurement is OK, so accumulate the differences
*                     from the good networks
                     if( eps.eq.0 ) eps = ep
                     do k = 1, num_net
                         if( sat_OK(j,k) ) then
                            dc = sat_clk(j,i,ep) - sat_clk(j,k,ep)
                            W = 1.d0/(sat_sgc(j,i,ep)**2+
     .                                sat_sgc(j,k,ep)**2)
                            n = n + 1
                            stats(1) = stats(1)+W
                            stats(2) = stats(2)+dc*W
                            stats(3) = stats(3)+dc**2*W
                            epe = ep
                         end if
                      end do
                   else
*                     Ok, we have reached a break in the clock.  
*                     Determine the offset to be applied if we have enough
*                     data
                      if ( n.gt.1 ) then 
                         dcm = stats(2)/stats(1)
                         rms  = sqrt((stats(3)-dcm*stats(2))/
     .                  stats(1))
                         sig = sqrt(rms**2/(n-1))
                         write(*,250) i,offcon(prn_list(j)),
     .                                 mod(prn_list(j),100), eps, epe, 
     .                                dcm*1.d9, sig*1.d9, rms*1.d9, n
 250                     format('For Net ',i2,1x,a1,i2.2,' Epochs ',i4,
     .                          ' to ',i4,' Mean offset ',F8.3,' +- ',
     .                          F8.3,' RMS ',F8.3,' ns, Num ',i5)
                         do k = eps, epe
                            sat_clk(j,i,k) = sat_clk(j,i,k)-dcm
                         end do
                         n = 0
                         do k = 1, 3
                             stats(k) = 0.d0
                         end do
                         eps = 0
                      end if
                   end if
*                  See if we are done yet
                   if( ep.eq.num_ep ) done = .true.
               end do
*              Finally see if we have some data left over
               if ( n.gt.1 ) then 
                  dcm = stats(2)/stats(1)
                  rms  = sqrt((stats(3)-dcm*stats(2))/
     .                   stats(1))
                  sig = sqrt(rms**2/(n-1))
                  write(*,260) i,offcon(prn_list(j)),
     .                        mod(prn_list(j),100), eps, epe,  
     .                        dcm*1.d9, sig*1.d9, rms*1.d9, n
 260              format('End Net ',i2,1x,a1,i2.2,' Epochs ',i4,
     .                   ' to ',i4,' Mean offset ',F8.3,' +- ',
     .                    F8.3,' RMS ',F8.3,' ns, Num ',i5)
                  do k = eps, epe
                     sat_clk(j,i,k) = sat_clk(j,i,k)-dcm
                  end do
               end if
            end if       ! Sat not OK
         end do          ! Looping over satellites
      end do             ! Looping over networks

* MOD TAH 031130: Now check those sites and satellites that have
*     no contiguous data at any site
      write(*,'(a)') 'OFFSET CHECKING PROCESS: Non-contiguos data'
*
*     Loop over sites first
      do j = 1, num_site
*
*        See if non-contiguous on all networks
         done = .false.
         m = 0
         do i = 1, num_net
            if( rcv_avail(j,i) ) m = m + 1
         enddo
         do i = 1, num_net
            if( rcv_OK(j,i) ) done = .true.
         end do
*        See if we need adjustment.  Only do if station appears in
*        more than one network.
         if( .not.done .and. m.gt.1 ) then
             write(*,320) site_codes(j)
 320         format('Data not contiguous at Site ',a4)
*
*            Now get the intervals of contiguous data
             eps = 0
             epe = 0
             ep = 0
             n = 0
             done = .false.
             do while ( .not.done )
                ep = ep + 1
*               See if any nets have data at this time
                OK = .false.
                do i = 1, num_net
                   if ( rcv_flgc(j,i,ep).eq.1 .and.
     .                  rcv_sgc(j,i,ep).lt. 10.e-9 ) OK = .true.
                end do
*               OK, good data found.  See if start or end of block
                if ( n.eq.0 .and. OK ) then
*                   No data so far and good found
                    eps = ep
                    n = 1
                endif
*               If we have reached the end of data, set OK to
*               false to force final calculation
                if( ep .eq. num_ep ) then
                    done = .true.
                    if( OK ) epe = ep
                    OK = .false.
                endif
                
                if ( n.gt.0 ) then  ! Already in block
                    if( OK ) then
                        n = n + 1
                        epe = ep
                    elseif ( n.gt.0 ) then  ! OK, we have come to the
*                                   !  end of good block.  Calculate 
*                                   ! offsets between networks
*                       Compute average offsets between the networks
*                       relative to network 1 and then apply
                        do i = 2, num_net
                           m = 0
                           do k = 1,3
                              stats(k) = 0.d0
                           end do
                           do k = eps, epe
                              if(rcv_flgc(j,1,k).eq.1 .and.
     .                           rcv_sgc(j,1,k).lt. 10.e-9 .and. 
     .                           rcv_flgc(j,i,k).eq.1 .and.
     .                           rcv_sgc(j,i,k).lt. 10.e-9  ) then
*                                 Data at site from Nets 1 and Net i OK
                                  dc = rcv_clk(j,i,k) - rcv_clk(j,1,k)
                                  W = 1.d0/(rcv_sgc(j,i,k)**2+
     .                                      rcv_sgc(j,1,k)**2)
                                  m = m + 1
                                  stats(1) = stats(1)+W
                                  stats(2) = stats(2)+dc*W
                                  stats(3) = stats(3)+dc**2*W
                              endif
                           end do
                           if( m.gt.1 ) then 
                              dcm = stats(2)/stats(1)
                              rms  = sqrt((stats(3)-dcm*stats(2))/
     .                                stats(1))
                              sig = sqrt(rms**2/(m-1))
                              write(*,350) i,site_codes(j), eps, epe, 
     .                                dcm*1.d9, sig*1.d9, rms*1.d9, m
 350                          format('For Net ',i2,' ',A4,' Epochs ',i4,
     .                            ' to ',i4,' Mean offset ',F8.3,' +- ',
     .                            F8.3,' RMS ',F8.3,' ns, Num ',i5)
                              do k = eps, epe
                                 rcv_clk(j,i,k) = rcv_clk(j,i,k)-dcm
                              end do
                           else
                              write(*,355) i, site_codes(j), eps,epe
 355                          format('For Net ',i2,' ', A4,' Epochs ',
     .                            i4,' to ',i4,
     .                               ' No overlap with Net #1')
                           endif
                        enddo
                        n = 0
                    endif
                endif
             enddo
         endif
      end do 

****  Now do the satellites                             
*
      do j = 1, num_sat
*
*        See if non-contiguous on all networks
         done = .false.
         do i = 1, num_net
            if( sat_OK(j,i) ) done = .true.
         end do
*        See if we need adjustment
         if( .not.done ) then
             write(*,420) offcon(prn_list(j)),mod(prn_list(j),100)
 420         format('Data not contiguous Satellite ',a1,I2.2)
*
*            Now get the intervals of contiguous data
             eps = 0
             epe = 0
             ep = 0
             n = 0
             done = .false.
             do while ( .not.done )
                ep = ep + 1
*               See if any nets have data at this time
                OK = .false.
                do i = 1, num_net
                   if ( sat_flgc(j,i,ep).ne.0 .and.
     .                  sat_sgc(j,i,ep).lt. 10.e-9 ) OK = .true.
                end do
*               OK, good data found.  See if start or end of block
                if ( n.eq.0 .and. OK ) then
*                   No data so far and good found
                    eps = ep
                    n = 1
                endif
*               If we have reached the end of data, set OK to
*               false to force final calculation
                if( ep .eq. num_ep ) then
                    done = .true.
                    if( OK ) epe = ep
                    OK = .false.
                endif
                
                if ( n.gt.0 ) then  ! Already in block
                    if( OK ) then
                        n = n + 1
                        epe = ep
                    elseif ( n.gt.0 ) then  ! OK, we have come to the
*                                   !  end of good block.  Calculate 
*                                   ! offsets between networks
*                       Compute average offsets between the networks
*                       relative to network 1 and then apply
                        do i = 2, num_net
                           m = 0
                           do k = 1,3
                              stats(k) = 0.d0
                           end do
                           do k = eps, epe
                              if(sat_flgc(j,1,k).ne.0 .and.
     .                           sat_sgc(j,1,k).lt. 10.e-9 .and. 
     .                           sat_flgc(j,i,k).ne.0 .and.
     .                           sat_sgc(j,i,k).lt. 10.e-9  ) then
*                                 Data at site from Nets 1 and Net i OK
                                  dc = sat_clk(j,i,k) - sat_clk(j,1,k)
                                  W = 1.d0/(sat_sgc(j,i,k)**2+
     .                                      sat_sgc(j,1,k)**2)
                                  m = m + 1
                                  stats(1) = stats(1)+W
                                  stats(2) = stats(2)+dc*W
                                  stats(3) = stats(3)+dc**2*W
                              endif
                           end do
                           if( m.gt.1 ) then 
                              dcm = stats(2)/stats(1)
                              rms  = sqrt((stats(3)-dcm*stats(2))/
     .                       stats(1))
                              sig = sqrt(rms**2/(m-1))
                              write(*,450) i,offcon(prn_list(j)), 
     .                                mod(prn_list(j),100), eps, epe, 
     .                                dcm*1.d9, sig*1.d9, rms*1.d9, m
 450                          format('For Net ',i2,1x,a1,I2.2,' Epochs ',
     .                            i4,' to ',i4,' Mean offset ',F8.3,
     .                            ' +- ',F8.3,' RMS ',F8.3,
     .                            ' ns, Num ',i5)
                              do k = eps, epe
                                 sat_clk(j,i,k) = sat_clk(j,i,k)-dcm
                              end do
                           else
                              write(*,455) i,offcon(prn_list(j)), 
     .                                mod(prn_list(j),100), eps,epe
 455                          format('For Net ',i2,1x,a1,I2.2,' Epochs ',
     .                            i4,' to ',i4,
     .                               ' No overlap with Net #1')
                           endif
                        enddo 
                        n = 0
                    endif
                endif
             enddo
         endif
      end do 

* MOD TAH 200624: After iteration 2, offset and align clocks
      if( it.lt.2 ) RETURN

      write(*,'(a,1x,I2)') 'Satellite net align: Iter ',it
*     Now offset the satellites: Loop over each satellite and
*     get offset to Network 1.
      do j = 1, num_sat
         do i = 2, num_net

*           Clear statistics
            do k = 1, 3
               stats(k) = 0.d0
            end do
            n = 0 
            do ep = 1, num_ep

*               See data good in both networks
                if( sat_flgc(j,1,ep).eq. 1 .and.
     .              sat_flgc(j,i,ep).eq. 1 ) then 
*                   This measurement is OK, so accumulate the differences
*                   from the good networks
                    dc = sat_clk(j,i,ep) - sat_clk(j,1,ep)
                    W = 1.d0/(sat_sgc(j,i,ep)**2+
     .                        sat_sgc(j,1,ep)**2)
                    stats(1) = stats(1)+W
                    stats(2) = stats(2)+dc*W
                    stats(3) = stats(3)+dc**2*W
                    n = n + 1
                end if
            end do
*           Finally see if we have some data left over
            if ( n.gt.1 ) then 
               dcm = stats(2)/stats(1)
               rms  = sqrt((stats(3)-dcm*stats(2))/
     .                stats(1))
               sig = sqrt(rms**2/(n-1))
               write(*,520) i,offcon(prn_list(j)), mod(prn_list(j),100),
     .                     dcm*1.d9, sig*1.d9, rms*1.d9, n
 520           format('SVS Net ',i2,1x,a1,i2.2,
     .               ' Mean offset to Net #1',F8.3,' +- ',
     .                 F8.3,' RMS ',F8.3,' ns, Num ',i5)
               do ep = 1, num_ep
                  sat_clk(j,i,ep) = sat_clk(j,i,ep)-dcm
               end do
            end if
         end do          ! Looping over satellites
      end do             ! Looping over networks
                              
      end

CTITLE GET_DOMES

      subroutine get_domes( unitc )

      implicit none
      
*     Routine to update the domes information based on the 
*     four character code of the site

      include '../includes/kalman_param.h'
      include 'merge_igs_clk.h'
 
* unitc  -- Comments file unit numbers

      integer*4 unitc
      
* LOCAL VARIABLES      
*   ierr        - IOSTAT error
*   trimlen     - Length of string
 
      integer*4 ierr, trimlen, ns, indx
 
*   start_found - Set true when start line found
*   stop_found  - Set true when stop line found
 
      logical start_found, stop_found
 
*   line        - Line read from comments
 
      character*128  line
            
*   domes       - DOMES number for site (9 characters needed)
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]
*   globk_extent - Last four characters of globk name

      character*12 domes
      character*4 type, pt, globk_extent, code
      character*8 full_name
      character*128 snx_comfile   ! Name of sinex comments file that contains
                                  ! domes numbers
      
****  Check to see if comment unit
      snx_comfile = 'head.snx'
      open(unitc, file=snx_comfile, iostat=ierr, status='old')
      
      if( ierr.ne.0 ) then

*         Try to open version in HELP_DIR (use line for name)
          call getenv('HELP_DIR',line)
          line(trimlen(line)+1:) = '/' // snx_comfile
          open(unitc, file=line, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open',line,1,'GET_DOMES')
          if ( ierr.ne.0 ) RETURN
      end if

***** See if +DOMES can be found

      line = ' '
      ierr = 0
      start_found = .false.
      stop_found = .false.
      do while ( ierr.eq.0 .and. .not.stop_found )
 
          read(unitc, '(a)', iostat=ierr) line
          if( index(line,'-DOMES').eq.1 ) stop_found = .true.
          if( ierr.eq.0 .and. .not.stop_found .and. start_found .and.
     .        trimlen(line).gt.0 ) then
              
              read(line, 120) code, pt, domes, type, globk_extent
 120          format(1x,a4,1x,a2,1x,a9,1x,a1,1x,a4)
              full_name = code(1:4) // globk_extent
              indx = 1
              call get_cmd(code(1:4), site_codes, num_site,
     .               ns, indx)
              if( ns.gt.0 ) then
                 site_names(ns) = domes(1:9)
              end if
          end if
*         Check for start afterwards so that the start line
*         is not echoed.
          if( index(line,'+DOMES').eq.1 ) start_found = .true.
      end do

****  Thats all
      return 
      end

CTITLE RESET_REFS

      subroutine reset_refs

      implicit none

*     Routine to compute the average residual at the reference clocks
*     and remove this value (but epoch) from all clocks.  This will
*     flatten the reference clocks

      include 'merge_igs_clk.h'

* LOCAL Variables
      integer*4 i,j,k   ! Loop counters
     .,         nr      ! Site number of reference site
     .,         n       ! Number of values in offset by epoch
     .,         ep      ! Epoch number

      real*8 stats(3)   ! Statistics summation
     .,      mean, rms  ! Mean and RMS of reference residuals at each epoch
     .,      res, W     ! Residual at reference and weight


***** OK: Based on the fits to the reference clocks, compute
*     average residal at each epoch
      do ep = 1, num_ep
*        Compute average residual at reference sites
         do k = 1, 3
            stats(k) = 0.d0
         end do
         n = 0
         do i = 1, num_ref
            nr = ref_sites(i)
            do j = 1, num_net
               if( rcv_flgc(nr,j,ep).eq. 1 .and.
     .             rcv_sgc(nr,j,ep).lt.10.d-9 ) then
*                  OK: Value is good.  Compute Residuals and add
*                  to stats
                   res = rcv_clk(nr,j,ep) - (rcv_fit(2,nr,j)+
     .                                       rcv_fit(3,nr,j)*(ep-1440))
                   
                   W = 1.d0/rcv_sgc(nr,j,ep)**2
                   n = n + 1
                   stats(1) = stats(1) + W
                   stats(2) = stats(2) + res*W
                   stats(3) = stats(3) + res**2*w
               end if
            end do
         end do

*        Complete the stats
         if( n.gt.1 ) then
             mean = stats(2)/stats(1)
             rms = sqrt((stats(3)-mean*stats(2))/stats(1))
C             write(*,120) ep, n, mean*1.d9, rms*1.d9
C120         format('RR Ep ',i4,' Num ',i4,' Mean ',F10.2, 
C     .              ' ns; RMS ',F8.3,' ns')
*            Now remove this mean from all clocks
             do i = 1, num_site
                do j = 1, num_net
                   if( rcv_flgc(i,j,ep).eq.1 ) then
                       rcv_clk(i,j,ep) =  rcv_clk(i,j,ep)-mean
                   end if
                end do
             end do
         
             do i = 1, num_sat
                do j = 1, num_net
                   if( sat_flgc(i,j,ep).eq.1 ) then
                       sat_clk(i,j,ep) =  sat_clk(i,j,ep)-mean
                   end if
                end do
             end do
         else
             write(*,140) ep
 140         format('RR Ep ',i4,' No data')
         endif
      end do

***** Thats all
      return 
      end

CTITLE EDIT_CLK

      subroutine edit_clk(it)

      implicit none

*     Routine to use the linear fits to the good clocks to remove
*     any outlier epochs. 


      include 'merge_igs_clk.h'

* PASSED VARIABLES
      integer*4 it  ! Iteration
 
* LOCAL VARIABLES

      integer*4 ep    ! Epoch number
      integer*4 nr    ! Receiver number
     .,         nt    ! Network number
     .,         i, j     ! Loop counter
     .,         njmp     ! NUmber of points in jump


      real*8 res ! residual seconds
     .,      cres   ! Check residual
     .,      avres  ! Average of residuals after a jump
     .,      offset ! Offset to be applied for jumps

      real*8 tol      ! Tolerance of fit of data 

      logical jump    ! Set true for a clock jump


*     Now loop reference sites
      write(*,120) it
 120  format('EDITING SITE CLOCKS: Iteration ',i2)

      do i = 1, num_site
*         nr = ref_sites(i)
         nr = i
         do nt = 1, num_net
*           Get the tolerance for editing this site
            if( rcv_fit(1,nr,nt).gt.0 ) then
                tol = 5*rcv_fit(4,nr,nt)
            else
                tol = 100.d-9
            end if
            if( tol.lt.0.25d-9 ) tol = 0.25d-9
            if( tol.lt.10.d-9 ) then
               offset = 0.d0 
               do ep = 1, num_ep
                  res = rcv_clk(nr,nt,ep) - (rcv_fit(2,nr,nt)+
     .                          rcv_fit(3,nr,nt)*(ep-1440)) - offset
                  if( abs(res).gt.tol .and. 
     .                rcv_flgc(nr,nt,ep).eq. 1 ) then

*                     Check the next epochs to see if this is a jump
*                     in the clock.
                      jump = .false.
                      avres = 0
                      njmp = 0
                      do j = ep+1, min(num_ep,ep+10)
                         cres = rcv_clk(nr,nt,j) - (rcv_fit(2,nr,nt)+
     .                            rcv_fit(3,nr,nt)*(j-1440)) - offset
                         if( abs(res-cres).lt.tol .and.
     .                       rcv_flgc(nr,nt,j).eq. 1 ) then
                             avres = avres + cres
                             njmp = njmp + 1
                         endif
                      end do

*                     See how many points in jump and if we should offset
                      if( njmp.gt.4 ) then
                          offset = offset + avres/njmp
                          jump = .true.
                          write(*,220) ep, site_codes(nr), nt, 
     .                                 offset*1.d9, it
 220                      format('CLK JUMP Epoch ',i4,1x,a4,' Net ',i2,
     .                           ' Jump ',F8.2,' ns, Iter ',i2)
                      endif

                      if( .not. jump ) then
                          write(*,250) ep, site_codes(nr), nt, res*1.d9,
     .                                 tol*1.d9, it
 250                      format('DELETING Epoch ',i4,1x,a4,' Net ',i2,
     .                           ' Res  ',F8.2,' ns, Tol ',f8.2,
     .                           ' ns, Iter ',I2)
                          rcv_flgc(nr,nt,ep) = 3
                      endif
                  end if
                  if( abs(res).lt.tol .and. 
     .                rcv_flgc(nr,nt,ep).eq. 3 ) then
                      write(*,260) ep, site_codes(nr), nt, res*1.d9,
     .                             tol*1.d9, it
 260                  format('RESTORING Epoch ',i4,1x,a4,' Net ',i2,
     .                       ' Res ',F8.2,' ns, Tol ',f8.2,
     .                       ' ns, Iter ',I2)
                      rcv_flgc(nr,nt,ep) = 1
                  end if
               end do
            end if
          end do
      end do


****  Thats all
      return
      end

CITTLE COPY_FIN

      subroutine copy_fin

      implicit none

*     Routine to copy values for get final values (used when num_net.eq.1

      include 'merge_igs_clk.h'

* LOCAL VARIABLES
      integer*4 ep  ! Epoch counter
     .,         i,j,k   ! Loop counters
     .,         num_fin  ! Final number of reference clocks
     .,         num      ! Number counter

* MOD TAH 200528: Add GNSS separation
      integer*4 prn        ! PRN with GNSS offset removed
      character*1 offcon   ! function gets GNSS type from PRN
     .,         gnss       ! GNSS type


*     OK Final adjustment
      num_fin = num_ref
      best_ref = ref_sites(1)

      write(*,320) num_fin
 320  format('COPY WITH ',i2,' CLOCKS')
c      debug = .true.
      do ep = 1, num_ep
         do j = 1, num_site
            rcv_fin(j,ep) = rcv_clk(j,1,ep) 
            rcv_sig(j,ep) = rcv_sgc(j,1,ep)
            rcv_flag(j,ep)= rcv_flgc(j,1,ep)
         end do
         do j = 1, num_sat
            sat_fin(j,ep) = sat_clk(j,1,ep)
            sat_sig(j,ep) = sat_sgc(j,1,ep)
            sat_flag(j,ep)= sat_flgc(j,1,ep)
         end do
      end do

      end


