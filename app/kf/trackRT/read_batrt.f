ctitle read_bat

      subroutine read_batrt( rerun, pass_debug )

      implicit none

*     Routine to read through the batch file and save the
*     information about the parameters and the Markov elements 
*     Each set of runs of the kalman filter
*     ends with a -1 starting in column 1, except for the last
*     run which ends in a -2 denoting end of file.
*     For each additional kalman run only the parameters to be changed
*     need to be specified.
*     MOD G. CHEN 971028 change buffer to 126 
*                        add read_neu
* MODded from trackRT

      include 'trackRT_cmds.h'
      include 'trackRT.h'

* PASSED 
      integer*4 rerun  ! Set 1 when this is a re-run and zero when
                       ! first run (from trackTcomm.cpp)
      integer*4 pass_debug(10)  ! Debug array (passed so that test_trackRT
                       ! can use debug options)

*
*     Local Variables
*     ---------
*     buffer   -- 126 character buffer for reading the Batch file
*     site(max_site) -- casefolded version of site names
*     add_site -- Site name to be added when none passed in command line
*     end_of_run -- logical variable which tells when a -1 is encountered
*    
      character buffer*256, site(max_site)*4, add_site*4
      integer*4 itype
 
*   i,j,k   - Loop counters.
*   ierr    - IOSTAT error during reading.
 
      integer*4 i, j, k, ierr, indx,jndx, icon, sit_num, trimlen
      integer*4 jel    ! Function to return baseline number
      integer*4 itest  ! Test of see if command in sub_set mode.
  
      real*8 bl        ! Baseline length from coordinates

 
c
      logical end_of_run, sub_sets
 

*     eof_file  -- Indicates the end of the batch file.  This variable
*            is set true if a -2 is encountered or the physical end
*            of file is reached.

****  See if this is an update call or first call
      if( rerun.eq.0 ) then 

*        Initialize defaults for program
         call init_trackRT

c....    Initialize site_apr input
         do i=1, max_site
            read_site_apr(i) = .false. 
         enddo
            
      
c....    If neccessary open the batch file
         open(100,file=bat_file,status='old',iostat=ierr)
         call report_error('IOSTAT',ierr,'open',bat_file,1,'TrackRT')

****     Set the default names of files if the -n option was used to 
*        specify the prt_root name
         if( len_trim(prt_root).gt.0 ) then
             sum_file = prt_root
             posit_root = prt_root
         end if

      endif

c.... now start reading file
 
      end_of_run = .false.
      sub_sets = .false.
      icon = 1
      
      do while (.not. end_of_run)
c
*        We may have already read line when icon = 0 ; so only
*        read if not zero.
         if( icon.ne.0 ) read(100,'(a)',iostat=ierr) buffer
         if ( ierr.ne.0 ) end_of_run = .true.
c
c....    see if comment line or end of file or run marker
*                                       ! not a comment line
         if( buffer(1:1).eq.' ' .and. ierr.eq.0 ) then
c
c....       Look for a control  command
            if(.not.sub_sets) then

               indx = 0
               call get_cmd(buffer,commands, max_xk_cmds,itype,indx)
               if( itype.eq.-2 ) stop
               if( itype.eq.-1 ) then
                  write(*,200) buffer(1:trimlen(buffer))
 200              format('TRACKRT: Warning Unknown command: ',a)
               end if 
c....          if itype > 0 then this is a control card
*                                 ! find out type of control(1-4)
               if( itype.gt.0.and.itype.lt.999999 ) then
                   call assign_conrt(itype,icon)
		
c                  If icon = 1, 2 or 3, search the following lines.
                   if( icon.le.3.and.icon.ge.1 ) sub_sets = .true.

c....              if icon = 4, the rest of parameters are on this line.
c                  then decode the rest of this line
                   if( icon.eq.4 )  call get_missrt(buffer,itype,indx)

               end if

*              readin the subset commands (site or GPS sat. related)            
            else
          
*              check the next line until a blank string is found
*                                       ! if reach a blank string 
*              See if valid command given
               indx = 0
               call get_cmd(buffer,commands, max_xk_cmds,itest,indx)
               if( itest.gt.0 .and. itest.ne.999999) then
*                  valid command so do not try to interprett as site command
                   itype = -3
                   icon = 0
               end if

               if( itype.eq.-3 ) sub_sets = .false.
               
*                                     ! Observation files expected
               if( icon.eq.1) then
C                  call get_obsfile(buffer, itype)

*                                     ! if this is a blank string
                  write(*,320) trim(buffer)
 320              format('**WARNING** No obs_file commands in trackRT.',
     .                   1x,a,' not valid')
c
               endif
*                                     ! site-related files expected
               if( icon.eq.2) then
                  indx = 0
                
c ....            change commands to the upper case
                  do i = 1,num_site
                      site(i) = site_names(i)
                      call casefold(site(i))
                  enddo
                  call get_cmd(buffer,site,num_site, sit_num,indx)

****              See if we need to build a list of sites because no
*                 site names have been passed in command line (num_prc is
*                 number of sites given in the command line).
                  if( num_prc.le.1 .and. sit_num.eq.-1 ) then
*                     Need to add this site
                      num_site = num_site + 1
                      if( num_site.gt.max_site ) then
                          call report_stat('FATAL','TRACKRT',
     .                       'read_batrt',' ',
     .                       'Too many sites: Max allowed',
     .                        max_site)
                      endif
                      jndx = 0
                      call GetWord(buffer, add_site, jndx) 
                      site_names(num_site) = add_site
                      call casefold(add_site)
                      site(num_site) = add_site
                      sit_num = num_site
                  endif
                
*                                     ! if this is not a blank string
                  if( sit_num.eq.-3) sub_sets = .false.
*                                                   ! get site information
                  if( sit_num. gt.0 ) then
c
c....             see if sit_num=999999, if it is do all sites
*                                                   ! do all sites
                      if ( sit_num.eq.999999 ) then
                         if( num_prc.gt.0 ) then
                            do i = 1,num_site
                               jndx = indx
                               call get_sitert(buffer,i,itype,jndx)
                            end do
                         else  ! No sites yet, so fill max allowed
                            do i = 1,max_site
                               jndx = indx
                               call get_sitert(buffer,i,itype,jndx)
                           end do
                         end if
*                                                   ! do just this site
                      else
                         call get_sitert(buffer,sit_num,itype,indx)
                      end if
                  end if
               endif

c
            end if
c---
         endif
      end do

****  Check the number of kinematic sites we have
      k = 0
      do i = 1, num_site
         if( site_type(i).eq.1 ) k = k + 1
      end do

      if( k.gt.max_kine ) then
          call report_stat('FATAL','TRACKRT','read_batrt',' ',
     .                     'Too many kinematic sites: Max allowed',
     .                     max_kine)
      end if

* MOD TAH 090121: Update the process noise values for the time_unit
      if( tu_to_ep.ne.1.0 ) then
          do i = 1,num_site
             do j = 1,3
                 mar_site(j,i) = mar_site(j,i)/tu_to_ep
             enddo
             mar_atm(i)  = mar_atm(i)/tu_to_ep
          end do
      end if
*
****  Output the baseline lengths based on apriori coordimate
      write(*,610)
 610  format('Initial site positions and baseline lengths')
      do i = 1, num_site 
         write(*,620) i, site_names(i), site_apr(:,i), 
     .                site_vel(:,i), site_ep(i)
 620     format('Sites ',i3,1x,a,1x,3(F15.4,1x),3(F8.4,1x),F15.3)
*        Could add a test here to see if aprioris have been read.
      end do
      write(*,'(a,1x,a4)') 'Baseline lengths from',site_names(1)
      do i = 1,num_site-1
         do j = 2, num_site
            baselens(jel(i,j)) = sqrt((site_apr(1,j)-site_apr(1,i))**2
     .            +(site_apr(2,j)-site_apr(2,i))**2 
     .            +(site_apr(3,j)-site_apr(3,i))**2)
         end do
      end do
      do i = 2,num_site
         bl = baselens(jel(i,1))
         write(*,640) site_names(i),site_names(1), bl
 640     format('Distance ',a4,'-',a4,1x,F12.3,' m')
      end do

****  See if we will output reference site results
      do i = 1,num_site
         if( ref_code.eq.site_names(i) .and.
     .       mar_atm(i).eq. 0.0 ) then
             write_pos(i) = .false.
         else 
             write_pos(i) = .true.
         endif
      enddo

*     Save number of processed sites 
      if( num_prc.eq.0 ) num_prc = num_site

***** Pass back the debug array
      pass_debug = debug

      return
      end

* Routine below not currently used (maybe in future versions)

CTITLE SUBSRDWRT
       
      subroutine subsrdwrt(buffer, runday, runweek, runstr)

      implicit none

****  Routine to scan buffer for <day> or <week> and replace with the
*     runday and runweek strings

      character*(*) buffer, runday, runweek, runstr(10)

* LOCAL
      character*256 new_buffer
      integer*4 trimlen, indx, lenrd, lenrw, i

      character*5 strlabl(10), strlabu(10)

      data strlabl / '<s01>', '<s02>', '<s03>', '<s04>', '<s05>',
     .               '<s06>', '<s07>', '<s08>', '<s09>', '<s10>' /
      data strlabu / '<S01>', '<S02>', '<S03>', '<S04>', '<S05>',
     .               '<S06>', '<S07>', '<S08>', '<S09>', '<S10>' /


****  See if we can find <day> or <DAY> in string
      if( trimlen(buffer).eq.0 ) RETURN
      new_buffer = buffer
      indx = index(buffer,'<day>')
      if( indx.eq.0 ) indx = index(buffer,'<DAY>')
      if( indx.gt.0 .and. buffer(1:1).eq.' ' ) then
          lenrd = trimlen(runday)
          if( lenrd.eq.0 ) then
              write(*,120) buffer(1:trimlen(buffer))
 120          format('**FATAL** Track command ',a,/,
     .               'used <day> option but not passed with -d option')
              stop 'TRACKRT: -d option not in runstring'
          endif

*         OK now "splice" in the new string
          new_buffer = buffer
          new_buffer(indx:) = runday(1:lenrd) // buffer(indx+5:)
          buffer = new_buffer
      endif
         
      indx = index(new_buffer,'<week>')
      if( indx.eq.0 ) indx = index(buffer,'<WEEK>')
      if( indx.gt.0 .and. buffer(1:1).eq.' ' ) then
          lenrw = trimlen(runweek)
          if( lenrw.eq.0 ) then
              write(*,140) buffer(1:trimlen(buffer))
 140          format('**FATAL** Track command ',a,/,
     .               'used <week> option but not passed with -w option')
              stop 'TRACKRT: -w option not in runstring'
          endif

*         OK now "splice" in the new string
          new_buffer(indx:) = runweek(1:lenrw) // buffer(indx+6:)
          buffer = new_buffer
      endif
      i = 1
      do while ( i.le.10 .and. trimlen(runstr(i)).gt.0) 
         if( trimlen(runstr(i)).gt.0 ) then
            indx = index(buffer,strlabl(i))
            if( indx.eq.0 ) indx = index(buffer,strlabu(i))
            if( indx.gt. 0 ) then
                lenrw = trimlen(runstr(i))
*               OK now "splice" in the new string
                if( runstr(i)(1:lenrw).eq.'space' ) then
                    new_buffer(indx:) = ' ' // buffer(indx+5:)
                else
                     new_buffer(indx:) = runstr(i)(1:lenrw) // 
     .                                   buffer(indx+5:)
                endif
                buffer = new_buffer
            else 
                i = i + 1   ! Only move to next one after all cases have been found
            endif
         endif
      end do

           
****  Thats all
      close(100)

      return
      end

  
