ctitle read_bat
      subroutine read_bat

      implicit none

*     G. CHEN 
*
*     Routine to read through the batch file and save the
*     information about the parameters and the Markov elements 
*     Each set of runs of the kalman filter
*     ends with a -1 starting in column 1, except for the last
*     run which ends in a -2 denoting end of file.
*     For each additional kalman run only the parameters to be changed
*     need to be specified.
*     MOD G. CHEN 971028 change buffer to 126 
*                        add read_neu

      include 'makexk_cmds_bd.h'

      include '../includes/xfile_def.h'
      include 'track_com.h'

*
*     Local Variables
*     ---------
*     buffer   -- 126 character buffer for reading the Batch file
*     batf_open -- logical variable which indicates that the batch file
*            is open
*     end_of_run -- logical variable which tells when a -1 is encountered
*    
      character*256 buffer
      integer*4 itype
 
*   i,j,k   - Loop counters.
*   ierr    - IOSTAT error during reading.
 
      integer*4 i, j, k, ierr, indx,indx1, icon, sit_num, trimlen
      real*8 bl  ! Baseline length (m)
 
c
      logical bat_open, end_of_run, sub_sets
 

      data bat_open / .false. /


*     eof_file  -- Indicates the end of the batch file.  This variable
*            is set true if a -2 is encountered or the physical end
*            of file is reached.

c.... Initialize site_apr input
      do i=1, max_site
         read_site_apr(i) = .false. 
      enddo
      
c.... If neccessary open the batch file
      if( .not. bat_open ) then
           open(100,file=bat_file,status='old',iostat=ierr)
           call report_error('IOSTAT',ierr,'open',bat_file,0,'MAKEXK')
           if( ierr.ne.0 ) then
             stop
           else
             bat_open = .true.
           end if
      end if
c....    end of file opening

c.... now start reading file
 
      end_of_run = .false.
      sub_sets = .false.
      
      do while (.not. end_of_run)
c
         read(100,'(a)', end=1000) buffer

*        Replace any variable strings
         call subsrdw(buffer, runday, runweek, runstr)

c        call casefold(buffer)
c
c....    see if comment line or end of file or run marker
*                                       ! not a comment line
         if( buffer(1:1).eq.' ' ) then
c
c....      Look for a control  command
           if(.not.sub_sets) then

             indx = 0
             call get_cmd(buffer,commands, max_xk_cmds,itype,indx)
             if( itype.eq.-2 ) stop
             if( itype.eq.-1 ) then
                write(*,200) buffer(1:trimlen(buffer))
 200            format('TRACK: Warning Unknown command: ',a)
             end if 
c             write(*,*) "--",itype,"--",buffer
c
c....        if itype > 0 then this is a control card
*                                 ! find out type of control(1-4)
             if( itype.gt.0.and.itype.lt.999999 ) then
               call assign_con(itype,icon)
		
c              If icon = 1, 2 or 3, search the following lines.
               if( icon.le.3.and.icon.ge.1 )  sub_sets = .true.

c....          if icon = 4, the rest of parameters are on this line.
c              then decode the rest of this line
               if( icon.eq.4 )  call get_miss(buffer,itype,indx)

             end if

*              readin the subset commands (site or GPS sat. related)            
           else
          
*            check the next line until a blank string is found
*                                     ! if reach a blank string 
             if( itype.eq.-3 ) then
                sub_sets = .false.
c ....          change commands to the upper case
                do i = 1,num_site
                   site_NUC(i) = site_names(i)
                   call casefold(site_NUC(i))
                enddo
             end if
               
*                                     ! Observation files expected
             if( icon.eq.1) then
                call get_obsfile(buffer, itype)

*                                     ! if this is a blank string
                if( itype.eq.-3) then
                    sub_sets = .false.
c ....              change commands to the upper case
                    do i = 1,num_site
                       site_NUC(i) = site_names(i)
                       call casefold(site_NUC(i))
                    enddo
                end if

c
             endif
                  
*                                     ! site-related files expected
             if( icon.eq.2) then
                indx = 0
                
                call get_cmd(buffer,site_NUC,num_site,
     .                                  sit_num,indx)
                
*                                     ! if this is not a blank string
                if( sit_num.eq.-3) sub_sets = .false.
*                                                   ! get site information
                if( sit_num. gt.0 ) then
c
c....              see if sit_num=999999, if it is do all sites
*                                                   ! do all sites
*                  itype = 6 is TIMEDEP_PROCNS where ALL is saved
*                  as site -1
                   if ( sit_num.eq.999999 .and. itype.ne.6 ) then
                       do i = 1,num_site
                         indx1= indx
                         call get_site(buffer,i,itype,indx1)
                      end do
*                                                   ! do just this site
                   else
                      call get_site(buffer,sit_num,itype,indx)
                   end if
c
                end if
             endif

c
            end if
c---
         endif
      end do

c.... link the parameters        
 1000 continue

****  Now update any un-resolved entries
      if( float_type(1:3).eq.'def' ) float_type = search_type

****  Check the number of kinematic sites we have
      k = 0
      do i = 1, num_site
         if( site_type(i).eq.1 ) k = k + 1
      end do
      if( k.gt.max_kine ) then
          call report_stat('FATAL','TRACK','read_bat',' ',
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
          do i = 1, num_tdep
             do j = 1,3
                mar_tdep(j,i) = mar_tdep(j,i)/tu_to_ep
             enddo
          end do
      end if
*
****  Output the baseline lengths based on apriori coordimate
      if( read_site_apr(1) ) then 
         write(*,610)
 610     format('Initial site positions and baseline lengths')
         do i = 1, num_site 
            write(*,620) i, site_names(i), site_apr(:,i) 
C    .,                site_vel(:,i), site_ep(i)
 620        format('Sites ',i3,1x,a,1x,3(F15.4,1x),3(F8.4,1x),F15.3)
*           Could add a test here to see if aprioris have been read.
         end do
         write(*,'(a,1x,a4)') 'Baseline lengths from',site_names(1)
         do i = 2,num_site
            bl = sqrt((site_apr(1,i)-site_apr(1,1))**2
     .            +(site_apr(2,i)-site_apr(2,1))**2 
     .            +(site_apr(3,i)-site_apr(3,1))**2)
            write(*,640) site_names(i),site_names(1), bl
 640        format(a4,'-',a4,1x,F12.3,' m')
         end do
      end if
*
      return
      end

CTITLE SUBSRDW
       
      subroutine subsrdw(buffer, runday, runweek, runstr)

      implicit none

****  Routine to scan buffer for <day> or <week> and replace with the
*     runday and runweek strings

      integer*4 max_snn   ! Max number of <sNN> strings allowed
      parameter ( max_snn = 20 ) 

      character*(*) buffer, runday, runweek, runstr(max_snn)

* LOCAL
      character*256 new_buffer
      integer*4 trimlen, indx, lenrd, lenrw, i

      character*5 strlabl(20), strlabu(20)

      data strlabl / '<s01>', '<s02>', '<s03>', '<s04>', '<s05>',
     .               '<s06>', '<s07>', '<s08>', '<s09>', '<s10>',
     .               '<s11>', '<s12>', '<s13>', '<s14>', '<s15>',
     .               '<s16>', '<s17>', '<s18>', '<s19>', '<s20>' /
      data strlabu / '<S01>', '<S02>', '<S03>', '<S04>', '<S05>',
     .               '<S06>', '<S07>', '<S08>', '<S09>', '<S10>',
     .               '<S11>', '<S12>', '<S13>', '<S14>', '<S15>',
     .               '<S16>', '<S17>', '<S18>', '<S19>', '<S20>' /


****  See if we can find <day> or <DAY> in string
      if( trimlen(buffer).eq.0 ) RETURN
      new_buffer = buffer
      if( buffer(1:1).eq.' ' .or.  buffer(1:1).eq.'<') then
         i = 1
         do while ( i.le.max_snn .and. trimlen(runstr(i)).gt.0) 
            if( trimlen(runstr(i)).gt.0 ) then
               indx = index(buffer,strlabl(i))
               if( indx.eq.0 ) indx = index(buffer,strlabu(i))
               if( indx.gt. 0 ) then
                   lenrw = trimlen(runstr(i))
*                  OK now "splice" in the new string
                   if( runstr(i)(1:lenrw).eq.'space' ) then
                       new_buffer(indx:) = ' ' // buffer(indx+5:)
                   else
                        new_buffer(indx:) = runstr(i)(1:lenrw) // 
     .                                      buffer(indx+5:)
                   endif
                   buffer = new_buffer
               else 
                   i = i + 1   ! Only move to next one after all 
                               ! cases have been found
               endif
            endif
         end do
      end if

      if( buffer(1:1).eq.' ' ) then 
         indx = index(buffer,'<day>')
         if( indx.eq.0 ) indx = index(buffer,'<DAY>')
         if( indx.gt.0 ) then
             lenrd = trimlen(runday)
             if( lenrd.eq.0 ) then
                 write(*,120) buffer(1:trimlen(buffer))
 120             format('**FATAL** Track command ',a,/,
     .              'used <day> option but not passed with -d option')
                 stop 'TRACK: -d option not in runstring'
             endif

*            OK now "splice" in the new string
             new_buffer = buffer
             new_buffer(indx:) = runday(1:lenrd) // buffer(indx+5:)
             buffer = new_buffer
         endif
         
         indx = index(new_buffer,'<week>')
         if( indx.eq.0 ) indx = index(buffer,'<WEEK>')
         if( indx.gt.0 ) then
             lenrw = trimlen(runweek)
             if( lenrw.eq.0 ) then
                 write(*,140) buffer(1:trimlen(buffer))
 140             format('**FATAL** Track command ',a,/,
     .               'used <week> option but not passed with -w option')
                 stop 'TRACK: -w option not in runstring'
             endif

*            OK now "splice" in the new string
             new_buffer(indx:) = runweek(1:lenrw) // buffer(indx+6:)
             buffer = new_buffer
         endif
      endif 

           
****  Thats all
      return
      end

  
