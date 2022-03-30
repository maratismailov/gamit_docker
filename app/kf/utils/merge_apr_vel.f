      program merge_apr_vel

      implicit none 

*     Program to merge an apriori coordinate file with velocity field file
*     by replacing the velocities in the apriori file with those from the
*     velocity file. 

      include 'merge_apr_vel.h'

* LOCAL VARIABLES

      integer*4 ierr, jerr  ! IOSTAT errors
     .,         trimlen     ! Length of line
     .,         js          ! Velocity site number that matches current site
     .,         indx        ! Position in string
     .,         i           ! Loop counter


      real*8 aprs(7)        ! Position, velocity and epoch from
*                           ! apriori file
     .,      mind           ! Mimimum distance from apr site to velocity
                            ! site (m)
     .,      dist           ! Distance bwteen two points
     .,      calcd          ! computes distance bwteen sites (m)

      character*8 site_name ! Name of site from apriori file

      character*256 line    ! Line read from input file


****  Get the runstring for the program and initialize.
      call get_run_mav

****  Read the complete velocity file (apriori is read and written
*     simulataneously so that comments can be passed through)/
      call read_vel_mav 

****  OK now read apriori file
      open(50,file=in_apr_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',in_apr_file,
     .                   1,'input Apriori file')
      open(70,file=out_apr_file,status='unknown',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',out_apr_file,
     .                   1,'output Apriori file open')

*     Write header for new file
      write(70,120) in_apr_file(1:trimlen(in_apr_file)),
     .              in_vel_file(1:trimlen(in_vel_file))
 120  format('* MERGE_APR_VEL apriori file',/,
     .       '* Input apr file : ',a,/,
     .       '* Input vel file : ',a)
      if( zero_uvel ) then
          write(70,140)
 140      format('* Height velocities are forced to zero')
      end if

****  Loop over the apriori file
      do while ( ierr.eq.0 )
         read(50,'(a)',iostat=ierr) line
         if( ierr.eq.0 ) then
             if( trimlen(line).gt.0 .and.
     .           line(1:1).ne.' ' ) then
                 write(70,'(a)',iostat=jerr) line(1:trimlen(line))
             elseif (trimlen(line).gt.0 .and.
     .           index(line,'EXTENDED').gt.0 ) then  ! Just echo extended lines
                 write(70,'(a)',iostat=jerr) line(1:trimlen(line))
             else if( line(1:1).eq.' ' ) then 

*                Get the site name and position
                 num_apr_site = num_apr_site + 1
                 indx = 0
                 call GetWord(line, site_name, indx)
                 read(line(indx:),*,iostat=jerr) aprs
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                             'Main/read apriori')

*  
*                Now see if we can match a name
                 indx = 0
                 js = -1
                 call get_cmd(site_name, vel_name, num_vel_site, 
     .                            js, indx)
                 if( js.gt.0 ) then
                     num_apr_nch = num_apr_nch + 1
                     mind = 0.d0
                 else
*                    See if we should check position
                     mind = 1.d9
                     js = -1
                     do i = 1, num_vel_site
                        dist = calcd(aprs,i)
                        if( dist.lt.mind ) then
                            mind = dist
                            js = i
                        end if
                     end do
                     if( js.gt.0 .and. mind.le.dist_tol ) then
                         num_apr_pch = num_apr_pch + 1
                     else
                         js = -1
                     end if
                 end if

****             We we have a comparison site save now
                 if( js.gt.0 ) then
                     call change_vel(aprs,js)
                 end if

*                Write out the apriori line
                 if( js.gt.0 ) then
                    write(70,220,iostat=jerr) site_name, aprs, 
     .                           vel_name(js), mind/1000.d0
                 else
                    write(70,220,iostat=jerr) site_name, aprs 
                 end if
 220             format(1x,a8,1x,3(F13.4,1x),3(F7.4,1x),f8.3,
     .                  2x,a8,1x,F7.3,' km')
                 call report_error('IOSTAT',ierr,'writ',out_apr_file,1,
     .                             'Main')
             end if
         end if
      end do

      write(*,410) num_apr_site, num_vel_site, num_apr_nch, num_apr_pch
 410  format('* There were ',i4,' Apriori sites and ',i4,
     .       ' Velocity sites',/,
     .       ' Site changed due to names: ',i4,' due to position ',i4)

****  Thats all
      end
              
CTITLE GET_RUN_MAV

      subroutine get_run_mav

      implicit none 

*     Routine to read the runstring for merge_apr_vel

      include 'merge_apr_vel.h'

* LOCAL VARIABLES

      integer*4 jerr  ! IOSTAT errors
     .,         lenr  ! Length of runstring
     .,         arg   ! Argument for runstring number
     .,         rcpar ! Function to read runstring

      character*256 line    ! Line read from input file


****  Initialize the variables in the program
      num_vel_site = 0
      num_apr_site = 0
      num_apr_nch  = 0
      num_apr_pch  = 0

      zero_uvel = .false.
      dist_tol  = -1

****  See what options were passed
      lenr = 1
      arg = 0
      options = 0
      do while ( lenr.gt.0 )
          arg = arg + 1
          lenr = rcpar(arg,line)
          if( lenr.gt.0 ) then
*             See which option
              if( line(1:2).eq.'-f' ) then
                  arg = arg + 1
                  lenr = rcpar(arg,in_apr_file)
                  if( lenr.gt.0 ) call sbit(options,1,1)
              else if ( line(1:2).eq.'-o' ) then
                  arg = arg + 1
                  lenr = rcpar(arg,out_apr_file)
                  if( lenr.gt.0 ) call sbit(options,2,1)
              else if ( line(1:2).eq.'-v' ) then
                  arg = arg + 1
                  lenr = rcpar(arg,in_vel_file)
                  if( lenr.gt.0 ) call sbit(options,3,1)
              else if ( line(1:2).eq.'-z') then
                  zero_uvel = .true.
              else if ( line(1:2).eq.'-m') then
                  arg = arg + 1
                  lenr = rcpar(arg,line)
                  if( lenr.gt.0 ) then
                      read(line,*,iostat=jerr) dist_tol
                      call report_error('IOSTAT',jerr,'decod',
     .                     line,0,'GET_RUN_MAV/DIST_TOL')
                      call sbit(options,5,1)
                  end if
              end if
          end if
       end do

*****  See if all needed options were passed.
       if( iand(options,7).ne.7 ) then
           write(*,220) options
 220       format('MERGE_APR_VEL: Incomplete runstring ',
     .            'options -f -o and -v all needed. Options ',o8)
           call proper_runstring('merge_apr_vel.hlp','merge_apr_vel',1)
       end if

****   Thats all 
       return
       end

CTITLE READ_VEL_MAV

      subroutine read_vel_mav

      implicit none 

*     Routine to read the velocity file and save the velocities.  Duplicate
*     names are replaced with later values.

      include 'merge_apr_vel.h'

* LOCAL VARIABLES

      integer*4 ierr, jerr  ! IOSTAT errors
     .,         indx  ! Position in string
     .,         js    ! Site number
     .,         i     ! Loop counter

      real*8 vel_args(12)  ! Velocity data line

      character*256 line    ! Line read from input file
      character*4 cd        ! Dummy characeters
      character*8 name      ! site name

****  Initialize the variables in the program
      num_vel_site = 0

****  Open the velocity file
      open(51,file=in_vel_file,status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',in_vel_file,1,
     .     'READ_VEL_MAV')

****  Now loop over the file
      do while ( ierr.eq.0 )
         read(51,'(a)',iostat=ierr) line
         if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
            indx = 0
            do i = 1,12
               call read_line(line,indx,'R8',jerr,vel_args(i),cd)
            end do
            if( jerr.eq.0 ) then
                call GetWord(line, name, indx)
                indx = 0
                call get_cmd(name,vel_name,num_vel_site, js,indx)
                if( js.gt.0 ) then
*                   Duplicate
                    write(*,120) vel_name(js)
 120                format('**DUPLICATE Velocity for ',a,
     .                     ' using later value')
                else
                    num_vel_site = num_vel_site + 1
                    js = num_vel_site
                    vel_name(js) = name
                endif

*               Now copy values
                vel_ll(1,js) = vel_args(1)
                vel_ll(2,js) = vel_args(2)
                vel_neu(1,js) = vel_args(4)/1000.d0
                vel_neu(2,js) = vel_args(3)/1000.d0
                if( zero_uvel ) then
                   vel_neu(3,js) = 0.d0
                else
                   vel_neu(3,js) = vel_args(10)/1000.d0
                end if
            end if
         end if
      end do

****  Thats all 
      return
      end

CTITLE CALCD

      real*8 function calcd(aprs,js)

      implicit none 

*     Routine to calculate distance between current site and 
*     velocity site js

      include 'merge_apr_vel.h'
      include '../includes/const_param.h'

* PASSED VARIABLES
      real*8 aprs(3)  ! XYZ of site
      integer*4 js    ! Velocity site number


* LOCAL VARIABLES

      real*8 rot_mat(3,3)  ! Rotation from XYZ to NEU
     .,      geod_pos(3)   ! Co-lat long (rads) Height (m)
     .,      dlat, dlng    ! Difference in lat and long (rads)

****  Convert XYZ to Geodetic Lat and log
*     (Geod_pos is colat and long in radians)
      call XYZ_to_GEOD(rot_mat,aprs,geod_pos)
      dlat = (pi/2-geod_pos(1))-vel_ll(2,js)*pi/180
      dlng = geod_pos(2) - vel_ll(1,js)*pi/180

      calcd = 6378153.d0*sqrt(dlat**2+(dlng*sin(geod_pos(1)))**2)

****  Thats all 
      return
      end

CTITLE CHANGE_VEL

      subroutine change_vel(aprs,js)

      implicit none 

*     Routine to calculate distance between current site and 
*     velocity site js

      include 'merge_apr_vel.h'

* PASSED VARIABLES
      real*8 aprs(6)  ! XYZ of site and velocity
      integer*4 js    ! Velocity site number


* LOCAL VARIABLES
      integer*4 i

      real*8 rot_mat(3,3)  ! Rotation from XYZ to NEU
     .,      geod_pos(3)   ! Co-lat, long and height
     .,      dXYZ(3)       ! Velocity in XYZ frame


****  OK, Rotate the velocity to XYZ
      call rotate_geod(vel_neu(1,js),dXYZ,'NEU','XYZ',
     .   aprs, geod_pos, rot_mat)

*     Now save the result
      do i = 1,3
         aprs(3+i) = dXYZ(i)
      end do
 
****  Thats all 
      return
      end
 




