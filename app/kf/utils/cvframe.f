      program cvframe

      implicit none 

*     Program to change the velocity field in a .vel file
*     into a new frame.

      include '../includes/const_param.h'

      integer*4 ierr, jerr, rcpar, trimlen, len_run, jndx, i
      real*8 lat, long, ve, vn, vu, de, dn,du, se, sn, su, rho
      real*8 rot_prime(3), geod_pos(3), loc_coord(3), site_pos(3),
     .       site_vel(3), neu_vel(3), rot_matrix(3,3), sang

      character*256 invel, outvel, line, test
      character*8 inframe, outframe, cd
      character*9 site
      logical user_rotp

      len_run = rcpar(1,invel)
      if( len_run.eq.0 ) then
           call proper_runstring('cvframe.hlp','cvframe/invel',1)
      endif
      len_run = rcpar(2,outvel)
      if( len_run.eq.0 ) then
           call proper_runstring('cvframe.hlp','cvframe/outvel',1)
      endif
      len_run = rcpar(3,inframe)
      call casefold(inframe)
      if( len_run.eq.0 ) then
           call proper_runstring('cvframe.hlp','cvframe/inframe',1)
      endif
      len_run = rcpar(4,line )
      if( len_run.eq.0 ) then
           call proper_runstring('cvframe.hlp','cvframe/outframe',1)
      endif
      
* MOD TAH 981204: See if vector passed.
      jndx = 0
      call sub_char( line,':',' ')
      call getword( line, test, jndx)
      call check_num( test, ierr )
      if( ierr.eq.0 ) then
          jndx = 0
          call multiread(line,jndx, 'R8',jerr, rot_prime, cd, 3)
          if( jerr.ne.0 ) then
             write(*,105) jerr,line(1:trimlen(line))
 105         format('IOSTAT Error ',i5,' reading pole ',a,/,
     .              'Exiting')
          endif

          do i = 1,3
             rot_prime(i) = rot_prime(i)*pi/180.d0/1.d6
          end do
          outframe = 'USER'
          user_rotp = .true.
      else
          user_rotp = .false.
          call casefold(line)
          outframe = line
      endif      

****  OPen up the input and output
      open(100,file=invel, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',invel,1,'cvframe')
      open(200,file=outvel, status='unknown', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',outvel,1,'cvframe')

*     Get the rotation rate between the two frames
      if( .not.user_rotp ) then
          call frame_to_frame(outframe, inframe, rot_prime)
          write(*,110) (rot_prime(i)*180/pi*1.d6, i=1,3)
 110      format('CVFRAME: Rotation Pole ',3F12.6,' deg/Myr')
      end if
      
      if( abs(rot_prime(1))+abs(rot_prime(2))+abs(rot_prime(3)).lt.
     .    1.d-15 ) then
          write(*,120) inframe, outframe
 120      format('No rotation rate between ',a8,' and ',a8)
          stop 'CVFRAME: No rotation rate difference'
      end if

*     Write header:
      write(200,140) invel(1:trimlen(invel)), inframe, outframe
 140  format('* CVFRAME: Vel file ',a,' rotated from ',a,' to ', a)
      write(200,145) (rot_prime(i)*180/pi*1.d6, i=1,3)
 145  format('* Rotation Pole ',3F12.6,' deg/Myr')
      write(200,155)
 155  format('*  Long.       Lat. ',7x,'E & N Rate ',3x,
     .       ' E & N Adj. ',2x,' E & N +-',1x,
     .       ' RHO ',5x,' H Rate  H adj.   +-',1x,'SITE',/,
     .       '*  (deg)      (deg) ',2x,3(6x,'(mm/yr)'),15x,
     .          '(mm/yr)' )


*     Start reading the input velocity file
      do while (ierr.eq.0 )
         read(100,'(a)',iostat=ierr) line
         if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .       trimlen(line).gt.0  ) then
             jndx = 0
             call read_line(line, jndx, 'R8', jerr, long, cd)
             call read_line(line, jndx, 'R8', jerr, lat , cd)
             call read_line(line, jndx, 'R8', jerr, ve  , cd)
             call read_line(line, jndx, 'R8', jerr, vn  , cd)
             call read_line(line, jndx, 'R8', jerr, de  , cd)
             call read_line(line, jndx, 'R8', jerr, dn  , cd)
             call read_line(line, jndx, 'R8', jerr, se  , cd)
             call read_line(line, jndx, 'R8', jerr, sn  , cd)
             call read_line(line, jndx, 'R8', jerr, rho , cd)
             call read_line(line, jndx, 'R8', jerr, vu  , cd)
             call read_line(line, jndx, 'R8', jerr, du  , cd)
             call read_line(line, jndx, 'R8', jerr, su  , cd)
             call read_line(line, jndx, 'CH', jerr, su  , site)
C            read(line,220,iostat=jerr) long, lat, ve, vn, de,dn, 
C    .                se, sn, rho, vu, du, su, site
             if( jerr.eq.0 ) then

*                OK, so convert
*                convert lat and long to XYZ
                 geod_pos(1) = (90.d0-lat)*pi/180.d0
                 geod_pos(2) = long*pi/180.d0
                 geod_pos(3) = 0.d0

                 call geod_to_XYZ( geod_pos, site_pos )

*                Now get the velocity in the new frame
                 call cross_prod(rot_prime, site_pos, site_vel, sang)
                 call rotate_geod(site_vel, neu_vel, 'XYZ', 'NEU',
     .                            site_pos, loc_coord, rot_matrix )

*                Now substract the frame velocity and output results
                 ve = ve - neu_vel(2)*1000.d0
                 vn = vn - neu_vel(1)*1000.d0

                 write(200,220) long, lat, ve, vn, de,dn,
     .                se, sn, rho, vu, du, su, site      
 220             format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .                3(1x,f7.2), 1x,a9)
              end if
          else
              if( ierr.eq.0 ) write(200,'(a)' ) 
     .                        line(1:max(1,trimlen(line)))
          end if
      end do
      if( jerr.ne.0 ) then
          write(*,300) jerr
 300      format('IOSTAT error ',i5,' occurred decoding some lines')
      end if


****  Close files
      close(100)
      close(200)
      end



