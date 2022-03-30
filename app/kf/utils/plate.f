      program plate

      implicit none 

 
*     This program read an apriori global apr file and a file containing
*     site names and velocity vectors and ouput a new file with the
*     velocities included.
* MOD TAH 120120: Added list of available frames to help
 
      include '../includes/const_param.h'
 
*   max_update_sites    - Maxiumun number of sites which can be
*                       - updated
 
      integer*4 max_update_sites
 
      parameter ( max_update_sites = 2048 )
 
*   len_run     - Length of runstring
*   i           - Loop counter
*   iel         - Site number from get_cmd
*   ierr        - IOSTAT error
*   jerr        - IOSTAT error reading input
*   indx        - pointer for readline
*   len_buffer  - Length of buffer
*   num_sites   - Number of sites in the plate file
*   rcpar       - Runstring reader
*   trimlen     - String length function
 
      integer*4 len_run, i, iel, ierr, jerr, indx, jndx, len_buffer,
     .    num_sites, rcpar, trimlen
 
*   plate_site(max_update_sites)   - Names of the sites from the
*               - plate file
      character*8 plate_site(max_update_sites)
 
*   plate_file  - Name of plate file
*   in_file     - Name of input file
*   out_file    - Name of the output file
*   vel_file    - Velocity file name
 
      character*64 plate_file, in_file, out_file, vel_file

*   ref_sys     - Reference system for velocities (default
*                 NUV-NNR)
*   site_prime  - Primary plate for site
*   site_sec    - Seconary plate for site

      character*10 ref_sys, site_prime, site_sec
 
*   buffer      - Line read from input
 
      character*132 buffer
 
*   plate_motion(3,max_update_sites) - Three rotation rates for
*               - for site
*   sang        - Sin of angle from cross product
*   values(7)   - Values for XYZ XYZ dot and epoch

*   rot_prime(3)  - Rotation vector for primary plate
*   rot_sec(3)    - Rotation vector for secondat plate
*   scale_sec     - Secondary scale factor
*   next_num      - Dummy numeric value to see if name or number
*                   given
 
 
      real*8 plate_motion(3,max_update_sites), sang, values(7),
     .       rot_prime(3), rot_sec(3), scale_sec, next_num 

*   longitude, latitude - of site
*   rat_neu(3)  - RAte in local coords
*   rot_mat(3,3) - Rotation from XYZ to NEU
*   loc_coord(3) - Local coordinates of site

      real*8 longitude, latitude, rat_neu(3), rot_matrix(3,3), 
     .       loc_coord(3)

*   next_word  - Next string pulled from buffer

      character*20 next_word
 
***** Decode the runstring
 
      len_run = rcpar(1, plate_file )
      if( len_run.eq.0 ) then
          call proper_runstring('plate.hlp','plate',0)
          call frame_to_frame('LIST', ref_sys, rot_prime)
          stop 
      end if
 
      len_run = rcpar(2, in_file )
      if( len_run.eq.0 ) call proper_runstring('plate.hlp','plate',1)
 
      len_run = rcpar(3, out_file)
      if( len_run.eq.0 ) call proper_runstring('plate.hlp','plate',1)
 
      len_run = rcpar(4, vel_file)
      if( len_run.eq.0 ) vel_file = out_file(1:trimlen(out_file)) 
     .                              // '.vel'

      len_run = rcpar(5, ref_sys )
      if( len_run.eq.0 ) ref_sys = 'NUV-NNR'
      call casefold(ref_sys)

*     Now read the plate file
      open(100, file=plate_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',plate_file,1,'PLATE')
 
      num_sites = 0
 
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) buffer
          len_buffer = trimlen(buffer)
          if( buffer(1:1).ne.' ' ) len_buffer = 0
*                                                       ! Extract info
          if( ierr.eq.0 .and. len_buffer.gt.0 ) then
              indx = 1
              num_sites = num_sites + 1
              call read_line(buffer,indx,'CH',ierr,i,
     .                       plate_site(num_sites) )
              call casefold(plate_site(num_sites))
 
*             Now get the velocity.  See if names of plates are given
              jndx = indx
              call getword(buffer, next_word, jndx )

*             see if we can read a numeric value:
              read(next_word,'(F20.10)', iostat=jerr ) next_num
              if( jerr.eq.0 ) then
                  call multiread(buffer, indx, 'R8', ierr,
     .                            plate_motion(1,num_sites),' ',3 )
 
*                         ! Convert from rads/my to rads/year
                  do i = 1,3
                      plate_motion(i,num_sites) =
     .                    plate_motion(i,num_sites)*1.d-6
                  end do
              else

*                 We seem to be passing names so decode the plate
*                 names
                  site_prime = next_word
                  call casefold(site_prime)
                  if( index(ref_sys,':A').gt.0 ) then
                      site_prime(trimlen(site_prime)+1:) = ':A'
                  end if
                  if( index(ref_sys,':O').gt.0 ) then
                      site_prime(trimlen(site_prime)+1:) = ':O'
                  end if
                  call frame_to_frame(site_prime, ref_sys, rot_prime)

*                 See if we have a scale factor to another plate
                  next_word = ' '
                  call getword(buffer, next_word, jndx )
                  if( trimlen(next_word).gt.0 .and. 
     .                next_word(1:1).ne.'!' ) then
                      read(next_word,*, iostat=jerr) scale_sec

*                     Now get the names of secondary plate
                      call getword(buffer, site_sec , jndx)
                      call casefold(site_sec)
                      if( index(ref_sys,':A').gt.0 ) then
                          site_sec(trimlen(site_sec)+1:) = ':A'
                      end if
                      if( index(ref_sys,':O').gt.0 ) then
                          site_sec(trimlen(site_sec)+1:) = ':O'
                      end if
                      call frame_to_frame(site_sec, ref_sys, rot_sec)
                  else
                      call frame_to_frame(site_prime, ref_sys, rot_sec)
                      scale_sec = 0.d0
                  end if

*****             Now create the composite pole of rotation
                  do i = 1,3 
                     plate_motion(i,num_sites) = rot_prime(i) +
     .                      scale_sec*(rot_sec(i)-rot_prime(i))
                  end do
              end if

          end if
      end do
 
*     Now open input and output files
      write(*,100)
 100  format(/' PLATE:  Local coordinate system velocities')

      open(100, file=in_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_file,1,'PLATE')
 
      open(200, file=out_file, iostat=ierr, status='unknown')
      call report_error('IOSAT',ierr,'open',out_file,1,'PLATE')

      open(201, file=vel_file, iostat=ierr, status='unknown')

*     write some header information to the output file
      write(200,120) plate_file(1:trimlen(plate_file)), ref_sys,
     .               in_file(1:trimlen(in_file)),
     .               out_file(1:trimlen(out_file)),
     .               vel_file(1:trimlen(vel_file))
 120  format('* PLATE: plate motion velocities of sites',/,
     .       '* Plate file         : ',a,' Reference system: ',a,/,
     .       '* Input apriori file : ',a,/,
     .       '* Output file        : ',a,/,
     .       '* Velocity file      : ',a,/,'*')
      
      write(201,120) plate_file(1:trimlen(plate_file)), ref_sys,
     .               in_file(1:trimlen(in_file)),
     .               out_file(1:trimlen(out_file)),
     .               vel_file(1:trimlen(vel_file))

      write(200,140)
 140  format('*  Name',11x,' X ',12x,' Y ',12x,' Z ',8x,
     .       ' dX/dt',3x,' dY/dt',3x,' dZ/dt',3x,
     .       'Ref Epoch',/,
     .       '*      ',11x,'(m)',12x,'(m)',12x,'(m)',8x,
     .       '(m/yr)',3x,'(m/yr)',3x,'(m/yr)')

      write(201,145)
 145  format('* Velocity field',/,
     .       '*  Long      Lat       E & N vel ',t80,' Up vel',/,
     .       '*  deg       deg         mm/yr  ',t80,'(mm/yr)' )

      jerr = 0
      do while ( jerr.eq.0 )
          read(100,'(a)', iostat=jerr ) buffer
          call casefold(buffer)
 
*         See if we can find the site name
          indx = 1
          call get_cmd( buffer, plate_site, num_sites, iel, indx)

          if( buffer(1:1).ne.' ' .or. jerr.ne. 0 ) iel = 0

*         Decode the line and add new information if found
          if( iel.gt.0 .and. jerr.eq.0 .and. buffer(1:1).eq.' ') then
              call multiread( buffer, indx, 'R8', ierr, values,
     .                        ' ',7)
 
              call cross_prod(plate_motion(1,iel), values, values(4),
     .                        sang)
 
              write(buffer,200, iostat=ierr) plate_site(iel), values
  200         format(1x,a8,3(1x,f15.5),1x,3(f9.5,1x),1x,f9.4)

*             Now output the horizontal velocity field estimates
              call rotate_geod(values(4), rat_neu, 'XYZ','NEU',
     .                         values(1), loc_coord, rot_matrix)
              latitude = 90.d0-loc_coord(1)*180/pi
              longitude = loc_coord(2)*180/pi

              write(201,220)   longitude, latitude,
     .            rat_neu(2)*1.d3, rat_neu(1)*1.d3,
     .            0.d0, 0.d0, 0.d0, 0.d0, 0.001d0,
     .            rat_neu(3)*1.d3, 0.d0, 0.d0,
     .            plate_site(iel)
 220          format(2(1x,f9.4),1x,6(1x,f7.2),1x,f6.3,2x,
     .               3(1x,f7.2), 1x,a8)

*
****      Write out error message
          else
              if( iel.eq.-1 ) write(*,250) buffer(1:indx)
 250          format('* No entry in plate file for        ',a)
              if( iel.eq.-2 ) write(*,260) buffer(1:indx)
 260          format('* Duplicate entry in plate file for ',a)
          end if
 
          if( jerr.eq.0 ) then
              len_buffer = max(1,trimlen(buffer))
              write(200,'(a)', iostat=ierr) buffer(1:len_buffer)
          end if
      end do
 
***** Thats all
      close(100)
      close(200)
      close(201)
      end
 
 
