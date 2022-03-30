
      program apr_to_vel

      implicit none 
 
*     This program read an apriori global apr file and a file containing
*     site names and velocity vectors and ouput a new file with the
*     velocities included.
 
      include '../includes/const_param.h'
 
*   max_update_sites    - Maxiumun number of sites which can be
*                       - updated
 
      integer*4 max_update_sites
 
      parameter ( max_update_sites = 1 )
 
*   len_run     - Length of runstring
*   i,j,k       - Loop counters
*   iel         - Site number from get_cmd
*   ierr        - IOSTAT error
*   jerr        - IOSTAT error reading input
*   indx        - pointer for readline
*   len_buffer  - Length of buffer (not used)
*   num_sites   - Number of sites in the apr_to_vel file
*   rcpar       - Runstring reader
*   trimlen     - String length function
 
      integer*4 len_run, iel, ierr, jerr, indx,
     .    rcpar, trimlen
 
*   site_name(max_update_sites)   - Names of the sites from the
*               - apr_to_vel file
      character*8 site_name(max_update_sites)
 
*   apr_to_vel_file  - Name of apr_to_vel file  (not used)
*   in_file     - Name of input file            
*   out_file    - Name of the output file      
 
      character*64 in_file, out_file  
 
*   buffer      - Line read from input
 
      character*132 buffer
 
*   apr_to_vel_motion(3,max_update_sites) - Three rotation rates for 
*               - for site  (not used)
*   sang        - Sin of angle from cross product  (not used)
*   values(7)   - Values for XYZ XYZ dot and epoch
 
 
c      real*8 apr_to_vel_motion(3,max_update_sites), sang (not used)
      real*8 values(7)

*   longitude, latitude - of site
*   rat_neu(3)  - RAte in local coords
*   rot_mat(3,3) - Rotation from XYZ to NEU
*   loc_coord(3) - Local coordinates of site

      real*8 longitude, latitude, rat_neu(3), rot_matrix(3,3), 
     .       loc_coord(3)
 
***** Decode the runstring
 
      len_run = rcpar(1, in_file )
      if( len_run.eq.0 ) call proper_runstring('apr_to_vel.hlp',
     .                   'apr_to_vel',1)
 
      len_run = rcpar(2, out_file)
      if( len_run.eq.0 ) call proper_runstring('apr_to_vel.hlp',
     .                   'apr_to_vel',1)
 
      open(100, file=in_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_file,1,'apr_to_vel')
 
      open(200, file=out_file, iostat=ierr, status='unknown')
      call report_error('IOSAT',ierr,'open',out_file,1,'apr_to_vel')
      jerr = 0

      write(*,100)
 100  format(' apr_to_vel:  Local coordinate system velocities')

      write(200,110) in_file(1:trimlen(in_file))
 110  format('* Velocity field contained in ',a,/,
     .       '*  Long      Lat       E & N vel ',t80,' Up vel',/,
     .       '*  deg       deg         mm/yr  ',t80,'(mm/yr)' )         
 
      do while ( jerr.eq.0 )
          read(100,'(a)', iostat=jerr ) buffer
 
*         See if we can find the site name
          indx = 1
          if( trimlen(buffer).gt.0 .and. buffer(1:1).eq.' ' ) then
              call getword( buffer, site_name, indx)
              iel = 1
          else
              iel = 0
          end if
 
*         Decode the line and add new information if found
          if( iel.gt.0 ) then
              call multiread( buffer, indx, 'R8', ierr, values,
     .                        ' ',7)
 
*             Now output the horizontal velocity field estimates
              if( abs(values(1)).gt.25.d0 ) then
                call rotate_geod(values(4), rat_neu, 'XYZ','NEU',
     .                           values(1), loc_coord, rot_matrix)
                latitude = 90.d0-loc_coord(1)*180/pi
                longitude = loc_coord(2)*180/pi

                write(200,220)   longitude, latitude,
     .              rat_neu(2)*1.d3, rat_neu(1)*1.d3,
     .              0.d0, 0.d0, 0.d0, 0.d0, 0.001d0,
     .              rat_neu(3)*1.d3, 0.d0, 0.d0,
     .              site_name(iel)
 220            format(2(1x,f9.4),1x,6(1x,f8.2),1x,f6.3,2x,
     .                 3(1x,f8.2), 1x,a8)
              end if

          end if
      end do
 
***** Thats all
      close(100)
      close(200)
      end
 
 
