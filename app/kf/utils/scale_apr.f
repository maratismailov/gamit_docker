      program scale_apr
 
      implicit none 

*     This program will read one of my apriori solution files
*     and rescale it for the 2*10-8 scale difference between the
*     old relativity correction and the new.
 
*   ierr, jerr  - IOSTAT errors
*   rcpar   - Reads the runstring
*   len_run - Length of runstring
*   trimlen - Length of string
*   i,j     - Loop counters
*   indx    - Position of Site or Station to see that
*           - we have reached station data.
 
      integer*4 ierr, jerr, rcpar, len_run, trimlen, i,j, indx
 
*   x,y,z       - X,Y and Z from input file
*   vx,vy,vz        - velocities in x and y
*   epoch       - Epoch of postion.
 
      real*8 x,y,z, vx,vy,vz, epoch
 
*   not_site        - Indicates this is not a site coordinate
 
      logical not_site
 
*   site_name   - Name of site
 
      character*8 site_name
 
*   input_file  - Read from runstring
*   output_file - Read from runstring
 
      character*128 input_file, output_file, help_file
 
*   line        - Line read from input.
 
      character*256 line

      data help_file / 'scale_apr' / 

****  Decode the runstring
      len_run = rcpar(1, input_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( help_file, 6, 1)
      end if
      open(100, file=input_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',input_file,1,
     .                'scale_apr')
 
****  Get the output and open
      len_run = rcpar(2,output_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( help_file, 6, 1)
      end if
      open(200, file=output_file, iostat=ierr, status='new')
      call report_error('IOSTAT',ierr,'creat',output_file,1,
     .                'scale_apr')
 
 
****  Now start reading file looking for stations
 
      ierr = 0
      not_site = .true.
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
 
*         Try to find Station  or Site
          if( not_site .and. ierr.eq.0 ) then
              indx = index(line,'Station')
              if( indx.eq.0 ) indx = index(line,'Site')
*                                         ! We now look for site
              if( indx.gt.0 ) then
*                                     ! coordinates
                  not_site = .false.
              end if
          end if
 
*         See if comment
          if( ierr.eq.0 ) then
              if( line(1:1).ne.' ' .or. trimlen(line).eq.0 .or.
     .            not_site                                     ) then
 
*                Write out line
                 write(200,'(a)', iostat=ierr)
     .               line(1:max(1,trimlen(line)))
                 call report_error('IOSTAT',ierr,'writ',line,1,
     .                   'scale_apr')
 
*                         ! Decode as station position
              else
                 indx = 1
                 call Getword( line, site_name, indx )
                 read(line(indx:),*, iostat=jerr) x,y,z,vx,vy,vz, epoch
                 call report_error('IOSTAT',jerr,'decod',line,1,
     .                    'scale_apr')
 
*                Rescale
                 x = (1.d0 + 2.d-8) * x
                 y = (1.d0 + 2.d-8) * y
                 z = (1.d0 + 2.d-8) * z
 
*                Now write the line
                 write(200,200, iostat=ierr) site_name, x,y,z, vx,vy,vz,
     .                                    epoch
 200             format(1x,a8,3(1x,f13.4),1x,3(1x,f7.4),1x,f10.4)
*                     ! No error reading file
             end if
         end if
*                     ! Looping over input
      end do
 
****  Thats all
      close(100)
      close(200)
      end
 
 
 
