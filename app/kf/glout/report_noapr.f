CTITLE report_noapr
 
      subroutine report_noapr(iout)

      implicit none 
 
*     Routine to report the names of sites with no aprori coordinate
*     updates/
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'
 
* PASSED variables
 
*   iout        - Output unit
 
 
      integer*4 iout
 
* LOCAL variables
 
*   i,j     - Loop counters
*     no      - NUmber of entris in line
 
 
      integer*4 i, trimlen, no
 
 
      character*100 line
 
 
      logical header_written, kbit
 
***** Start, see if we have any stations to report:

      if( num_apr_files.eq.0 ) RETURN
 
      header_written = .false.
      no = 0
      line = ' '
      do i = 1, gnum_sites
          if( kbit(guse_site,i) .and. .not. kbit(gapr_updated,i) ) then
              no = no + 1
              if( no.gt.10 ) then
                  if( .not.header_written ) then
                      write(iout,120)
 120                  format(/' SITES WITH NO UPDATED APRIORI',
     .                      ' COORDINATES:')
                      header_written = .true.
                  end if
                  write(iout,'(a)') line
                  line = ' '
                  no = 1
              end if
              line((no-1)*10+2:) = gsite_names(i)
          end if
      end do
 
*     See if we still have lines to write and header not written yet.
      if( no.gt.0 .and. .not.header_written ) then
        write(iout,120)
      end if
 
      if( trimlen(line).gt.0 ) then
          write(iout,'(a)') line
      end if
 
****  Thats all
 
 
****  Thats all
      return
      end
 
 
 
 
 
