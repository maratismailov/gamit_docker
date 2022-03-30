CTITLE RW_NAMES_BLOCK
 
      subroutine rw_names_block( option )

      implicit none 
 
 
*     Routine to read the names blocks from a global file and
*     match them to the names in the GLOBK common (or add them
*     if they are not already present)
 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/globk_common.h'
      include '../glinit/qsvi_rec.h'
 
*   i,j         - Loop counters
*   ierr        - FmpRead error.
*   iname(128)  - Integer equivalence to names_buffer (used
*               - for FmpRead.
*   len_read    - Length of record read from file (Not used)
*   lnum        - Number of current site or source in global
*               - list
*   site_boundary   - Record at which we shift from site names
*               - source names
 
      integer*4 i,j, ierr, iname(128), len_read, lnum, site_boundary
      integer*4 trimlen  ! length of string
 
*   option      - Option passed to routine (either R or W)
*   opt         - Upper case version of option.
 
      character*1 option, opt
 
*   names_buffer(32)    - Buffer containing the names to be
*               - read or written to the file.
*   next_name   - Next name read from the names buffer
 
      character*8 names_buffer(64), next_name, temp_name
      character*1 ssys      ! GNSS system.

      equivalence ( names_buffer, iname )
 
***** Get the option in upper case
      opt = option
      call casefold( opt )
 
***** Now act on option.
*                                 ! Read names.
      if( opt.eq.'R' ) then
 
*****     Do site names and source names together.
          do i = 1,cnum_names
              call readd(cglb_dcb,ierr, iname, 128, len_read,
     .                   crec_names+i-1)
 
*****         Now loop over the names in the names buffer.
              do j = 1, cnum_sites
                  call get_name_from_buffer( names_buffer, j, i, 0,
     .                    next_name)

* MOD TAH 130221: See if eq_reset is set and it is, reset the the
*                 exent name 
                  if( eq_reset .and. next_name(5:6).ne.'_X' .and.
     .                trimlen(next_name).gt.0 ) then
                      next_name = next_name(1:5) // 'GPS'
                  end if
                  call add_name_to_list( next_name, lnum,
     .                    gsite_names, gnum_sites )

*                 Check to see if number of sites is not greater
*                 limit.
                  if( gnum_sites.gt.max_glb_sites ) then
                      write(*,900) gnum_sites,  max_glb_sites
 900                  format('***DISASTER*** Too many sites.',
     .                        i5,' so far with ',i5,' maximum',/,
     .                       '               Modify  max_glb_sites',
     .                       ' in kalman_param.h to increase number',
     .                       ' of sites',/)
          
                      call report_stat('FATAL','GLOBK','rw_name_block',
     .                               ' ','Too many sites',0)
                  end if     
                  if( lnum.gt.0 ) then
                      times_used(lnum) = times_used(lnum) + 1
                      ltog_sites(j) = lnum
		           end if
              end do
 
*****         Get the source names
              do j = 1, cnum_sources
                  call get_name_from_buffer( names_buffer, j, i,
     .                    cnum_sites, next_name)
                  call add_name_to_list( next_name, lnum,
     .                    gsource_names, gnum_sources )
                  if( lnum.gt.0 ) ltog_sources(j) = lnum
              end do

*****         Get the SVS names
              do j = 1, cnum_svs
                  call get_name_from_buffer( names_buffer, j, i,
     .                    cnum_sites+cnum_sources, next_name)
                  if( trimlen(next_name).gt.0 .and. 
     .                trimlen(next_name).le.8 ) then
* MOD TAH 180402: Test to see if GNSS or PRN names
                      if( next_name(1:1).eq.'P' .and. use_prnn ) then  ! Old PRN Name 
*                         if GNSS name then svn number alreddy set                     
*                         Match the PRN number 
                          write(next_name(7:8),'(I2.2)') qsvi_svn(j)
* MOD TAH 180402: Update name to new satellite name type.  If PRN form must
*                 be GPS
                      elseif( next_name(1:1).eq.'P' ) then
                          ssys = 'G'
                          write(next_name,'(a1,I3.3,"_",I2.2)') ssys,
     .                      qsvi_svn(j), qsvi_prn(j) 
                      endif
* MOD TAH 180402: See if PRN names will be forced for GPS 
                      if( next_name(1:1).eq.'G' .and. use_prnn ) then
*                         Use temp_name to avoid writing and reading from 
*                         string at same time
                          write(temp_name,'("PRN_",a2,a2)') 
     .                        next_name(6:7), next_name(3:4) 
                          next_name = temp_name
!                         next_name = 'PRN_' // next_name(6:7) //
!    .                                 next_name(3:4) 
                      end if
                  endif

                  call add_name_to_list( next_name, lnum,
     .                    gsvs_names, gnum_svs )
                  if( lnum.gt.0 ) ltog_svs(j) = lnum
              end do
*                     ! Looping over records
          end do
*                     ! Read names.
      end if
 
***** Write names block.
 
*                                 ! Write out the list of names
      if( opt.eq.'W' ) then
 
*****     Write each record out with the names
          site_boundary = (cnum_sites - 1)/32 + 1
 
*                                     ! Loop over all of the records we
          do i = 1, cnum_names
*                                     ! need. NOTE: CNUM_NAMES was
*                                     ! computed when the glg_header was
*                                     ! written.
              call add_names_to_buffer( gsite_names, cnum_sites, i,
     .                gsource_names, cnum_sources, names_buffer)
 
              call writd(cglb_dcb, ierr, iname, 128, crec_names+i-1)
              call report_error('FmpWrite',ierr,'writ',names_buffer,
     .                    0,'RW_NAMES_BLOCK')
          end do
      end if
 
***** Thats all
      return
      end
 
