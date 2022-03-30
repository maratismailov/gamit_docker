CTITLE SORT_APRIORI
 
      subroutine sort_apriori

      implicit none 
 
 
*     This routine will sort the sites by increasing west-longitude
*     and the sources by increasing RA. For the sites, only the
*     names, positions and axis offsets rate rearranged (All of the
*     tidal parameters are assumed to be site independent and therefore
*     are not changed).  For sources only the name and position
*     are changed.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'
 
*   i,j,k       - Loop counters
*   temp_i4     - I*4 variaable
 
      integer*4 i,j,k, temp_i4
 
*   longs(max_glb_sites)   - Longitudes used for sorting
*   temp_4      - Temporary Real*4 for switching
 
      real*4 longs(max_glb_sites), temp_4
 
*   temp_8      - Temporary Real*8 for switching
 
      real*8 temp_8
 
*   temp_ch     - Temporary character for switching.
 
      character*8 temp_ch
      character*8 prnj, prnjp  ! Strings with PRN numbers j, and j+1
 
***** Sort the sites.  First compute longtiudes
 
      do i = 1, gnum_sites
          if( apr_val_site(2,1,i).ne.0 .and.
*                                                     ! Check to avoid ATAN
     .        apr_val_site(1,1,i).ne.0       ) then
*                                                     ! error
              longs(i) = atan2(apr_val_site(2,1,i),apr_val_site(1,1,i))
          else
              longs(i) = 0.d0
          end if
 
          if( longs(i).lt.0 ) then
              longs(i) = longs(i) + 2*pi
          end if
      end do
 
***** Now do sort
      do i = 1, gnum_sites-1
          do j = 1, gnum_sites-i
*                                                ! West
* MOD TAH 031202: Placed a tolerance on how much less so that sites
*             that have been renamed will keep their ordder (ie.
*             oldest site name first).  Tolerance is 1d-5 which
*             is about 60m.
* MOD TAH 090911: Dropped tolerance 6 meters to reduce interleaved
*             sites.
* MOD TAH 161117: Dropped tolerance 1 meters to reduce interleaved
*             sites.

*             if( longs(j).lt.longs(j+1)-1.d-6 ) then ! 090911 code
              if( longs(j).lt.longs(j+1)-0.16d-6 ) then
*                                                ! longitude sort
                  call switch_ch(gsite_names(j), gsite_names(j+1),
     .                           temp_ch )
                  call switch_4(longs(j), longs(j+1), temp_4)
                  call switch_8(apr_val_axo(1,j), apr_val_axo(1,j+1),
     .                             temp_8)
                  call switch_8(apr_val_axo(2,j), apr_val_axo(2,j+1),
     .                             temp_8)
                  call switch_8(apr_val_atm(j),    apr_val_atm(j+1),
     .                             temp_8)
                  do k = 1,3
*                                                              ! Value
                      call switch_8(apr_val_site(k,1,j),
     .                              apr_val_site(k,1,j+1),temp_8)
*                                                              ! Rate
                      call switch_8(apr_val_site(k,2,j),
     .                              apr_val_site(k,2,j+1),temp_8)
                  end do
                  call switch_i4(times_used(j),times_used(j+1),temp_i4)
              end if
          end do
      end do
 
 
***** Now do sort for source positions
 
      do i = 1, gnum_sources-1
          do j = 1, gnum_sources-i
 
              if( apr_val_source(1,1,j).gt.
*                                                ! Switch locations
     .            apr_val_source(1,1,j+1) ) then
 
                  call switch_ch(gsource_names(j), gsource_names(j+1),
     .                           temp_ch )
 
*                             ! loop over RA and dec
                  do k = 1,2
*                                                                ! Value
                      call switch_8(apr_val_source(k,1,j),
     .                                 apr_val_source(k,1,j+1),temp_8)
*                                                                ! Rate
                      call switch_8(apr_val_source(k,2,j),
     .                                 apr_val_source(k,2,j+1),temp_8)
                  end do
              end if
          end do
      end do
 
***** Now do sort SVS in increasing order of PRN number

      do i = 1, gnum_svs-1
          do j = 1, gnum_svs-i

* MOD TAH Sort by PRN even if using G034_18 rather than PRN_1834 form.
*             Create test strings
              if( use_prnn ) then
                 prnj  = gsvs_names(j)
                 prnjp = gsvs_names(j+1)
              else
                 prnj  = gsvs_names(j)(6:)
                 prnjp = gsvs_names(j+1)(6:)
              endif

C             if( gsvs_names(j).gt. gsvs_names(j+1) ) then
              if( prnj         .gt. prnjp           ) then

                  call switch_ch(gsvs_names(j), gsvs_names(j+1),
     .                           temp_ch )

*                             ! loop over max_svs_elem orbital elements
                  do k = 1,max_svs_elem
*                                                                ! Value
                      call switch_8(apr_val_svs(k,j),
     .                              apr_val_svs(k,j+1),temp_8)
                  end do
              end if
          end do
      end do

***** Thats all
      return
      end

CTITLE CHECK_ALL_CRD

      subroutine check_all_crd

      implicit none 
 
 
*     This routine will check that all the sites have some apriori
*     coordinate.  With earthquake renames for sites that have no
*     measurements before the earthquake, we can be left in a situation
*     where the original _GPS site name has no coordinates.  Fix that
*     problem here.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'
 
*   i,j       - Loop counters
*   copy_site - Number of site to take copy of corrdinates from.
 
      integer*4 i,j, copy_site
 
***** Loop over the sites and if any have zero apriori
 
      do i = 1, gnum_sites
          if( sqrt(apr_val_site(1,1,i)**2+
     .             apr_val_site(2,1,i)**2+
     .             apr_val_site(3,1,i)**2).lt.6.d6 ) then
*             Values too small (< than 6000km) so find alternative
              copy_site = 0
              do j = 1, gnum_sites
                 if( j.ne.i ) then     ! See if first 4 characters match
                    if( gsite_names(j)(1:4).eq.gsite_names(i)(1:4) .and.
     .                  sqrt(apr_val_site(1,1,i)**2+
     .                  apr_val_site(2,1,i)**2+
     .                  apr_val_site(3,1,i)**2).lt.6.d6 ) then
                        copy_site = j
                    end if 
                 end if
              end do
*             Copy site coordinates over if we found one
              if( copy_site.gt.0 ) then
                  do j = 1,3
                     apr_val_site(j,1,i) = apr_val_site(j,1,copy_site)
                     apr_val_site(j,2,i) = apr_val_site(j,2,copy_site)
                  end do
C                 write(*,120) gsite_names(i), gsite_names(copy_site)
 120              format('No initial coordinate for ',a8,' copying ',
     .                   a8,' coordinates')
              else
                  write(*,140) gsite_names(i)
 140              format('**WARNING** No coordinate for ',a8,
     .                   ' or any site with same first 4 characters',/,
     .                   'If this is causes problems, put coordinates ',
     .                   'in apr_file')
              end if
          end if
      end do

****  Thats all
      return
      end
 
