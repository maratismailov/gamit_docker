CTITLE TRY_UPD_APR
 
      subroutine try_upd_apr( ns )

      implicit none
 
*     This routine will try to update the apriori station
*     coordinates by finding the original from either an
*     earthquake or a rename operation.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED VARIABLES
 
*         ns        - Number of site with no new apriori
 
      integer*4 ns
 
* LOCAL VARIABLES
 
*   i,j,k           - Loop counters
*   os              - Old site number to be copied
*   indx            - Position in string
 
      integer*4 i,j,k, l, os, indx
 
*   kbit            - Checks if bit is set
 
      logical kbit
 
*   new_name        - Trial new name for a site.
 
 
 
      character*8 new_name
 
***** Start by looping over the renames
      do i = 1, num_renames
 
          if( gsite_names(ns).eq.rn_codes(2,i) ) then
 
*             See if the original site name has an updated
*             coordinate.
              indx = 1
              call get_cmd(rn_codes(1,i),gsite_names, gnum_sites,
     .                    os, indx)

* MOD TAH 140722: if os returns as non-unique, try to find the last site
*             before this one with the same 4-char ID
              if( os.lt.0 ) then
                  do j = ns-1,1,-1
                     if( gsite_names(j)(1:4).eq.rn_codes(1,i)(1:4) .and.
     .                   kbit(gapr_updated,j) ) then
                         os = j
                         exit
                     endif
                  enddo
              endif
* MOD TAH 140722: Only report a problem if the renamed site is used
*             (Often EQ will change the name from the rename leaving
*              the renamed site un-used (e.g., _2PS -> _2HT with no
*              _2PS data remaining).
              if( os.le.0 .and. times_used(ns).gt.0 ) then
                  write(*,110) rn_codes(1,i), rn_codes(2,i), ns 
  110             format('** No site match with updated aprioris for ',
     .                   ' rename ',a,' to ',a,' Site # ',i5)
                  RETURN
              end if

              if( kbit(gapr_updated,os) ) then
 
*                 Use the old coordinate
                  do j = 1,3
                      apr_val_site(j,1,ns) = apr_val_site(j,1,os)
                      apr_val_site(j,2,ns) = apr_val_site(j,2,os)
                  end do
                  site_epoch(ns) = site_epoch(os)

* MOD TAH 140722: See if extended need to be copied too.  Only
*                 do this if we have used the site
                  if( times_used(ns).gt.0 ) then 
                    k = num_nonsec
                    do j = 1, k
                       if( param_nonsec(1,j).eq.os ) then
*                          Copy values to site ns
                           num_nonsec = num_nonsec+1
                           if( num_nonsec.gt.max_nonsec ) then
                               call report_stat('FATAL','GLOBK',
     .                             'try_upd_apr',' ',
     .                             'Too many non-secular extended mode', 
     .                              num_nonsec)
                           endif 

                           param_nonsec(1, num_nonsec) = ns
                           param_nonsec(2, num_nonsec) = 
     .                               param_nonsec(2,j)
                           apr_val_nonsec(:,num_nonsec) =
     .                         apr_val_nonsec(:,j)
                           write(*,115) gsite_names(ns), 
     .                                  j, gsite_names(os)
                       endif
                    end do
                  endif

                  call sbit(gapr_updated,ns,1)
                  if( times_used(ns).gt.0 )
     .            write(*,120) gsite_names(ns), rn_codes(1,i),
     .                   times_used(ns)
                  RETURN
              end if
          end if
      end do

 115  format('Updating RN non-sec terms at ',a,' from ',
     .                          I4,' at ',a)
 
 120   format('** No apriori coordinates for ',a8,
     .        ' using values for renamed site ',a8,
     .        ' Used ',i4,' times')

****  Now try earthquake site names
      do i = 1, num_eq
 
****      See if last two letters match.  If they try to find
*         original site
          if( gsite_names(ns)(7:8).eq.eq_codes(i)(1:2) ) then
 
              do j = ns-1,1,-1
 
*                 Form new name
                  new_name = gsite_names(j)
                  new_name(7:8) = eq_codes(i)(1:2)
                  call sub_char(new_name,' ','_')
                  
                  if( new_name.eq.gsite_names(ns) ) then
 
*                     Found, update position
                      if( kbit(gapr_updated,j) ) then
     
*                         Use the old coordinate
                          do k = 1,3
                              apr_val_site(k,1,ns) = 
     .                                 apr_val_site(k,1,j)
                              apr_val_site(k,2,ns) = 
     .                                 apr_val_site(k,2,j)
                          end do
                          site_epoch(ns) = site_epoch(j)

                          call sbit(gapr_updated,ns,1)
* MOD TAH 140722:         See if extended need to be copied too
*                         Only test if site is used.
                          if( times_used(ns).gt.0 ) then
                          k = num_nonsec
                            do l = 1, k
                                 if( param_nonsec(1,l).eq.j ) then
*                                    Copy values to site ns
                                     num_nonsec = num_nonsec+1
                                     param_nonsec(1, num_nonsec) = 
     .                                  ns 
                                     param_nonsec(2, num_nonsec) = 
     .                                      param_nonsec(2,l)
                                     apr_val_nonsec(:,num_nonsec) =
     .                                   apr_val_nonsec(:,l)
                                     write(*,150) gsite_names(ns), 
     .                                           l, gsite_names(j)
                                 endif
                            end do
                          endif
                          if( num_nonsec.gt.max_nonsec ) then
                              call report_stat('FATAL','GLOBK',
     .                           'try_upd_apr',' ',
     .                        'Too many non-secular extended mode', 
     .                            num_nonsec)
                          endif 
 
                          RETURN
                      end if
                  end if
              end do
          end if
      end do

150   format('Updating EQ non-sec terms at ',a,' from ',
     .                          I4,' at ',a)
 
****  Could not find a match, warning user
      if( num_apr_files.gt.0 .and. kbit(guse_site,ns) ) then
          write(*,300) gsite_names(ns)
 300      format('**WARNING** ',a8,' not in apr_files')
      end if
 
****  Thats all
      return
      end
