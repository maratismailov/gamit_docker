CTITLE FIND_PREQ
 
      subroutine find_preq(dir, ns, ne, ps)

      implicit none 
 
*     Routine to search the list of sites and earthquakes to
*     find the one which last effectd this site.i
* NOTE: This code assumes earthquakes are sorted in time.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED variables
 
*   ne      - Curremt earthquake number
*   ns      - Site that we are trying to match with the
*           - previous version of the site
*   ps      - The site to be found.  This is the site
*           - that was renamed to get the current site (ns)
 
      integer*4 ne, ns, ps
 
*   dir     - Direction (FO if solution is running forward
*           - and we want previous site name; or bakward (BA)
*           - and we want next site.
 
      character*(*) dir
 
* LOCAL variables
 
*   i,j,k   - Loop counters
*   indx    - Position in string
 
      integer*4 i,j,k, indx
 
*   found       - Indicates site found
*   orig_site   - remains true while this could still be the
*               - site we are lokking for.
 
      logical found, orig_site, kbit
 
*   new_name    - Generated name of site.
 
      character*8 new_name   

      orig_site = .false.

***** Check to see if we are using this site (i.e. it may not
*     exist due to renames

      ps = 0
      if( .not. kbit(guse_site,ns) ) RETURN
 
****  DO the search backwards in time when we are running the
*     forward solution.
      if( dir(1:2).eq.'FO' ) then
          found = .false.
          i = ne - 1
          do while ( .not.found .and. i.gt.0 )
 
*             Generate name for the site if this earthquake renmaed
              if( eq_rename(i) ) then
                  new_name = gsite_names(ns)
                  new_name(7:8) = eq_codes(i)
                  call sub_char(new_name,' ','_')
                  indx = 1
                  call get_cmd(new_name, gsite_names, gnum_sites,
     .                        ps, indx)
                  if( ps.gt.0 ) found = .true.
              end if
              i = i - 1
          end do
 
****      If we did not find a site.  It must have its original name
*         so find the site that generates the current name of the
*         site.
          if( .not.found ) then
              j = 0
              do while ( .not. found .and. j.lt. gnum_sites )
                  j = j + 1
                  if( j.ne.ns ) then
*                     Make sure not another earthquake site.
                      orig_site = .false.
                      if( gsite_names(j)(1:6).eq.
     .                    gsite_names(ns)(1:6) ) orig_site = .true.
                      do k = 1, num_eq
                          if( eq_rename(k) ) then
                              if( gsite_names(j)(7:8) .eq.
     .                            eq_codes(k) ) orig_site = .false.
                          end if
                      end do
                  end if
 
*                 See if this site matches the one we want if this
*                 could be the original site
                  if( orig_site ) then 
                      ps = j 
                      found = .true.
                  end if
              end do
          end if
 
*         Final check to make sure ps is non-zero
          if( ps.eq.0 ) then
              write(*,100) gsite_names(ns), ne, eq_codes(ne), dir
 100          format('**NOTE** Could not find original site for ',
     .                a,' for earthquake # ',i3,' Code ',a2,' for',
     .                ' a ',a,' direction')
          end if
 
*     Else we are check a back ward solution
      else
 
*         In this case start looking forward.  This is much simpler
          new_name = gsite_names(ns)
          new_name(7:8) = eq_codes(ne)
          call sub_char(new_name,' ','_')
          indx = 1
          call get_cmd(new_name, gsite_names, gnum_sites, ps, indx)
          if( ps.le.0 ) then
              write(*,100) gsite_names(ns), ne, eq_codes(ne), dir
              ps = 0
          end if
      end if
 
****  Thats all
      return
      end
 
CTITLE COPY_COVEQ
 
      subroutine copy_coveq( ps, ns, cov_parm, sol_parm )

      implicit none 
 
*     Routine to copy the covarinace rows and cols from site
*     ps to site ns.  Nothing happens if ns or ps is zero.
 
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED variables
 
*   ns      - Site that we are copying to
*   ps      - Site that we are copying from
 
      integer*4 ns, ps
 
*   cov_parm(num_glb_parn,num_glb_parn)     - Solution covariance
*           - matrix
*   sol_parm(num_glb_parn)      - Solution vector
 
      real*8 cov_parm(num_glb_parn,num_glb_parn),
     .    sol_parm(num_glb_parn)
 
* LOCAL VARIABLES
*   np, nn  - Parameter numbers of X coordinate of ps and ns
*   i       -  loop counter
 
      integer*4 np, nn, i
 
****  If either ps or ns is zero return
      if( ps.le.0 .or. ns.le.0 ) RETURN
 
      np = parn_site(1,1,ps)
      nn = parn_site(1,1,ns)
 
      if( np.eq.0 .or. nn.eq. 0 ) then
          write(*,100) gsite_names(ps), gsite_names(ns)
 100      format('***WARNING*** Zero parameter number copying ',
     .            a,' to ',a)
          RETURN
      end if
 
****  Now do the copy, column first
      write(*,9000) gsite_names(ps), gsite_names(ns)
 9000 format('Copying ',a,' covariances to ',a)
      do i = 0,2
          call dwmov(cov_parm(1,np+i),1,cov_parm(1,nn+i),1,
     .               num_glb_parn)
c           call dwmov(cov_parm(np,np+i),1,cov_parm(nn,nn+i),1,
c    .                3)
      end do
*     copy row
      do i = 0,2
          call dwmov(cov_parm(np+i,1),num_glb_parn,
     .               cov_parm(nn+i,1), num_glb_parn, num_glb_parn)
      end do
 
****  do solution
      do i = 0,2
          sol_parm(nn+i) = sol_parm(np+i)
      end do
 
****  Thats all
      return
      end
 
