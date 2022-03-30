CTITLE SUMM_EQ

      subroutine summ_eq(iout, options, cov_parm, sol_parm )

      implicit none 

*     Routine to check sites that have names which suggest that they are the
*     site and see how well there estimates match.  For example, GOLD_GHT would be
*     compared to GOLD_GLA etc.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
*                                           ! To get nutation names
      include '../includes/globk_cmds.h'
      include '../includes/sd_common.h'
      include '../includes/const_param.h'

* PASSED VARIABLES
      integer*4 iout   ! Output unit number
      integer*4 options   ! Options for output (nothing specific yet).

      real*8   cov_parm(num_glb_parn,num_glb_parn), 
     .         sol_parm(num_glb_parn)     ! Solution and covariance function

* LOCAL VARIABLES
      integer*4 np(2)   ! Two parameter numbers associated with X position of matched sites
      integer*4 nv(2)   ! Two parameter numbers associated with X velocity of matched sites
      integer*4 i,j,k   ! Loop counters
      integer*4 eq_num  ! Number of earthquake
      integer*4 ieu, ierr     ! Unit number for equate line output 
      integer*4 iel, jel      ! Parameter numbers
      real*8 comp_jd    ! Date for comparing positions at
      real*8 dlog(3), slog(3), dchi(3)   ! Difference and sigma of log differences

      character*8 nsi, nsj  ! Names of two sites being considered
      character*128 eqs_file ! Name of file to which equate lines are written)

      character*4 vlog(3)

      logical logout   ! Set true to output log estimates
      logical kbit, norename     ! test bit, set false once a rename match has been found.

      data  vlog / 'nlog','elog','ulog' /
      
*     Open the output for the equates lines that can be used by sh_exeqs
      call gen_org_name( eqs_file, newoutfile, 'eqs' )

      ieu = 203
      call open_lu(ieu, eqs_file, ierr,'unknown')
      write(ieu,120)
 120  format('* Equate lines: Use sh_exeqs to extract lines')

      call report_error('IOSTAT',ierr,'open',eqs_file,0,'SUMM_EQ')
      if( ierr.ne.0 ) then
          ieu = 6
          write(ieu,120)
      end if

*     Loop over all the site names to see which ones seem to match.  Use any rename
*     entries to guide the date the comparison should be made at.
      write(iout,210)
 210  format(/,'* RENAME REPORT (RNRP option)',/,
     .         '*   Sites             dN        sN       dE       sE',
     .         '     dU         sU  Units Compare date',
     .         '  EQ Dist EQ Name') 

      
      do i = 1, gnum_sites-1

*        See which parameters have been estimated for this site
         if( KBIT(guse_site,i) ) THEN
         np(1) = parn_site(1,1,i)
         nv(1) = parn_site(1,2,i)
         nsi = gsite_names(i)

*        Check for name matches
         norename = .true.
         j = i
         do while ( j.lt.min(i+5,gnum_sites) .and. norename )
            j = j + 1
            if( KBIT(guse_site,j) ) THEN
            nsj = gsite_names(j)
            if( nsi(1:4).eq.nsj(1:4) ) then
*               Site names match.  See which parameters have been
*               estimated for each site
                np(2) = parn_site(1,1,j)
                nv(2) = parn_site(1,2,j)
*
*               For the position comparison, see if we can find a 
*               rename entry that will give us a date.
                comp_jd = 0
                do k = 1, num_renames
                  if(  rn_codes(2,k)(1:6).eq.nsj(1:6) )  then

*                      Now check that the last two characters match
                       if( nsi(7:8).eq.nsj(7:8) ) then 
                           comp_jd = rn_times(1,k)
                           norename = .false.
                       end if
                   end if
                end do
*
*               Now see if there is an earthquake relationship 
*               between the sites
                eq_num = 0
                if( comp_jd.eq.0 ) then
                    do k = 1, num_eq
                       if( nsj(7:8).eq.eq_codes(k)(1:2)    ) then
                           comp_jd = eq_epoch(k)
                           eq_num = k
                           norename = .false.
                       endif
                    end do
                end if
*
*               If we still have not found, check again to see if 
*               user renamed with just the last two characters.  Done
*               at this time to avoid conflict with earthquake renames
                if( comp_jd.eq.0 ) then
                   do k = 1, num_renames
                      if( rn_codes(2,k)(1:6).eq.nsj(1:6) )  then
                          comp_jd = rn_times(1,k)
                          norename = .false.
                      end if
                   end do
                end if

*               If we still do not have a date, compare at the end
*               epoch for the moment
C                if ( comp_jd.eq.0 ) then
C                     comp_jd = gepoch_end
C                endif

*               Now do the comparision
                if( np(1).gt.0 .and. np(2).gt.0 .and.comp_jd.gt.0 ) then
                    call compare_est(iout, ieu, options,i,j, eq_num,
     .                   cov_parm, sol_parm,np,nv,comp_jd,'POS')
                end if
                if( nv(1).gt.0 .and. nv(2).gt.0 .and.comp_jd.gt.0 ) then
                    call compare_est(iout, ieu, options, i,j, eq_num,
     .                   cov_parm, sol_parm,np,nv,comp_jd,'VEL')
                endif
            endif
            ENDIF
         end do
         ENDIF
      end do

***** Now report on log estimates that overlap in time at a site.
*     Search for all sites common to an earthquake
      write(iout,310)
 310  format('* EARTHQUAKE LOG FITS',/,
     .       '*  Sites             dN        sN       dE       sE',
     .       '     dU         sU   mm') 
      do i = 1, gnum_sites-1
          do j = i+1, gnum_sites
*            See if both used
             logout = .false.
             if( kbit(guse_site,i) .and. kbit(guse_site,j) ) then
*                if beginning and ends of names match
                 if( gsite_names(i)(1:4).eq.gsite_names(j)(1:4).and.
     .               gsite_names(i)(7:8).eq.gsite_names(j)(7:8) ) then
*
*                    Get the parmater numbers for the log estimates
                     do k = 1,3
                        iel = parn_log(k,i)
                        jel = parn_log(k,j)
                        if( iel.gt.0 .and. jel.gt.0 ) then
                             dlog(k) = sol_parm(iel)-sol_parm(jel)
                             slog(k) = (cov_parm(iel,iel) +
     .                                  cov_parm(jel,jel) -
     .                                2*cov_parm(iel,jel))
                             if( slog(k).gt.0 ) then
                                 slog(k) = sqrt(slog(k))
                             else
                                 slog(k) = 0.d0
                             end if
                             logout = .true.
                             if( slog(k).gt.1.d-6) then
                                 dchi(k) = abs(dlog(k)/slog(k))
                                 if( dchi(k).gt.99.99d0 ) 
     .                                           dchi(k) = 99.99d0 
                             else 
                                 dchi(k) = 99.99d0
                             endif
                        else
                             dlog(k) = 0.d0
                             slog(k) = 0.d0
                             dchi(k) = 99.99d0
                        end if
                     end do
*                    Output the results
                     if( logout ) then
                        dchi(1) = max(dchi(1),dchi(2))
                        dchi(2) = dchi(1)
                        write(iout,320) gsite_names(i), gsite_names(j),
     .                        (dlog(k)*1.d3,slog(k)*1.d3,k=1,3)

 320                    format(A8,'-',A8,3(F8.2,' +- ',f6.2),' mm DLOG')

*                       Write out the equate commands
                        do k = 1,3
                            write(ieu,340) 'Lo',dchi(k), gsite_names(i), 
     .                          vlog(k),  gsite_names(j), vlog(k)
 340                        format('EQ',A2,1x,F5.2,'  equate ',
     .                            a8,1x,a4,1x,a8,1x,a4)
                        end do
                     end if
                 end if
             end if
          end do
      end do 


*     Thats all
      return
      end 
 
CTITLE COMPARE_EST

      subroutine compare_est(iout, ieu, options,s1,s2,eq_num, 
     .           cov_parm, sol_parm,np,nv, comp_jd,type)

      implicit none 

*      Routine to compare either POS or VEL type estimates at time comp_jd
*      (time only needed for POS estimate when velocities are estimated).

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/const_param.h'

* PASSED VARIABLES
      integer*4 iout, ieu  ! Output unit numbers
      integer*4 options   ! Options for output (nothing specific yet).
      integer*4 s1, s2    ! Two site numbers
      integer*4 np(2), nv(2)  ! Parameter numbers of X component of 
                              ! position and velocity
      integer*4 eq_num        ! number of earthquake for this rename
                              ! 0 if not an earthquake

      real*8   cov_parm(num_glb_parn,num_glb_parn), 
     .         sol_parm(num_glb_parn),     ! Solution and covariance function
     .         comp_jd        ! JD at which positions should be compared

      character*(*) type      ! Either POS or VEL comparison type

* LOCAL VARIABLES
      real*8 rot_matrix(3,3)  ! Rotation matrix from XYZ to NEU
      real*8 loc_coord(3)     ! Geodetic co-lat, long and height
      real*8 covar(12,12)     ! Covariance matrix between XYZ pos and velocity
                              ! at the two sites
      real*8 tempc(12,12)     ! Scratch intermediate matrix
      real*8 cov_dNEU(3,3)    ! Covariance matrix between Delta NEU
      real*8 map(3,12)   ! Matrix to map XYZ value and velocity to dNEU
      real*8 dt               ! Time difference from expt_epoch and comp_jd (years)
      real*8 dNEU(3)          ! Difference in North, East and Up at comp_jd
      real*8 dXYZ(12)         ! Differnce in XYZ value and velocity
      real*8 eq_dist          ! Distance to eq hypocenter
      real*8 dchi(3)          ! Chi of adjustment (forced to be same for NE)

      integer*4 i,j,k, date(5)
      real*8 sectag
      real*8 dts1, dts2       ! Time difference (yrs) between site epoch and
                              ! experiment epoch.

      character*8 code_eq 
      character*4 pet(3), vet(3)

      data pet / 'npos','epos','upos' /, vet / 'ndot','edot','udot' /

***** Get the elements we need from the covaraince matrix and the solution vector
*     Get the transformation into NEU
      call XYZ_to_NEU(rot_matrix, apr_val_site(1,1,s1),loc_coord)

      call jd_to_ymdhms(comp_jd, date, sectag)
      dt = (comp_jd - gepoch_end)/365.25d0
      dts1 = (gepoch_end-site_epoch(s1))/365.25d0
      dts2 = (gepoch_end-site_epoch(s2))/365.25d0

****  MOD TAH 070126: If the site epoch is zero (constant apriori
*     make the dts time difference be zero
      if( site_epoch(s1).eq.0 ) dts1 = 0.d0
      if( site_epoch(s2).eq.0 ) dts2 = 0.d0

      if( eq_num.ne.0 ) then
          call eval_dist( eq_pos(1,eq_num), apr_val_site(1,1,s1), 
     .                    eq_dist)
          code_eq = eq_codes(eq_num)
      else
          eq_dist = 0.d0
          code_eq = 'None'
      endif

*     Clear the covariance matrix before filling (this way we do not 
*     need to worry about the entries that are not filled
      call dwint(0.d0,covar,1,144)
      call dwint(0.d0,map,1,36)

****  Get the difference in velocity if type is vel
      if( type.eq.'VEL') then
*         Get the covariance entries for the covar matrix 
*         (only the 6x6 part used)
          call cp_cov(nv(1),nv(1),3,covar(1,1),12,
     .                 cov_parm,num_glb_parn)
          call cp_cov(nv(2),nv(2),3,covar(4,4),12,
     .                 cov_parm,num_glb_parn)
          call cp_cov(nv(1),nv(2),3,covar(1,4),12,
     .                 cov_parm,num_glb_parn)
          call cp_cov(nv(2),nv(1),3,covar(4,1),12,
     .                 cov_parm,num_glb_parn)

*         Get the velocity diffences (estimate+any difference in the apriori)
          do i = 1,3
             dXYZ(i)   = sol_parm(nv(1)+i-1)+apr_val_site(i,2,s1)
             dXYZ(i+3) = sol_parm(nv(2)+i-1)+apr_val_site(i,2,s2)
          end do
*
*         Create the mapping matrix
          do i=1,3
             do j=1,3
                map(i,j)   =  rot_matrix(i,j)
                map(i,j+3) = -rot_matrix(i,j)
             end do
          end do

*         Generate the estimates
          do i = 1,3
             dNEU(i) = 0.0d0
             do j = 1,6
                dNEU(i) = dNEU(i)+map(i,j)*dXYZ(j)
             end do
          end do

*         Now compute the covariance matrix of the difference
          call var_comp(map,covar,cov_dNEU,tempc,3,12,1)

          call comp_dchi(dNEU,cov_dNEU, dchi)

*         Output the results
* MOD TAH 040309: change the order so that result will be 
*         coseismic velocity change (similar mod made the position
*         output)
          write(iout,220) gsite_names(s2), gsite_names(s1),
     .         (-dNEU(i)*1.d3,sqrt(abs(cov_dNEU(i,i)))*1.d3,i=1,3),
     .         (date(k),k=1,3), eq_dist/1.d3, code_eq
 220      format(A8,'-',A8,3(F8.2,' +- ',f6.2),' mm/yr ',
     .          1x,I4,1x,i2,1x,i2,1x,F7.1,' km, Eq ',a4)

*         Write out the equate commands
          do i = 1,3
             write(ieu,240) code_eq(1:2),dchi(i), gsite_names(s1), 
     .                       vet(i),  gsite_names(s2), vet(i)
 240         format('EQ',A2,1x,F5.2,'  equate ',a8,1x,a4,1x,a8,1x,a4) 
                                    
          end do
      end if

****  Now do position
      if( type.eq.'POS') then
*         Get the covariance entries for the covar matrix (only the 6x6 part used)
*         (Check that velocities are estimated before copying matrix.  Since covar
*         is cleared before the start, we can just not copy of velocities are not
*         estimated.)
          call cp_cov(np(1),np(1),3,covar(1,1),12,
     .                 cov_parm,num_glb_parn)
          if( nv(1).gt.0 ) 
     .    call cp_cov(nv(1),nv(1),3,covar(4,4),12,
     .                 cov_parm,num_glb_parn)
          call cp_cov(np(2),np(2),3,covar(7,7),12,
     .                 cov_parm,num_glb_parn)
          if( nv(2).gt.0 )
     .    call cp_cov(nv(2),nv(2),3,covar(10,10),12,
     .                 cov_parm,num_glb_parn)
*         Off-diagonal terms
          if( nv(1).gt.0 )
     .    call cp_cov(np(1),nv(1),3,covar(1,4),12,
     .                 cov_parm,num_glb_parn)
          call cp_cov(np(1),np(2),3,covar(1,7),12,
     .                 cov_parm,num_glb_parn)
          if( nv(2).gt.0 )
     .    call cp_cov(np(1),nv(2),3,covar(1,10),12,
     .                 cov_parm,num_glb_parn)

          if( nv(1).gt.0 )
     .    call cp_cov(nv(1),np(1),3,covar(4,1),12,
     .                 cov_parm,num_glb_parn)
          call cp_cov(np(2),np(1),3,covar(7,1),12,
     .                 cov_parm,num_glb_parn)
          if( nv(2).gt.0 )
     .    call cp_cov(nv(2),np(1),3,covar(10,1),12,
     .                 cov_parm,num_glb_parn)

*         next term in
          if( nv(1).gt.0 )
     .    call cp_cov(nv(1),np(2),3,covar(4,7),12,
     .                 cov_parm,num_glb_parn)
          if( nv(1).gt.0 .and. nv(2).gt.0 )
     .    call cp_cov(nv(1),nv(2),3,covar(4,10),12,
     .                 cov_parm,num_glb_parn)

          if( nv(1).gt.0 )
     .    call cp_cov(np(2),nv(1),3,covar(7,4),12,
     .                 cov_parm,num_glb_parn)
          if( nv(1).gt.0 .and. nv(2).gt.0 )
     .    call cp_cov(nv(2),nv(1),3,covar(10,4),12,
     .                 cov_parm,num_glb_parn)
 
*         Last off diagonal block
          if( nv(2).gt.0 )
     .    call cp_cov(np(2),nv(2),3,covar(7,10),12,
     .                 cov_parm,num_glb_parn)

          if( nv(2).gt.0 )
     .    call cp_cov(nv(2),np(2),3,covar(10,7),12,
     .                 cov_parm,num_glb_parn)


*         Get the velocity diffences (estimate+any difference in the apriori)
          do i = 1,3
* MOD TAH 060111: Updated velocity use for case when velocities are estimated
*            and maybe different between sites
             if( nv(1).eq.0 ) then
                dXYZ(i)   = sol_parm(np(1)+i-1)+
     .              (apr_val_site(i,1,s1)+apr_val_site(i,2,s1)*dts1)
             else
                dXYZ(i)   = sol_parm(np(1)+i-1)+
     .              (apr_val_site(i,1,s1)+
     .              (sol_parm(nv(1)+i-1)+apr_val_site(i,2,s1))*dts1)

                dXYZ(i+3) = sol_parm(nv(1)+i-1)+apr_val_site(i,2,s1)
             end if
* MOD TAH 060111: Same mode for velocities
             if( nv(2).eq.0 ) then
                dXYZ(i+6) = sol_parm(np(2)+i-1)+
     .               (apr_val_site(i,1,s2)+apr_val_site(i,2,s2)*dts2)
             else
                dXYZ(i+6) = sol_parm(np(2)+i-1)+
     .               (apr_val_site(i,1,s2)+
     .               (sol_parm(nv(2)+i-1)+apr_val_site(i,2,s2))*dts2)

                dXYZ(i+9) = sol_parm(nv(2)+i-1)+apr_val_site(i,2,s2)
             endif
          end do
*
*         Create the mapping matrix

          do i=1,3
             do j=1,3
                map(i,j)   =  rot_matrix(i,j)
                map(i,j+3) =  rot_matrix(i,j)*dt

                map(i,j+6)  = -rot_matrix(i,j)
                map(i,j+9)  = -rot_matrix(i,j)*dt
  
             end do
          end do

*         Generate the estimates
          do i = 1,3
             dNEU(i) = 0.0d0
             do j = 1,12
                dNEU(i) = dNEU(i)+map(i,j)*dXYZ(j)
             end do
          end do

*         Now compute the covariance matrix of the difference
          call var_comp(map,covar,cov_dNEU,tempc,3,12,1) 
          call comp_dchi(dNEU,cov_dNEU, dchi)

*         Output the results. Switch the order to get the coseismic
*         offsets
          write(iout,320) gsite_names(s2), gsite_names(s1),
     .         (-dNEU(i)*1.d3,sqrt(abs(cov_dNEU(i,i)))*1.d3,i=1,3),
     .         (date(k),k=1,3), eq_dist/1.d3, code_eq
 320      format(A8,'-',A8,3(F8.2,' +- ',f6.2),' mm    ',
     .          1x,I4,1x,i2,1x,i2,1x,F7.1,' km, Eq ',a4)

*         Write out the equate commands
          do i = 1,3
             write(ieu,340) code_eq(1:2),dchi(i), gsite_names(s1), 
     .                       pet(i),  gsite_names(s2), pet(i)
 340         format('EQ',A2,1x,F5.2,'  equate ',a8,1x,a4,1x,a8,1x,a4) 
                                   
          end do

      end if

      return
      end

      subroutine cp_cov(nv1,nv2,nel,covar,ndim,
     .                 cov_parm,ncov)

      implicit none 

*     Routine to move a nelxnel block of the covariance matrix
*     from cov_parm (num_covxnum_cov) to covar (ndimxndim)

* PASSED VARIABLES
      integer*4 nv1, nv2 ! Row and column of start to move
      integer*4 nel      ! number of rows and columns to copy
      integer*4 ndim     ! Dimension of matrix being moved into
      integer*4 ncov     ! Dimension of main covaitance matrix

      real*8 covar(ndim,ndim)      ! Copy to matrix
      real*8 cov_parm(ncov,ncov)   ! Source matrix

      integer*4 i,j

      do i = 1, nel
         do j = 1, nel
            covar(i,j) = cov_parm(nv1+i-1,nv2+j-1)
         end do
      end do

*
      return
      end
          
CTITLE COMP_DCHI

      subroutine  comp_dchi(dNEU,cov_dNEU, dchi)

      implicit none 

*     Compute the change is chi**2 for the match (force north and
*     east to have the same value)
      real*8 dNEU(3)   ! Change in NEU
      real*8 cov_dNEU(3,3)   ! Covariance matrix
      real*8 dchi(3) 

      integer*4 i

*     Compute the sqrt of chi**2
      do i = 1,3
         if( cov_dNEU(i,i).gt.1.d-12 ) then
             dchi(i) = sqrt(dNEU(i)**2/cov_dNEU(i,i))
             if( dchi(i).gt.99.99d0 ) dchi(i) = 99.99d0
         else
             dchi(i) = 99.99d0
         endif
      end do
      dchi(1) = max(dchi(1),dchi(2))
      dchi(2) = dchi(1)

      return
      end


     







                 



               




