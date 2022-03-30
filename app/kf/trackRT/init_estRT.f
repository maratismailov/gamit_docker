      subroutine init_estRT

      implicit none

*     Routine to initialize the estimation covarinace and solution
*     vectors for positions and atmosphere delays (ambiquity parameters
*     are added and removed as needed).

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* LOCAL
      integer*4 i,j,np  ! Counters

*     Clear full space for covariance matrix
      do i = 1,max_parm
         sol_vecp(i) = 0
         sol_vecm(i) = 0
         do j = 1, max_parm
             cov_parmp(j,i) = 0
             cov_parmm(j,i) = 0
         end do
      end do

*     Starting counting and asigning paramters
      np = 0

*     Site parameters
      do i = 1, num_site
         do j = 1, 3
            site_parn(j,i) = 0
            if( apr_site(j,i).ne.0 ) then
                np = np + 1
                site_parn(j,i) = np
                cov_parmp(np,np) = apr_site(j,i)
            endif
         end do

*        See if atmospheric delay parameters are to be estimated
         atm_parn(i) = 0
         if( apr_atm(i).ne.0 ) then
             np = np + 1
             atm_parn(i) = np
             cov_parmp(np,np) = apr_atm(i)
         end if
      end do

****  Save the number of parameters
      num_parm = np
      non_amb_parm = np

      avbad_tot = 0   ! Total number of bad residuals being tracked

      return
      end

      subroutine reset_state(ep, mjd)

      implicit none

*     Routine to reset the state of some or all sites.  The reset
*     re-initializes the parameter estimates and the state of the
*     bias fixing.

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 ep
      real*8 mjd   ! Current time

* LOCAL
      integer*4 i, j, k, m, np
      integer*4 cnt(2)  ! Number of ambiguities reset:
                        ! 1 -- Fixed, 2 -- Free
      integer*4 date(5)
      real*8 sectag
      logical kbit


****  Loop over sites to see which ones are to be reset
      do i = 1,num_site
          if( reset(i) ) then
*            Reset the parameter estimates: Site position
             do j = 1,3
                if( site_parn(j,i).gt.0 ) then
                    np = site_parn(j,i)
                    do k = 1, num_parm
                       cov_parmp(k,np) = 0.d0
                       cov_parmp(np,k) = 0.d0
                    enddo
                    cov_parmp(np,np) = apr_site(j,i)
                    sol_vecp(np) = 0 
                 end if
             end do
*            Reset atmospheric estimate
             if( atm_parn(i).gt.0 ) then
                 np = atm_parn(i)
                 do k = 1, num_parm
                    cov_parmp(k,np) = 0.d0
                    cov_parmp(np,k) = 0.d0
                 enddo
                 cov_parmp(np,np) = apr_atm(i)
                 sol_vecp(np) = 0
             end if

*****        For all expect first, reset ambiguity state (all
*            ambigities at first site are set fixed
             cnt(:) = 0
             if( i.gt.1 ) then
                 do j = 1, num_ambs
                    if( bf_ents(1,j).eq.i ) then
*                      This ambiquity is for this site
                       if( kbit(bf_ents(5,j),2) ) then
                           cnt(1) = cnt(1)+1  ! Fixed
                       else
                           cnt(2) = cnt(2)+1  ! Free
                       end if
                       call reset_cslip( j, ep )

                    end if
                 end do
             end if
             call mjd_to_ymdhms(mjd,date,sectag)
             write(lus,120) ep, date, sectag, site_names(i), cnt
 120         format('RESET Epoch ',i6,' Time ',i4,4i3.2,1x,F6.3,
     .              ' Site ',a,' BFlags Fxd ',i3,' Free ',i3)
         end if
         reset(i) = .false.   ! Only reset once
      end do

      if( ep.ge.debug(5) .and. ep.le.debug(6) )
     .call print_parRT(6, ep, 'P', 'After Reset')

****  Thats all
      end

CTITLE RESET_CSLIP

      subroutine reset_cslip( na, ep )
                  
      implicit none

*     Routine to reset the state of cycle slips.  This includes
*     marking the ambiguity as not-fixed, WL averages back to
*     zero and re-fresh parameter estimate

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 na   ! Ambiguity being reset
     .,         ep   ! Epoch numbe of reset

* LOCAL
      integer*4 k, m    ! Loop counter
     .,         np    ! Parameter number

      call sbit(bf_ents(5,na),3,1)
      call sbit(bf_ents(5,na),2,0)  ! Set resolved bit to 0

      bf_ents(3,na) = ep
      bf_ents(4,na) = ep

*     Reset the parameter estimates for these parameter
      do m = 1,2
         if( amb_parn(m,na).gt.0 ) then
            np = amb_parn(m,na)
            do k = 1, num_parm
               cov_parmp(np,num_parm) = 0
               cov_parmp(num_parm,np) = 0
            end do
            sol_vecp(np) = 0
            cov_parmp(np,np) = 100.d0
         elseif( amb_parn(m,na).eq.0 ) then
            amb_parn(m,na) = -1  ! Cause new parameter to be added
            amb_save(m,na) = -1  ! Cause new parameter to be added
         end if
      end do

****  Now reset the WL-sums
      do m = 1,2
         WLS_sum(m,na) = 0
         WLS_sqr(m,na) = 0
       end do
       WLS_num(na) = 0

      return
      end 


      
