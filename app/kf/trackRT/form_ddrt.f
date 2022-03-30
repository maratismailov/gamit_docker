      subroutine form_ddrt(na, comb_coeff, dd_obs)

      implicit none

*     Routine to form double differenc o-minus-c for the current
*     ambiguity (na)

      include '../includes/const_param.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 na    ! Current double difference combination being
                      ! processed.

      real*8 comb_coeff(4)  ! Combination of one-way data 
      real*8 dd_obs         ! OMC double difference to be used in KF

* LOCAL
      integer*4 i, j   ! Loop counters
     .,   ob, ns, pn   ! Obs number, site and PRN

      real*8 elev      ! Elevation angle (rad)
     .,   ow_dat(4)    ! One-way data combination that is then differenced
     .,   obs_var(4)   ! Variance of L1, L2 P1 and P2 observable for
                       ! each one-way
      
****  OK form the four one-way combinations.  First we form the
*     information for each one-way (i.e. station/satellite 
*     combination)
      do i = 1, 4
         ob = wls_obn(i,na)              ! Observation number
         ns = RT_sitenum(wls_obn(i,na))  ! Site
         pn = RT_satnum(wls_obn(i,na))   ! Satellite number
*        NOTE: Ambiquities were removed in get_bflag routine and
*        do not need to be removed from the phase here.
         ow_dat(i) = comb_coeff(1)*RT_omcraw(1,wls_obn(i,na))
     .             + comb_coeff(2)*RT_omcraw(2,wls_obn(i,na))
     .             + comb_coeff(3)*RT_omcraw(3,wls_obn(i,na))
     .             + comb_coeff(4)*RT_omcraw(4,wls_obn(i,na)) 

*        This code could be in its own routine looped over the
*        one ways (calcs are repeated here).
         do j = 1,4  ! Loop over L1, L2, P1, P2 to get noise
            elev = RT_azel(2,wls_obn(i,na))*pi/180 ! Rads
            obs_var(j) = data_var(j,pn)*
     .         (1.d0+(data_var(5,pn)/sin(elev)**2)) 
         end do

*        Now from one-way variance (saved in common)
         ow_var(ob) = comb_coeff(1)**2*obs_var(1)
     .              + comb_coeff(2)**2*obs_var(2)
     .              + comb_coeff(3)**2*obs_var(3)
     .              + comb_coeff(4)**2*obs_var(4)
C         print *,'VARS ',na, ob, ns, pn, elev, data_var(:,pn),
C     .       obs_var(:),ow_var(ob)      
      end do

****  Now form the double difference observable
      dd_obs = (ow_dat(1)-ow_dat(2))-(ow_dat(3)-ow_dat(4))

****  Thats all
      return
      end


      subroutine form_apart(na, part_fact, nd )

      implicit none

*     Routine to form the double difference site and atmospheric
*     delay paritials.

      include '../includes/const_param.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 na    ! Current double difference combination being
                      ! processed.
     .,         nd    ! Double difference being processed.

      real*8 part_fact      ! Factor to convert range partial to 
                      ! observable type (1/lambda for phase)

* LOCAL
      integer*4 i     ! Loop counter

****  Form the partial by differencing the partials of all 
*     estimated parameters
      do i = 1, non_amb_parm
         apart(i,nd) = ((ow_part(i,wls_obn(1,na))
     .                  -ow_part(i,wls_obn(2,na)))
     .                 -(ow_part(i,wls_obn(3,na))
     .                  -ow_part(i,wls_obn(4,na))))*part_fact
      end do

*     Clear the rest of the partials array
      do i = non_amb_parm+1, num_parm
          apart(i,nd) = 0
      end do

****  Thats all
      return 
      end

      subroutine form_covdd(ep, nd )

      implicit none

*     Routine to form the double difference covarinace matrix
*     for the data

      include '../includes/const_param.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 ep         ! Epoch
      integer*4 nd         ! Data type counter.  Whnen nd=1 we
                           ! clear cov_obs for num_dd*num_anal_type

* LOCAL
      integer*4 i,j, k, l  ! Loop counters

      integer*4 dd(4)      ! Double difference operator

      data  dd  / 1, -1, -1, 1 /
 
***** Clear the off-diaongal parts of matrix
      if( nd.eq.1 ) then   ! Clear cov obs
         do i = 1, num_anal_type*num_dd
            do j = 1, num_anal_type*num_dd
               cov_obs(i,j) = 0
            end do
         end do
      end if

***** Loop over dd's
      do i = 1, num_ambs
        if( amb_dd(i).ne.0 ) then
           do j = i, num_ambs
              if( amb_dd(j).ne. 0 ) then
c                cov_obs(amb_dd(i),amb_dd(j)) = 0
                 do k = 1,4
                    do l = 1,4 
                        if( wls_obn(k,i).eq.wls_obn(l,j) ) then
                            cov_obs(amb_dd(i),amb_dd(j)) =
     .                          cov_obs(amb_dd(i),amb_dd(j)) +
     .                          dd(k)*dd(l)*ow_var(wls_obn(k,i))
C                           if( k.ne.l ) then
C                               write(*,220) ep, i,j, k,l, 
C    .                             amb_dd(i), amb_dd(j),
C    .                             wls_obn(:,i), wls_obn(:,j) 
C220                            format('DD_NEG Ep ',i6,' Ambs ',2i4,
C    .                                 ' IDX ',2i2,' DD ',2i4,
C    .                                 ' OBN ',2(4i4,2x))
C                           end if
                        end if
                    end do
                 end do
                 cov_obs(amb_dd(j),amb_dd(i)) =
     .                cov_obs(amb_dd(i),amb_dd(j)) 
              end if
           end do
        end if
      end do

***** Thats all
      return
      end

 
      subroutine reset_ambest(ep, na)          

      implicit none

*     Routine to reste the estimate and variance of a bias parameter
*     with a cycle slip

      include '../includes/const_param.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 ep   ! Epoch number
     .,     na       ! Ambiquity number

* LOCAL
      integer*4 np   ! Parameter number
     .,     i, j     ! Loop counters


      do j = 1, neam
         np = amb_parn(j,na)
         if( np.gt.0 ) then
            do i = 1, num_parm
               cov_parmp(i,np) = 0
               cov_parmp(np,i) = 0
            end do
            cov_parmp(np,np) = 100.d0
            sol_vecp(np) = 0.d0
         endif
      end do

      return
      end


        
      



