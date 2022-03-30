      subroutine remap_ambparn( ep, old_sv, new_sv, ref_am ) 

      implicit none

*     Subroutine to update WLS ambiquity and L1/L2 ambqiquity estimates

      include '../includes/const_param.h'
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 ep  ! Current epoch
      integer*4 old_sv(max_site)  ! The old reference satellite
     .,         new_sv(max_site)  ! New reference satellite (ideally the
                        ! ambiguity has been resolved and we don't
                        ! need to do anything
     .,         ref_am(max_site)  ! Ambiguity number for the new
                        ! reference satellite

* LOCAL
      integer*4 i,j,k   ! Loop counters
     .,    np_ref       ! parameter number of new reference SV
                        ! (if fixed then value is zero and nothing done)
     .,    ns, ne       ! Start and end parameter numbers for biases
                        ! at the site being processed
      logical rm_parm(max_parm)  ! Logical set true for each parameter
                        ! that needs to be remapped (ambiquity at site)

      real*8 csave(max_parm)  ! Saved column of covariance matrix to allow
                        ! fix in place
     .,      ssave      ! Solution vector save



****  Loop over all of the stations, updating the reference
*     SV at each station
      do i = 2,num_site
         if( old_sv(i).ne.new_sv(i) .and.
     .       ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .   write(*,120) ep, site_names(i), old_sv(i), new_sv(i), 
     .                ref_am(i), amb_parn(1,ref_am(i))
 120     format('REMAP_AMB Ep ',i6,1x,a,' Old/New G ',I2.2,' G ',I2.2,
     .       ' Ref Amb ',i3,' Amb param number ',I5)
      end do
      j = 0
      do i = 2,num_site
         if( amb_parn(1,ref_am(i)).gt.0 ) j = j + 1
      end do

      if( j.eq.0 ) RETURN   ! We don't need to remap covariances
                            ! Old and new are all fixed

      if( ep.ge.debug(3) .and. ep.le.debug(4) ) then
         call print_parRT(6, ep,'P','Before remap')
         do i = 1,num_parm
            write(*,900) 'CP', i,rm_parm(i),
     .                   (cov_parmp(j,i),j=1,num_parm)
         end do
      endif

      do i = 2, num_site

****     First see if the new sv is being estimated 
*        ******** (Should loop over neam
*        when L1 and L2 are separate). ***********
         if( ref_am(i).gt.0 ) then
            np_ref = amb_parn(1,ref_am(i))  ! Most likely no data
         else
            np_ref = 0
         end if
*        Set the array for parameters to mapped to false.  (Bias
*        parameters for current site will be set true below.)
         do j = 1, num_parm
            rm_parm(j) = .false.
         end do 

         if( np_ref.gt.0 ) then
*            This satellite at this site still has float ambiquity
*            Update the covariance matix and solution vector
*            Find range of parameters here
             ns = 10000
             ne = 0
             do k = 1, num_ambs
                if( bf_ents(1,k).eq.i .and. 
     .              amb_parn(1,k).gt.0 ) then  ! Correct site and estimate
                    rm_parm(amb_parn(1,k)) = .true.  ! Set to be remapped
                    if ( amb_parn(1,k).lt.ns  ) then
                       ns = amb_parn(1,k)
                    endif
                    if( amb_parn(1,k).gt.ne ) then
                       ne = amb_parn(1,k)
                    endif
                end if
             end do
             if( ep.ge.debug(3) .and. ep.le.debug(4) ) then
                 print *,'RM_PARM Site ',i,' L ',rm_parm(1:num_parm)
                 print *,'Range Site ',i,ns,ne, ' L ', rm_parm(ns:ne)
             end if
             do k = 1, num_ambs
                if( bf_ents(1,k).eq.i .and. 
     .              bf_ents(2,k).eq.old_sv(i) ) then
*                   This ambiguity should now point to parameter
*                   that was previously the new reference
                    amb_parn(1,k) = np_ref
                endif
                if( bf_ents(1,k).eq.i .and. 
     .              bf_ents(2,k).eq.new_sv(i) ) then
*                   This ambiguity should now point to parameter
*                   that was previously the new reference
                    amb_parn(1,k) = 0
                endif
             end do

*            ns and ne are the start and end parameter numbers of this
*            sites ambiquity parameters.
*            Save the reference value column
             do k = 1, num_parm
                csave(k) = cov_parmp(k,np_ref)
             end do
             ssave = sol_vecp(np_ref)
*            Update the sol_vec entries
             do k = ns,ne
                if( rm_parm(k) ) then  ! only change values to be remapped
                   if( k.ne.np_ref ) then
                      sol_vecp(k) = sol_vecp(k)-ssave
                   else
                      sol_vecp(k) = -ssave
                   endif
                endif
             end do

****         Now update the covariance matrix  
             do k = 1, num_parm
                do j = ns,ne
*                   Only do the j elements that are part of current
*                   ambiguities being remapped
                    if( rm_parm(j) ) then 
*                      Do parameters in the C D' block (ie. not
*                      the bias parameters at this site)
                       if( .not.rm_parm(k) ) then
                           if( j.ne.np_ref ) then
                              cov_parmp(k,j) = cov_parmp(k,j)-csave(k)
                              cov_parmp(j,k) = cov_parmp(k,j)
                           else
                              cov_parmp(k,np_ref) = -csave(k)
                              cov_parmp(np_ref,k) = cov_parmp(k,np_ref)
                           end if
                       else
*                          Now do the bias parameter block (D C D')
                           if( k.eq.np_ref .and. j.eq.np_ref ) then
                              cov_parmp(k,k) = csave(np_ref)
                           elseif ( k.eq.np_ref ) then
                              cov_parmp(k,j) = csave(np_ref) -
     .                                         cov_parmp(k,j)
                           elseif ( j.eq.np_ref ) then
                              cov_parmp(k,j) = csave(np_ref) -
     .                                         cov_parmp(k,j)
                           else
                              cov_parmp(k,j) =(cov_parmp(k,j)-csave(j))-
     .                                        (csave(k)-csave(np_ref))
                           endif
                       endif
                    end if
                end do
             end do
         end if
      end do
      if( ep.ge.debug(3) .and. ep.le.debug(4) ) then
         call print_parRT(6, ep,'P','After remap')
         do i = 1,num_parm
            write(*,900) 'CA', i,rm_parm(i),
     .                   (cov_parmp(j,i),j=1,num_parm)
 900        format(a,1x,i3,1x,L1,1x,100(F10.6,1x))
         end do
      end if


****  Thats all
      return
      end

   


                 


