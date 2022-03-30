      subroutine NLx_formneq(lv, dexwl, dmwwl, sig_awl, sig_bwl, 
     .                       NLx_neq, NLx_bvc)

      implicit none

*     Routine to form normal equations for extimating cycle slips from
*     the MW-WL and EX-WK

      include 'track_com.h'   

* PASSED
* INPUT VALUES 
      integer*4 lv  ! Satellite number.  Needed to get frequency

      real*8 dexwl(1), dmwwl  ! Differences in EX and MW widelines (after-before)
     .,      sig_awl(2), sig_bwl(2)   ! Sigma of after and and before EX(1) and 
                                      ! MW (2) widelanes
* RETURNED VALUES
      real*8 NLx_neq(2,2)  ! Normal equations for L1/L2 cycle estimates 
     .,      NLx_bvc(2)    ! B-vector for estimate


* LOCAL
      real*8 Wght         ! Weight for each estimate (1/sigma^2)
     .,      WLap(1,2)    !  Parials array for EX-WL


*     Start with MW-WL
      Wght = 1.d0/(sig_bwl(2)**2+sig_awl(2)**2)
      NLx_neq = reshape((/ Wght , -Wght , -Wght, Wght /),shape(NLx_neq))
      NLx_bvc = (/ -Wght*dmwwl , Wght*dmwwl /)

*     Add the EX-WL
      Wght = 1.d0/(sig_bwl(1)**2+sig_awl(1)**2)
* MOD AZ 190305: alter the signs of these two terms,
*                which I believe are not right
      WLap = reshape((/ -fL1(lv)/fR1, fL2(lv)/fR2*fL1(lv)/fL2(lv) /),
     .       shape(WLap))
      NLX_neq = NLx_neq + matmul(transpose(WLap),WLap*Wght)
      NLx_bvc = NLx_bvc + matmul(transpose(WLap),dexwl*Wght)

      write(*,320) NLX_neq, NLx_bvc, dexwl(1), dmwwl, sig_awl, sig_bwl
 320  format('NLX Neq ',4F10.3,' Bvc ',2F10.3,' DWL ',2F10.3,
     .       ' Asig ',2F10.3,' Bsig ',2F10.3)
 
****  Thats all
      return 
      end

CTILE NLx_est   
    
      subroutine NLx_est(lv, NLx_neq, NLx_bvc, damb, csfix)

      implicit none

*     Routine to finish estimate and judge robustness based
*     on chi^2 calculations with full covariance matrix.

      include 'track_com.h'   

* PASSED
* INPUT VALUES 
      integer*4 lv  ! Satellite number.  Needed to get frequency

      real*8 NLx_neq(2,2)  ! Normal equations for L1/L2 cycle estimates 
     .,      NLx_bvc(2,1)    ! B-vector for estimate

      real*8 dexwl(1), dmwwl  ! Differences in EX and MW widelines (after-before)
     .,      sig_awl(2), sig_bwl(2)   ! Sigma of after and and before EX(1) and 
                                      ! MW (2) widelanes
* RETURNED VALUES
      real*8 damb(2)    ! L1 / L2 ambiquity.  Will be non-integer for non-GPS
                        ! GNSS processing.
      logical csfix     ! Set true if chi^2 ratio greater than relrank_limit

* LOCAL
      real*8 NLx_inv(2,2)  ! Inverse of normal equations
     .,      dsol(2,1)     ! Solution vector with float estimates of L1 and L2 cycles.
     .,      dsig(2)       ! Sigma for estimates of dsol 
     .,      dL2s          ! L1-L2 sigma (used to set range)
     .,      dchi_min, dchi_nxt  ! Smalltest ad next smallest chi**2
     .,      dchi          ! Change in chi**2 with test T1 and T2 cycles
     .,      det           ! Determinate for 2x2 inversion
     .,      dtest(2)    ! Test values for difference between estimate and test
                           ! ambguities.

      integer*4 dNLx(2)  ! Integer ambuguities with smallest chi**2
     .,     t1, t2         ! Test integer values to find smallest chi**2
     .,     NL1, NL2       ! First estimates of ambiguities.  
     .,     i, j          ! Loop counters
* MOD AZ 190305: additional variables for 0.5 step of search
      real*8 t1_start, t1_end   ! outer loop counter
     .,     t2_start, t2_end    ! inner loop counter

****  Invert the normal equations
      det = 1.d0/(NLx_neq(1,1)*NLx_neq(2,2)- NLx_neq(2,1)*NLx_neq(1,2))
      NLx_inv = reshape((/ NLx_neq(2,2) ,-NLx_neq(2,1) , 
     .     -NLx_neq(1,2) , NLx_neq(1,1) /),shape(NLx_neq))*det
      dsol = matmul(NLx_inv, NLx_bvc)

      NL1 = nint(dsol(1,1)*fL1(lv)/fR1) 
      NL2 = nint(dsol(2,1)*fL2(lv)/fR2)
      dsig = (/ sqrt(NLx_inv(1,1)),sqrt(NLx_inv(2,2)) /) 
      dL2s = sqrt(NLx_inv(1,1) + NLx_inv(2,2) -2*NLx_inv(1,2))
      write(*,220) dsol, dsig,  dL2s,  NL1, NL2
 220  format('NLX DSOL dL1/2 ',2F10.3,' +- ',3F6.3,' NL1/2 ',2I3,' cyc')

      dchi_min = 1.d20
      dchi_nxt = 1.d20
      
      t1_start = NL1-2*int(dsig(1)+1)
      t1_end = NL1+2*int(dsig(1)+1)
* MOD AZ 190305: Shorten the search step from 1 to 0.5
      do while (t1_start.le.t1_end)
         t2_start = t1_start-2*int(dL2s+1)
         t2_end = t1_start+2*int(dL2s+1)
         do while (t2_start.le.t2_end)
  
            dtest = (/ dsol(1,1) - t1_start*fR1/fL1(lv) , 
     .                 dsol(2,1) - t2_start*fR2/fL2(lv) /)
            dchi = 0.d0
            do i = 1, 2
               do j = 1,2 
                  dchi = dchi + dtest(i)*NLx_neq(i,j)*dtest(j)
               end do
            end do

            if( dchi.lt.dchi_min ) then
                dchi_nxt = dchi_min
                dchi_min = dchi
                dNLx = (/ t1_start, t2_start /)
            else if( dchi.lt.dchi_nxt ) then 
                dchi_nxt = dchi 
            endif

            t2_start = t2_start+0.5

         end do
        
         t1_start = t1_start+0.5  

      end do
*     Now scan to see minimum chi**2
!      do t1 = NL1-2*int(dsig(1)+1), NL1+2*int(dsig(1)+1)
!         do t2 =  t1-2*int(dL2s+1), t1+2*int(dL2s+1)
!            dtest = (/ dsol(1,1) - t1*fR1/fL1(lv) , 
!     .                 dsol(2,1) - t2*fR2/fL2(lv) /)
!            dchi = 0.d0
!            do i = 1, 2
!               do j = 1,2 
!                  dchi = dchi + dtest(i)*NLx_neq(i,j)*dtest(j)
!               end do
!            end do
!
!            if( dchi.lt.dchi_min ) then
!                dchi_nxt = dchi_min
!                dchi_min = dchi
!                dNLx = (/ t1, t2 /)
!            else if( dchi.lt.dchi_nxt ) then 
!                dchi_nxt = dchi 
!            endif
C           write(*,320) t1,t2, dtest, dchi, dchi_min, dNLx
C320        format('NLX L1 L2 ',2I3,' Err ',2F10.3,' dChi ',2F10.3,ne 
C    .             ' Best ',2I3)
!         end do
!      end do

      write(*,320) dchi_min, dchi_nxt, NLx_RR, NLx_addChi,  NLx_minChi
 320  format('NLx Chis Min ',F10.2,' Nxt ', F10.2,' Tols ',3F8.2)
      if( dchi_nxt/max(dchi_min,NLx_addChi).gt.NLx_RR   .and. 
     .    dchi_min.lt. NLx_minChi ) then
          csfix = .true.
      else
          csfix = .false.
      endif

      damb = (/ dNLx(1)*fR1/fL1(lv) , dNLx(2)*fR2/fL2(lv) /)

      write(*,340) csfix, dchi_nxt/max(dchi_min,1.d0), NLx_RR, 
     .             NLx_addChi,  NLx_minChi, damb 
 340  format('NLX Fix ',L1,' Ratio ', F10.2, ' RR Limit ',F10.2,
     .       ' Add/Min Chi ',2F6.2,' dL1/2 ',2F10.5)

****  Thats all
      return 
      end
 
 

