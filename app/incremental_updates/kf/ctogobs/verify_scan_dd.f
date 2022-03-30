CTITLE VERIFY_SCAN_DD 
 
      subroutine verify_scan_dd( bep, ns, lv, js, kv, 
     .        L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, verified)

      implicit none
 
*     This routine will check whether a slip really seems to 
*     have happened by fitting linear trend to L1 and L2 with 
*     4 points before and after.  Verified is returned true if:
*     (1) There is not enough data for the fit
*     (2) RMS of the fits is two large (happens when there are
*         slips in rapid succession.
*     (3) If the difference in the L1 and L2 fit at the break
*         epoch exceeds dd_lc_tol(2)*2 (normally 0.7 cycles)

 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv, js, kv  - First station, first satellite, second station
*                     second statellite in double difference 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   bep             - The epoch at which the bias flag should be added
*                     if one needs to be added.
 
      integer*4 ns, lv, js, kv, bep,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2

 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep)

*   verified        - Set true if there is a break or we can't
*                     determine if there is

      logical verified
 
* LOCAL PARAMETERS
 
 
* LOCAL VARIABLES

*     >> Autcln_fit
*     >> Chi =
*         0.3000   -0.4000   -0.1000    0.2000
*        -0.4000    0.7000   -0.2000   -0.1000
*        -0.1000   -0.2000    0.7000   -0.4000
*         0.2000   -0.1000   -0.4000    0.3000
*     Left estimate: [ -4 -3 -2 -1 ] 
*     >> EstL(1,:)
*        -0.5000         0    0.5000    1.0000
*     Right Estimate [  0  1  2  3 ]
*     >> EstR(1,:)
*         0.7000    0.4000    0.1000   -0.2000

      real*8 EstL(4), EstR(4)  ! Vector to muliply dL1 and dL2 by
                               ! get offset estimate at bias epoch
      real*8 ChiMat(4,4)       ! Matrix to muliply dL1 and dL2 by
                               ! to generate sum of residuals squared.

*   c11, c12, c21, c22 - Channel numbers for each one-way in
*                 double difference
*   ltoc        - Function to return channel number for a 
*                 site, satellite list and epoch.
*   ep          - Epoch for data
*   i           - Loop counter

      integer*4 c11, c12, c21, c22, ltoc, ep,  i

      logical kbit   ! Test bit.
 
 
*   dd  - Double dofferece (tmp)
*   dd_L1(4), dd_L2(4) - Double differences for L1 and L2 at the four
*         epoch fit (-4,-3,-2,1) or (0,1,2,3)

*   tol         - Value that the tested changes in the double differnces
*                 must be less than or will be flagged bad.
 
      real*8 dd, dd_L1(4), dd_L2(4), tol

* Estimated values
      real*8 L1L, L2L  ! :Left estimates at bias epoch
      real*8 V1L, V2L  !  sum residual**2/2 for left L1L and L2L

      real*8 L1R, L2R  ! :Right estimates at bias epoch
      real*8 V1R, V2R  ! sum residual**2/2 for right L1L and L2L

      real*8 Ref_L1, Ref_L2  ! Offset values at L1 and L2 to keep values
                       ! near zero.


****  Data statements with matrices needed. Generated with Autcln_fit.m
      data EstR   /  0.7000,  0.4000,  0.1000,  -0.2000 /
      data EstL   / -0.5000,  0.0000,  0.5000,   1.0000 / 
      data ChiMat /  0.3000, -0.4000, -0.1000,   0.2000, 
     .              -0.4000,  0.7000, -0.2000,  -0.1000,
     .              -0.1000, -0.2000,  0.7000,  -0.4000, 
     .               0.2000, -0.1000, -0.4000,   0.3000 /


****  Start: Make sure we can do this
      verified = .true.
      if( bep.le.5 ) RETURN
      if( bep.gt.num_ep-2 ) RETURN

****  Start with right hand side and save the offset at bep
*     for L1 and L2 to remove from other values (to keep numerical
*     errors small):
      do ep = bep,bep+3
         c11 = ltoc(ctol_cse(1,ns,ep), lv, actual_max_chan)
         c12 = ltoc(ctol_cse(1,ns,ep), kv, actual_max_chan)
         c21 = ltoc(ctol_cse(1,js,ep), lv, actual_max_chan)
         c22 = ltoc(ctol_cse(1,js,ep), kv, actual_max_chan)
 
         If( c11.gt.0 .and. c12.gt.0 .and.
     .       c21.gt.0 .and. c22.gt.0   ) then
            dd   = (L1r_phs_cse(c11,ns,ep)+
     .                 L1_cyc_cse(c11,ns,ep)) -
     .             (L1r_phs_cse(c12,ns,ep)+
     .                 L1_cyc_cse(c12,ns,ep)) -
     .             (L1r_phs_cse(c21,js,ep)+
     .                 L1_cyc_cse(c21,js,ep)) +
     .             (L1r_phs_cse(c22,js,ep)+
     .                 L1_cyc_cse(c22,js,ep))
           if( ep.eq.bep)  Ref_L1 = dd
           dd_L1(ep-bep+1) = dd - Ref_L1
*          L2 Phase double difference
* MOD TAH 990512: Make sure we have L2 everywhere:
           if( L2r_phs_cse(c11,ns,ep).ne.0.d0 .and.
     .         L2r_phs_cse(c12,ns,ep).ne.0.d0 .and.
     .         L2r_phs_cse(c21,js,ep).ne.0.d0 .and.
     .         L2r_phs_cse(c22,js,ep).ne.0.d0  ) then
               dd =   (L2r_phs_cse(c11,ns,ep)+
     .                 L2_cyc_cse(c11,ns,ep)) -
     .                (L2r_phs_cse(c12,ns,ep)+
     .                 L2_cyc_cse(c12,ns,ep)) -
     .                (L2r_phs_cse(c21,js,ep)+
     .                 L2_cyc_cse(c21,js,ep)) +
     .                (L2r_phs_cse(c22,js,ep)+
     .                 L2_cyc_cse(c22,js,ep))
               if( ep.eq.bep ) Ref_L2 = dd
               dd_L2(ep-bep+1) = dd - Ref_L2
           else
               dd_L2(ep-bep+1) = 0    ! For L1 only data case
           end if
         else
*****      Missing data so let the bias flag get set 
*          (gap-py data that will probably be deleted)
           RETURN
         end if
      end do

****  OK; Get Left estimate of L1L and L2L at bias epoch
      L1R = dot_product(EstR,dd_L1)
      L2R = dot_product(EstR,dd_L2)
      V1R = dot_product(dd_L1,matmul(ChiMat,dd_L1))/2
      V2R = dot_product(dd_L2,matmul(ChiMat,dd_L2))/2


****  OK see of we can fit data before: Generate left hand data set
      do ep = bep-4,bep-1
         c11 = ltoc(ctol_cse(1,ns,ep), lv, actual_max_chan)
         c12 = ltoc(ctol_cse(1,ns,ep), kv, actual_max_chan)
         c21 = ltoc(ctol_cse(1,js,ep), lv, actual_max_chan)
         c22 = ltoc(ctol_cse(1,js,ep), kv, actual_max_chan)
 
         If( c11.gt.0 .and. c12.gt.0 .and.
     .      c21.gt.0 .and. c22.gt.0   ) then
            dd   = (L1r_phs_cse(c11,ns,ep)+
     .                 L1_cyc_cse(c11,ns,ep)) -
     .             (L1r_phs_cse(c12,ns,ep)+
     .                 L1_cyc_cse(c12,ns,ep)) -
     .             (L1r_phs_cse(c21,js,ep)+
     .                 L1_cyc_cse(c21,js,ep)) +
     .             (L1r_phs_cse(c22,js,ep)+
     .                 L1_cyc_cse(c22,js,ep))
           dd_L1(ep-bep+5) = dd - Ref_L1
*          L2 Phase double difference
* MOD TAH 990512: Make sure we have L2 everywhere:
           if( L2r_phs_cse(c11,ns,ep).ne.0.d0 .and.
     .         L2r_phs_cse(c12,ns,ep).ne.0.d0 .and.
     .         L2r_phs_cse(c21,js,ep).ne.0.d0 .and.
     .         L2r_phs_cse(c22,js,ep).ne.0.d0  ) then
               dd =   (L2r_phs_cse(c11,ns,ep)+
     .                 L2_cyc_cse(c11,ns,ep)) -
     .                (L2r_phs_cse(c12,ns,ep)+
     .                 L2_cyc_cse(c12,ns,ep)) -
     .                (L2r_phs_cse(c21,js,ep)+
     .                 L2_cyc_cse(c21,js,ep)) +
     .                (L2r_phs_cse(c22,js,ep)+
     .                 L2_cyc_cse(c22,js,ep))
               dd_L2(ep-bep+5) = dd - Ref_L2
           else
               dd_L2(ep-bep+5) = 0    ! For L1 only data case
           end if
         else
*****      Missing data so let the bias flag get set 
*          (gap-py data that will probably be deleted)
           RETURN
         end if
      end do

****  OK; Get Left estimate of L1L and L2L at bias epoch
      L1L = dot_product(EstL,dd_L1)
      L2L = dot_product(EstL,dd_L2)
      V1L = dot_product(dd_L1,matmul(ChiMat,dd_L1))/2
      V2L = dot_product(dd_L2,matmul(ChiMat,dd_L2))/2

****  Now see if break is real
      if( abs(L1R-L1L).lt.2*dd_lc_tol(2) .and.
     .    abs(L2R-L2L).lt.2*dd_lc_tol(2) .and.
     .    sqrt(V1R+V1L).lt. dd_lc_tol(2) .and.
     .    sqrt(V2R+V2L).lt. dd_lc_tol(2) ) then
          verified = .false.
      endif
      if( kbit(status_rep,10) ) 
     .write(*,220) verified, bep, L1R, L1L, L2R, L2L, 
     .    sqrt(V1R), sqrt(V1L), sqrt(V2R), sqrt(V2L)
 220  format(' VERIFIED ',L1,' EP ',I5,' RL: L1 ',2F6.2,
     .       ' L2 ',2F6.2,' S1 ',2F6.2,' S2 ',2F6.2) 

      
      end

