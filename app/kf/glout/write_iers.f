CTITLE WRITE_IERS

      subroutine write_iers( iout, cov_parm, sol_parm )

      implicit none 

*     Routine to write out entries in the IERS/IGS polar motion/UT1
*     submission format.  The lines are preceeded by IERS for easy
*     grep'ing from the output files 

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'

* PASSED VARIABLES 

* iout  -- Output unit number
* cov_parm(num_glb_parn,num_glb_parn) -- Covariance matrix
* sol_parm(num_glb_parn) -- Solution adjustments
 
      integer*4 iout
      real*8 cov_parm(num_glb_parn,num_glb_parn), 
     .       sol_parm(num_glb_parn) 

* LOCAL VARIABLES

* nf -- Number of fiducial sites
* ns -- Number of sites
* i,j  -- Loop counters
* nxp(4) -- Pole position parameter numbers
* nut(2) -- UT1 and LOD parameter numbers 

      integer*4 ns, nf, nt, i,j,k,l, nxp(4), nut(2)

* xpr(4) -- Pole position and rates (uas and uas/d)
* spr(4) -- Sigmas for pole position and rates
* utr(2) -- UT1 and rate (sign changed to LOD at output)
* sut(2) -- Sigmas for UT1 and LOD
* rho(3) -- correlations between xy, xut1 and yut1

* sp     -- Temp variable for computing position constraint
* mjd    -- Modified Julian date for PMU value

      real*8 xpr(4), spr(4), utr(2), sut(2), rho(3), sp,
     .       mjd

* kbit   -- Checks is bit set

      logical kbit 

* dq    -- Double quote character

      character*1 dq

      data dq / '"' /

* Start counting out the information that we need
* Get the number of stations constrained
* MOD TAH 120913: Change compuation of NF to number used in 
*     defining the reference frame.
      nf = 0
      ns = 0
      do i = 1, gnum_sites
!         sp = sqrt(apr_neu(1,1,i)**2  + apr_neu(2,1,i)**2  + 
!     .             apr_neu(3,1,i)**2  +
!     .             apr_site(1,1,i)**2 + apr_site(2,1,i)**2 +
!     .             apr_site(3,1,i)**2   ) 
!         if( sp.lt.0.030d0 .and. kbit(guse_site,i) ) nf = nf + 1
          if( kbit(use_pos,i)   )  nf = nf + 1
          if( kbit(guse_site,i) )  ns = ns + 1
      end do
      write(iout,120) dq, dq, dq, dq
 120  format('IERS  MJD      Xpole   Ypole  UT1-UTC    LOD  Xsig',
     .       '  Ysig   UTsig  LODsig  Nr  Nf  Nt    Xrt    Yrt',
     .       '  Xrtsig  Yrtsig XYcorr XUTcor YUTcor',/,
     .       'IERS             (10**-6',a1,')       (0.1 usec)',
     .       '    (10**-6',a1,')     (0.1 usec)           ',
     .       '   (10**-6',a1,'/d)    (10**-6',a1,'/d)' )
*     Now do the output
      if( num_mul_pmu.eq.0 ) then

* NOTE: If midpt output, then num_mul_pmu will not be zero and 
*       gepoch_expt will equal gepoch_out when num_mul_pmu==0 

*         Extract the pole position and rates
          do j = 1,4
              nxp(j) = parn_wob(j)
              if( nxp(j).gt.0 ) then
                  xpr(j) = apr_val_wob(j) + gwob_apr(j) +
     .                                      sol_parm(nxp(j))
                  spr(j) = sqrt(cov_parm(nxp(j),nxp(j)))*1000.d0 
              else
                  xpr(j) = apr_val_wob(j) + gwob_apr(j)
                  spr(j) = 0.d0
              end if
*             Convert results to mico-arc-seconds
              xpr(j) = xpr(j)*1000.d0
          end do

*         Now do UT1 and LOD
          do j = 1,2
             nut(j) = parn_ut1(j)
             if( nut(j).gt.0 ) then
                 utr(j) = apr_val_ut1(j) + gut1_apr(j) +
     .                                     sol_parm(nut(j))
                 sut(j) = sqrt(cov_parm(nut(j),nut(j)))/15.d0*1.d4
             else
                 utr(j) = apr_val_ut1(j) + gut1_apr(j)
                 sut(j) = 0.d0
             endif

*            Remove leap seconds for UT1 if needed and convert the
*            rate to LOD
             if( j.eq.1 ) then
                 call get_leapsec( gepoch_expt, utr(1), gtai_utc ) 
                 utr(1) = (utr(1)/15.d0 - gtai_utc*1000)*1.d4
             else
                 utr(2) = -utr(2)/15.d0*1.d4
             endif
          end do 

*         Now get the correlations
          if( nxp(1).ne.0    .and. nxp(2).ne.0 .and. 
     .        spr(1).gt.0.d0 .and. spr(2).gt.0.d0  ) then
              rho(1) = cov_parm(nxp(1),nxp(2))/
     .            sqrt(cov_parm(nxp(1),nxp(1))*cov_parm(nxp(2),nxp(2)))
          else
              rho(1) = 0.d0
          endif
          if( nxp(1).ne.0    .and. nut(1).ne.0 .and.
     .        spr(1).gt.0.d0 .and. sut(1).gt.0.d0 ) then
              rho(2) = cov_parm(nxp(1),nut(1))/
     .            sqrt(cov_parm(nxp(1),nxp(1))*cov_parm(nut(1),nut(1)))
          else
              rho(2) = 0.d0
          endif
          if( nxp(2).ne.0    .and. nut(1).ne.0 .and.
     .        spr(2).gt.0.d0 .and. sut(1).gt.0.d0  ) then
              rho(3) = cov_parm(nxp(2),nut(1))/
     .            sqrt(cov_parm(nxp(2),nxp(2))*cov_parm(nut(1),nut(1)))
          else
              rho(3) = 0.d0
          endif 

****      Now write out the line
          mjd = gepoch_expt - 2400000.5d0
          write(iout,150) mjd, nint(xpr(1)),nint(xpr(2)),nint(utr(1)),
     .                    nint(utr(2)),nint(spr(1)),nint(spr(2)),
     .                    nint(sut(1)), nint(sut(2)), 
     .                    ns,nf,gnum_svs,
     .                    nint(xpr(3)),nint(xpr(4)),
     .                    nint(spr(3)),nint(spr(4)), rho
 150      format('IERS',F8.2,2I8,I9,I7,2I6,2I8,I4,2I4,4I7,2x,3F7.3)

      else 

*         Loop over the multi-day PMU values
          do k = 1, num_mul_pmu
*             Extract the pole position and rates
              do j = 1,2
                  do i = 1, 2
                     l = 2*(i-1) + j
                     nxp(l) = parn_mul_pmu(i,j,k)
                     if( nxp(l).gt.0 ) then
                         xpr(l) = apr_val_mul_pmu(i,j,k)+
     .                                          sol_parm(nxp(l))
                         spr(l) = sqrt(cov_parm(nxp(l),nxp(l)))*1000.d0 
                     else
                         xpr(l) = apr_val_mul_pmu(i,j,k)
                         spr(l) = 0.d0
                     end if
*                    Convert results to mico-arc-seconds
                     xpr(l) = xpr(l)*1000.d0
                  end do
              end do

*             Now do UT1 and LOD
              do i = 1,2
                 nut(i) = parn_mul_pmu(i,3,k)
                 if( nut(i).gt.0 ) then
                     utr(i) = apr_val_mul_pmu(i,3,k) +
     .                                     sol_parm(nut(i))
                     sut(i) = sqrt(cov_parm(nut(i),nut(i)))/15.d0*1.d4
                 else
                     utr(i) = apr_val_mul_pmu(i,3,k)
                     sut(i) = 0.d0
                 endif

*                Remove leap seconds for UT1 if needed and convert the
*                rate to LOD
                 if( i.eq.1 ) then
                     call get_leapsec( gmul_pmu_ep(k), utr(1), gtai_utc) 
                     utr(1) = (utr(1)/15.d0 - gtai_utc*1000)*1.d4
                 else
                     utr(2) = -utr(2)/15.d0*1.d4
                 endif
              end do 

*             Now get the correlations
              if( nxp(1).ne.0    .and. nxp(2).ne.0 .and. 
     .            spr(1).gt.0.d0 .and. spr(2).gt.0.d0  ) then
                  rho(1) = cov_parm(nxp(1),nxp(2))/
     .                sqrt(cov_parm(nxp(1),nxp(1))*
     .                     cov_parm(nxp(2),nxp(2)))
              else
                  rho(1) = 0.d0
              endif
              if( nxp(1).ne.0    .and. nut(1).ne.0 .and.
     .            spr(1).gt.0.d0 .and. sut(1).gt.0.d0 ) then
                  rho(2) = cov_parm(nxp(1),nut(1))/
     .                sqrt(cov_parm(nxp(1),nxp(1))*
     .                     cov_parm(nut(1),nut(1)))
              else
                  rho(2) = 0.d0
              endif
              if( nxp(2).ne.0    .and. nut(1).ne.0 .and.
     .            spr(2).gt.0.d0 .and. sut(1).gt.0.d0  ) then
                  rho(3) = cov_parm(nxp(2),nut(1))/
     .                sqrt(cov_parm(nxp(2),nxp(2))*
     .                     cov_parm(nut(1),nut(1)))
              else
                  rho(3) = 0.d0
              endif

****          Now count up number of stations and satellites
              ns = 0
              nf = 0
              do i = 1, gnum_sites
                 if( kbit(mul_pmu_site(1,k),i) ) ns = ns + 1 
                 if( kbit(mul_pmu_fidu(1,k),i) ) nf = nf + 1
              end do
              nt = 0
              do i = 1, gnum_svs
                 if( kbit(mul_pmu_svs(k),i) ) nt = nt + 1
              end do
  
****          Now write out the line
              mjd = gmul_pmu_ep(k) - 2400000.5d0
              write(iout,150) mjd, nint(xpr(1)),nint(xpr(2)),
     .                        nint(utr(1)),
     .                        nint(utr(2)),nint(spr(1)),nint(spr(2)),
     .                        nint(sut(1)), nint(sut(2)), 
     .                        ns,nf,nt,
     .                        nint(xpr(3)),nint(xpr(4)),
     .                        nint(spr(3)),nint(spr(4)), rho
          enddo
          write(iout,'(1x)')

      endif

            
      end

      



         





