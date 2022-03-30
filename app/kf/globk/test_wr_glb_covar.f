      program test_wr_glb_covar

      implicit none

*     Test writing >32767 parameter matric
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/globk_markov.h'

      integer*4 mnp

      parameter ( mnp = 36000 )

      real*8 cov_parm(mnp,mnp), sol_parm(mnp)

      common / big / cov_parm, sol_parm

*   dcb(1)  - Dcb buffer

      integer*4 dcb(16), ierr, i,j, ic
      integer*8 len_run, rcpar

      real*8 errc, errs

      integer*8 I8

      character*(1) rw    ! Passed set to W to write, R to read

      data I8 / 1 /

*     See if reading or writing
      len_run = rcpar(1,rw)
      if( len_run.eq.0 ) then
          print *,'% test_wr_glb_covar <W/R>'
          stop 'No option'
      endif
*     Set name of sol file 
      glb_sol_file = 'test_wr_glb_covar.sol' 
      num_glb_parn = mnp 
      call dwint8( 0.d0, cov_parm, 1, (num_glb_parn+I8)*num_glb_parn )
      if ( rw.eq.'W' ) then
         print *,'Creating ',trim(glb_sol_file)
         do i = 1,mnp
            do j = -2,2
               ic = max(i+j,1)
               if( ic.gt.mnp ) ic = mnp
               cov_parm(i,ic) = (i-1)*mnp + j
               sol_parm(i) = -i
C              print *,'Set ',i,j,ic, cov_parm(i,ic),sol_parm(i), 
C    .                                cov_parm(i,mnp+1)
            end do
         end do
         print *,'Opening file'
         call rw_glb_covar('C', dcb, cov_parm)
         print *,'Writing file'
         call rw_glb_covar('W', dcb, cov_parm)
      else
         print *,'Reading ',trim(glb_sol_file)
         print *,'Opening file'
         call rw_glb_covar('O', dcb, cov_parm)
         print *,'Reading file'
         call rw_glb_covar('L', dcb, cov_parm)
         print *,'Checking file'
         do i = 1,mnp
            do j = -2,2
               ic = max(i+j,1)
               if( ic.gt.mnp ) ic = mnp
               if( ic.eq.i+j) then 
                  errc = abs(cov_parm(i,ic) - ((i-1)*mnp + j))
                  errs = abs(sol_parm(i) + i)
                  if( errc.gt.1.d-6 .or. errs.gt.1.d-6 ) then
                      
                      print *,'END RC ',i,ic,' Diff ',errc,errs
                      stop 'Matrix bad'
                  endif
                endif
            end do
         end do
      end if

      print *,'Done'
      end
 



 
                  

