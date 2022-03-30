CTITLE MOV_COV
 
      subroutine mov_cov(parn_1,parn_2, nel, covar, ndim, cov_parm,
     .   par_num)

      implicit none 
 
 
c
c     This routine will move the covariances from ema storage
c     into main memory in packed form.  The routine is used
c     to compute variances of quanities derived from the estimated
c     parameter.
c
c Variables
c ---------
c parn_1 -- the parameter numbers of the first entry
c parn_2 -- the parameter numbers of the second entry
c covar  -- the main memory covariance matrix
c ndim   -- the dimensioned size of covar
c nel    -- the number of parameters to moved for each of parn_1 and
c     parn_2
c cov_parm -- the ema stored covariance matrix
c par_num -- number of parameters estimated
c
      integer*4 parn_1(1), parn_2(1), ndim, nel, par_num
 
c
      real*8 covar(ndim,1)
 
c
      real*8 cov_parm(par_num,1)
 
c
c
c Local parameters
c ----------------
c jel and kel -- local index variables
c
 
      integer*4 jel, kel
 
*   i,j,k   - Loop counters
      integer*4 i,j,k
 
c
c.... First clear the copy of the covariance matrix in main memory
      do i = 1, nel
C        call dwsmy(0.d0,covar(1,i),1,covar(1,i),1, nel)
         call dwint(0.d0,covar(1,i),1, nel)
      end do
c
c.... Now loop over the elements
      do j = 1,nel
*                                        ! get the first entry
         if( parn_1(j).ne.0 ) then
            jel = parn_1(j)
c
            do k = 1, nel
*                                          ! get the second entry
              if( parn_2(k).ne.0 ) then
                 kel = parn_2(k)
                 covar(j,k) = cov_parm(jel,kel)
               end if
            end do
c
         end if
      end do
c
      return
      end
 
