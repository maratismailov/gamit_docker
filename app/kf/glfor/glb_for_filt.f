CTITLE GLB_FOR_FILTER
 
      subroutine glb_for_filter( cov_parm, sol_parm, row_copy,
     .                           part_pnt, a_part,
     .                           temp_gain, kgain, cov_obs, sol_obs,
     .                           glb_used )

      implicit none 
  
*     Routine to do the forward Kalman filtering.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j,k   - Loop counter
*   part_pnt(2, max_glb_deriv, cnum_parn) - pointers to the
*   nc          - counter for position in a_part
*   num         - number of contiguous partials
*   pn          - number of parameters for this observation
*   start       - First partial number
*   trimlen     - HP function for length of string
 
 
      integer*4 part_pnt(2, max_glb_deriv, cnum_parn)
 
*   a_part(max_glb_deriv,cnum_parn) - compressed form of the
*           - partials matrix
 
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix for
*           - the global parameters
*   cov_obs(cnum_parn, cnum_parn)   - covariance matrxo for the
*           - parameters read from disk
 
*   kgain(cnum_parn,num_glb_parn)   - Kalman filter Gain matrix
 
*   row_copy(num_glb_parn)          - Copy of one row of the
*                           - covariance matrix.  Used when column
*                           - in state transsion matrix is less than
*                           - row
*   sol_parm(num_glb_parn)  - The solution vector for the global
*           - parameters
 
*   sol_obs(cnum_parn)  - Solution vector read from disk
 
*   temp_gain(cnum_parn,num_glb_parn) - Working space for computing
*           - the Kalman gain.
 
      real*8 a_part(max_glb_deriv,cnum_parn),
     .    cov_parm(num_glb_parn,num_glb_parn),
     .    cov_obs(cnum_parn, cnum_parn), kgain(cnum_parn,num_glb_parn),
     .    row_copy(num_glb_parn), sol_parm(num_glb_parn),
     .    sol_obs(cnum_parn), temp_gain(cnum_parn,num_glb_parn)

* glb_used - Logical to denote if global file used

      logical glb_used 
 
***** Increment the filter for the next experiment being added.
!     print *,'DEBUG: gwob_apr ',gwob_apr(1:4), gut1_apr(1:2)
 
*     Predict the values at the next step 
      call glb_predict( cov_parm, sol_parm, row_copy )
  
!     if( cnum_used.le.16 ) then 
!        do i = 1, cnum_used
!           write(*,120) i, cov_obs(i,1:cnum_used) 
!120        format('COV_OBS ',i3,100(1x,E15.6))
!        end do

!        do i = 1, cnum_used
!           write(*,125) i, 
!    .            (cov_obs(i,j)/sqrt(cov_obs(i,i)*cov_obs(j,j)),
!    .             j=1,cnum_used)
!125        format('COR_OBS ',i3,100(1x,F8.5))
!        end do

!     endif

*     Now compute the filter gain for this new data
      call glb_kalman_gain( cov_parm, a_part, part_pnt, cov_obs,
     .                      temp_gain, kgain )
 

*     Now add the information from this experiment
* MOD TAH 960905: Moved the pre-fit chi**2 calc into glb_update
*     so that bad experiments could be detected.

      call glb_update( cov_parm, sol_parm, a_part, part_pnt,
     .                 temp_gain, kgain, sol_obs, cov_obs, 'Y',
     .                 glb_used  )
!     if( num_glb_parn.lt.16 ) then 
!         write(*,130) (i, sol_parm(i), 
!    .             sqrt(cov_parm(i,i)),i=1,num_glb_parn)
!130      format(16('PARAM UPDATE ',i4,2(1x,E15.6),/))
!     endif

*     Increment the "pre-fit" statistics.  Moved to inside update
 
C     call glb_inc_prefit( cov_obs, sol_obs, dchi )

***** Thats all
      return
      end
 
CTITLE CHECK_DC
 
      subroutine check_dc( cov_parm, dir_copy )

      implicit none 
 
*     Routine to check if we can do a direct copy of the initial
*     covariance matrix.  This is done by checking if apriori
*     variance on all the parameters is greater than or equal
*     to 100.d0 (equivalent to 10 meter for position and
*     10 m/yr for velocity/
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   i   - Loop counter
      integer*4  i
 
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix for
*           - the global parameters
 
      real*8 cov_parm(num_glb_parn,num_glb_parn)

*   dir_copy - Logical set true to indicate that we can copy

      logical dir_copy

***** Set dir_copy true and see if we have any violation.  See
*     apriori sigmas are large or that there are no derived
*     parameters (for example, polar motion or translations)
      dir_copy = .true.
      do i = 1, num_glb_parn
         if( cov_parm(i,i).lt.100.d0 .or.
     .       indx_pnt(i).gt. 1      )  dir_copy = .false.
      end do

***** Thats all
      return
      end


CTITLE COPY_DIR
 
      subroutine copy_dir( cov_parm, sol_parm, row_copy,
     .                     part_pnt, a_part, cov_obs, sol_obs )

      implicit none 
 
*     Routine to directly copy the initial covariance matrix 
*     from cov_obs into cov_parm
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   part_pnt(2, max_glb_deriv, cnum_parn) - pointers to the
 
      integer*4 part_pnt(2, max_glb_deriv, cnum_parn)
 
*   a_part(max_glb_deriv,cnum_parn) - compressed form of the
*           - partials matrix
 
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix for
*           - the global parameters
*   cov_obs(cnum_parn, cnum_parn)   - covariance matrxo for the
*           - parameters read from disk
*   row_copy(num_glb_parn)          - Copy of one row of the
*                           - covariance matrix.  Used when column
*                           - in state transsion matrix is less than
*                           - row
 
*   sol_parm(num_glb_parn)  - The solution vector for the global
*           - parameters
 
*   sol_obs(cnum_parn)  - Solution vector read from disk
 
      real*8 a_part(max_glb_deriv,cnum_parn),
     .    cov_parm(num_glb_parn,num_glb_parn),
     .    cov_obs(cnum_parn, cnum_parn),  row_copy(num_glb_parn),
     .    sol_parm(num_glb_parn), sol_obs(cnum_parn)

*  i, j, - Loop counters overs data
*  np,mp - Loop counters over parameters

      integer*4 i,j, np, mp

****  Execute the glb_predict routine so that we will acount for
*     any earthquakes.

      call glb_predict( cov_parm, sol_parm, row_copy )

***** Run down the data covariance maxtrix and solution
      do i = 1, cnum_parn

*****    See if we should copy
         if( indx_pnt(i).gt.0 .and. 
     .       part_pnt(1,1,i).gt.0 ) then

             np = part_pnt(1,1,i)
             sol_parm(np) = sol_obs(i)

*            now check the rows that should be copied
             do j = 1, cnum_parn
                if( indx_pnt(j).gt.0 .and. 
     .              part_pnt(1,1,j).gt.0 ) then
                    mp = part_pnt(1,1,j)
                    cov_parm(np,mp) = cov_obs(i,j)
                endif
             end do
          end if
      end do

***** Thats all
      return
      end 

