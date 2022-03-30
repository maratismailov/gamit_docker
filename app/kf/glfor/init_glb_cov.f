CTITLE INIT_GLB_COV
 
      subroutine init_glb_cov( cov_parm, sol_parm )

      implicit none  
 
*     Routine to the apriori variances in the covariance matrix and
*     to clear the solution vector
*
* MOD TAH 911015: Added markov neu computation.
*     This routine will also compute the markov step covariance matrix
*     for Markov neu on the sites.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   i       - Loop counter
*   np      - Parameter number of site
 
      integer*4 i, j, k, l, np
 
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix
*   sol_parm(num_glb_parn)              - Solution vector
 
*   init_val        - DEBUG value
 
      real*8 cov_parm(num_glb_parn,num_glb_parn),
     .    sol_parm(num_glb_parn), init_val

*   rot_mat(3,3)     - Transformation from NEU to XYZ (transposed at
*         at first
*   loc_coord(3)     - Dummy lat long and height
*   cov_neu(3,3)     - Covariance matrix of neu and up
*   cov_xyz(3,3)     - cov_neu projected to XYZ
*   temp_cov(3,3)    - Work space during comps

      real*8 rot_mat(3,3), loc_coord(3), cov_neu(3,3), cov_xyz(3,3), 
     .       temp_cov(3,3)

*  pa(2)  -- Parameter numbers for j offset and rate EOP
*  pb(2)  -- Parameter numbers for k offser and rate EOP

      integer*4 pa(2), pb(2)

*  cov_init(2,2) -- Initial variance covariance of offset and rate
*  cov_irw(2,3)  -- Variance of offset and rate and covariance between
*                   them
*  mar_rw, mar_irw -- Markov process noise to use.

      real*8 cov_init(2,2), cov_irw(2,2)

      real*4 mar_rw, mar_irw
 
      logical kbit 

****  Clear sol_parm and add apriori variances
 
      init_val = 0.d0
      do i = 1, num_glb_parn
 
          sol_parm(i) = init_val
 
*         Clear current columm
          call DWINT( 0.d0, cov_parm(1,i),1, num_glb_parn)
          
          cov_parm(i,i) = cov_apr(i)
      end do

***** Now add the NEU station constrains.  Also count the number of
*     3x3 constraint matrices put on the solution.  This is used when
*     write the apriori constraints to combined global files.

      gnum_off_diag = 0
      gcons_type = 2
      do i = 1, gnum_sites

*        Loop over value and rate.
         do j = 1,2
            np = parn_site(1,j,i)
            if( np.ne.0 ) then

                gnum_off_diag = gnum_off_diag + 1
*               Check the constraint level here.  If any constraint
*               is small then set the class based on the smallest
*               value.
                do k = 1,3
                   if( apr_neu(k,j,i).lt.0.050d0 .and.
     .                 gcons_type.gt.1 ) gcons_type = 1                           
                   if( apr_neu(k,j,i).lt.0.005d0 .and.
     .                 gcons_type.gt.0 ) gcons_type = 0 
                end do                          
                

*               Get the transformation from XYZ to NEU
                call xyz_to_neu( rot_mat, apr_val_site(1,1,i), 
     .                           loc_coord)
* MOD TAH 990802: Changed to xyz_to_geod call
* MOD TAH 030116: Changed back to NEU
C               call xyz_to_geod( rot_mat, apr_val_site(1,1,i), 
C    .                           loc_coord)

*               Transpose the matrix
                do k = 1,2
                    call dvswp(rot_mat(k,k+1),3,
     .                         rot_mat(k+1,k),1, 3-k)
                end do


*               Set up the diagonal NEU matrix
                do k = 1,3
                   do l = 1,3
                       cov_neu(k,l) = 0.d0
                   end do
                   cov_neu(k,k) = apr_neu(k,j,i)**2
                end do

                call var_comp(rot_mat, cov_neu, cov_xyz,
     .                        temp_cov, 3,3,1 )

                do k = 0,2
                   do l = 0,2
                      cov_parm(np+k,np+l) = 
     .                      cov_parm(np+k,np+l) + cov_xyz(k+1,l+1)
                   end do
                end do

*               Now do the NEU Markov processing.  Do only for the
*               position value (i.e., random walk markov)
                if( j.eq.1 ) then
                   do k = 1,3
                      do l = 1,3
                          cov_neu(k,l) = 0.d0
                      end do
                      cov_neu(k,k) = mar_neu(k,1,i)
                   end do

                   call var_comp(rot_mat, cov_neu, cov_xyz,
     .                           temp_cov, 3,3,1 )

                   do k = 1,3
                      do l = 1,3
                         cov_mar_neu(k,l,i) = cov_xyz(k,l)
                      end do
                   end do
                end if

*                       ! Parameter estimated
             end if
*                       ! looping on value and rate
         end do
*                       ! Looping over the sites.
      end do 

*     Now see if we need to modify the enters for the multi-day
*     polar motion.  Here we set up the apriori covariance for
*     multiday PMU.  If user has said IND estimates then we are
*     all set. 
      if( .not. kbit(mul_pmu_opt,1) .and. num_mul_pmu.gt.0 ) then
*
*         First check to see if any apriori's are -1 (indicating
*         F used as option).
          do i = 1, 4
             if( apr_wob(i).lt.0 ) apr_wob(i) = 0.0
          end do
          do i = 1, 2
             if( apr_ut1(i).lt.0 ) apr_ut1(i) = 0.0
          end do

*         Loop over all the pole and UT1 components
          do i = 1, 3   ! X, Y and UT1
*            Set up the apriori covariance matrix
             cov_init(1,2) = 0.d0
             cov_init(2,1) = 0.d0

*            Based on type of EOP get the statistics that we need.
             if( i.eq.1 ) then
                 cov_init(1,1) = apr_wob(1)**2
                 cov_init(2,2) = apr_wob(3)**2
                 mar_rw        = mar_wob(1)
                 mar_irw       = mar_wob(3)
             else if( i.eq.2 ) then
                 cov_init(1,1) = apr_wob(2)**2
                 cov_init(2,2) = apr_wob(4)**2
                 mar_rw        = mar_wob(2)
                 mar_irw       = mar_wob(4)
             else if ( i.eq.3 ) then
                 cov_init(1,1) = apr_ut1(1)**2
                 cov_init(2,2) = apr_ut1(2)**2
                 mar_rw        = mar_ut1(1)
                 mar_irw       = mar_ut1(2)
             end if
*
*            Now loop over the multi epochs
             do j = 1, num_mul_pmu
                do k = j, num_mul_pmu
*                  Get the parameter numbers that we need for the
*                  offset and rate terms
                   pa(1) = parn_mul_pmu(1,i,j)
                   pa(2) = parn_mul_pmu(2,i,j)  
                   pb(1) = parn_mul_pmu(1,i,k)
                   pb(2) = parn_mul_pmu(2,i,k) 

*                  Now get the covariance matrix contribution for these
*                  two values separated in time. 
                   call irw_cov(cov_irw, cov_init,  mar_rw, 
     .                          mar_irw, j, k,  
     .                          spacing_mul_pmu) 

*                  Now add to covariance matrix
                   if( pa(1).gt.0 .and. pb(1).gt.0 )
     .                cov_parm(pa(1),pb(1)) = cov_irw(1,1)
                   if( pa(1).gt.0 .and. pb(2).gt.0 )
     .                cov_parm(pa(1),pb(2)) = cov_irw(1,2)
                   if( pa(2).gt.0 .and. pb(1).gt.0 )
     .                cov_parm(pa(2),pb(1)) = cov_irw(2,1)
                   if( pa(2).gt.0 .and. pb(2).gt.0 )
     .                cov_parm(pa(2),pb(2)) = cov_irw(2,2)

*                  Save into the lower block as well
                   if( j.ne.k ) then
                       if( pa(1).gt.0 .and. pb(1).gt.0 )
     .                    cov_parm(pb(1),pa(1)) = cov_irw(1,1)
                       if( pa(1).gt.0 .and. pb(2).gt.0 )
     .                    cov_parm(pb(2),pa(1)) = cov_irw(1,2)
                       if( pa(2).gt.0 .and. pb(1).gt.0 )
     .                    cov_parm(pb(1),pa(2)) = cov_irw(2,1)
                       if( pa(2).gt.0 .and. pb(2).gt.0 )
     .                    cov_parm(pb(2),pa(2)) = cov_irw(2,2)
                   end if 
                end do 
             end do
          end do
      end if
            
***** Thats all
      return
      end

CTITLE FULL_IRW

      subroutine full_irw( num, mar, dt,  cov_full ) 

      implicit none 

      integer*4 num
      real*4 mar
      real*8 dt, cov_full(42,42), st(42,42), cov_rw(42,42),
     .       phi_irw, cov_int(42,42)

      integer*4 i,j,k, l

****  Generate the random walk covariance
      phi_irw = mar*dt/365.d0
      do i = 1, num+1
         do j = i, num+1
            cov_rw(i,j) = (i-1)*phi_irw
            cov_rw(j,i) = (i-1)*phi_irw
         end do
      end do

*     Generate state transition
      do i = 1, num+1
         st(1,i) = 0.d0
         do j = 2,i
            st(i,j) = dt
         end do
         do j = i+1,num+1
            st(i,j) = 0.d0
         end do

*        Now generate the rate part
         do j = 1, num+1
            if( i.eq.j ) then
                st(num+1+i,j) = 1.d0
            else
                st(num+1+i,j) = 0.d0
            end if
         end do
      end do
C     st(1,1) = 0.d0
      write(*,*) 'State Transition'
      do i = 1, 2*num+2
         write(*,270) i,(st(i,j),j=1,2*num+2)
 270     format(i2,24F7.1)
      end do

****  Now multiple to get the covariance matrix
      do i = 1, 2*num+2
         do j = 1, 2*num+2
            cov_int(i,j) = 0.d0
            do k = 1, num+1
               do l = 1, num+1
                  cov_int(i,j) = cov_int(i,j) +
     .                   st(i,k)*cov_rw(k,l)*st(j,l)
               end do
            end do
         end do
      end do
      write(*,*) 'COV INT'
      do i = 1, 2*num+2
         write(*,280) i,(cov_int(i,j),j=1,2*num+2)
 280     format(i2,24F7.2)
      end do

*     Now take the average of value pairs
      do i = 1, 2*num
         do j = 1, 2*num+2
            st(i,j) = 0.d0
         end do
      end do

      do i = 1, num
         st(i,i) = 0.5d0
         st(i,i+1) = 0.5d0
      end do
      do i = 1, num
          st(num+i,num+2+i) = 1.d0
      end do
      write(*,*) 'State Transition'
      do i = 1, 2*num
         write(*,290) i,(st(i,j),j=1,2*num+2)
 290     format(i2,24F7.1)
      end do
     

      do i = 1, 2*num
         do j = 1, 2*num
            cov_full(i,j) = 0.d0
            do k = 1, 2*num+2
               do l = 1, 2*num+2
                  cov_full(i,j) = cov_full(i,j) +
     .                   st(i,k)*cov_int(k,l)*st(j,l)
               end do
            end do
         end do
      end do

*     Remove the first variance 
C     varf = cov_full(1,1)
C     do i = 1, 2*num 
C        varf = cov_full(1,i)
C        do j = 1, i
C           cov_full(i,j) = cov_full(i,j) - varf
C           if( j.ne.i ) cov_full(j,i) =  cov_full(i,j)
C        end do
C     end do


      write(*,*) 'IRW Noise'
      do i = 1, 2*num
         write(*,300) i,(cov_full(i,j),j=1,2*num)
 300     format(i2,24F8.2)
      end do
     
      return
      end
             

CTITLE IRW_COV

      subroutine irw_cov(cov_irw, cov_init, mar_rw, mar_irw,
     .                   i1, i2, dt)

      implicit none 

*     Routine to compute the variance and covariance of the random
*     walk and integrated random walk processs.

* PASSED variables
*----------------- 
*  cov_init(2,2) -- Initial variance covariance of offset and rate
*  cov_irw(2,2)  -- Variance of offset and rate and covariance between
*                   them
*  dt            -- Spacing between values in day

      real*8 cov_init(2,2), cov_irw(2,2), dt

*  mar_rw        -- RW process noise mas**2/yr
*  mar_irw       -- IRW process noise (mas/day)**2/yr
      real*4 mar_rw, mar_irw
      real*8 phi_rw, phi_irw

*  i1, i2  -- epochs to be correlated (indices starting at 1; j<=k)

      integer*4 i1, i2

* LOCAL VARIABLES
* ---------------
* s1(2,2)   -- State transition used to get apriori variance effects
* s2(2,2)   -- State transition for second epoch
* dvar(2,2) -- Process noise contribution
* inc       -- Amount to increment sums by
* mul       -- Multiplier for covariance matrix
* add       -- additive constant
* coeff     -- Final value for coefficient 

      real*8 s1(2,2), s2(2,2), dvar(2,2), inc, mul, add, coeff 

      integer*4 i,j,k,l

****  Start the calculation.  First compute the affects of the
*     apriori variances 
      call gen_st(i1, dt, s1)
      call gen_st(i2, dt, s2)
      do i = 1,2
         do j = 1,2
            cov_irw(i,j) = 0.d0
            do k = 1, 2
               do l = 1, 2
                  cov_irw(i,j) = cov_irw(i,j) + 
     .                           s1(i,k)*cov_init(k,l)*s2(j,l)
               end do
            end do
         end do
      end do 

****  Now add in the contribution from the process noise
      Phi_rw = mar_rw*dt/365.d0
      Phi_irw = mar_irw*dt/365.d0
      dvar(1,1) = 0
      dvar(2,2) = 0
      dvar(1,2) = 0
      dvar(2,1) = 0

*     Save the random walk offset part first
      dvar(1,1) = i1*Phi_rw
*     Save the Rate IRW term (ie., the random walk for the rate)
      dvar(2,2) = i1*Phi_irw

*     Now compute the Offset to Offset covariance term
      inc = 0.25d0
      coeff = (2.d0*(i1-1)+1.d0)/4.d0 
      mul   = 2.d0
      add   = 1.d0
      do j = 2, i1
         mul = mul + 4
         inc = (mul*(i1-j+1)+add)/4.d0
         coeff = coeff + inc
         add   = add   + mul
      end do
      inc = inc + 0.25d0
      do j = i1+1, I2
         coeff = coeff + inc
      end do
      
      dvar(1,1) = dvar(1,1) + coeff*Phi_irw

*     Now get the upper Offset to Rate covariance term
      coeff = (i1-0.5d0)
      inc = coeff
      do j = 2, i1
         inc = inc - 1.d0
         coeff = coeff + inc
      end do
      dvar(1,2) = dvar(1,2) + coeff*Phi_irw

*     Now get the Rate to offset covariance term
      coeff = (i2-0.5d0)
      inc = coeff
      do j = 2, i1
         inc = inc - 1.d0
         coeff = coeff + inc
      end do
      dvar(2,1) = dvar(2,1) + coeff*Phi_irw


*     Now add the two contributions together
      do i = 1,2
         do j = 1,2 
            cov_irw(i,j) = cov_irw(i,j) + dvar(i,j)
         end do
      end do

****  Thats all
      return
      end 

CTITLE GEN_ST

      subroutine gen_st( ep, dt, st )

      implicit none 

*     Generate the state transition matrix
*     ep is assumed to start at 1 for the first epoch.

      integer*4 ep
      real*8 dt, st(2,2)

      st(1,1) = 1.d0
      st(1,2) = (ep-1)*dt
      st(2,1) = 0.d0
      st(2,2) = 1.d0 

      return
      end

