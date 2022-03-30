CTITLE GLB_PREDICT
 
      subroutine glb_predict( cov_parm, sol_parm, row_copy )

      implicit none 
 
*     Routine to step the parameters and there covariance matrix
*     forward one experiment.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   chk_start   - copy of col_start for checking if next transission
*           - uses the current row
*   col     - Column affecting current transission
*   col_num - position in trans_col for the column number
*   col_start   - Start position in trans_col for current  transission
*           - row
*   i,j,k,l - Loop counters
*   iel     - element number for covariance matrix (used for
*           - convenience)
*   mat     - Index to parameter (row number) which is affected by
*           - a transission with another markov parmeter
*   nt      - number of columns in current transission
*   row     - Affected by current transsion
*   nj,nk   - Column numbers for adding cov_mar_neu
*   nel     - Number of satellite elements to reset. (Added
*             980518 with rad_rese command.
 
      integer*4 chk_start, col, col_start, i,j,k, mat,
     .    nt, row, nj, nk, nel
 
*   coeff   - Transission coefficient
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix
*                   - of the global parameters
*   sol_parm(num_glb_parn)      - Estimates of the global
*                   - parameters
*   sol_copy        - Copy of solution value if we need to use
*                   - it later
*   row_copy(num_glb_parn)      - Copy of the row of the
*                   - covariance matrix if state transission
*                   - affects parameters before itself
 
      real*8 coeff, cov_parm(num_glb_parn,num_glb_parn),
     .    sol_parm(num_glb_parn), sol_copy, row_copy(num_glb_parn)
 
      sol_copy = 0.d0
 
***** Start updating for the state transissions
*     First do:
*      m      m  T       m      m
*     C   =  C  S   and x  = S x
*      m      m          m      m
*
 
C      call out_cov( cov_parm, num_glb_parn, 'PREDICT START')
  
      do i = 1, num_trans_row
 
*         Get the basic information about the rowa and columns
*         affected by current transission
          row       = trans_row(1,i)
          col_start = trans_row(2,i)
          nt        = trans_row(3,i)
 
*         Loop doing the diagonal terms first
          do j = 1, nt
              col = trans_col(col_start+j-1)
*                                     ! Process
              if( col.eq.row ) then
                  call get_trans_coeff( col,row, coeff, deltat,
     .                 trans_coeff(col_start+j-1))
 
*                 Use VIS scalar multiply
                  call dwsmy( coeff, cov_parm(1,row),1,
     .                        cov_parm(1,row),1, num_glb_parn)
                  sol_parm(row) = coeff*sol_parm(row)
              end if
          end do
 
*         Check ahead and make sure parameter affected here is not
*         to be used later
 
*                                         ! Dont look past last entry
          if( i.lt.num_trans_row ) then
              chk_start = trans_col(trans_row(2,i+1))
*                                         ! Next transission affects
              if( chk_start.eq.row ) then
*                                         ! current column, so copy first
                  sol_copy = sol_parm(row)
 
                  call DWMOV(cov_parm(1,row),1, row_copy,1,
     .                       num_glb_parn )
 
C                 do j = 1, num_glb_parn
C                     row_copy(j) = cov_parm(j,row)
C                 end do
 
              end if
          end if
 
 
*         Now loop over the remaining columns in the current transission
          do j = 1, nt
 
              col = trans_col(col_start+j-1)
              call get_trans_coeff( col, row, coeff, deltat,
     .                              trans_coeff(col_start+j-1) )
 
 
*             Check to see if we should use the copy of the transission
*             Diagonal has already been done so we do not do here
*                                         ! use the copy (because the
              if( col.lt.row ) then
*                                         ! row value has already
*                                         ! changed)
                  sol_parm(row) = sol_parm(row) +
     .                           sol_copy*coeff
 
*                 Now do column covariance matrix
                  call dwpiv(coeff, row_copy,1,
     .                 cov_parm(1,row),1, cov_parm(1,row),1,
     .                 num_glb_parn)
 
              end if
 
*                                     ! Use the matrix itself
              if( col.gt.row ) then
                  sol_parm(row) = sol_parm(row) +
     .                           sol_parm(col)*coeff
 
*                 Do column of covariance matrix
                  call dwpiv(coeff,cov_parm(1,col),1,
     .                 cov_parm(1,row),1, cov_parm(1,row),1,
     .                 num_glb_parn)
 
              end if
*                         ! Looping over the columns for this row
          end do
*                         ! Looping over transission rows
      end do
 
 
C      call out_cov( cov_parm, num_glb_parn, 'END OF COLUMNS')
 
 
***** Start updating for the state transissions
*     Now do:
*      m        m
*     C   =  S C
*      m        m
* 
      do i = 1, num_trans_row
 
*         Get the basic information about the rowa and columns
*         affected by current transission
          row       = trans_row(1,i)
          col_start = trans_row(2,i)
          nt        = trans_row(3,i)
 
*         Loop doing the diagonal terms first
          do j = 1, nt
              col = trans_col(col_start+j-1)
*                                     ! Process
              if( col.eq.row ) then
                  call get_trans_coeff( col,row, coeff, deltat,
     .                 trans_coeff(col_start+j-1))
 
*                 Use VIS scalar multiply
                  call dwsmy( coeff, cov_parm(row,1),num_glb_parn,
     .                 cov_parm(row,1),num_glb_parn, num_glb_parn)
              end if
          end do
 
*         Check ahead and make sure parameter affected here is not
*         to be used later
 
*                                         ! Dont look past last entry
          if( i.lt.num_trans_row ) then
              chk_start = trans_col(trans_row(2,i+1))
*                                         ! Next transission affects
              if( chk_start.eq.row ) then
*                                         ! current row, so copy first
 
                  call DWMOV(cov_parm(row,1),num_glb_parn, row_copy,1,
     .                       num_glb_parn )
 
              end if
          end if
 
 
*         Now loop over the remaining columns in the current transission
          do j = 1, nt
 
              col = trans_col(col_start+j-1)
              call get_trans_coeff( col, row, coeff, deltat,
     .                              trans_coeff(col_start+j-1) )
 
 
*             Check to see if we should use the copy of the transission
*             Diagonal has already been done so we do not do here
*                                         ! use the copy (because the
              if( col.lt.row ) then
*                                         ! row value has already
*                                         ! changed)
 
*                 Now do row covariance matrix
                  call dwpiv(coeff, row_copy,1,
     .                 cov_parm(row,1),num_glb_parn,
     .                 cov_parm(row,1),num_glb_parn, num_glb_parn)
 
              end if
 
*                                     ! Use the matrix itself
              if( col.gt.row ) then
 
*                 Do row of covariance matrix
                  call dwpiv(coeff,cov_parm(col,1),num_glb_parn,
     .                 cov_parm(row,1),num_glb_parn,
     .                 cov_parm(row,1),num_glb_parn,  num_glb_parn)
 
              end if
*                         ! Looping over the rows for this column
          end do
*                         ! Looping over transission rows
      end do
 
 
C      call out_cov( cov_parm, num_glb_parn, 'END OF ROWS')
  
 
***** Now we need to update markov elements
 
*     First loop through adding the markov statistics directly
      do i = 1, num_glb_mar
          row = ind_mar(i)
* MOD TAH 070216: Allow negative process noise which is independent
*         of time.
          if( cov_mar(i).gt.0 ) then
             cov_parm(row,row) = cov_parm(row,row) +
     .                           cov_mar(i)*abs(deltat)
          else
             cov_parm(row,row) = cov_parm(row,row) +
     .                           abs(cov_mar(i))
          endif
      end do
 
* MOD TAH 981215: See if we should increment translation process noise
*     due to white-noise translation process.
      do i = 1,3
         if( mar_tran(i,3).gt.0 .and. parn_tran(i,1).gt.0 ) then
            row = parn_tran(i,1)
            cov_parm(row,row) = cov_parm(row,row) +  mar_tran(i,3)
         end if
      end do

* MOD TAH 981215: See if we should increment rotation process noise
*     due to white-noise rotation process.
      do i = 1,3
         if( mar_rot(i,3).gt.0 .and. parn_rot(i,1).gt.0 ) then
            row = parn_rot(i,1)
            cov_parm(row,row) = cov_parm(row,row) +  mar_rot(i,3)
         end if
      end do

* MOD TAH 130716: See if white noise scale term invoke
      if( mar_scale(3).gt.0 .and. parn_scale(1).gt.0 ) then
          row = parn_scale(1)
          cov_parm(row,row) = cov_parm(row,row) +  mar_scale(3)
      end if

*     Now add the effects of transissions on the markov statistics
      do i = 1, num_glb_mar
 
          row = ind_mar(i)
 
*         Scan the transission rows to see if any transissions
*         with this parameter
          do j = 1, num_trans_row
 
*             Add direct effect and the covariance terms
              col_start = trans_row(2,j)
              nt        = trans_row(3,j)
 
*             Now see how many of the transissions are markov
              do k = 1, nt
                  col = trans_col(col_start+k-1)
 
*                 Now see if this col from the transissions matches
*                 the markov parameter we are considering
*                                         ! Yes, it does
                  if( col.eq.row ) then
*                     Get the row affected
                      mat = trans_row(1,j)
 
*                     Get transission coeff and add
                      call get_trans_coeff( mat,col, coeff,
     .                    deltat, trans_coeff(col_start+k-1))
 
*                     Add contribution
                      if( mat.eq.col ) coeff = 0.

* MOD TAH 990226: Check to see which IRW model to use:  The
*                     old_irw is not correct but is included
*                     for backward compatability.
* MOD TAH 070216: Only allow transition if not white noise process
                      if( .not. old_irw .and. cov_mar(i).gt.0 ) then
                          cov_parm(mat,mat) = cov_parm(mat,mat) +
     .                        cov_mar(i)*abs(deltat)*
     .                        (coeff*deltat)**2/3.d0
                          cov_parm(mat,col) = cov_parm(mat,col) +
     .                        cov_mar(i)*abs(deltat)*
     .                        (coeff*deltat)/2.d0
                          cov_parm(col,mat) = cov_parm(mat,col)
                      elseif( cov_mar(i).gt.0  ) then
                          cov_parm(mat,mat) = cov_parm(mat,mat) +
     .                        cov_mar(i)*abs(deltat)*
     .                        (coeff)**2
                          cov_parm(mat,col) = cov_parm(mat,col) +
     .                        cov_mar(i)*abs(deltat)*
     .                        (coeff)
                          cov_parm(col,mat) = cov_parm(mat,col)
                      end if
*                                     ! Current transission matchs a
                  end if
*                                     ! markov element
*                                     ! this row involved in
              end do
*                                     ! transissions
*                                     ! Looping over transissions
          end do
*                                     ! looping over all markov
      end do
*                                     ! parameters
***** See if we habe additional markov parameters for the satellite
*     orbits
      if( gnum_svs.gt.0 ) then 
          do i = 1, num_mar_svs
             nt = parn_mar_svs(i)
             cov_parm(nt,nt) = cov_parm(nt,nt) + cov_mar_svs(i)*
     .                                           abs(deltat)
          end do
      end if

***** Now treat satellite orbits.  Here if the ephemeris time has
*     changed, re-add the aproirio sigmas to covariance matrix.

* MOD TAH 980517: Changed update so that only the six element 
*     IC is updated (ie., Markoov process needed to get the
*      radiation parameters to change.)

      if( gnum_svs.gt.0 .and. deltaephem.ne.0.d0 ) then

* MOD TAH 980517: Changed max_svs_elem to 6, for position and
*         velocity of satellite.
          if( rad_reset ) then
* MOD TAH 981020: Never reset the antenna offset values (use mar_svan instead)
              nel = max_svs_elem - 3
          else
              nel = 6
          endif

          do i = 1, gnum_svs
* Change to number of elements here:
C            do j = 1, max_svs_elem
             do j = 1, nel
                 nt = parn_svs(j,i)
                 if( nt.gt.0 .and. apr_svs(j,i).gt.0.d0 ) then

*                    MOD TAH 950113: Clear the row and column and
*                    put the variance on the diagonal
                     call dwint(0.d0, cov_parm(1,nt),1,num_glb_parn)
                     call dwint(0.d0, cov_parm(nt,1),num_glb_parn,
     .                                num_glb_parn)

                     cov_parm(nt,nt) = apr_svs(j,i)**2
                     sol_parm(nt) = 0.d0
                 end if
             end do
* MOD TAH 190625: Added option to reset the satellite antenna offsets
*            Used in back solutions with GLX file saved (Last 3 values)
             do j = max_svs_elem - 2, max_svs_elem
                if( mar_svs(j,i).lt.0 ) then
*                    Clear row and column and re-assign aprrio sigma
* MOD TAH 200608: Added computation of nt parameter number
* MOD TAH 200731: For white-noise constraint gets applied a second
*                    during back solution so add large sigma constraint
*                    on the back
                     nt = parn_svs(j,i)
                     call dwint(0.d0, cov_parm(1,nt),1,num_glb_parn)
                     call dwint(0.d0, cov_parm(nt,1),num_glb_parn,
     .                                num_glb_parn)

* MOD TAH 200731: See if this is back solution
                     if( sort_direction*deltat.ge.0.0 ) then
*                        Apply apriori sigma
                         cov_parm(nt,nt) = apr_svs(j,i)**2
                     elseif( sort_direction*deltat.lt.0.0 ) then
*                        Make large variance for back solution.
*                        If zero leave un-changed
                         cov_parm(nt,nt) = 10.d0**2
                     endif
                     sol_parm(nt) = 0.d0
                 endif
              enddo
          end do
      end if

****  Now add the NUE markov statistics for the stations positions
      do i = 1, gnum_sites
         do j = 1,3
            nj = parn_site(j,1,i)
            do k = 1,3
               nk = parn_site(k,1,i)
               if( nj.gt.0 .and. nk.gt.0 ) then
                   cov_parm(nj,nk) = cov_parm(nj,nk) + 
     .                               cov_mar_neu(j,k,i)*abs(deltat)
               end if
            end do
         end do
      end do

****  Now update the covariance matrix for any earthquake related
*     events.  This can include pre-seismic, co-seismic and post-
*     seismic process noise.

      call eq_glb_pred(cov_parm, sol_parm)

****  Now update the multi-day PMU Process
      call mul_pmu_pred(cov_parm, sol_parm) 
 
*     DEBUG
C      call out_cov( cov_parm, num_glb_parn, 'END OF MARKOV')
 
****  Thats all
      return
      end
 
