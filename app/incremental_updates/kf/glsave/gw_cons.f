CTITLE GW_CONS

      subroutine gw_cons

      implicit none  
 
*     Routine to write the apriori covariance matrix used in the
*     analysis into the global file.
 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
 
      include '../includes/globk_common.h'
 
*   rot_mat(3,3)     - Transformation from NEU to XYZ (transposed at
*         at first
*   loc_coord(3)     - Dummy lat long and height
*   cov_neu(3,3)     - Covariance matrix of neu and up
*   cov_xyz(3,3)     - cov_neu projected to XYZ
*   temp_cov(3,3)    - Work space during comps
*   cov              - General real*8 covariance value (used since
*                      cov_apr is real*4).

      real*8 rot_mat(3,3), loc_coord(3), cov_neu(3,3), cov_xyz(3,3), 
     .       temp_cov(3,3), cov

*   record(128) -  Record to be output (paired as p1,p2, covariance)
*                  Therefore 32 triplets per record.
*   rec_acvc    -  Current record number of the apriori covariance
*                  matrix
*   num_in_rec  -  Number of triplets in current record.  When this
*                  equals 32 it is written to the file.
*   pir          - Position in record of first parameter number

      integer*4 record(128),rec_acvc, num_in_rec, pir, i,j,k,l,
     .          ierr, np
 
*   param_done(max_glb_parn)  - Logical to indicate that the
*            apriori constraint has been output for a parameter
*            (mainly because stations are done by themselves)

      logical  param_done(max_glb_parn)
      
****  Start      
      do i = 1, num_glb_parn
         param_done(i) = .false.          
      end do
      
****  Initialize the record counters
      rec_acvc =  crec_acvc 
      num_in_rec = 0
      cnum_acvc  = 0   
      
***** Now add the NEU station constrains.  Also count the number of
*     3x3 constraint matrices put on the solution.  This is used when
*     write the apriori constraints to combined global files.

      do i = 1, gnum_sites

*        Loop over value and rate.
         do j = 1,2
            np = parn_site(1,j,i)
            if( np.ne.0 ) then

*               Get the transformation from XYZ to NEU
                call xyz_to_neu( rot_mat, apr_val_site(1,1,i), 
     .                           loc_coord)

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

****            Now loop over XYZ, and add any apr_site variance
*               as well.
                do k = 1,3
                   cov_xyz(k,k) = cov_xyz(k,k) + cov_apr(np+k-1)
                   param_done(np+k-1) = .true.
                                   
*                  Now put the values into the output record
                   do l = 1,k
                      num_in_rec = num_in_rec + 1
                      pir = (num_in_rec-1)*4 + 1
                      record(pir) = np+k-1
                      record(pir+1) = np+l-1
*                     Move two values for real*8 variable                      
                      call dwmov(cov_xyz(k,l), 1, record(pir+2),1,1)
                      

*****                 Now see if we need to write the record
                      if( num_in_rec.eq.32 ) then
                          call writd(cglb_dcb, ierr, record, 128,
     .                               rec_acvc)
                          rec_acvc = rec_acvc + 1
                          cnum_acvc = cnum_acvc + num_in_rec
                          num_in_rec = 0
                      end if
                   end do
                end do          
*                       ! Parameter estimated
             end if
*                       ! looping on value and rate
         end do
*                       ! Looping over the sites.
      end do

****  Now do the polar motion UT1 parameters. 
* MOD TAH 200806 Modified delta calc to account for 24-hr duration
*     and gepoch_start 12 hours for expt epoch (i.e, +12 to start) 
      if( sort_direction.eq. 1 ) then
          deltat = (gepoch_expt-(gepoch_start+0.5d0))/365.25d0
      else
          deltat = (gepoch_expt-(gepoch_end-0.5d0))/365.25d0
      end if

      do i = 1, 4
         np = parn_wob(i)
         if( np.gt.0 ) then
             param_done(np) = .true.
             cov = apr_wob(i)**2 + mar_wob(i)*abs(deltat)
             if( i.le.2 ) then
*                See of we have a rate contribution
                 cov = cov + apr_wob(i+2)**2*(deltat*365.25d0)**2 +
     .                       mar_wob(i+2)*abs(deltat)*
     .                                       (deltat*365.25d0)**2
             end if
             num_in_rec = num_in_rec + 1
             pir = (num_in_rec-1)*4 + 1
             record(pir) = np
             record(pir+1) = np
             call dwmov(cov, 1, record(pir+2),1,1)

*****        Now see if we need to write the record
             if( num_in_rec.eq.32 ) then
                 call writd(cglb_dcb, ierr, record, 128,
     .                      rec_acvc)
                 rec_acvc = rec_acvc + 1
                 cnum_acvc = cnum_acvc + num_in_rec
                 num_in_rec = 0
             end if
          end if
      end do

*     Now do UT1   
      do i = 1, 2   
         np = parn_ut1(i)
         if( np.gt.0 ) then
             param_done(np) = .true.
             cov = apr_ut1(i)**2 + mar_ut1(i)*abs(deltat)
             if( i.eq.1 ) then
*                See of we have a rate contribution
                 cov = cov + apr_ut1(i+1)**2*(deltat*365.25d0)**2 +
     .                       mar_ut1(i+1)*abs(deltat)*
     .                                       (deltat*365.25d0)**2
             end if
             num_in_rec = num_in_rec + 1
             pir = (num_in_rec-1)*4 + 1
             record(pir) = np
             record(pir+1) = np
             call dwmov(cov, 1, record(pir+2),1,1)

*****        Now see if we need to write the record
             if( num_in_rec.eq.32 ) then
                 call writd(cglb_dcb, ierr, record, 128,
     .                      rec_acvc)
                 rec_acvc = rec_acvc + 1
                 cnum_acvc = cnum_acvc + num_in_rec
                 num_in_rec = 0
             end if
          end if
      end do

****  Now write out the values for the other parameters
      do np = 1, num_glb_parn
         if( .not.param_done(np) ) then 
             num_in_rec = num_in_rec + 1
             pir = (num_in_rec-1)*4 + 1
             record(pir) = np
             record(pir+1) = np

*            Move value to real*8 and then copy
             cov = cov_apr(np) 

*            see if Markov process
             do i = 1, num_glb_mar
                if( np.eq.ind_mar(i) ) then
                   cov = cov + 
     .                   cov_mar(i)*abs(deltat)
                end if
             end do
             call dwmov(cov, 1, record(pir+2),1,1)
                      

*****        Now see if we need to write the record
             if( num_in_rec.eq.32 ) then
                 call writd(cglb_dcb, ierr, record, 128,
     .                      rec_acvc)
                 rec_acvc = rec_acvc + 1
                 cnum_acvc = cnum_acvc + num_in_rec
                 num_in_rec = 0
             end if
          end if
      end do           

****  See if we have a residual amount of record to write to
*     file
      if( num_in_rec.gt.0 ) then

*         Clear the trailing part of the record      
          do i = num_in_rec*4+1, 128
              record(i) = 0
          end do         
          cnum_acvc = cnum_acvc + num_in_rec
          call writd(cglb_dcb, ierr, record, 128,
     .               rec_acvc)
          call report_error('VWRIT',ierr,'writ','COV_PARM',
     .                      0,'GW_CONS')
      end if
 
 
***** Thats all
      return
      end
 
