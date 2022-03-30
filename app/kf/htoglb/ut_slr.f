
CTITLE GET_UT_SLR_PARN 

      subroutine get_ut_slr_parn( unit, line, np)

      implicit none

*     This routine will read the SLR(UT) solution part
*     of the covariance matrix to find out how many parameters
*     have been estimated.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'


* unit   - Unit for reading the file
* np     - Number of parameters in this solution

      integer*4 unit, np

* line   - Line read from file

      character*(*) line 

* LOCAL VARIABLES

* i,j      - Loop counters
* ierr     - IOSTAT error
* trimlen  - Length of string
* ns       - temp site counter.

      integer*4 date(5), ns, i, ierr,
     .          trimlen, sn, indx

* sectag   - Seconds part of time
* param_val(3) - Value of parameter (XYZ or NEU dot)
* param_adj(3) - Adjustment to apriori (XYZ of NEU dot)
* xyz_vel(3)   - XYZ velocity (m/yr)
* xyz_adj(3)   - Adjustment to XYZ velocity. 
* loc_coord(3) - Geod lat, long and height
* rot_matrix(3,3) - Rotation from XYZ to NEU

      real*8 sectag,  param_val(3), param_adj(3), xyz_vel(3), 
     .       xyz_adj(3), loc_coord(3), rot_matrix(3,3)

* temp_name  - Temprary names for stations when rate checked.

      character*8 temp_name

*
****  This is an SLR file. Skip two lines to get to start of station 
*     coordinates
      read(unit,'(a)', iostat=ierr ) line
      read(unit,'(a)', iostat=ierr ) line

*     Set the nominal epoch to 1994/1/1
      date(1) = 1994
      date(2) = 1
      date(3) = 1
      date(4) = 0
      date(5) = 0
      sectag = 0.d0

      call ymdhms_to_jd( date, sectag, sepoch )

      write(*,100) sdata_base, date
 100  format(/' Found SLR solution ',a,' Reference Epoch ',
     .        i4,'/',i2,'/',i2,1x,i2,':',i2)

****  Now get the liast of stations and read the parameter estimates
      np = 0
      ns = 0
      do while ( trimlen(line).gt.0 )
         read(unit,'(a)', iostat=ierr) line
         if( trimlen(line).gt.0 ) then
             read(line(21:),*) param_val, param_adj

*            See type of parameter
             ns = ns + 1
             qnum_sites = ns
             qsite_names(ns) = line(11:16) // 'UT'
             do i = 1,3
                 np = np + 1
                 qparn_sites(i,ns) = np
                 qscale(np)     = 1.d0
                 site_pos(i,ns) = param_val(i) - param_adj(i)
             end do
          end if
      end do

****  Add underscores to names
      call underscore( qsite_names, qnum_sites)

****  Now do the velocity field. Skip two lines
      read(unit,'(a)', iostat =ierr ) line
      read(unit,'(a)', iostat =ierr ) line

      do while ( trimlen(line).gt.0 .and. 
     .           line(1:2).ne.'--' .and. ierr.eq.0 )

          read(unit,'(a)', iostat =ierr ) line
          if( trimlen(line).gt.0 .and.
     .           line(1:2).ne.'--' .and. ierr.eq.0 ) then

              temp_name = line(2:7) // 'UT'
              call underscore( temp_name, 1)
              indx = 1
              call get_cmd(temp_name, qsite_names, qnum_sites, 
     .                         sn,indx )

*             Read the rates of the line
              read(line(13:),*,iostat=ierr) param_val, param_adj

*             Switch the North and East values
              call swapne(param_val)
              call swapne(param_adj)

*             Now convert NEU velocity to XYZ
              call rotate_geod(param_val, xyz_vel, 'NEU', 'XYZ',
     .                         site_pos(1,sn), loc_coord, rot_matrix)
              call rotate_geod(param_adj, xyz_adj, 'NEU', 'XYZ',
     .                         site_pos(1,sn), loc_coord, rot_matrix)

*             Convert the values to meters.
              do i = 1,3
                 xyz_vel(i) = xyz_vel(i)/1000.d0
                 xyz_adj(i) = xyz_adj(i)/1000.d0

*                Save parameter numbers
                 np = np + 1
                 qparn_vel(i,sn) = np
                 qscale(np) = 1.d0
                 site_vel(i,sn) = xyz_vel(i) - xyz_adj(i)
              end do
          end if
      end do

      qnum_parn = np
             

***** Thats all
      return
      end

CTITLE MAKE_UT_SLR

      subroutine make_ut_slr( unit, line, nt, cov_parm, sol_parm )

      implicit none

*     Routine to read the slr file and get the solution vector
*     and covariance matrix.  The file is assumed to be have
*     been rewound back to the beginning

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'


* unit   - Unit for reading the file
* nt     - Total Number of parameters in this solution

      integer*4 unit, nt

* cov_parm(nt,nt) - Covarince matrix
* sol_parm(nt)    - Solution vector

      real*8 cov_parm(nt,nt), sol_parm(nt)

* line   - Line read from file

      character*(*) line 

* LOCAL VARIABLES

* nq       - Generic parameter number
* i,j      - Loop counters
* num_in_file = number of solutions in file (1 only for JPL solns)
* ierr     - IOSTAT error
* trimlen  - Length of string

      integer*4 i,j,k,  num_in_file, ierr,
     .          trimlen, np, nq, mp, mq

* cov_mean(6,6) - Covarinace matrix of position and velocity at 
*                 mean epoch
* cov_ref(6,6)  - Site position covariance at reference epoch.
* dt            - Time difference between ref and mean epocjh (yrs)
* cov_neu(3,3)  - NEU Covariance matrix
* cov_xyz(3,3)  - XYZ Covariamce matrix
* loc_coord(3)  - Geod Lat/Long/hgt
* rot_matrix(3,3)  - Ratation from NEU to XYZ
* param_est(3)  - Estimate parameters
* param_adj(3)  - Adjustment to aprorio
* param_sig(3)  - Sigma of parameters
* param_rho(3)  - Correlations of the parameters

* xyz_est(3)    - XYZ velocity estimated
* xyz_adj(3)    - XYZ velcoity adjustment

* trans(6,6)    - Transformation for mean epoch to solution epoch
* cov_temp(6,6) - Temporary storage

* mean_epoch(max_glb_sites) - Mean epochs of the global sites

      real*8 cov_mean(6,6), cov_ref(6,6), dt, cov_neu(3,3),
     .       cov_xyz(3,3), loc_coord(3), rot_matrix(3,3), param_est(3),
     .       param_adj(3), param_sig(3), param_rho(3),
     .       mean_epoch(max_glb_sites), xyz_est(3), xyz_adj(3),

     .       trans(6,6), cov_temp(6,6)

****  Clear the entire covariance matrix before we start
      do i = 1, qnum_parn
         do j = 1, qnum_parn
            cov_parm(i,j) = 0.d0
         end do
      end do

*
****  This is an SLR(UT) file: Skip down to the start of the 
*     parameter estimates
      do while ( index(line,'NAME').eq.0 )
         read(unit,'(a)' ) line
      end do

*     Now read the solution vector
      do i = 1, qnum_sites
         read(unit,'(a)', iostat=ierr) line

****     Now decode the line
         read(line,*) mean_epoch(i)
         read(line(21:),*) param_est, param_adj, param_sig,
     .                    param_rho
         np = qparn_sites(1,i)

*        Save the estimate on sol_parm.  We find velocity later
*        we will update these values to final epoch.
         do j = 1,3
            sol_parm(np-1+j) = param_est(j)
         end do

*        Save the local covarince elements
         cov_parm(np  ,np  ) = param_sig(1)**2
         cov_parm(np  ,np+1) = param_sig(1)*param_sig(2)*param_rho(1)
         cov_parm(np  ,np+2) = param_sig(1)*param_sig(3)*param_rho(2)
        
         cov_parm(np+1,np  ) = param_sig(2)*param_sig(1)*param_rho(1)
         cov_parm(np+1,np+1) = param_sig(2)**2
         cov_parm(np+1,np+2) = param_sig(2)*param_sig(3)*param_rho(3)

         cov_parm(np+2,np  ) = param_sig(3)*param_sig(1)*param_rho(2)
         cov_parm(np+2,np+1) = param_sig(3)*param_sig(2)*param_rho(3)
         cov_parm(np+2,np+2) = param_sig(3)**2
      end do

*     Skip three lines
      read(unit,'(a)') line
      read(unit,'(a)') line
      read(unit,'(a)') line

*     Now read the the velocities   
      do i = 1, qnum_sites

*        If there is a parameter number assigned then it means there
*        if a velocity
         if( qparn_vel(1,i).ne.0 ) then
             read(unit,'(a)') line
             read(line(13:),*) param_est, param_adj, param_sig, 
     .                         param_rho

*****        Convert all estimates to m/yr
             do j = 1,3
                 param_est(j) = param_est(j)/1000.d0
                 param_adj(j) = param_adj(j)/1000.d0
                 param_sig(j) = param_sig(j)/1000.d0
             end do

*            Swap the north and east values (so we have NEU)
             call swapne(param_est)
             call swapne(param_adj)
             call swapne(param_sig)

*            Form up the NEU covariamnce matrix so we can rotate
*            to XYZ. NOTE: RHO is EN, EU and NU so swap the last two
*            entries
             call swapne(param_rho(2))
*            Save the local covarince elements
             cov_neu(1  ,1  ) = param_sig(1)**2
             cov_neu(1  ,1+1) = param_sig(1)*param_sig(2)*param_rho(1)
             cov_neu(1  ,1+2) = param_sig(1)*param_sig(3)*param_rho(2)
        
             cov_neu(1+1,1  ) = param_sig(2)*param_sig(1)*param_rho(1)
             cov_neu(1+1,1+1) = param_sig(2)**2
             cov_neu(1+1,1+2) = param_sig(2)*param_sig(3)*param_rho(3)

             cov_neu(1+2,1  ) = param_sig(3)*param_sig(1)*param_rho(2)
             cov_neu(1+2,1+1) = param_sig(3)*param_sig(2)*param_rho(3)
             cov_neu(1+2,1+2) = param_sig(3)**2

****         Now rotate the estimated velocity into XYZ
             call rotate_geod(param_est, xyz_est, 'NEU', 'XYZ',
     .                        site_pos(1, i), loc_coord, rot_matrix)
             call rotate_geod(param_adj, xyz_adj, 'NEU', 'XYZ',
     .                        site_pos(1, i), loc_coord, rot_matrix)

****         Compute velocity covariance in XYZ
             call var_comp(rot_matrix, cov_neu, cov_xyz, cov_temp,
     .                     3,3,1)

*            Transferr into cov_mean
             np = qparn_sites(1,i)
             do j = 1,3
                do k = 1,3
                   cov_mean(j  ,k  ) = cov_parm(np-1+j,np-1+k)
*                  Since we are at mean epoch there should be no
*                  velocity, position correlation
                   cov_mean(j  ,k+3) = 0.d0
                   cov_mean(k+3,j  ) = 0.d0
                   cov_mean(j+3,k+3) = cov_xyz(j,k)
                end do
             end do
       
****         Now move the position to the solution epoch
             dt = ((sepoch - 2 400 000.5)-mean_epoch(i))/365.25d0

****         Form the transformation
             do j = 1, 3
                do k = 1,6
                   trans(j  ,k) = 0.d0
                   trans(j+3,k) = 0.d0
                end do
                trans(j  ,j  )  = 1.d0
                trans(j  ,j+3)  = dt
                trans(j+3,j+3)  = 1.d0

****            Update solution vector
                nq = qparn_vel(j,i)
                np = qparn_sites(j,i)
                sol_parm(nq) = xyz_est(j)
                sol_parm(np) = sol_parm(np) + xyz_est(j)*dt
                site_pos(j,i) = site_pos(j,i) + site_vel(j,i)*dt
             end do

****         Compute covariance matrix at the reference epoch:
             call var_comp(trans,cov_mean, cov_ref, cov_temp, 6,6,1)

****         Now copy back into covar
             do j = 1,3
                do k = 1,3
                   np = qparn_sites(j,i)
                   mp = qparn_sites(k,i)
                   nq = qparn_vel(j,i)
                   mq = qparn_vel(k,i)
                   cov_parm(np,mp) = cov_ref(j  ,k  )
                   cov_parm(np,mq) = cov_ref(j  ,k+3)
                   cov_parm(mq,np) = cov_ref(k+3,j  )
                   cov_parm(nq,mq) = cov_ref(j+3,k+3)
               end do
             end do
          end if
      end do
   
*     Write out some information
      write(*,300) qnum_sites, hfile(1:trimlen(hfile)) 
 300  format(i5,' sites found in ',a)
      do i = 1, qnum_sites
         qfull_names(i) = qsite_names(i) // ' from SLR(UT) file'
         write(*,310) i, qsite_names(i)
 310     format(i4,'. ',a8)
      end do

*     Generate the Fake KalObs file name
      sKalobs_file = hfile(1:trimlen(hfile))
 
*     Save some of the epochs:
      qstart_epoch = sepoch - 0.5d0
      qend_epoch   = sepoch + 0.5d0
      ssvs_epoch   = 0.d0  
     
 
****  Now generate the apriori and parameter codes
 
      call qgen_codes
 
***** Now write out the global file in GLOBK format.
      hf_ext = 'glu'
 
      call gen_gname( hfile, sepoch, num_in_file, glb_dir,
     .                glb_file, hf_ext, name_format, expt_code)
      write(*,600) glb_file(1:trimlen(glb_file))
 600  Format(' Creating Binary file ',a)
 
***** Write out the header (create file if need be)
      call qw_glb_header
 
*     Write out the names
      call qw_glb_names
 
*     Write out the full names
      call qw_glb_full
 
*     Write out the solution description
      call qw_description
 
*     Write out the codes
      call qw_codes
 
*     Write out the apriori values
      call qw_aprioris
 
*     Write out the solution and covariance matrix
      call qw_soln( cov_parm )
 
***** Thats all for this solution.  Now see if there is another
*     one
      return
      end

CTITLE SWAPNE

      subroutine swapne( param )

      implicit none

*     Swaps the first two elements of param
      real*8 param(2), temp

      temp = param(1)
      param(1) = param(2)
      param(2) = temp

      return
      end
