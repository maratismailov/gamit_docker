CTITLE MAKE_GSFC_VLBI

      subroutine make_gsfc_vlbi( unit, line, np, cov_parm, sol_parm )

      implicit none

*     Routine to read and make a binary hfile for GSFC VLBI files:
*     Lots of entries need to be dummied because information is
*     not in these files.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'


* unit   - Unit for reading the file
* np     - Number of parameters in this solution

      integer*4 unit, np

* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector

      real*8 cov_parm(np,np), sol_parm(np)

* line   - Line read from file

      character*(*) line 

* LOCAL VARIABLES

* date(5)  - Calender date
* num_month - Function to return numeric month given ascii month
* nq       - Generic parameter number
* ns       - Simple form of number of sites
* i,j      - Loop counters
* num_in_file = number of solutions in file (1 only for JPL solns)
* ierr     - IOSTAT error
* trimlen  - Length of string

      integer*4 date(5), nq, ns, i,j, num_in_file, ierr,
     .          trimlen, k

* sectag   - Seconds part of time
* param_var - Parmeter variance
* param_adj - Adjustment to parameter
* param_est - Estimate of parameter

       real*8 sectag, param_var, param_adj, param_est

* name - Generic site name

       character*8 name

* ayr  - ASCII Year for site position (when coesmics are
*        estimated.

       character*2 ayr

      ns = 0
*
****  This is VLBI GSFC File:  Read the number of paramters and date from
*     the first line.
      qnum_parn = np
      sdata_base = line(trimlen(line)-9:)
      hfile_type = 'VLB'

****  Skip the next line to get to reference date
      read(unit,'(a)' ) line
      read(unit,'(a)' ) line

*     Read the reference date:
      read(line,100, iostat=ierr) date(1), date(2), date(3)
 100  format(11x,3i2)
      date(4) =  0
      date(5) =  0
      sectag  = 0.0
      call ymdhms_to_jd( date, sectag, sepoch )

*     Skip blank line
      read(unit,'(a)') line

****  Now get the liast of stations and read the parameter estimates
      qnum_sites = 0
      do i = 1, qnum_parn
         read(unit,'(a)', iostat=ierr) line

*        Now see what the line is:
         if( index(line,'X COMP').gt.0 ) then
             qnum_sites = qnum_sites + 1
             ns = qnum_sites
             qsite_names(ns) = ' '
             read(line,200) np, qsite_names(ns), param_var, param_adj,
     .                          param_est
 200         format(i6,1x,A8,15x,d22.16,d23.16,d23.16)
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_sites(1,ns) = np 
             site_pos(1,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else if( index(line,'Y COMP').gt.0 ) then
             read(line,200) np, name, param_var, param_adj,
     .                          param_est
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_sites(2,ns) = np 
             site_pos(2,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else if( index(line,'Z COMP').gt.0 ) then
             read(line,200) np, name, param_var, param_adj,
     .                          param_est
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_sites(3,ns) = np 
             site_pos(3,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else if( index(line,'X VELO').gt.0 ) then
             read(line,200) np, name, param_var, param_adj,
     .                          param_est
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_vel(1,ns) = np 
             site_vel(1,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else if( index(line,'Y VELO').gt.0 ) then
             read(line,200) np, name, param_var, param_adj,
     .                          param_est
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_vel(2,ns) = np 
             site_vel(2,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else if( index(line,'Z VELO').gt.0 ) then
             read(line,200) np, name, param_var, param_adj,
     .                          param_est
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_vel(3,ns) = np 
             site_vel(3,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else if( index(line,'AXIS O').gt.0 ) then
             read(line,200) np, name, param_var, param_adj,
     .                          param_est

*            If parameter has zero variance set to finite value
             if( param_var.eq.0 ) param_var = 1.d0
             cov_parm(np,np) = param_var
*            DUE to Bug in file, save only the adjusemt to axo
C            sol_parm(np)    = param_est
C            site_axo(ns) = param_est - param_adj
             sol_parm(np) = param_adj
             site_axo(ns) = 0.d0
             qparn_axo(ns) = np 
             qscale(np) = 1.d0
         else if( index(line(10:),' X').gt.0 ) then
             qnum_sites = qnum_sites + 1
             ns = qnum_sites
             qsite_names(ns) = ' '
             read(line,220) np, qsite_names(ns), ayr, param_var, 
     .                          param_adj, param_est
 220         format(i6,1x,A8,2x,a2,11x,d22.16,d23.16,d23.16)
             qsite_names(ns)(7:8) = ayr
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_sites(1,ns) = np 
             site_pos(1,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else if( index(line(10:),' Y').gt.0 ) then
             read(line,200) np, name, param_var, param_adj,
     .                          param_est
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_sites(2,ns) = np 
             site_pos(2,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else if( index(line(10:),' Z').gt.0 ) then
             read(line,200) np, name, param_var, param_adj,
     .                          param_est
             cov_parm(np,np) = param_var
	     sol_parm(np)    = param_est
             qparn_sites(3,ns) = np 
             site_pos(3,ns) = param_est - param_adj
             qscale(np) = 1.d0
         else 

*             We don;t know what line is
             write(*,290) line(1:max(1,trimlen(line)))
 290         format('Unkown line in GSFC VLBI file: ',a)
         end if
      end do

*     Write out some information
      write(*,300) qnum_sites, hfile(1:trimlen(hfile)) 
 300  format(i5,' sites found in ',a)
      call underscore(qsite_names, qnum_sites)
      do i = 1, qnum_sites
         qfull_names(i) = qsite_names(i) // ' from GSFC VLBI file'
         write(*,310) i, qsite_names(i)
 310     format(i4,'. ',a8,1x,a32)
      end do

*     Now read the correlation maxtrix.
      do i = 2, qnum_parn
         read(unit,410) np, (cov_parm(np,j),j=1,i-1)
 410     format(i6,1x, 5d23.16,:,/,1000(7x,5d23.16,:,/))

*        Copy row
         do j = 1, i-1
            cov_parm(j,np) = cov_parm(np,j)
         end do
      end do

****  Now add in the translation variation in covariance matrix.
      do i = 1, qnum_sites
         do j = 1, qnum_sites
            do k = 1,3

*              Position (adds 10 m uncertain to all coordinates 
*              including covariances)
               np = qparn_sites(k,i)
               nq = qparn_sites(k,j)
               if( np.gt.0 .and. nq.gt.0 ) then
                   cov_parm(np,nq) = cov_parm(np,nq) + 100.d0
               end if

*              Rates (adds 10 m/yr to variances and covariances)
               np = qparn_vel(k,i)
               nq = qparn_vel(k,j)
               if( np.gt.0 .and. nq.gt.0 ) then
                   cov_parm(np,nq) = cov_parm(np,nq) + 100.d0
               end if
            end do
         end do
      end do


****  Now rescale the matrix (meter everywhere except solar radiation
*     which will be ratio to direct effect)
      do i = 1, qnum_parn
          call dwsmy(qscale(i), cov_parm(1,i),1, 
     .                          cov_parm(1,i),1,  qnum_parn)
          call dwsmy(qscale(i), cov_parm(i,1),qnum_parn,
     .                          cov_parm(i,1),qnum_parn, qnum_parn)
 
          sol_parm(i) = sol_parm(i)*qscale(i)
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
      hf_ext = 'glg'
 
      call gen_gname( hfile, sepoch, num_in_file, glb_dir,
     .                glb_file, hf_ext, name_format,expt_code)
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

