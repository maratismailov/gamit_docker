CTITLE MAKE_CTOG

      subroutine make_ctog( unit, line, np, cov_parm, sol_parm )

      implicit none

*     Routine to read and make a binary hfile for ctog files:
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
* nq       - Generic parameter number
* ns       - Simple form of number of sites
* i,j      - Loop counters
* num_in_file = number of solutions in file (1 only for ctog solns)
* ierr     - IOSTAT error
* trimlen  - Length of string

      integer*4 date(5), nq, ns, i,j, num_in_file, ierr,
     .          trimlen

* sectag   - Seconds part of time
* param_adj - Parmeter adjustment read from file
* corr      - Parameter corelation read from file

       real*8 sectag, param_adj, corr

       ns = 0

****  This is ctog file:  Read the number of paramters and date from
*     the first line.
      qnum_parn = np

*     Skip the next two lines
      read(unit,'(a)') line 
      read(unit,'(a)') line 

*     Convert character month to numeric month
      read(unit,'(a)') line 
      read(line(20:),* ) date, sectag
      call ymdhms_to_jd( date, sectag, sepoch )

*     Skip the next two lines
      read(unit,'(a)') line 
      read(unit,'(a)') line 

****  Now get the liast of stations and read the parameter estimates
      qnum_sites = 0
      qnum_svs   = 0
      do i = 1, qnum_parn
         read(unit,'(a)', iostat=ierr) line

*        Now see what the line is:
         if( index(line,'X-coordinate').gt.0 ) then
             qnum_sites = qnum_sites + 1
             ns = qnum_sites
             read(line,200) np, qsite_names(ns), sol_parm(np),param_adj
 200         format(i5,2x,A8,1x,t45,F24.10,f12.4)
             qparn_sites(1,ns) = np 
             qscale(np) = 1.d0
             site_pos(1,ns) = sol_parm(np) - param_adj
         else if( index(line,'Y-coordinate').gt.0 ) then
             read(line,200) np, qsite_names(ns), sol_parm(np),param_adj
             qparn_sites(2,ns) = np 
             qscale(np) = 1.d0
             site_pos(2,ns) = sol_parm(np) - param_adj
         else if( index(line,'Z-coordinate').gt.0 ) then
             read(line,200) np, qsite_names(ns), sol_parm(np),param_adj
             qparn_sites(3,ns) = np 
             qscale(np) = 1.d0
             site_pos(3,ns) = sol_parm(np) - param_adj

****     See if satellite orbital elements
         else if( index(line,'Inert X').gt.0 ) then
             qnum_svs = qnum_svs + 1
             ns = qnum_svs
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(1,ns) = np 
             qscale(np) = 1.d0
             svs_pos(1,ns) = sol_parm(np) - param_adj
         else if( index(line,'Inert Y').gt.0 ) then
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(2,ns) = np 
             qscale(np) = 1.d0
             svs_pos(2,ns) = sol_parm(np) - param_adj
         else if( index(line,'Inert Z').gt.0 ) then
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(3,ns) = np 
             qscale(np) = 1.d0
             svs_pos(3,ns) = sol_parm(np) - param_adj
         else if( index(line,'In dX/dT').gt.0 ) then
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(4,ns) = np 
             qscale(np) = 1.d0
             svs_pos(4,ns) = sol_parm(np) - param_adj
         else if( index(line,'In dY/dT').gt.0 ) then
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(5,ns) = np 
             qscale(np) = 1.d0
             svs_pos(5,ns) = sol_parm(np) - param_adj
         else if( index(line,'In dZ/dT').gt.0 ) then
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(6,ns) = np 
             qscale(np) = 1.d0
             svs_pos(6,ns) = sol_parm(np) - param_adj
         else if( index(line,'Direct R').gt.0 ) then
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(7,ns) = np 
             qscale(np) = 1.d0
             svs_pos(7,ns) = sol_parm(np) - param_adj
         else if( index(line,'Y-Bias  ').gt.0 ) then
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(8,ns) = np 
             qscale(np) = 1.d0
             svs_pos(8,ns) = sol_parm(np) - param_adj
         else if( index(line,'Z-Bias  ').gt.0 ) then
             read(line, 200) np, qsvs_names(ns), sol_parm(np),param_adj
             qparn_svs(9,ns) = np 
             qscale(np) = 1.d0
             svs_pos(9,ns) = sol_parm(np) - param_adj
         else

             write(*,290) line(1:max(1,trimlen(line)))
 290         format('Unkown line in ctog file: ',a)
         end if
      end do

*     Write out some information
      write(*,300) qnum_sites, hfile(1:trimlen(hfile)) 
 300  format(i5,' sites found in ',a)
      do i = 1, qnum_sites
         qfull_names(i) = qsite_names(i) // ' from ctog Stacov file'
         write(*,310) i, qsite_names(i)
 310     format(i4,'. ',a8)
      end do

      write(*,320) qnum_svs
 320  format(i5,' Satellites found')

*     Now read the correlation maxtrix.
*     Skip two lines (1 blank)
      read(unit,'(a)') line 
      read(unit,'(a)') line 
      do i = 1, qnum_parn
         do j = 1,i
            read(unit,*) np, nq, corr
            cov_parm(np,nq) = corr
            cov_parm(nq,np) = cov_parm(np,nq)
         end do
      end do

*     Generate the Fake KalObs file name
      sKalobs_file = hfile(1:trimlen(hfile))
 
*     Get Mfile name (claim this is database)
      sdata_base = hfile(trimlen(hfile)-13:trimlen(hfile)-4)

*     Save some of the epochs:
      qstart_epoch = sepoch - 0.5d0
      qend_epoch   = sepoch + 0.5d0
      ssvs_epoch   = sepoch
     
 
****  Now generate the apriori and parameter codes
 
      call qgen_codes
 
***** Now write out the global file in GLOBK format.
      hf_ext = 'glm'
 
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

