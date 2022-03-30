CTITLE CRESNX_SOLN
 
      subroutine cresnx_soln(unit, unitc, np, cov_parm, sol_parm )
 
      implicit none 

*     Routine to create the solution  blocks from SINEX:
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   unitc       - Unit number of comments file
*   np      - number of parameters expected
 
      integer*4 unit, unitc, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
 
* LOCAL VARIABLES
  
      logical  kbit
 
****  Start looking for blocks
      call cre_soln_epochs(unit,unitc)
      call cre_soln_stats(unit,unitc)
      call cre_soln_est(unit,unitc, np, cov_parm, sol_parm)
      call cre_soln_apr(unit,unitc, np, cov_parm, sol_parm)
      if( .not.kbit(glbtosnx_opt,1) ) then
          call cre_soln_mat(unit,unitc, np, cov_parm, sol_parm)
          call cre_soln_con(unit,unitc, np, cov_parm, sol_parm)
      end if
 
****  Thats all
      return
      end
 
 
CTITLE CRE_SOLN_epochs
 
      subroutine cre_soln_epochs(unit,unitc )
 
      implicit none 

*     Routine to decode the solution epochs to get the start and
*     stop times for each station
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit, unitc
 
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   syr, sdoy, ssec - start yr, doy and ses
*   eyr, edoy, esec - stop yr, doy and see
*   indx        - Position in string
*   occ         - Occupation number
 
      integer*4  sn, syr, sdoy, ssec,  occ,
     .          eyr, edoy, esec, myr, mdoy, msec

 
*   code    - Site code read from Sinex file
*   pt      - Point characters
*   full_name  -  Full site name with occ and pt added
 
      character*8 code
      character*4 pt, type
 
****  Start reading the block
      write(unit,'(a)') '+SOLUTION/EPOCHS'
      call cp_comments(unit, unitc, 'SOLUTION/EPOCHS')
      do sn = 1, cnum_sites

* MOD TAH 1906127: Only output entry if site actually used.
         if( gtol_sites(sn).ge.0 ) then 

            call name_to_id( qsite_names, sn, code, pt, occ)
            type  = qfull_names(sn)(32:32)

            call jd_to_yds( qdata_st(sn),syr, sdoy, ssec )
            call jd_to_yds( qdata_en(sn),eyr, edoy, esec )
            call jd_to_yds( (qdata_st(sn)+qdata_en(sn))/2,myr, 
     .                                                     mdoy, msec )
            write(unit,210) code,pt,occ, type, mod(syr,100), sdoy, ssec,
     .               mod(eyr,100), edoy, esec, mod(myr,100), mdoy, msec
 210        format(1x,a4,2x,a1,1x,i4,1x,a1,1x,i2.2,':',i3.3,':',i5.5,
     .          1x,i2.2,':',i3.3,':',i5.5,1x,i2.2,':',i3.3,':',i5.5)
         endif 
      end do
      write(unit,'(a)') '-SOLUTION/EPOCHS'

****  Thats all
      return
      end
      
CTITLE CRE_SOLN_stats
 
      subroutine cre_soln_stats(unit,unitc )
 
      implicit none 

*     Routine to output the solution stats block to get number of
*     data and chiqsq 
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit, unitc
 
* LOCAL VARIABLES
 
 
****  Start writing the block
      write(unit,'(a)') '+SOLUTION/STATISTICS'
      call cp_comments(unit, unitc, 'SOLUTION/STATISTICS')
 
      write(unit,210) 'NUMBER OF OBSERVATIONS', cnum_obs
 210  format(1x,a,t32,i10)
      write(unit,210) 'NUMBER OF UNKNOWNS', cnum_parn
      
      write(unit,220) 'VARIANCE FACTOR', cchisq
 220  format(1x,a,t32,f22.15)

      write(unit,'(a)') '-SOLUTION/STATISTICS'             
 
****  Thats all
      return
      end

CTITLE CRE_SOLN_EST
 
      subroutine cre_soln_est(unit,unitc, np, cov_parm, sol_parm)
 
      implicit none 

*     Routine to decode the solution estimates part of the sinex
*     files.  These values are written directly into the sol_parm
*     vector.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, unitc, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
  
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   cyr, cdoy, csec - Reference yr, doy and sec
*   indx        - Position in string
*   nr          - "real" parameter number in final covariance matrix
*   occ         - Occupation number
*   cons_type   - Constraint type on station (0,1 or 2)
*   ne          - Pointer to epoch number for multi-pmu
*   nf          - Sequential pointer to each multipmu type.
 
      integer*4 ierr,  cyr, cdoy, csec, occ, i, k,
     .          cons_type, j, ni, no, ne, nf
 
  
*   code    - Site code read from Sinex file
*   stype   - Parameter type read from sinex file
*   type    - System type
*   estu        - Units for the estimate
*   pt      - Point characters
*   full_name  -  Full site name with occ and pt added
*   stype_pos(3), stype_vel(3), stype_pmu(2,3) -- Names of parameter types
*   sunit_pos(3), sunit_vel(3), sunit_pmu(2,3) -- Names of parameter types

      character*8 code
      character*6 pt, stype_pos(3), stype_vel(3), stype_trn(3,2),
     .                stype_scl(2), stype_pmu(2,3), stype_sva(3)
      character*4 type
      character*4     sunit_pos(3), sunit_vel(3), sunit_trn(2),
     .                sunit_scl(2), sunit_pmu(2,3)
      character*1 offcon    ! Function to return constellation letter 
                            ! based on offset PRN (original PRN for 
                            ! mod(prn,100))
      character*1 sys       ! GREC system code

 
      data stype_pos / 'STAX  ', 'STAY  ', 'STAZ  ' /
      data stype_vel / 'VELX  ', 'VELY  ', 'VELZ  ' /
      data stype_sva / 'SATA_X', 'SATA_Y', 'SATA_Z' /

      data stype_trn / 'XGC   ', 'YGC   ', 'ZGC   ',
     .                 'XGCR  ', 'YGCR  ', 'ZGCR  ' /
      
      data stype_scl / 'SBIAS ', 'SBIASR' /
      
      data stype_pmu / 'XPO   ', 'XPOR  ', 'YPO   ', 'YPOR  ', 
     .                 'UT    ', 'LOD   ' /

      data sunit_pos / 'm   ', 'm   ', 'm   ' /
      data sunit_vel / 'm/yr', 'm/yr', 'm/yr' /
      data sunit_trn / 'm   ', 'm/yr' /
      data sunit_scl / 'ppb ', 'pb/y' /

      data sunit_pmu / 'mas ', 'ma/d', 'mas ', 'ma/d', 'ms  ', 'ms  ' /


****  Start reading the block
      write(unit,'(a)') '+SOLUTION/ESTIMATE'
      call cp_comments(unit, unitc, 'SOLUTION/ESTIMATE')

****  Now loop over each of the types of parameters
      do i = 1, cnum_sites

         call name_to_id( qsite_names, i, code, pt, occ)
         type  = qfull_names(i)(32:32)

*****    Loop over XYZ and then Xdot, Ydot, Zdot
         do j = 1, 3

*****       Get the input and output parameter number
            ni = qparn_sites(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)               
            if( no.gt.0 ) then

               call jd_to_yds( cepoch_expt, cyr, cdoy, csec)
               call get_ctype( 0.01d0, qapr_cov(j,j,i),
     .                         cons_type )
               write(unit,210,iostat=ierr) 
     .             no, stype_pos(j), code,pt,occ,mod(cyr,100),cdoy,csec,
     .             sunit_pos(j), cons_type, sol_parm(ni), 
     .             sqrt(cov_parm(ni,ni))
 210          format(i6,1x,a6,1x,a4,2x,a1,1x,i4,1x,i2.2,':',i3.3,':',
     .               i5.5,1x,a4,i2,1x,e21.15,1x,e11.6)
            end if
         end do

*****    Now do velocity
         do j = 1, 3

*****       Get the input and output parameter number
            ni = qparn_vel(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)               
            if( no.gt.0 ) then

               call jd_to_yds( cepoch_expt, cyr, cdoy, csec)

               call get_ctype( 0.01d0, qapr_vel(j,j,i),
     .                         cons_type )
               write(unit,210,iostat=ierr) 
     .             no, stype_vel(j),code,pt,occ,mod(cyr,100),cdoy,csec,
     .             sunit_vel(j), cons_type, sol_parm(ni), 
     .             sqrt(cov_parm(ni,ni))
            end if
         end do

      end do

* MOD TAH 050622: Output any satellite antenna offsets that have been 
*     estimated.
      ne = 0
      do i = 1, cnum_svs
         do j = 1, 3
            ni = qparn_svs(max_svs_elem-3+j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 ) then
* MOD TAH 180401: Output for GNSS
               sys =  offcon( qsvi_prn(i))

               write(code,320) sys, qsvi_svn(i)
 320           format(a1,I3.3)
               call get_ctype( 0.01d0, qapr_svs_ant(j,i),
     .                         cons_type )
                
*              
                write(unit,330,iostat=ierr) 
     .             no, stype_sva(j), code,'LC','1',
     .             mod(cyr,100),cdoy,csec,
     .             sunit_pos(j), cons_type, sol_parm(ni), 
     .             sqrt(cov_parm(ni,ni))
 330           format(i6,1x,a6,1x,a4,1x,a2,1x,a4,1x,i2.2,':',i3.3,':',
     .               i5.5,1x,a4,i2,1x,e21.15,1x,e11.6)
             end if
          end do
      end do

***** Now do the translation and rate parameters (Loop through XYZ
*     and then XYZ dot)
      do j = 1,2  
         do i = 1,3
            ni = qparn_tran(i,j)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 ) then
                call jd_to_yds( cepoch_expt, cyr, cdoy, csec)
                call get_ctype( 0.01d0, qapr_tran(i,j),
     .                          cons_type )
   
                occ = 1
                write(unit,210,iostat=ierr) 
     .              no, stype_trn(i,j), '----','-',occ,mod(cyr,100),
     .              cdoy, csec, sunit_trn(j), cons_type, sol_parm(ni), 
     .              sqrt(cov_parm(ni,ni))
            end if
         end do
      end do
      
***** Now do the scale and rate parameters 
      do i = 1,2
         ni = qparn_scale(i)
         no = 0
         if( ni.gt.0 ) no = itoo(ni)
         if( no.gt.0 ) then
             call jd_to_yds( cepoch_expt, cyr, cdoy, csec)
             call get_ctype( 0.01d0, qapr_scale(i),
     .                       cons_type )
   
             occ = 1
             write(unit,210,iostat=ierr) 
     .           no,stype_scl(i),'----','-',occ,mod(cyr,100),cdoy,csec,
     .           sunit_scl(j), cons_type, sol_parm(ni), 
     .           sqrt(cov_parm(ni,ni))
         end if
      end do

***** Now do the Earth rotation parameters
      do i = 1, 3
         do j = 1, 2
            ni = qparn_pmu(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 ) then

               call jd_to_yds( cepoch_expt, cyr, cdoy, csec)
               if( i.lt.3 ) then
                   call get_ctype( 0.10d0, qapr_pmu(j,i),
     .                             cons_type )
               else
                   call get_ctype( 0.007d0, qapr_pmu(j,i),
     .                             cons_type )
               end if

               occ = 1
               write(unit,210,iostat=ierr) 
     .             no, stype_pmu(j,i),'----','-',occ,mod(cyr,100),
     .             cdoy,csec,sunit_pmu(j,i), cons_type, sol_parm(ni), 
     .             sqrt(cov_parm(ni,ni))
            end if
         end do
      end do
      
***** Now do the multi-epoch Earth rotation parameters
      ne = 0
      do i = 1, 3
         do j = 1, 2
            nf = 0
            do k = 1, qnum_mul_pmu(j,i)
            
               ni = qparn_mul_pmu(j,i,k)
               no = 0
               if( ni.gt.0 ) no = itoo(ni)
               if( no.gt.0 ) then
                   ne = ne + 1
                   call jd_to_yds( qmul_par_ep(ne), cyr, cdoy, csec)
                   if( i.lt.3 ) then
                       call get_ctype( 0.10d0, qapr_pmu(j,i),
     .                                 cons_type )
                   else
                       call get_ctype( 0.007d0, qapr_pmu(j,i),
     .                                cons_type )
                   end if
                        
                   nf = nf + 1
                   write(unit,210,iostat=ierr) 
     .                no, stype_pmu(j,i), '----','-',nf, mod(cyr,100), 
     .                cdoy, csec,
     .                sunit_pmu(j,i), cons_type, sol_parm(ni), 
     .                sqrt(cov_parm(ni,ni))
               end if
            end do
         end do
      end do

      
      write(unit,'(a)') '-SOLUTION/ESTIMATE'
  
****  Thats all
      return
      end
 
 
CTITLE CRE_SOLN_APR
 
      subroutine cre_soln_apr(unit,unitc, np, cov_parm, sol_parm)

      implicit none 
 
*     Routine to decode the solution apriori part of the sinex
*     files.  These values are written directly into the site_pos
*     vector.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, unitc, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   cyr, cdoy, csec - Reference yr, doy and sec
*   indx        - Position in string
*   nr          - "real" parameter number in final covariance matrix
*   occ         - Occupation number
*   cons_type   - Constraint type on station (0,1 or 2)
*   na          - Counter for the apriori parameter number.  Some
*                 of these could ve skipped if no apriori.  Here
*                 we output all values but they may not appear
*                 in the covariance matrix
 
      integer*4 ierr, cyr, cdoy, csec,  occ, i, k,
     .          cons_type, j, ni, no, na, ne, nf
 

*   code    - Site code read from Sinex file
*   stype   - Parameter type read from sinex file
*   type    - System type
*   estu        - Units for the estimate
*   pt      - Point characters
*   full_name  -  Full site name with occ and pt added
*   stype_pos(3), stype_vel(3), stype_pmu(2,3) -- Names of parameter types
*   sunit_pos(3), sunit_vel(3), sunit_pmu(2,3) -- Names of parameter types

      character*8 code
      character*6 pt, stype_pos(3), stype_vel(3), stype_trn(3,2),
     .                stype_scl(2), stype_pmu(2,3), stype_sva(3)
      character*4 type
      character*4     sunit_pos(3), sunit_vel(3), sunit_trn(2),
     .                sunit_scl(2), sunit_pmu(2,3)

      character*1 offcon    ! Function to return constellation letter 
                            ! based on offset PRN (original PRN for 
                            ! mod(prn,100))
      character*1 sys       ! GREC system code

 
      data stype_pos / 'STAX  ', 'STAY  ', 'STAZ  ' /
      data stype_vel / 'VELX  ', 'VELY  ', 'VELZ  ' /
      data stype_sva / 'SATA_X', 'SATA_Y', 'SATA_Z' /
      data stype_trn / 'XGC   ', 'YGC   ', 'ZGC   ',
     .                 'XGCR  ', 'YGCR  ', 'ZGCR  ' /
      
      data stype_scl / 'SBIAS ', 'SBIASR' /
      
      data stype_pmu / 'XPO   ', 'XPOR  ', 'YPO ', 'YPOR  ', 
     .                 'UT    ', 'LOD   ' /

      data sunit_pos / 'm   ', 'm   ', 'm   ' /
      data sunit_vel / 'm/yr', 'm/yr', 'm/yr' /
      data sunit_trn / 'm   ', 'm/yr' /
      data sunit_scl / 'ppb ', 'pb/y' /
      
      data sunit_pmu / 'mas ', 'ma/d', 'mas ', 'ma/d', 'ms  ', 'ms  ' /


****  Start reading the block
      write(unit,'(a)') '+SOLUTION/APRIORI'
      call cp_comments(unit, unitc, 'SOLUTION/APRIORI')

****  Copy over the wobble and UT1 parameters to pmu_pos.  Change
*     sign on ut1 rate to get LOD. and convert rates to ma/s for
*     pole position.
      do i = 1,2
         pmu_pos(1,i) = cwob_apr(i)
* MOD TAH 160510: No unit conversion need (mas/day already)
!        pmu_pos(2,i) = cwob_apr(i+2)/86400.d0
         pmu_pos(2,i) = cwob_apr(i+2)
      end do
      pmu_pos(1,3) =  cut1_apr(1)/15.d0
      pmu_pos(2,3) = -cut1_apr(2)/15.d0

****  Now loop over each of the types of parameters
      na = 0
      do i = 1, cnum_sites

         call name_to_id( qsite_names, i, code, pt, occ)
         type  = qfull_names(i)(32:32)

*****    Loop over XYZ and then Xdot, Ydot, Zdot
         do j = 1, 3

*****       Get the input and output parameter number
            ni = qparn_sites(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)               
            if( no.gt.0 ) then

               call jd_to_yds( cepoch_expt, cyr, cdoy, csec)
               cons_type = 0
               na = na + 1
               write(unit,210,iostat=ierr) 
     .             no, stype_pos(j),code,pt,occ,mod(cyr,100),cdoy,csec,
     .             sunit_pos(j), cons_type, site_pos(j,i), 
     .             sqrt(qapr_cov(j,j,i))
 210          format(i6,1x,a6,1x,a4,2x,a1,1x,i4,1x,i2.2,':',i3.3,':',
     .               i5.5,1x,a4,i2,1x,e21.15,1x,e11.6)
            end if
         end do

*****    Now do velocity
         do j = 1, 3

*****       Get the input and output parameter number
            ni = qparn_vel(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)               
            if( no.gt.0 .and. qapr_vel(j,j,i).ne.0.d0 ) then

               call jd_to_yds( cepoch_expt, cyr, cdoy, csec)

               cons_type = 0
               na = na + 1
               write(unit,210,iostat=ierr) 
     .             no, stype_vel(j),code,pt,occ,mod(cyr,100),cdoy,csec,
     .             sunit_vel(j), cons_type, site_vel(j,i),
     .             sqrt(qapr_vel(j,j,i))
            end if
         end do

      end do

* MOD TAH 050622: Output any satellite antenna offsets that have been 
*     estimated.
      ne = 0
      do i = 1, cnum_svs
         do j = 1, 3
            ni = qparn_svs(max_svs_elem-3+j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 ) then
* MOD TAH 180401: Output for GNSS
               sys =  offcon( qsvi_prn(i))
               write(code,320) sys,qsvi_svn(i)
 320           format(a1,I3.3)
               cons_type = 0
               na = na + 1  
*              print *,'OUT ',i,ni,no,na, qsvi_svn(i), qsvi_prn(i)              
*
               if( svs_pos(max_svs_elem-3+j,i).le.1.d-12 ) 
     .             svs_pos(max_svs_elem-3+j,i) = 0.d0
              
               write(unit,330,iostat=ierr) 
     .             no, stype_sva(j), code,'LC','1',
     .             mod(cyr,100),cdoy,csec,
     .             sunit_pos(j), cons_type, 
     .             svs_pos(max_svs_elem-3+j,i),
     .             sqrt(qapr_svs_ant(j,i))
 330          format(i6,1x,a6,1x,a4,1x,a2,1x,a4,1x,i2.2,':',i3.3,':',
     .               i5.5,1x,a4,i2,1x,e21.15,1x,e11.6)
              end if
          end do
      end do

 
***** Now do the translation and rate of change parmateres
      do j = 1,2 
         do i = 1,3
            ni = qparn_tran(i,j)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 ) then
                call jd_to_yds( cepoch_expt, cyr, cdoy, csec)
                cons_type = 0
                occ = 1
                na = na + 1
                write(unit,210,iostat=ierr) 
     .              no, stype_trn(i,j),'----','-',occ,mod(cyr,100),
     .              cdoy, csec, sunit_trn(j), cons_type, 0.d0, 
     .              sqrt(qapr_tran(i,j))
            end if
         end do
      end do
      
***** Now do the scale and rate of change parmateres
      do i = 1,2
         ni = qparn_scale(i)
         no = 0
         if( ni.gt.0 ) no = itoo(ni)
         if( no.gt.0 ) then
             call jd_to_yds( cepoch_expt, cyr, cdoy, csec)
             cons_type = 0
             occ = 1
             na = na + 1
             write(unit,210,iostat=ierr) 
     .           no, stype_scl(i), '----','-',occ,mod(cyr,100),
     .           cdoy, csec, sunit_scl(i), cons_type, 0.d0, 
     .           sqrt(qapr_scale(i))
         end if
      end do

***** Now do the Earth rotation parameters
      do i = 1, 3
         do j = 1, 2
            ni = qparn_pmu(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 .and. qapr_pmu(j,i).gt.0 ) then

               call jd_to_yds( cepoch_expt, cyr, cdoy, csec)
               cons_type = 0
               occ = 1

               na = na + 1
               write(unit,210,iostat=ierr) 
     .             no, stype_pmu(j,i), '----','-',occ,mod(cyr,100),
     .             cdoy, csec, sunit_pmu(j,i), cons_type, pmu_pos(j,i), 
     .             sqrt(qapr_pmu(j,i))
            end if
         end do
      end do
 
***** Now do the multi-parameter Earth rotation parameters
      ne = 0
      do i = 1, 3
         do j = 1, 2
            nf = 0
            do k = 1, qnum_mul_pmu(j,i)
               ni = qparn_mul_pmu(j,i,k)
               no = 0
               if( ni.gt.0 ) no = itoo(ni)
               if( no.gt.0  ) then
                  ne = ne + 1
                 
                  call jd_to_yds( qmul_apr_ep(ne), cyr, cdoy, csec)
                  cons_type = 0
                  nf = nf + 1
                  na = na + 1
                  write(unit,210,iostat=ierr) 
     .                no, stype_pmu(j,i), '----','-',nf, 
     .                mod(cyr,100), cdoy, csec, sunit_pmu(j,i), 
     .                cons_type, qpmu_mul_apr(j,i,nf), 
     .                sqrt(qpmu_mul_asig(j,i,nf))
               end if
            end do
         end do
      end do
     
      write(unit,'(a)') '-SOLUTION/APRIORI'
  
****  Thats all
      return
      end
 
CTITLE CRE_SOLN_MAT
 
      subroutine cre_soln_mat(unit,unitc, np, cov_parm, sol_parm )
 
      implicit none 

*     Routine to decode the covariance matrix of the estimates from
*     the sinex files.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, unitc, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 

* LOCAL VARAIBLES
 
 
*   ierr        - IOSTAT error
*   np1, np2    - Parameter numbers
*   nr1, nr2    - Parameter numbers in the output solution
*   i,j     - Loop counters
*   lim     - Number of elements to copy with each read
*   io, ko  - Input parameter numbers corresponding to output
*             parameter i,k
 
      integer*4 i,j,k, lastk, io, ko
      real*8 cov(3)
 
****  Start reading the block
      write(unit,'(a)') '+SOLUTION/MATRIX_ESTIMATE L COVA'
      call cp_comments(unit, unitc, 'SOLUTION/MATRIX_ESTIMATE')

      do i = 1, qnum_parn
         io = otoi(i)
         do j = 1,i,3
            if( i-j.lt.3 ) then
                lastk = i-j+1
            else
                lastk = 3
            end if
            do k = 1, lastk
               ko = otoi(j+k-1)
               cov(k) = cov_parm(ko,io)
            end do
            write(unit,210) i,j, (cov(k),k=1,lastk)
 210        format(2i6,3(1x,e21.14))
         end do
      end do
 
      write(unit,'(a)') '-SOLUTION/MATRIX_ESTIMATE L COVA'
****  Thats all
      return
      end
 
 
CTITLE CRE_SOLN_CON
 
      subroutine cre_soln_con(unit,unitc, np, cov_parm, sol_parm )
 
      implicit none 

*     Routine to decode the apriori contraint matrix.  If the 
*     contraints are tight then they are removed from the solution.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, unitc, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
  
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   cyr, cdoy, csec - Reference yr, doy and sec
*   indx        - Position in string
*   nr          - "real" parameter number in final covariance matrix
*   occ         - Occupation number
*   cons_type   - Constraint type on station (0,1 or 2)
*   na          - Counter for the apriori parameter number.  Some
*                 of these could ve skipped if no apriori.  Here
*                 we output all values but they may not appear
*                 in the covariance matrix
 
      integer*4  i, j, k,  na

* MOD TAH 200102: Added to make covariance output consistent with
*     apriori value output
      integer*4 ni, no  ! Input and output parameter numbers
 

****  Start reading the block
      write(unit,'(a)') '+SOLUTION/MATRIX_APRIORI L COVA'
      call cp_comments(unit, unitc, 'SOLUTION/MATRIX_APRIORI')

****  Now loop over each of the types of parameters
      na = 0
      do i = 1, cnum_sites

*****    Loop over XYZ and then Xdot, Ydot, Zdot
         do j = 1, 3

* MOD TAH 200102: Check to see if site output.
            ni = qparn_sites(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)               
            if( no.gt.0 ) then
               if( qapr_cov(j,j,i).ne.0 ) then
                   write(unit,210) no, no-j+1, (qapr_cov(j,k,i), k=1,j)
 210               format(2i6,3(1x,e21.14))
               end if
            end if
         end do

*****    Now do velocity
         do j = 1, 3
            ni = qparn_vel(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)               
            if( no.gt.0 .and. qapr_vel(j,j,i).ne.0 ) then
                write(unit,210) no, no-j+1, (qapr_vel(j,k,i), k=1,j)
            end if
         end do

      end do

* MOD TAH 050622: Output any satellite antenna offsets that have been 
*     estimated.
      do i = 1, cnum_svs
         do j = 1, 3
            ni = qparn_svs(max_svs_elem-3+j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 .and. qapr_svs_ant(j,i).ne.0 ) then
                write(unit,210) no, no, qapr_svs_ant(j,i)
            end if
         end do
      end do
 
***** Now do the translation parmateres
      do j = 1,2
         do i = 1, 3
            ni = qparn_tran(i,j)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 .and. qapr_tran(i,j).ne.0 ) then
                write(unit,210) no, no, qapr_tran(i,j)
            end if
         end do
      end do
      
***** Now do the scale parmateres
      do i = 1, 2
         ni = qparn_scale(i)
         no = 0
         if( ni.gt.0 ) no = itoo(ni)
         if( no.gt.0 .and. qapr_scale(i).ne.0 ) then
             write(unit,210) no, no, qapr_scale(i)
         end if
      end do

***** Now do the Earth rotation parameters
      do i = 1, 3
         do j = 1, 2
            ni = qparn_pmu(j,i)
            no = 0
            if( ni.gt.0 ) no = itoo(ni)
            if( no.gt.0 .and. qapr_pmu(j,i).gt.0 ) then
                write(unit,210) no, no, qapr_pmu(j,i)
            end if
         end do
      end do
      
      write(unit,'(a)') '-SOLUTION/MATRIX_APRIORI L COVA'
  
****  Thats all
      return
      end
 
CTITLE GET_CTYPE

      subroutine get_ctype( typ, apr, cons_type )

      implicit none 

*     Routine to return the constraint type
*
* typ  - Type sigma for this type os parameter
*        0.01 m for position and velocity
*        0.10 mas/mts for UT1 and pole
* apr  - Apriori varaince for parameter
* cons_type  - Constraint type 0 - tight, 1 - constrained,
*              2 - loose

      real*8 typ, apr
      integer*4 cons_type

****  If the constraint is zero, then assume loose
      if( apr.lt.typ**2 ) then
          cons_type = 0
      else if( apr.lt.5*typ**2 ) then
          cons_type = 1
      else
          cons_type = 2
      end if

****  Thats all
      return
      end
            
