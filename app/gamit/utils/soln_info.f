      subroutine soln_info(ih,isnx,qnum_sites,qnum_sats,qsite_code
     .                    ,apr,adj,est_type,ref_time,num_orb,qnum_parn
     .                    ,qparn_sites,qparn_svs,qparn_eop,cov_parm
     .                    ,qnum_eop,fix_eop,constr,cov_apr,sol_type)

c  Purpose: write out the soln estimate,std dev and the soln vcv blocks to
c           a sinex file
c
c IN:
c     in, isnx        - unit numbers                                     I*4
c     qnum_sites      - number of sites in hfile                         I*4
c     qnum_sats       - number of sats in hfile                          I*4
c     qsite_code      - 4 char site code                                 C*4(qnum_sites)
c     apr             - apriories of all parameters                      R*8(qnum_parn)
c     adj             - adjustments of all parameters                    R*8(qnum_parn)
c     est_type        - type of constraint on param (0 - tight           I*4(maxsit*3)
c                             1 - medium; 2 - loose)
c     ref_time        - referemce epoch of solution                      C*12
c     num_orb         - number of orbital parameters estimated in hfile  I*4
c     qnum_parn        - total number of parameters in hfile             I*4
c     fix_eop         - apriori values of any fixed eop parameters       R*8(6)
c     constr          - apriori constraints of parameters                R*8(qnum_parn)
c     cov_apr         - apriori constraints on parameters                R*8(qnum_parn,qnum_parn)
c                       (passed in so array is dimensioned properly)
c     sol_type    - solution type indicators ( constraints,type of
c                   params, eop estimated, bias handling)                C*4
c
c
c  P Tregoning
c  11th August, 1995
c
c  NOTES: this subroutine uses the TAH routine read_cov to read in the vcv, transform
c         it from spherical to cartesian and rescale the units of the parameters

      implicit none

      include '../includes/dimpar.h'

      integer ih,isnx,qnum_sites,qnum_sats,qnum_parn
     .      , est_type(qnum_parn),qparn_sites(3,qnum_sites)
     .      , qparn_svs(maxorb,qnum_sats),qparn_eop(2,3)
     .      , qnum_eop,sat_par_num,np,count,num_orb,i,j

      real*8 apr(qnum_parn),adj(qnum_parn),constr(qnum_parn)
     .      ,par_sdv(maxsit*3+maxsat*maxorb+6)
     .      ,cov_parm(qnum_parn,qnum_parn)
     .      ,soln_parm(maxsit*3+maxsat*maxorb+6),fix_eop(6),scale_sv(15)
     .      ,scale_pmu(6),site_pos(3,maxsit)
     .      ,llr_xyz(3*maxsit,3),cov_apr(qnum_parn,qnum_parn)
     .      ,qscale(maxsit*3+maxsat*maxorb+6)
      character*4 qsite_code(qnum_sites),p_type(3),eop_type(6)
     .            ,sol_type(4),eop_units(6)
      character*12 ref_time


      data p_type/'STAX','STAY','STAZ'/
     .     ,eop_type/'XPO ','XPOR','YPO ','YPOR','UT  ','LOD '/
     .     ,eop_units/'mas ','ma/s','mas ','ma/s','ms  ','ms  '/

      data scale_sv / 1000.d0, 1000.d0, 1000.d0,
     .                1000.d3, 1000.d3, 1000.d3,
     .                   1.d0,    1.d0,    1.d0,
     .                   1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0  /

      data scale_pmu / 1.d3, 86400.d3, 1.d3, 86400.d3, 1.d3, 86400.d3 /

c use TAH routine to read the vcv and rotate it to cartesian space
      call rd_cov(ih,cov_parm,soln_parm,qnum_parn,qnum_sites
     .             ,qnum_sats,apr,adj,qparn_sites,qparn_svs
     .             ,qparn_eop,num_orb,scale_sv,scale_pmu
     .             ,site_pos,llr_xyz,qscale )


c now read in the apriori constraints
      call read_apr(ih,qnum_sites,qnum_sats,qnum_eop,num_orb,qnum_parn
     .              ,qscale,qparn_eop,cov_apr)

c assign the estimated standarv deviations
       do i=1,qnum_parn
         par_sdv(i) = dsqrt(cov_parm(i,i))
       enddo

c assign the estimation type for the station coords - arbitrarily set tight
c constraint to be less than 10mm, otherwise loose!
       do i=1,qnum_sites*3
         if(dsqrt(cov_apr(i,i)).lt.0.010)then
           est_type(i) = 0
         else
           est_type(i) = 2
         endif
       enddo

c assign the estimation types for the eop parameters
       do i=1,3
         do j=1,2
           np = qparn_eop(j,i)
c           if(qparn_eop(j,i).ne.0)then
c             est_type(qparn_eop(j,i)) = 2
           if(np.ne.0)then
             est_type(np) = 2

           endif
         enddo
       enddo

      write(*,'(/,a,/)')' Writing the solution ......'
c  open the SOLUTION/ESTIMATE block
      write(isnx,100)
100   format('+SOLUTION/ESTIMATE',/
     .      ,'*INDEX TYPE CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED'
     .      ,' VALUE___ _STD_DEV______')


c now loop through the sites
      count = 0
      do i=1,qnum_sites
        do j=1,3
          count = count + 1
          write(isnx,110)count,p_type(j),qsite_code(i),ref_time
     .               ,est_type(3*i-3+j),soln_parm(qparn_sites(j,i))
     .               ,par_sdv(qparn_sites(j,i))
110       format(i6,' ',a4,' ',a4,'  A    1 ',a12,' m   ',i2,e21.13
     .           ,1x,e14.8)
        enddo
      enddo

c  now write out the eop estimates

      do i=1,3
        do j=1,2
          if(qparn_eop(j,i).ne.0)then
            count = count + 1
            write(isnx,120)count,eop_type(2*i-2+j),ref_time
     .                 ,eop_units(2*i-2+j),est_type(qparn_eop(j,i))
     .                ,soln_parm(qparn_eop(j,i)),par_sdv(qparn_eop(j,i))
120         format(i6,' ',a4,' ---- --    1 ',a12,' ',a4,i2,e21.13
     .           ,1x,e14.8)
          endif
        enddo
      enddo

c  close the block
      write(isnx,'(a)')'-SOLUTION/ESTIMATE'
      call write_dash(isnx)

c now save the a priori constraints in a vector
       do i=1,qnum_sites*3
          constr(i) = dsqrt(cov_apr(i,i))
       enddo

c now the eop
      do i=1,3
        do j=1,2
          np = qparn_eop(j,i)
          if(np.ne.0)then
            constr(np) = dsqrt(cov_apr(np,np))
          endif
        enddo
      enddo

c      call printmat(constr,qnum_parn,1,'apriori constraints ')

c  write the soln_apriori block
      call wrtsol_apr(isnx,qnum_sites,qnum_eop,qnum_parn,qsite_code
     .                    ,qparn_eop,qparn_sites,apr,fix_eop,est_type
     .                    ,constr,ref_time,eop_type,eop_units,scale_pmu
     .                    ,site_pos,sol_type,llr_xyz)

c  write out the vcv
      call wrtcov(isnx,cov_parm,qnum_parn,qnum_sites,qparn_eop,qnum_eop)

c  compute the number of satellite parameters estimated
       sat_par_num = qnum_sats*num_orb

      call wrtcov_apr(isnx,cov_apr,scale_pmu,qnum_sites,qnum_parn
     .               ,qnum_eop,qparn_sites,qparn_eop,site_pos,llr_xyz
     .               ,sat_par_num)

      return
      end
