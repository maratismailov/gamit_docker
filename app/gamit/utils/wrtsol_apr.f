      subroutine wrtsol_apr(isnx,qnum_sites,qnum_eop,qnum_parn
     .                    ,qsite_code,qparn_eop,qparn_sites,apr,fix_eop
     .                    ,est_type,constr,ref_time,eop_type,eop_units
     .                    ,scale_pmu,site_pos,sol_type,llr_xyz)

c Purpose: write out the solution/apriori block
c
c IN:
c     isnx        - unit number of sinex file                          I*4
c     qnum_sites  - number of estimated sites in hfile                 I*4
c     qnum_eop    - number of eop parameters estimated in soln         I*4
c     qnum_parn    - total number of parameters estimated in soln       I*4
c     qsite_code  - 4 char site codes of estimated sites               C*4(qnum_sites)
c     qparn_eop   - matrix of eop parameter indices                    I*4(2,3)
c     qparn_sites - matrix of site parameter indices                   I*4(3,qnum_sites)
c     apr         - apriori values for estimated parameters            R*8(qnum_parnn)
c     fix_eop     - apriori values of fixed eop parameters             R*8(6)
c     eop_type    - eop parameter names                                C*4(6)
c     eop_units   - eop parameter units                                C*4(6)
c     est_type    - (0=tight,1=looser,2=loose)                         I*4(qnum_parn)
c     constr      - apriori constraints of parameters                  R*8(qnum_parn)
c     ref_time    - reference eopch                                    C*12
c     eop_type    - name of eop parameters                             C*4(6)
c     eop_units   - units of eop parameters                            C*4(6)
c     scale_pmu   - scaling matrix for eop parameters                  R*8(6)
c     site_pos    - xyz apriori site coords of estimated sites         R*8(3,maxsit)
c     sol_type    - solution type indicators ( constraints,type of
c                   params, eop estimated, bias handling)              C*4
c     llr_xyz     - rotation matrices for all sites                    R*8(3*maxsit,3)
c
c P Tregoning
c 17th August, 1995
c
c
c FIXED PARAMETERS: Any fixed parameters in the solution (eg UT1 etc) must be output
c                   in this block, and will be assigned a zero standard deviation to
c                   indicate that they are fixed.

      implicit none

      include '../includes/dimpar.h'

      integer isnx,qnum_sites,qnum_eop,qnum_parn,i,j,count
     .      ,qparn_eop(2,3),qparn_sites(3,qnum_sites)
     .      ,est_type(qnum_parn)
      real*8 fix_eop(6),apr(qnum_parn),constr(qnum_parn),zero
     .       ,scale_pmu(6),site_pos(3,maxsit),llr_xyz(3*maxsit,3)
     .       ,wght(3,3),tmp(3,3),jtrans(3,3),wght_xyz(3,3)
      character*4 qsite_code(qnum_sites),eop_type(6),eop_units(6)
     .            ,p_type(3),sol_type
      character*12 ref_time

      data p_type/'STAX','STAY','STAZ'/
      data zero/0.d0/

c  open the apriori block
      write(isnx,100)
100   format('+SOLUTION/APRIORI',/
     .      ,'*ATTENTION: STD_DEV IS GIVEN IN LOCAL N,E,UP SYSTEM! USE'
     .      ,' FULL COVARIANCE MATRIX',/
     .      ,'*INDEX TYPE CODE PT SOLN _REF_EPOCH__ UNIT S __APRIORI '
     .      ,'VALUE_____ _STD_DEV______')

c now loop through the sites
      count = 0
      do i=1,qnum_sites

c  form a matrix of the apriori weights and the neu_xyz rotation matrix
        do j = 1,3
            wght(j,j) = constr(3*i-3+j)**2
        enddo

c  now rotate the apriori llr weights to xyz
c                            t
c  covar (xyz) = J cov(llr) J

        call matmpy(llr_xyz(3*i-3+1,1),wght,tmp,3,3,3)
        call transp(llr_xyz(3*i-3+1,1),jtrans,3,3)
        call matmpy(tmp,jtrans,wght_xyz,3,3,3)

        do j=1,3
          count = count + 1
          write(isnx,110)count,p_type(j),qsite_code(i),ref_time
     .             ,est_type(3*i-3+j),site_pos(j,i),dsqrt(wght_xyz(j,j))
110       format(i6,' ',a4,' ',a4,'  A    1 ',a12,' m   ',i2,e21.13
     .           ,1x,e14.8)
        enddo
      enddo

c now write out any eop parameters that were estimated in the solution
      do i=1,3
        do j=1,2
          if(qparn_eop(j,i).ne.0)then
            count = count + 1
            write(isnx,120)count,eop_type(2*i-2+j),ref_time
     .                 ,eop_units(2*i-2+j),est_type(qparn_eop(j,i))
     .                ,apr(qparn_eop(j,i))*scale_pmu(2*i-2+j)
     .                ,constr(qparn_eop(j,i))
120         format(i6,' ',a4,' ---- --    1 ',a12,' ',a4,i2,e21.13
     .           ,1x,e14.8)
          endif
        enddo
      enddo


c  now write out any eop parameters that were fixed in the solution
      do i=1,3
        do j=1,2
          if(qparn_eop(j,i).eq.0)then
            count = count + 1
            write(isnx,120)count,eop_type(2*i-2+j),ref_time
     .                 ,eop_units(2*i-2+j),int(zero)
     .                ,fix_eop(2*i-2+j)*scale_pmu(2*i-2+j),zero
          endif
        enddo
      enddo

c close the block
      write(isnx,'(a)')'-SOLUTION/APRIORI'
      call write_dash(isnx)

      return
      end
