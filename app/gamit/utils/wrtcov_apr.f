      subroutine wrtcov_apr(isnx,cov_apr,scale_pmu,qnum_sites,qnum_parn
     .               ,qnum_eop,qparn_sites,qparn_eop,site_pos,llr_xyz
     .               ,sat_par_num)

c Purpose: write out the apriori covariance matrix of the solution.
c
c IN:
c      isnx         - unit number of sinex file                       I*4
c      cov_apr      - vcv of parameter apriori constraints            R*8(qnum_parn)
c      scale_pmu    - scaling factors for eop parameters              R*8(6)
c      qnum_sites   - number of estimated sites in hfile              I*4
c      qnum_eop     - number of estimated eop parameters in soln      I*4
c      qnum_parn    - total number of parameters                      I*4
c      qparn_sites  - matrix of site parameter indices                I*4(3,qnum_sites)
c      qparn_eop    - matrix of eop parameters                        I*4(2,3)
c      site_pos     - apriori xyz coords of sites                     R*8(3,maxsit)
c      llr_xyz      - llr to xyz jacobian matrices for all sites      R*8(3*maxsit,3)
c      sat_par_num  - number of satellite parameters                  I*4
c
c P Tregoning
c 20th August, 1995

      implicit none

      include '../includes/dimpar.h'

      integer isnx,qnum_sites,qnum_eop,qnum_parn,i,j,k,np
     .        ,qparn_sites(3,qnum_sites),qparn_eop(2,3)
     .      ,sat_par_num

      real*8 cov_apr(qnum_parn,qnum_parn),scale_pmu(6),tmp(3,3)
     .       ,site_pos(3,maxsit),wght(3,3),wght_xyz(3,3) ,jtrans(3,3)
     .       ,llr_xyz(3*maxsit,3),llrxyz_mat(3,3)




c  open block in sinex file
      write(isnx,100)
100   format('+SOLUTION/MATRIX_APRIORI L COVA',/
     .      ,'*PARA1 PARA2 ____PARA2+0_________ ____PARA2+1_________ '
     .      ,'____PARA2+2_________')


c  loop through all the sites
      do i=1,qnum_sites

c  form a matrix of the apriori weights and the llr_xyz rotation matrix
        do j = 1,3
            wght(j,j) = cov_apr(3*i-3+j,3*i-3+j)
            do k = 1,3
              llrxyz_mat(j,k) = llr_xyz(3*i-3+j,k)
            enddo
        enddo

c  now rotate the apriori neu weights to xyz
c                            t
c  covar (xyz) = J cov(llr) J

        call matmpy(llrxyz_mat,wght,tmp,3,3,3)
        call transp(llrxyz_mat,jtrans,3,3)
        call matmpy(tmp,jtrans,wght_xyz,3,3,3)

c  now write out the apriori constraint values - but only output the
c  3 x 3 lower triangular matrix for each site. There are no off diagonal
c  correlations between sites!
        do j=1,3
          if(j.eq.1)write(isnx,101)qparn_sites(j,i)
     .                  ,qparn_sites(1,i),(wght_xyz(j,k),k=1,j)
101       format(i6,i6,e21.13)
          if(j.eq.2)write(isnx,102)qparn_sites(j,i)
     .                  ,qparn_sites(1,i),(wght_xyz(j,k),k=1,j)
102       format(i6,i6,2e21.13)
          if(j.eq.3)write(isnx,103)qparn_sites(j,i)
     .                  ,qparn_sites(1,i),(wght_xyz(j,k),k=1,j)
103       format(i6,i6,3e21.13)
        enddo

      enddo


c  now write out the apriori eop constraints - only for the estimated parameters
      do i=1,3
        do j =1,2
          np = qparn_eop(j,i)
          if(np.ne.0)then
c  need to correct the eop parameter numbers in the cases where the satellite
c  parameters are not written out. Calculate the number of sat parameters by
c  subtracting the number of satellite parameters from the eop param number
            write(isnx,101)np-sat_par_num,np-sat_par_num,cov_apr(np,np)
          endif
        enddo
      enddo

c  close the block
      write(isnx,'(a)')'-SOLUTION/MATRIX_APRIORI L COVA'
      call write_dash(isnx)

      return
      end



