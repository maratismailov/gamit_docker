      subroutine read_apr(ih,qnum_sites,qnum_sats,qnum_eop,num_orb
     .                   ,qnum_parn,qscale,qparn_eop,cov_apr)


c  Purpose: read in the apriori vcv elements for the estimated parameters from an hfile
c           in order to output them to a sinex file.
c
c  IN:
c        ih            - unit number of hfile                         I*4
c        qnum_sites    - number of sites estimated in hfile           I*4
c        qnum_sats     - number of satellites in hfile                I*4
c        qnum_eop      - number of eop parameters estimated           I*4
c        num_orb       - number of orbital elements/sat estimated     I*4
c        qnum_parn     - total number of parameters estimated         I*4
c        qscale        - vector to scale apriori constraints          R*8(qnum_parn)
c        qparn_eop     - indices for eop parameter numbers            I*4(2,3)
c  OUT:
c        cov_apr       - a priori covariance matrix                   R*8(qnum_parn,qnum_parn)
c

      implicit none

      integer ih,qnum_sites,qnum_sats,qnum_eop,num_orb,qnum_parn,i,j
     .       ,par1,par2,par_diff,par_old,np,qparn_eop(2,3)

      real*8 cov_apr(qnum_parn,qnum_parn),qscale(qnum_parn),value
      character*80 line

      print*,' Reading apriori constraints ....'

c  skip 1 line
      read(ih,'(a)')line

c  first zero out the lower triangular matrix
      do i=1,qnum_parn
        do j=1,i
          cov_apr(i,j) = 0.d0
        enddo
      enddo

c now read in all the site apriori variances. The units will be radians, radians and km (for
c coords which are assumed to be spherical. There are no off diagonal elements for site
c apriori constraints
      i=1
      do while (i.le.qnum_sites*3)
         read(ih,'(14x,e22.16)')cov_apr(i,i)
        i=i+1
      enddo

c get the satellite apriori variances
      i=0
      par_old = 0
      par_diff = 0
      do while (i.lt.qnum_sats*num_orb)
        read(ih,'(a)')line
        read(line,'(2(1x,i4),d26.16)')par1,par2,value

c  need to compute the offset in parameter numbers due to the ambiguities. The
c  first satellite parameter is number qnum_sites*3+1, but it won't be known as
c  that in the hfile
        if(i.eq.0)then
          par_diff = par1 - (qnum_sites*3+1)
        endif

c  now correct the parameter numbers, and write the value to the appropriate
c  location in the cov_apr matrix
        cov_apr(par1-par_diff,par2-par_diff) = value

c  only increment i if this par1 is different to the last time
        if(par1-par_diff.ne.par_old)then
          i=i+1
          par_old = par1-par_diff
        endif

      enddo

c get the eop estimates - but only if eop has been estimated!
      if(qnum_eop.gt.0)then
        do i=1,3
          do j =1,2
            np = qparn_eop(j,i)
            if(np.ne.0)then
              read(ih,'(14x,e22.16)')cov_apr(np,np)
            endif
          enddo
        enddo
      endif

c  now rescale the apriori into the right units
      do i=1,qnum_parn
        do j=1,qnum_parn
           cov_apr(i,j) = qscale(i)*cov_apr(i,j)
           cov_apr(j,i) = qscale(i)*cov_apr(j,i)
        enddo
      enddo

      return
      end
