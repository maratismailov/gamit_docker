      subroutine wrtcov(isnx,cov_parm,qnum_parn,qnum_sites,qparn_eop
     .                 ,qnum_eop)

c  Purpose: write out the vcv lower triangular matrix in sinex format
c
c  IN:
c       isnx         - unit number of sinex file              I*4
c       cov          - full covariance matrix                 R*8(qnum_parn,qnum_parn)
c       qnum_parn    - total number of parameters in hfile    I*4
c       qnum_sites   - number of sites observed               I*4
c       qparn_eop    - indices of eop parameters in cov_parm  I*4(2,3)
c       qnum_eop     - number of eop parameters estimated     I*4
c
c P Tregoning
c 12th August, 1995

      implicit none

      include '../includes/dimpar.h'

      integer isnx,qnum_parn,i,j,k,l,para2,qparn_eop(2,3),qnum_sites,np
     .        ,qnum_eop,para1,no_eop
      real*8 cov_parm(qnum_parn,qnum_parn),tmp_cov(maxsit*3+
     .                                            maxorb*maxsat+6)


c write the SOLUTION/MATRIX_ESTIMATE L COVA block to sinex file
      write(isnx,10)
10    format('+SOLUTION/MATRIX_ESTIMATE L COVA',/
     .      ,'*PARA1 PARA2 _____PARA2+0________ _____PARA2+1________ '
     .      ,'_____PARA2+2________')


c  now need to write out lower trianlge vcv - only allowed three columns of values(cols 3-5). Col
c  1 is the parameter number (= row number), col 2 is the output col number that the
c  output is outputting.

c  first write out the station coord part
      do i = 1, qnum_sites*3

        j = i
        para2 = 1

        do while (j.ge.3)
          write(isnx,100)i,para2,(cov_parm(i,para2+k-1),k=1,3)
100       format(i6,i6,3e21.13)
          j = j - 3
          para2 = para2 + 3
        enddo

c  now write out the cases where there are less than three elements
        if(j.eq.1)write(isnx,101)i,para2,cov_parm(i,para2)
101       format(i6,i6,e21.13)
        if(j.eq.2)write(isnx,102)i,para2,(cov_parm(i,para2+k-1),k=1,2)
102       format(i6,i6,2e21.13)

      enddo

c  now skip all the satellite entries and just write out the eop entries wrt themselves
c  and to the site coords

c first write the site/eop elements and eop/eop elements to a temp variable

c  keep track of how many eop parameters were NOT estimated
        no_eop = 0

        do i=1,3
          do j=1,2
            np = qparn_eop(j,i)
            if(np.ne.0)then
              do k = 1, qnum_sites*3
                tmp_cov(i+k-1-no_eop) = cov_parm(i+qnum_sites*3,k)
              enddo
              tmp_cov(qnum_sites*3+2*i-2+j-no_eop) = cov_parm(np,np)

              para1 = qnum_sites*3+2*i-2+j-no_eop
              l = para1
              para2 = 1

              do while (l.ge.3)
                write(isnx,100)para1,para2,(tmp_cov(para2+k-1),k=1,3)
                l = l - 3
                para2 = para2 + 3
              enddo

c  now write out the cases where there are less than three elements
              if(l.eq.1)write(isnx,101)para1,para2,tmp_cov(para2)
              if(l.eq.2)write(isnx,102)para1,para2,(tmp_cov(para2+k-1)
     .                                                        ,k=1,2)
            else
c   parameter was NOT estimated
              no_eop = no_eop + 1
            endif
          enddo
        enddo

c now close the block
      write(isnx,'(a)')'-SOLUTION/MATRIX_ESTIMATE L COVA'
      call write_dash(isnx)

      return
      end
