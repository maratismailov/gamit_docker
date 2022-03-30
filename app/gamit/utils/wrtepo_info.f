      subroutine wrtepo_info(isnx,qnum_sites,qsite_code,times)

c Purpose: Compute the mean epoch time of the observations at each site
c          and write out the SOLUTION/EPOCHS block
c
c IN:  isnx        - sinex file unit #                                     I*4
c      qnum_sites  - # sites in hfile                                      I*4
c      qsite_code  - 4 char site codes                                     C*4(qnum_sites)
c       times        -  runtime,start,stop,ic,eop,mean times (snx form)    C*12(6)
c                       where mean is the mean time of obs used in soln
c
c P Tregoning
c 11th August, 1995

      implicit none

      integer isnx,qnum_sites,i
      character*4 qsite_code(qnum_sites)
      character*12 times(6)


c  open the block
      write(isnx,100)
100   format('+SOLUTION/EPOCHS',/
     .      ,'*CODE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_')

c  now loop through the sites
      do i=1,qnum_sites

c now write out the records
        write(isnx,110)qsite_code(i),times(2),times(3),times(6)
110     format(' ',a4,'  A    1 P',3(' ',a12))

      enddo

c now close the block
      write(isnx,'(a)')'-SOLUTION/EPOCHS'
      call write_dash(isnx)

      return
      end
