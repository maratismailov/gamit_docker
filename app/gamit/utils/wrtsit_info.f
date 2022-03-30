      subroutine wrtsit_info(isnx,istn,qnum_sites,qsite_code,rcvr
     .                     ,antenna,soln_time )

c  Purpose: write to sinex file the site/receiver information. This requires
c           the reading of a station.info file in order to get the start date
c           at which the receiver/firmware was started (required in sinex for
c           some reason!)
c
c IN:   isnx        -  sinex unit #                               I*4
c       istn        -  station.info unit number                   I*4
c       qnum_sites  -  number of sites in hfile                   I*4
c       qsite_code  -  4 char names of sites                      C*4(qnum_sites)
c       rcvr        -  receiver name, s/# and firmware            C*16(maxsit,3)
c       antenna     -  antenna type, s/#                          C*16(maxsit,2)
c       soln_time   -  epoch to which soln refers (EOP time)      C*12
c
c  P Tregoning
c  10th August, 1995

      implicit none

      include '../includes/dimpar.h'

      integer isnx,istn,qnum_sites,i
      character*4 qsite_code(qnum_sites)
      character*12 soln_time
      character*16 rcvr(maxsit,3),antenna(maxsit,2)

c  open the SITE/RECEIVER block
      write(isnx,100)
100   format('+SITE/RECEIVER',/
     . ,'*Code PT _Occ T _Data Start_ _Data End___ _receiver type__ '
     . ,'_S/N_ ___firmware____')


c  loop through the sites
      do i=1,qnum_sites

c  PT950810: Kluge while we don't have serial numbers
          if(rcvr(i,2)(1:5).eq.'     ')rcvr(i,2)(1:5)='-----'

c  write out the records
          write(isnx,110)qsite_code(i),soln_time,rcvr(i,1),rcvr(i,2)
     .                  ,rcvr(i,3)
110       format(' ',a4,'  A    1 P ',a12,' 00:000:00000 ',a16,1x,a5,1x
     .          ,a15)
      enddo

c  close the block
      write(isnx,'(a)')'-SITE/RECEIVER'
      call write_dash(isnx)


c  open the SITE/ANTENNA block
      write(isnx,120)
120   format('+SITE/ANTENNA',/
     . ,'*Code PT _Occ T _Data Start_ _Data End___ __antenna type__ '
     . ,'_S/N_')

c  loop through all the sites again
      do i=1,qnum_sites

c  PT950810: Kluge while we don't have serial numbers
          if(antenna(i,2)(1:5).eq.'     ')antenna(i,2)(1:5)='-----'

c  write out the records
          write(isnx,130)qsite_code(i),soln_time,antenna(i,1)
     .                   ,antenna(i,2)
130       format(' ',a4,'  A    1 P ',a12,' 00:000:00000 ',a16,1x,a5)

      enddo

c  close the block
      write(isnx,'(a)')'-SITE/ANTENNA'
      call write_dash(isnx)


      return
      end
