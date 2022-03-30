Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      Subroutine EPTIME( scr)
c
c     get the UTC time for every epoch
c                        
c     rwk 190516: input iepoch now in solve.h 

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      character*3 buf3
      character*6 buf6
      character*12 scr
c                     
      integer*4 i,k
      real*8 tmp
      real*8 tmptim(3)
c
      do 1 i=1,2
1     tmptim(i)=dble(t00(i))
c
      tmptim(3)=t00(3)+dble(iepoch-1)*inter
c
      do 100 k=3,2,-1
c
         tmp=tmptim(k)
         tmptim(k)=dmod(tmp+6000.d0,60.d0)
         tmp=tmp-tmptim(k)
         tmptim(k-1)=tmptim(k-1)+tmp/60.d0
100   continue
c
      tmptim(1)=dmod(tmptim(1)+2400.d0,24.d0)
c
      write(unit=buf3,fmt='(f3.0)') tmptim(1)
      read(unit=buf3,fmt='(a3)') scr(1:3)
      write(unit=buf3,fmt='(f3.0)') tmptim(2)
      read(unit=buf3,fmt='(a3)') scr(4:6)
      write(unit=buf6,fmt='(f6.3)') tmptim(3)
      read(unit=buf6,fmt='(a6)') scr(7:12)
      scr(3:3)=':'
      scr(6:6)=':'
c
      return
      end
