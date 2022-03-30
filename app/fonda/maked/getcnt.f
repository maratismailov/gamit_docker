      subroutine getcnt(i0,i1,xs,ys)
c
c     get geometric center coordinate. (2-D case)
c     All calculations are in the local geodetic coordinate.
c
c     mode = 1: without excluded sites
c     mode = 2: include only model sites
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
      integer i,i0,i1,is 
c
c     calculate geometric center's geodetic coordinate
      sumx = 0.0d0
      sumy = 0.0d0
      do 20 i = i0+1,i0+i1
      is = idgs(i)
      sumx = sumx+slo(is)
      sumy = sumy+sla(is)
 20   continue
      xs = sumx/dble(i1)
      ys = sumy/dble(i1)
c
 100  continue
      return
      end

