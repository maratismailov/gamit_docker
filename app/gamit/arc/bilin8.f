      subroutine bilin8(x,y,x0,y0,xstep,ystep,u1,u2,u3,u4,val)

c   computes a bilinear interpolation of the value of a single
c   point within a rectangle. I found this formula via google
c   and it seems to work!
c
c   INPUT:
c         x,y  :  the coordinates of the point to be interpolated
c        x0,y0 :  the coords of the lower left corner
c     xstep,ystep : the grid step sizes
c         u1   : value at lower left corner
c         u2   : value at lower right corner
c         u3   : value at upper right corner
c         u4   : value at upper left corner
c
c    OUTPUT:
c        val   : the interpolated value at the requested point
c
c   P. Tregoning
c   7 January 2004

      implicit none

c  argument variables
      real*8 x,y,x0,y0,xstep,ystep,u1,u2,u3,u4,val

c  local variables
      real*8 dy,dx,dy1,dx1

          dy  = abs(y-y0)/ystep
          dy1 = (ystep - dy*ystep)/ystep
          dx  = (x - x0)/xstep
          dx1 = (xstep - dx*xstep)/xstep
c then      
          val = dx1*dy1*u1+dx*dy1*u2+dx1*dy*u4+dx*dy*u3
      return
      end





