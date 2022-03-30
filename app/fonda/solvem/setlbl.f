      subroutine setlbl(idim,isit,mode)
c
c     set parameter label array
c        mode = 1: coordinates
c        mode = 2: velocity
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer idim,isit,mode,i1,j,j1
      character*16 temp(6),temp1(6)
      data temp/' geoc x ( m)    ',' geoc y ( m)    ',
     .           ' geoc z ( m)    ',' geod lat (dms) ',
     .           ' geod lon (dms) ',' height ( m)    '/
      data temp1/' Vx (mm/year)   ',' Vy (mm/year)   ',
     .           ' Vz (mm/year)   ',' Ve (mm/year)   ',
     .           ' Vn (mm/year)   ',' Vu (mm/year)   '/
c
c     every site has 12 parameter labels:
c        3 geocentric coordinates and 3 geocentric velocities
c        3 geodetic coordinates and 3 geodetic velocities
      i1 = 4
      if (idim.eq.3) i1 = 6
      do 10 j = 1,i1
         j1 = j
         if (idim.eq.2.and.j.gt.2) j1 = j+1
         label(j1)(1:8) = sname(isit)     
         if (mode.eq.1) label(j)(9:24) = temp(j1)
         if (mode.eq.2) label(j)(9:24) = temp1(j1)
 10   continue
c
      return
      end
