      subroutine getsit(idatf,iprif,idim)
c
c     get working site table 
c
c     basic logics:
c         1. Use SV6 coordinate as our reference frame.  All historical
c            geocentric coordinate should be transfered to SV6.
c         2. either lat-long-radius or x-y-z coordinate input.  but all
c            calculations are refered to xyz coordinate.
c   3. later consider append mode.!!!
c
c     coordinate system : (comode)
c         1. geodetic (ellipsoidal) coordinate (lat,lon,radius)
c         2. geocentric coordinate (x,y,z)
c         3. topocentric coordinate (x(east),y(north),z(height))
c         Usually, we name the geocentric coordinate as 'global' frame.
c         Meanwhile, the topocentric coordinate is named as 'local' frame.
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer idatf,iprif,idim,i
      character*8 sname2(maxsit)
      character*8 sitnam
c     character*12 fullnm(maxsit)
c
c     read site number
      read (idatf,*) nsit
c
c     read site name(id)
      do i = 1,nsit
         if (fcode.le.2) then
            read (idatf,'(a8)') sitnam
            sname(i) = sitnam
c           print*, 'sname(',i,')' ,sitnam
            if (fcode.eq.1) sname2(i) = sitnam
         endif
         if (fcode.eq.3) then
            read (idatf,'(a8,1x,a5)') sname(i),sname2(i)(1:5)
         endif
      enddo
c
c     get SV6 coordinate for these sites
      call getpri(iprif,sname2)
c     print*,(fullnm(ii),ii=1,maxsit)
c
      return
      end
