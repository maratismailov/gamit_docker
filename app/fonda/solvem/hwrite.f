      subroutine hwrite(ifile,nobs,version,drvfil,kobs,time1,time2)
c
c     write solution and covariance matrix to global h-file
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      character*1 sp
      character*16 uname
      character*64 drvfil
      character*80 line
      character*16 version
      integer ifile,i,lcmd,iyr,imon,iday,ihr,imn,isec,ijunk,i1,j,k
      integer ik,id1,i2,i3,isit,id,iid
      integer nobs
      integer iyr1, imo1, idy1, ihr1, min1
      integer iyr2, imo2, idy2, ihr2, min2
      integer iyr3, imo3, idy3, ihr3, min3
      integer idoy1,idoy2,idoy3
      integer kobs(maxsit)
      real*8 sec1,sec2,sec3

      dimension iid(12),temp(12)
      
      sp = ' '
c
c     header
      write (ifile,'(1x)')
      write (ifile,'(25x,a)') 'GLOBALT file'
      call getcmd(8,'output',line,lcmd,1)
      write (ifile,'(1x,a,2x,a)') 'Parallel output file:',line(1:lcmd)
      write (ifile,'(1x,a,1x,a)') 'SOLVEM version:',version
      call getdat(iyr,imon,iday)
      call gettim(ihr,imn,isec,ijunk)
      call getusr(uname)
      write (ifile,20) 'Running time:',iyr,
     .         '/',imon,'/',iday,ihr,':',imn,':',isec
 20   format(1x,a13,2x,i4,2(a1,i2),3x,2(i2,a1),i2)
      write (ifile,*) sp,'Operator:',sp,uname
      write (ifile,'(1x,a,1x,a8)') 'M-file name:',drvfil(1:8)
      write (ifile,'(1x,a)') 'Datum: cartesian coordinates'
      call getcmd(8,'loose',line,lcmd,1)
      write (ifile,'(1x,a,1x,a)') 
     .'keys: loose constraints (m, m/yr)',line(1:lcmd)
ckf   write (ifile,'(1x,a)') 'keys:'
ckf   write (ifile,'(1x,a)') 
ckf  .   'Assumed standard deviation of measurement error:'
ckf   write (ifile,'(1x,a)') 'Assumed ionosphere error:'

c     put these in for future use.
      write (ifile,'(1x,a)') 'N astronomic azimuth'
      write (ifile,'(1x,a)') 'N horizontal angle'
      write (ifile,'(1x,a)') 'N horizontal direction'
      write (ifile,'(1x,a)') 'N baseline length'
      write (ifile,'(1x,a)') 'N zenith height'
      write (ifile,'(1x,a)') 'N leveling'
      write (ifile,'(1x,a)') 'N astronomic azimuth rate'
      write (ifile,'(1x,a)') 'N horizontal angle rate'
      write (ifile,'(1x,a)') 'N horizontal direction rate'
      write (ifile,'(1x,a)') 'N baseline length rate'
      write (ifile,'(1x,a)') 'N zenith height rate'
      write (ifile,'(1x,a)') 'N leveling rate'
      write (ifile,'(1x,a)') 'N 3-D geocentric coordinate'
      write (ifile,'(1x,a)') 'N 3-D geocentric velocity'
      write (ifile,'(1x,a)') 'N 3-D geodetic coordinate'
      write (ifile,'(1x,a)') 'N 3-D geodetic velocity'
      write (ifile,'(1x,a)') 'N 3-D geocentric coord. diffs'
      write (ifile,'(1x,a)') 'N 3-D geocentric vel. diffs'
      write (ifile,'(1x,a)') 'N 3-D geodetic coord. diffs'
      write (ifile,'(1x,a)') 'N 3-D geodetic vel. diffs'
      write (ifile,'(1x,a)') 'N deflection'
      write (ifile,'(1x,a)') 'Assumed atmospheric model'
      write (ifile,'(1x,a)') 'Assumed scale factor'

      write (ifile,'(1x,a)') 'Number of sessions:   1'
      write (ifile,'(1x,a,1x,i4)') 'Number of stations:',gdsit
ckf   write (ifile,'(1x,a)') 'Number of svs:    0'
ckf   write (ifile,'(1x,a)') 'Ephemeris files: '
      write (ifile,'(1x,a)') 'Name of track stations'
      do i = 1,gdsit
         j = itoj(i)
         write(ifile,23) i,sname(j)(1:4),sname(j),kobs(j)
 23      format(1x,i3,2x,a4,2x,a12,3x,i6)
      enddo
      write (ifile,'(1x,a)') 'Data files:'
      call getcmd(8,'input',line,lcmd,1)
      write (ifile,'(1x,i3,2x,a)') 1,line(1:lcmd)
ckf   write (ifile,'(1x,a)') 'Satellite used:'
      write (ifile,'(1x,a,i10)') 'Number of observations:',nobs

c     start time 
      iyr1 = int(time1)
      idoy1 = int ((time1 - iyr1) * 365.25) + 1
      call cdate (iyr1,idoy1,imo1,idy1)
c     assume noon (why not?)
      ihr1 = 12
      min1 = 0
      sec1 = 0.d0

c     end time 
      iyr2 = int(time2)
      idoy2 = int ((time2 - iyr2) * 365.25) + 1
      call cdate (iyr2,idoy2,imo2,idy2)
c     assume noon (why not?)
      ihr2 = 12
      min2 = 0
      sec2 = 0.d0
 
c     IC time is reference time for solution.
      iyr3 = int(rtime)
      idoy3 = int ((rtime - iyr3) * 365.25) + 1
      call cdate (iyr3,idoy3,imo3,idy3)
c     assume noon (why not?)
      ihr3 = 12
      min3 = 0
      sec3 = 0.d0

      write(ifile,'(a,4i4,2x,i3,2x,f7.3)')
     .   ' Start time: ', iyr1, imo1, idy1, ihr1, min1, sec1

      write(ifile,'(a,4i4,2x,i3,2x,f7.3)')
     .   ' End time  : ', iyr2, imo2, idy2, ihr2, min2, sec2
      write(ifile,'(a,4i4,2x,i3,2x,f7.3)')
     .   ' Ref time  : ', iyr3, imo3, idy3, ihr3, min3, sec3

ckf   write (ifile,'(1x,a,,2x,f10.4)') 'ICs time  :',rtime
ckf   write (ifile,'(1x,a,2x,f10.4)') 'IC epoch(day & second) :'
ckf   write (ifile,'(1x,a)') 'E-rotation epoch(day & second) :'
ckf   write (ifile,'(1x,a)') 'UT1(sec) & rate(sec/day) :'
ckf   write (ifile,'(1x,a)') 'X pole(asec) & rate(asec/day) :'
ckf   write (ifile,'(1x,a)') 'Y pole(asec) & rate(asec/day) :'
ckf   write (ifile,'(1x,a)') 'Delta-psi(asec) & rate(asec/day) :'
ckf   write (ifile,'(1x,a)') 'Delta-eps(asec) & rate(asec/day) :'
      write (ifile,'(1x,a,2x,i6,5x,a,2x,i6)') 
     .   'Total parameters:',gdsit*6,'live parameters:',nlive
      write (ifile,'(1x,a)') 'Prefit:'
      write (ifile,'(2x)') 
      write (ifile,'(8x,a,17x,a,16x,a)') 'label (units)',
     .    'a priori','adjustment'
      i1 = 0
      do isit = 1,gdsit
         i = itoj(isit)
         call setlbl(3,i,1)
         do 25 j = 1,6
            temp(j) = 0.0d0
            iid(j) = (i-1)*6+j
            id = map(iid(j))
            iid(j+6) = id
            if (id.gt.0) temp(j) = bnorm(id)
 25      continue
         i1 = i1+1
         if (iid(7).gt.0) then
            write(ifile,40) i1,label(1),x(i),temp(1)
         else
            write(ifile,50) i1,label(1),x(i)
         endif
         i1 = i1+1
         if (iid(8).gt.0) then
            write(ifile,40) i1,label(2),y(i),temp(2)
         else
            write(ifile,50) i1,label(2),y(i)
         endif
         i1 = i1+1
         if (iid(9).gt.0) then
            write(ifile,40) i1,label(3),z(i),temp(3)
         else
            write(ifile,50) i1,label(3),z(i)
         endif
         i1 = i1+1

c        get velocity i1,label
c        omit velocities for the moment. Kurt
c        leave in m/yr Kurt
         call setlbl(3,i,2)
         do j = 1,3
            k = index(label(j),'mm')
            label(j)(k:k+1) = ' m'
         enddo
            
         if (iid(10).gt.0) then
            write(ifile,40) i1,label(1),vx(i) ,temp(4) 
         else
            write(ifile,50) i1,label(1),vx(i) 
         endif
         i1 = i1+1
         if (iid(11).gt.0) then
            write(ifile,40) i1,label(2),vy(i) ,temp(5) 
         else
            write(ifile,50) i1,label(2),vy(i) 
         endif
         i1 = i1+1
         if (iid(12).gt.0) then
            write(ifile,40) i1,label(3),vz(i) ,temp(6) 
         else
            write(ifile,50) i1,label(3),vz(i) 
         endif
      enddo
 40   format (1x,i4,1h*,a24,1x,d23.16,4x,d23.16)
 50   format (1x,i4,1x,a24,1x,d23.16)
c
c     output covariance matrix
      write (ifile,'(2x)') 
      write (ifile,'(1x,a)') 'Covariance matrix:'
      i1 = 0
      do isit = 1,gdsit
         i = itoj(isit)
         do 60 ik = 1,6
            id1 = (i-1)*6+ik
            id = map(id1)
            if (id.le.0) goto 60
            i1 = i1+1
            i2 = i1*(i1-1)/2+1
            i3 = i1*(i1+1)/2
            write(ifile,70) i1,(anorm(j),j=i2,i3)
 60      continue
      enddo
 70   format (1x,i3,'. ',50(5(1x,d23.16),:,/,6x))

c
 100  continue

c     end with a blank line
c     write (ifile,'(80x)')

      return
      end

