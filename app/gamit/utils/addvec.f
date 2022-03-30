      program addvec
c     add vectors and their uncertainties together



      real*8 x ,y ,cxx ,cyy ,cxy ,siga, sigr
      real*8 sigx ,sigy ,corxy,corra
      real*8 val(5),pi,caa,a,r,crr,cra,rmult
      real*8 c2g,s2g,aza,asq,bsq

      integer kode,i

      logical first

      first = .true.
      pi = 4.0d0*datan(1.0d0)
      x = 0.0d0
      y = 0.0d0
      r = 0.
      a = 0.
      cxx = 0.0d0
      cyy = 0.0d0
      cxy = 0.0d0
      cra = 0.0d0

 10   continue
c
      write(6,'(a)') 'Input format code and vector'
      write(6,'(a)') '  +/-1:  x sigx  y   sigy   cor'
      write(6,'(a)') '  +/-2:  x cxx   y   cyy    cxy'
      write(6,'(a)') '  +/-3:  r sigr  a   sigaz  cor'
      write(6,'(a)') '  +/-4:  r az    a   b      aza'
      write(6,'(a)') '  0 0 0 0 0 0 to stop'
      read(5,*) kode,(val(i),i=1,5)

      if( first .and. kode.eq.0 ) goto 999
      first = .false.

      if (kode .lt. 0) then
         rmult = -1.0d0
      else
         rmult = 1.0d0
      endif

      if (abs(kode) .eq. 1) then
         x    = x + rmult*val(1)
         y    = y + rmult*val(3)
         cxx  = cxx + val(2)**2
         cyy  = cyy + val(4)**2
         cxy  = cxy + val(2)*val(4)*val(5)
      else if (abs(kode) .eq. 2) then
         x    = x + rmult*val(1)
         y    = y + rmult*val(3)
         cxx  = cxx + val(2)
         cyy  = cyy + val(4)
         cxy  = cxy + val(5)
      else if (abs(kode) .eq. 3) then
c        assume azimuth is degrees clockwise from north
         r       = val(1)
         crr     = val(2)**2
         a       = val(3)*pi/180.0d0
         caa     = (val(4)*pi/180.0d0)**2
         cra     = val(2)*(val(4)*pi/180.0d0)*val(5)
         x       = x + rmult * r * dsin(a)
         y       = y + rmult * r * dcos(a)
         cxx     =  cxx
     .            + crr*(dsin(a)**2)
     .            + caa*(r*dcos(a))**2
     .            + 2.0d0*cra*r*dsin(a)*dcos(a)
         cyy     =  cyy
     .            + crr*(dcos(a)**2)
     .            + caa*(r*dsin(a))**2
     .            - 2.0d0*cra*r*dsin(a)*dcos(a)
         cxy     = cxy + crr*dsin(a)*dcos(a)
     .                           -caa*dsin(a)*dcos(a)*r*r
     .                           +cra*r*(dcos(a)**2 - dsin(a)**2)
      else if (abs(kode) .eq. 4) then
c        Clark et al. convention
c        assume azimuth is degrees clockwise from north
         r   = val(1)
         a   = val(2)*pi/180.0d0
         asq = val(3)**2
         bsq = val(4)**2
         aza = val(5)*pi/180.0d0
         x   = x + rmult * r * dsin(a)
         y   = y + rmult * r * dcos(a)
         c2g = dcos(2.d0*aza)
         s2g = dsin(2.d0*aza)
         cyy = cyy + (asq+bsq + (asq-bsq)*c2g) /2.d0
         cxx = cxx + (asq+bsq + (bsq-asq)*c2g) /2.d0
         cxy = cxy + (asq-bsq)*s2g /2.d0
      else if (kode .ne. 0) then
         write (6,'(a)') 'Invalid code'
      endif


      if (kode .ne. 0) then
         goto 10
      else
         sigx = dsqrt(cxx)
         sigy = dsqrt(cyy)
         corxy = cxy/(sigx*sigy)
         r    = dsqrt(x**2 + y**2)
         sigr = dsqrt(cxx*x*x + cyy*y*y + 2.0d0*cxy*x*y)/r
ckf      a=180.0d0*datan(x/y)/pi
         a=180.0d0*datan2(x,y)/pi
         caa  =(cxx + cyy*((x/y)**2) -2.0d0*cxy*x/y)/(y**2)
     .         / ((1.0d0+ (x/y)**2)**2)
         siga = dsqrt(caa)
         corra = cra/(sigr*siga)
         siga = 180.0d0 * siga/pi

         write(6,90) 'X','SIGX','Y','SIGY','CORRXY',
     .      'RATE','SIGRATE','AZ','SIGAZ','CORRRA'
  90     format (10(a7,1x))
         write(6,100) x,sigx,y,sigy,corxy,r,sigr,a,siga,corra
 100     format (4(f7.1,1x),f7.2,1x,4(f7.1,1x),f7.2,1x)
      endif

 999  stop
      end
