      program addvec
c     add vectors and their uncertainties together



      real*8 x ,y, z, cxx ,cyy ,czz, cxy ,cxz ,cyz,  siga, sigr
      real*8 sigx ,sigy ,sigz, corxy, corxz, coryz,corra
      real*8 val(9),pi,caa,a,r,crr,cra,rmult
      real*8 c2g,s2g,aza,asq,bsq

      integer kode,i

      pi = 4.0d0*datan(1.0d0)
      x = 0.0d0
      y = 0.0d0
      z = 0.0d0
      cxx = 0.0d0
      cxy = 0.0d0
      cyy = 0.0d0
      cxz = 0.0d0
      cyz = 0.0d0
      czz = 0.0d0

      write(6,'(a)') 'Input format code and vector'
      write(6,'(a)') '  +/-1:  x y z sigx sigy  sigz corxy corxz coryz'
      write(6,'(a)') '  +/-2:  x cxx   y   cyy    cxy'
      write(6,'(a)') '  +/-3:  r sigr  a   sigaz  cor'
      write(6,'(a)') '  +/-4:  r az    a   b      aza'
      write(6,'(a)') '  <CNTL> Z to stop'

 10   continue
      read(5,*,end=999) kode,(val(i),i=1,9)

      if (kode .lt. 0) then
         rmult = -1.0d0
      else
         rmult = 1.0d0
      endif

      if (abs(kode) .eq. 1) then
         x    = x + rmult*val(1)
         y    = y + rmult*val(2)
         z    = z + rmult*val(3)
         cxx  = cxx + val(4)**2
         cyy  = cyy + val(5)**2
         czz  = czz + val(6)**2
         cxy  = cxy + val(4)*val(5)*val(7)
         cxz  = cxz + val(4)*val(6)*val(8)
         cyz  = cyz + val(5)*val(6)*val(9)
         sigx = dsqrt(cxx)
         sigy = dsqrt(cyy)
         sigz = dsqrt(czz)
         if (sigx.ne.0.d0 .and. sigy.ne.0.d0) then
            corxy = cxy/(sigx*sigy)
         else
            corxy = 0.d0
         endif
         if (sigx.ne.0.d0 .and. sigy.ne.0.d0) then
            corxz = cxz/(sigx*sigz)
         else
            corxz = 0.d0
         endif
         if (sigx.ne.0.d0 .and. sigy.ne.0.d0) then
            coryz = cyz/(sigy*sigz)
         else
            coryz = 0.d0
         endif
         write(6,'(/,3f15.4,6f8.4)')
     .   x, y, z, sigx, sigy,  sigz, corxy, corxz, coryz
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
         cxx = cxx + (asq+bsq + (asq-bsq)*c2g) /2.d0
         cyy = cyy + (asq+bsq + (bsq-asq)*c2g) /2.d0
         cxy = cxy + (asq-bsq)*s2g /2.d0
      else if (kode .ne. 0) then
         write (6,'(a)') 'Invalid code'
      endif


      goto 10
 999  continue
C         sigx = dsqrt(cxx)
C         sigy = dsqrt(cyy)
C         corxy = cxy/(sigx*sigy)
C         r    = dsqrt(x**2 + y**2)
C         sigr = dsqrt(cxx*x*x + cyy*y*y + 2.0d0*cxy*x*y)/r
C         a=180.0d0*datan(x/y)/pi
C         caa  =(cxx + cyy*((x/y)**2) -2.0d0*cxy*x/y)/(y**2)
C     .         / ((1.0d0+ (x/y)**2)**2)
C         siga = dsqrt(caa)
C         corra = cra/(sigr*siga)
C         siga = 180.0d0 * siga/pi
C
C         write(6,90) 'X','SIGX','Y','SIGY','CORRXY',
C     .      'RATE','SIGRATE','AZ','SIGAZ','CORRRA'
C  90     format (10(a7,1x))
C         write(6,100) x,sigx,y,sigy,corxy,r,sigr,a,siga,corra
C 100     format (4(f7.1,1x),f7.2,1x,4(f7.1,1x),f7.2,1x)

      end
