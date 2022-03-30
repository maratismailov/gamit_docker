      subroutine read_ipgp_net(ifil1,ifil2,ifil3)
c
c     read IPGP site-table file and create FONDA format site-table file
c
c     assumes geocentric Cartesian coordinate (x,y,z)
c     All others are ignored.

c     unit:
c         x, y, z : km

      implicit real*8(a-h,o-z)

      include 'maked.fti'

      integer ifil1,ifil2,ifil3
      character*16 latstr,lonstr
      character*80 line
      character*8 code8
      integer*4 ierr
      integer icod,nsta,i
      real*8 p(7)
      real*8 xx,yy,zz
      real*8 xutm,yutm,zutm
      real*8 aa,bb,ep2,cc

c     grads to degrees
      g2r = 0.9d0

c     give JCR his constants from Dong's
      aa = radius
      bb = radius * (1.0d0 - 1.0d0/finv)
      ep2 = (aa*aa - bb*bb)/(bb*bb)
      cc = aa*aa/bb

c     there is usually a header.  Skip it
      read (ifil1,'(a)') line
c
c     write reference epoch
      write (ifil2,'(3x,a16,2x,f9.4)') 'reference epoch:',rtime

      ierr = 0

      write (ifil3,'(a)') ' Network distribution'
      icod = 1
      nsit = 0

      

      do 1000 while (icod .ne. 0) 
         read (ifil1,'(a)',iostat=ierr,end=9010) line
         read (line,*,iostat=ierr,err=9000) 
     .       ICOD,NSTA,xutm,yutm,zutm,ERPH,ERPZ
         if (icod .ne. 0) then
c           convert 8-digit code to 8-character code, padding
c           with zeros if necessary.
            write (code8,'(i8.8)') nsta

c           do this for UTM in Clark ellipsoid.
            if (.true.) then
               HE   = zutm
               CALL UTMGEO(XUTM,YUTM,XLON,YLAT,1)
c              convert from Geographic to geocentric XYZ
               CALL GEOCENTR(XLON,YLAT,HE,XX,YY,ZZ)

c              Transformation parameters from Clark 1880 Anglais
c              in Djibouti to WGS84
               TX     = 175.6081     
               TY     = 104.8617     
               TZ     = 215.3363     
               delta  = -5.5503d-5   
               RX     = 9.0e-9       
               RY     = 10.0d-9      
               RZ     = 11.0d-9      

               P(1)=TX
               P(2)=TY
               P(3)=TZ
               P(4)=DELTA
               P(5)=RX
               P(6)=RY
               P(7)=RZ

c              perform the transformation Clark to WGS84
               XT=XX+P(1)+P(4)*XX+P(5)*ZZ-P(6)*YY
               YT=YY+P(2)+P(4)*YY+P(6)*XX-P(7)*ZZ
               ZT=ZZ+P(3)+P(4)*ZZ+P(7)*YY-P(5)*XX
             endif

c           convert from cartesian X,Y,Z to lat,lon,rad and 
c           write to a priori coordinates file
            call SPHXYZ(ALAT,ALON,RADIAL,XT,YT,ZT,2)
            call wdms (1,alat,latstr)
            call wdms (2,alon,lonstr)
            WRITE(ifil2,57) code8,latstr,lonstr,radial
   57       FORMAT(A4,12x,1x,a16,1x,a16,1x,f12.4)

c           put the coordinates in storage
            nsit = nsit + 1
            x(nsit) = xt
            y(nsit) = yt
            z(nsit) = zt
            slo(nsit) = alon
            sla(nsit) = alat
            srad(nsit) = radial
            sname(nsit) = code8

c           write the map file in decimal degrees, zero velocity
c           first get ellipsoidal coordinates on current ellipsoid
            call CENTRGEO(cc,ep2,XT,YT,ZT,GLON,GLAT,GH)
            write (ifil3,40) glon*g2r,glat*g2r,0.,0.,code8
c40         format (1x,2f15.8,2f10.4,3x,a4,4x)
 40         format (1x,2f11.5,2f9.3,2x,a8)
         endif
 1000 continue 

c     quickly determine range of map
      slatmax = -500.0d0
      slonmax = -500.0d0
      slatmin =  500.0d0
      slonmin =  500.0d0
      do i = 1,nsit
         slatmax = max(slatmax,sla(i))
         slatmin = min(slatmin,sla(i))
         slonmax = max(slonmax,slo(i))
         slonmin = min(slonmin,slo(i))
      enddo

      print 2010,slonmin*rtod,slonmax*rtod,slatmin*rtod,slatmax*rtod
 2010 format (1x,'Map range: ',4f13.4) 
    

 9000 if (ierr .gt. 0) then
         print *,'READ_IPGP_NET: file error reading: '
         print *,line
         call ferror (ierr,6)
         stop 'READ_IPGP_NET: file error'
      endif

c     come here on end of file
 9010 continue

      return
      end

C-------------------------------------------------------------
      SUBROUTINE UTMGEO(XX,YY,LON,LAT,KH)
C --------------------------------------
C         TRANSFORMATION  UTM - GEOGRAPHIQUE
C           ENTREE :  COORD. UTM  EN  METRES
C           SORTIE :  COORD. GEOGRAPH. EN GRADES
C
C           HEMISPHERE NORD : KH=1
C                      SUD  : KH=-1
C
      
      IMPLICIT REAL*8 (A-H,L,O-Z)
      integer kh,ifu,mi
ckf      include 'maked.fti'

c     give JCR his constants from Dong's
c      a = radius
c      b = radius * (1.0d0 - 1.0d0/finv)
c      ep2 = (a*a - b*b)/(b*b)
c      cc = a*a/b

c     These constants are for Clarke 1880 (anglais)
      cc = 6400057.734886d0
      ep = 0.0827654283186422d0
      ep2 = ep*ep
      ifu = 38
      pi = 4.0d0 * datan(1.0d0)

      print*,'GEOCENTR: Clarke 1880 parameters Hardwired!'
c
      AK=0.9996d0
      KH=1

c     is this true everywhere?
      IFU = 38

      MI=IFU*6-183
      AMI=FLOAT(MI)/0.9D0
      RD=PI/200.D0
      AM0=AMI*RD
      YR=YY
      if(KH.EQ.1) goto 11
      YR=10.D6-YR
   11 YR=YR/AK
      P=KH*YR*PI*0.5D-7
      CO=DCOS(P)
      L=P*KH
      S=EP2*CO*CO
      R=CC/DSQRT(1.d0+S)
      PP=P*KH
      YR=YR-ARCME(L)
      U=(XX-5.D5)/AK/R
      V=S*U*U
      Q=YR/R*(1.D0-V/2.D0)
      U=U*(1.D0-V/6.D0)
      P=L+Q
      W=DEXP(U)
      AM=DATAN((W*W-1.D0)/(2.D0*W*DCOS(P)))
      V=DSIN(P)/DCOS(P)
      G=DATAN(V*COS(AM))
      V=1.D0+S-1.5D0*DSIN(L)*EP2*(G-L)*DCOS(L)
      LON=AM+AM0
      LAT=L+V*(G-L)
      LON=(AM+AM0)/RD
      LAT=LAT/RD*KH
c
      RETURN
      END
C-------------------------------------------------
       SUBROUTINE GEOCENTR(LON,LAT,H,XX,YY,ZZ)
C-----------------------
C---Passage coordonnees geographiques --> geocentriques
c         unites : grades et metres
c
      IMPLICIT REAL*8(A-H,L,O-Z)

ckf      include 'maked.fti'

c     give JCR his constants from Dong's
c      a = radius
c      b = radius * (1.0d0 - 1.0d0/finv)
c      ep2 = (a*a - b*b)/(b*b)
c      cc = a*a/b

c     These constants are for Clarke 1880 (anglais)
      cc = 6400057.734886d0
      ep = 0.0827654283186422d0
      ep2 = ep*ep
      pi = 4.0d0 * datan(1.0d0)
      rad = PI/200.d0
c
      print*,'GEOCENTR: Clarke 1880 parameters Hardwired!'

       LATR=LAT*RAD
       LONR=LON*RAD
       CB=DCOS(LATR)
       SB=DSIN(LATR)
       CL=DCOS(LONR)
       SL=DSIN(LONR)
C---Grande normale GN
c   Note: CC rayon de courbure au pole - EP2 = EPRIM **2
c         E = seconde excentricite - e2 =(a2-b2)/b2
       GN=CC/DSQRT(1.D0+EP2*CB*CB)
c
c---Coordonnees tridimensionnelles
c
       XX=(GN+H)*CB*CL
       YY=(GN+H)*CB*SL
       ZZ=(GN/(1.d0+EP2)+H)*SB
c
       RETURN
       END
C ------------------------------------------------
       REAL*8 FUNCTION ARCME(AL)
C--------------------------
C         Fonction Arc de mÇridien en fonction de la latitude :
c
      IMPLICIT REAL*8(A-H,L,O-Z)
c
ckf      include 'maked.fti'

c     give JCR his constants from Dong's
c      a = radius
c      b = radius * (1.0d0 - 1.0d0/finv)
c      ep2 = (a*a - b*b)/(b*b)
c      cc = a*a/b


c     These constants are for Clarke 1880 (anglais)
      cc = 6400057.734886d0
      ep = 0.0827654283186422d0
      ep2 = ep*ep

      print*,'GEOCENTR: Clarke 1880 parameters Hardwired!'


      AQ=EP2*3.d0/4.d0
      BQ=EP2*EP2*15.d0/16.d0
      CQ=(EP2**3)*35.d0/64.d0
      DQ=(EP2**4)*315.d0/512.d0
      U=DSIN(AL)*DCOS(AL)
      V=DCOS(AL)**2
      AJ2=AL+U
      AJ4=(AJ2*3.d0+U*V*2.d0)/4.d0
      AJ6=(AJ4*5.d0+2.d0*U*V**2)/3.d0
      AJ8=(AJ6*7.d0+4.d0*U*V**3)/8.d0
      ARCME=CC*(AL-AQ*AJ2+BQ*AJ4-CQ*AJ6+DQ*AJ8)
      RETURN
      END


C-------------------------------------------------
       SUBROUTINE CENTRGEO(cc,ep2,X,Y,Z,LON,LAT,H)
c---------------------------
c
c----Transformation de coordonnees geocentriques --> ellipsoidiques
c        unit‹s : m tres et grades
c
       IMPLICIT REAL*8 (A-H,L,O-Z)
ckf      COMMON/COMOB/CC,EP2,PI,RAD,IFU
C

       pi = 4.0d0 * datan(1.0d0)
       rad = PI/200.d0

       A=CC/DSQRT(1.D0+EP2)
       B=A*A/CC
       P= DSQRT(X**2+Y**2)
       PB=B*P
       ZA=Z*A
       TETA= DATAN2(ZA,PB)
       LBD= DATAN2(Y,X)
       ST=DSIN(TETA)
       CT=DCOS(TETA)
C
       FX=Z+EP2*B*ST*ST*ST
       FY=P-EP2*A*CT*CT*CT/(1.d0+EP2)
       FI= DATAN2(FX,FY)
       CFI=DCOS(FI)
       AN=CC/DSQRT(1.D0+EP2*CFI*CFI)
c
       LAT=FI*200.D0/PI
       LON=LBD*200.D0/PI
       H= P/CFI-AN
       H1=Z/DSIN(FI)-AN/(1.d0+EP2)
       RETURN
       END
