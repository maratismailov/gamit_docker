c Program confpol.f
c Calculate 1-sigma confidence ellipse from sige, sign, 
c and correlation cor.
c sige converted to degrees large circle
c
c Mikhail Kogan LDEO
c 20-JUN-99
c Ref: 
c Richard Reyment and K.G. Joreskog, Applied factor analysis in
c the natural sciences, 371 pp., Cambridge University Press,
c 2nd edition, 1993, ISBN 0-521-41242-0
c See pages: 54-68 
c
c Modified by S. McClusky to read multiple input lines
c
c      implicit none

      integer*4 ierr,plate_num,iplate,jerr,indx,iel,trimlen,num

      real*8 a(2,2), lam1,lam2,grad,rad 
      real*8 lat(30),lat_sig(30),lon(30),lon_sig(30),mag(30)
     .     , mag_sig(30),rho_latlon(30),rho_latmag(30),rho_lonmag(30)
     .     , values(9),coslat,sig1,sig2,cor,sig1s,sig2s,a11,a12,a21,a22
     .     , suma,det,pm,u11,u12,ulen,az,ae,be,ae95,be95,aekm,bekm
     .     , ae95km,be95km

      character*14 plate_name(30)
      character*80 in_file
      character*132 buffer

      grad=180./3.141592653
      rad=1./grad
c
c      a(1,1)=208.
c      a(1,2)=144.
c      a(2,1)=144.
c      a(2,2)=292.
c
c      open(1, file='confpol.in', status='old')
c
c-- Construct var/covar matrix from sig1, sig2, and correlation cor
c
c      read(1,*) a(1,1),a(1,2)
c      read(1,*) a(2,1),a(2,2)
c       read(1,*) sign, sige, cor, alt
c       read(1,*) sign, sige, cor, alt

c Read the command line
      call rcpar(1,in_file)

c check if user knows how to use program.
      if ( in_file.eq." ") then
        print*,'Get a brain, give me an input file'
        print*,'Format (free)'
        print*,' PLATE -  PLATE  Lat (deg)  +-  Long (deg)  +-  Mag (de\
     .g/My)  +-   RhoLaLg RhoLaMa RhoLgMa'
        print*,' NUVEL-1A-AF-EU  21.0       6.0  -20.6      0.7 0.12   \
     .     0.02     0.001   0       0   LLM'
        stop
      endif

      open(10, file=in_file, iostat=ierr, status='old')
      if ( ierr .ne. 0 ) then
        print*,'File does not exist.',in_file
        stop
      endif
       
      plate_num = 1
      iplate = 0 
      jerr = 0  

c     Read the input data file extract available plate names
c     and their associated euler vectors 
10    do while ( jerr.eq.0 )
        read(10,'(a)', iostat=jerr ) buffer
c        See if we can find the plate-motion vector names
        indx = 1
        if( trimlen(buffer).gt.0 .and. buffer(1:1).eq.' ' ) then
          iplate = iplate + 1
          call getword( buffer, plate_name(iplate), indx)  
cd          write(*,'(i4,a,a)') iplate,'. ', plate_name(iplate)
          iel = 1
        else
          iel = 0
        end if
c       Decode the line and add new information if found
        if( iel.gt.0 ) then
            call multiread( buffer, indx, 'R8', ierr, values,
     .                      ' ',9)
          lat(iplate)     =  values(1)
          lat_sig(iplate) =  values(2)
          lon(iplate)     =  values(3)
          lon_sig(iplate) =  values(4)
          mag(iplate)     =  values(5)
          mag_sig(iplate) =  values(6)
          rho_latlon(iplate)   =  values(7)
          rho_latmag(iplate)   =  values(8)
          rho_lonmag(iplate)   =  values(9)
cd          print*,'input ',plate_name(iplate),lat(iplate),lon(iplate),
cd     .    mag(iplate),lat_sig(iplate),lon_sig(iplate),mag(iplate),
cd     .    rho_latlon(iplate),rho_latmag(iplate),rho_lonmag(iplate)
        endif
      enddo
      iplate=iplate - 1

      num = 1
c Start loop on plate found
20    do while ( num .le. iplate )
        coslat=cos(rad*lat(num))
cd        print *, ' cos(lat): ', coslat
        sig1=lon_sig(num)*coslat
        sig2=lat_sig(num)
        cor=rho_latlon(num)
c      
        sig1s=sig1**2
        sig2s=sig2**2
        a(1,1)= sig1s
        a(2,2)= sig2s
        a(1,2)= cor*sqrt(sig1s*sig2s)
        a(2,1)= a(1,2)
c
        a11=a(1,1)
        a12=a(1,2)
        a21=a(2,1)
        a22=a(2,2)
cd        print *, '--- Var/covar matrix ---'
cd        print *, ' a(1,j): ', a(1,1),a(1,2)
cd        print *, ' a(2,j): ', a(2,1),a(2,2)
c
c-- Calculate eigenvalues lam1, lam2
c and first eigenvector u11, u12
c  
        suma=a11+a22
        det=a11*a22-a12*a12
        pm=sqrt(suma*suma -4.*det)
c        print *, ' pm: ', pm
        lam1=0.5*(suma+pm)
        lam2=0.5*(suma-pm)
c
        u11=1.
c        u12= -(a11-lam1)/a12
        u12=  -a21/(a22-lam1)
        ulen=sqrt(u11**2 + u12**2)
        u11=u11/ulen
        u12=u12/ulen
        az= grad*atan2(u11,u12)
        ae=sqrt(lam1)
        be=sqrt(lam2)
c
c--- Conversion 1-sigma to 95% confidence in 2-D case, see:
c Numerical Recipes, 2nd edition (1994), p. 692, table on top:
c ni=2, sqrt(delta chi**2) = 2.45 (after interpolation 
c between 90% - 95.4% - 99%
        ae95 = ae*2.45
        be95 = be*2.45
c Values in km are required to plot with psxy, option -SE
        aekm = ae*111.2
        bekm = be*111.2
        ae95km = ae95*111.2
        be95km = be95*111.2
c
cd        print *, ' max eigenvalue: ', lam1
cd        print *, ' min eigenvalue: ', lam2
cd        print *, ' First eigenvector:  ', u11,u12
cd        print *, ' az: '  , az
cd        print 1000, ae,be
cd        print 1010, ae95, be95
cd        print 1020, ae95km, be95km
          write(*,1030)lon(num),lat(num),az,aekm,bekm
c          print*,lon(num),lat(num),az,aekm,bekm

c End loop on plates
        num = num + 1
      enddo   

c
c Formats:
c
1000  format (' Semiaxes [deg, 1-sig] : ', 2x,f6.3, 2x, f6.3)
1010  format (' Semiaxes [deg, 95%]   : ', 2x,f6.3, 2x, f6.3)
1020  format (' Semiaxes [km, 95%]    : ', 2x,f6.3, 2x, f6.3)
1030  format (1x,2(f7.4,1x),1x,3(f7.4,1x))

      stop
      end   
