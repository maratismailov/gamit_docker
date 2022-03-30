      subroutine strain (idst,vcov,gamma1,gamma2,g1sig,g2sig
     .                   ,xv,cx,eps1,eps1sig,eps2,eps2sig
     .                   ,theta,thetasig
     .                   ,spin,spinsig,nstn,lerror,ch2)
c D. C. Cong 970103
 
c      subroutine strain (idst,vcov,gamma1,gamma2,g1sig,g2sig
c     .                   ,xv,eps1,eps1sig,eps2,eps2sig
c     .                   ,theta,thetasig
c     .                   ,spin,spinsig,nstn,lerror)
 
c     root : strain.f from Kurt's directory
c     modified to fit FONDA frame by Dong 910701
c
c     given the coordinates and velocities of N points,
c     estimate the components horizontal velocity gradient,
c     and derive the strain rate tensor from this
 
      include 'solvem.fti'
 
c     maximunm number of stations = jnet (in fti file)
c     actual number of stations
      integer nstn,jnet,jnsit,i1,i2
c     maximum number of stations = jnsit (local) (in fti file)
      parameter (jnet=maxsit)
      parameter (jnsit=maxsit)
 
c INPUT
c     site index array
      integer idst(jnet)
c     latitude, longitude in degrees
      real*8 dlats(jnet),dlons(jnet)
c     latitude, longitude of centroid in degrees
      real*8 dlat0,dlon0
c     velocities in mm/yr
c     slnve, slnvn (in solvem.fti)
c     covariances in (mm/yr)(e1,n1,e2,n2,e3,n3,... eN,nN)
      real*8 vcov(maxprm*6)
c OUTPUT
c     eigenvalues and their uncertainties
      real*8 eps1,eps1sig,eps2,eps2sig
c     EE,EN,NN components of strain rate
      real*8 edot11,edot12,edot22
c     clockwise spin rate and sigma in rad/yr
      real*8 spin,spinsig
c     components of shear
      real*8 gamma1,gamma2,g1sig,g2sig
c     azimuth of more compressive eigenvector (degrees) and uncertainty
      real*8 theta,thetasig
 
c D. C. Cong 970103:
      real*8 cx(21),ch2 
 
* LOCAL VARIABLES:
c     relative distances in kilometers
      real*8 delx,dely
c     design matrix and its inverse
 
c     real*8 aat(12*jnsit),cx(21),cx2(21)
 
c D. C. Cong 970103:
      real*8 aat(12*jnsit),cx2(21)
      real*8 ax(2*jnsit),vv(2*jnsit),vpv(1),temp(2*jnsit)
 
      real*8 ctest(21),bigg(12*jnsit)
c     data, model, and derived parameter vectors
      real*8 bv(2*jnsit),xv(6),x2v(6)    
c     second derivative for descriminatinon
      real*8 d2
c     .true. on disaster
      logical lerror
      equivalence (bigg,gvm)
      equivalence (aat,gmvm)
      equivalence (bv,strat)
     
c     others
      real*8 rat0,ron0,faz,baz,ds,dgdl11,dgdl21,dgdl12,dgdl22
     .,      term,term1,term2,term3,rat1,ron1
 
      integer i,j,irow
c
      lerror = .false.
c
c     get the site coordinate and velocity from index array
      do i = 1,nstn
         i1 = itoj(idst(i))
         dlats(i) = slat(i1)*rtod
         dlons(i) = slon(i1)*rtod
      enddo
c     find the centroid (dlat0,dlon0) in degrees
      call centroid (dlats,dlons,dlat0,dlon0,nstn)
c     convert to degrees
      rat0 = dlat0*dtor
      ron0 = dlon0*dtor
 
      call zero1d(1,12*nstn,aat)
      i2 = 2*nstn
      do i = 1,nstn
         i1 = idst(i)
c        get the delta x and y in kilometers
         rat1 = dlats(i)*dtor
         ron1 = dlons(i)*dtor
         call HELINV(rat0,ron0,rat1,ron1,faz,baz,ds)
         delx = ds * sin(faz)/1000.
         dely = ds * cos(faz)/1000.
 
c        form the A~ matrix with 6 rows
c        odd numbered row (for A~, it is column)
         irow = 2*i - 1
         aat(irow) = 1.0d0
         aat(2*i2+irow) = delx
         aat(3*i2+irow) = dely
c        E velocity in mm/yr
         bv(irow) = slnve(i1)*1.d3
 
c        even numbered row
         irow = 2*i 
         aat(i2+irow) = 1.0d0
         aat(4*i2+irow) = delx
         aat(5*i2+irow) = dely
c        N velocity in mm/yr
         bv(irow) = slnvn(i1)*1.d3
      enddo
 
c     inverse covariance matrix to get weighting matrix (vcov)       
      call cholsk(vcov,ctest,1,2*nstn,j)
      if (j.gt.100) then 
         lerror = .true.
         print*,' error in subnet: ',sname(itoj(idst(1))),'-',
     .      sname(itoj(idst(2))),'-',sname(itoj(idst(3)))
         goto 100
      endif
     
c     normal matrix a~wa = cx
      call atwa(6,2*nstn,aat,vcov,cx,bigg,1)
c     a~w
      call axb(6,2*nstn,2*nstn,aat,vcov,bigg,3,0)
c     right hand term a~wb
      call axb(6,2*nstn,1,bigg,bv,xv,1,0)
c     get solution
      call cholsk(cx,xv,3,6,j)
c     print*,'CD parameter: Veo,Vno,dVe/dxe,dVe/dxn,dVn/dxe,dVn/dxn'
c     write (*,'(6(1pe14.6,1x))') (xv(k),k=1,6)
 
c D. C. Cong 970103:
c     calculate chi-square: ch2
      call axb(2*nstn,6,1,aat,xv,ax,1,1)
      do i=1,2*nstn
         vv(i) = ax(i) - bv(i)
      enddo
      call atwa(1,2*nstn,vv,vcov,vpv,temp,1)
      ch2 = vpv(1)
 
c     calculate strain rate tensor and rotation     
c     propogation looks like:
c     0 0  1  0   0  0  Edot11
c     0 0  0 .5  .5  0  Edot12
c     0 0  0  0   0  1  Edot22
c     0 0  0 .5 -.5  0  omega
c     0 0  1  0   0 -1  gamma1
c     0 0  0  1   1  0  gamma2
      
      call zero1d(1,36,bigg)
c     row 1 corresponds to (1,1) component of strain rate
      bigg(3) = +1.00d0
c     row 2 corresponds to (1,2) component of strain rate
      bigg(10) = +0.50d0
      bigg(11) = +0.50d0
c     row 3 corresponds to (2,2) component of strain rate
      bigg(18) = +1.0d0
c     row 4 corresponds to spin rate (clockwise)
      bigg(22) = +0.50d0
      bigg(23) = -0.50d0
c     row 5 corresponds to gamma 1
      bigg(27) =  1.0d0   
      bigg(30) = -1.0d0
c     row 6 corresponds to gamma 2
      bigg(34) =  1.00d0
      bigg(35) =  1.00d0
 
      call axb(6,6,1,bigg,xv,x2v,1,0)
      call atwa(6,6,bigg,cx,ctest,aat,1)
 
c     spin rate and uncertainty
      spin   = x2v(4)/1.0d6
      spinsig= sqrt(ctest(10))/1.0d6
 
c     shear rates
      gamma1 = x2v(5)
      gamma2 = x2v(6)
      g1sig = sqrt(ctest(15))
      g2sig = sqrt(ctest(21))
 
c     components
      edot11 = x2v(1)
      edot12 = x2v(2)
      edot22 = x2v(3)
 
c     partials of gamma with respect to components of vel. grad. L
      term = .5/sqrt((xv(3)-xv(6))**2 + (xv(4)+xv(5))**2)
      dgdl11 = term * 2.0d0 *  (xv(3) - xv(6))
      dgdl12 = term * 2.0d0 *  (xv(4) + xv(5))
      dgdl21 = term * 2.0d0 *  (xv(4) + xv(5))
      dgdl22 =-term * 2.0d0 *  (xv(3) - xv(6))
 
c     partials of epsilon 1 with respect to comps of vel. grad. L
      call zero1d(1,36,bigg)
c     row 1 corresponds to epsilon 1
      bigg(3) = +1.00d0 * (.50d0 + .50d0 * dgdl11)
      bigg(4) = +1.00d0 * (        .50d0 * dgdl12)
      bigg(5) = +1.00d0 * (        .50d0 * dgdl21)
      bigg(6) = +1.00d0 * (.50d0 + .50d0 * dgdl22)
c     row 2 corresponds to epsilon 2
      bigg(9) = +1.00d0 * (.50d0 - .50d0 * dgdl11)
      bigg(10) = -1.00d0 * (        .50d0 * dgdl12)
      bigg(11) = -1.00d0 * (        .50d0 * dgdl21)
      bigg(12) = +1.00d0 * (.50d0 - .50d0 * dgdl22)     
c     row 3 corresponds to theta
      term1 = (xv(4)+xv(5))/(xv(6)-xv(3))
      term2 = 1.d0 + term1**2
      term3 = xv(6)-xv(3)
      bigg(15) = 0.50d0*term1/term3/term2
      bigg(16) = 0.50d0      /term3/term2
      bigg(17) = 0.50d0      /term3/term2
      bigg(18) =-0.50d0*term1/term3/term2
 
c     row 4 corresponds to maximum shear
      bigg(21) = dgdl11
      bigg(22) = dgdl12
      bigg(23) = dgdl21
      bigg(24) = dgdl22
 
      call atwa(4,6,bigg,cx,cx2,aat,1)
 
c     but of course, these parameters are not linear, so:
c     calculate epsilon 
      x2v(1) = .50d0 * (xv(3) + xv(6)
     .           + sqrt((xv(3)-xv(6))**2 + (xv(4)+xv(5))**2))
c     calculate epsilon 2
      x2v(2) = .50d0 * (xv(3) + xv(6) 
     .           - sqrt((xv(3)-xv(6))**2 + (xv(4)+xv(5))**2))
 
c     2) formula from L
      x2v(3) = .50d0 * datan((xv(4)+xv(5))/(xv(6)-xv(3)))
 
c     calculate gamma dot
      x2v(4) = sqrt((xv(3)-xv(6))**2 + (xv(4)+xv(5))**2)
c       print*,'2-D non-linear solution:'
c       write (*,'(6(1pe14.6,1x))') (x2v(k),k=1,6)
 
c     eigenvalues of strain rate (1/yr)
      eps1 = x2v(1)/1.0d6
      eps2 = x2v(2)/1.0d6
      eps1sig = dsqrt(cx2(1))/1.0d6
      eps2sig = dsqrt(cx2(3))/1.0d6
c     azimuth of more compressive eigenvector in degrees
      theta = x2v(3)
c     calculate second derivative
      d2  = 2.*(edot11-edot22)*cos(2.*theta)-4.*edot12*sin(2.*theta)
      if (d2 .lt. 0) theta = theta + 2.*atan(1.0)
      theta = theta / dtor
 
c      if (theta .lt. 360.) theta = theta + 360.
c      if (theta .gt. 360.) theta = theta - 360.
      thetasig = dsqrt(cx2(6))/dtor
 
c     unit = 1/yr
      edot11 = edot11*1.0d-6
      edot12 = edot12*1.0d-6
      edot22 = edot22*1.0d-6
 
 100  continue
      return
      end

