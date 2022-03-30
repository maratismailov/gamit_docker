Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      Subroutine COVBLN(istat,jstat,sigs,baseln,dsig,iflg,siglcl)
c
c     calculate covariance matrix of baseline
c
c Input:
c     istat     = index of first site
c     jstat     = index of second site
c     In parameters.h:
c     postvl    = postfit geocentric coordinates
c     sigma     = parameter uncertainties (unscaled)
c     free      = array indicating if parameters are adjusted or free
c Output:
c     sigs(*,1) = baseline
c     sigs(*,2) = uncertainties
c     sigs(*,3) = correlations (lat-lon,lat-rad,lon-rad)
c     baseln    = baseline length
c     dsig      = baseline length uncertainty (unscaled)
c     iflg      = flag if site is fixed or free in adjustment (logical)
c     siglcl(*,1) = neu baseline
c     siglcl(*,2) = neu uncertainties
c     siglcl(*,3) = neu correlations (n-e,n-u,e-u)
c
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'  
      include 'parameters.h'
      
      logical iflg

      integer istat,jstat,ix1,ix2,itemp,ind,kstat
     .      , ico,k1,k2,ilk,indk,ilj,indj,indxs,indxa,ico1,ico2
     .      , ilive1(6),ilive2(3),il11,il12,il21,il22,i,j,k

      real*8 sigs(3,3),g(3,6),h(1,3),siggeo(21),sigxyz(6)
     .     , temp1(3),temp2(6),sigcoo(2,3)
     .     , gg(3,3),siguvw(6),siglcl(3,3),uvw(3),xyz(3)
     .     , baseln,dsig,x1,y1,z1,x2,y2,z2,radius,sigd
     .     , sinlat,coslat,sinlon,coslon,alat,along,aradius

c.......................................................................
c
      ix1=3*(istat-1)
      ix2=3*(jstat-1)
c
      call zero2d(3,6,3,6,g)
      call zero1d(1,21,siggeo)

      call sphxyz( postvl(ix1+1),postvl(ix1+2),postvl(ix1+3)
     .           , x1,y1,z1,1)
      call sphxyz( postvl(ix2+1),postvl(ix2+2),postvl(ix2+3)
     .           , x2,y2,z2,1)
c
c Compute baseline vector
      xyz(1)=x2-x1
      xyz(2)=y2-y1
      xyz(3)=z2-z1
      baseln=dsqrt((xyz(1))*(xyz(1))+(xyz(2))*(xyz(2))+
     $ (xyz(3))*(xyz(3)))      
c      write(*,'(a,2(/,2f20.13,f20.7,3f20.4))'),'COVBLN ictyp postvl x '
c     .,postvl(ix1+1),postvl(ix1+2),postvl(ix1+3),x1*1.d3,y1*1.d3,z1*1.d3    
c     .,postvl(ix2+1),postvl(ix2+2),postvl(ix2+3),x2*1.d3,y2*1.d3,z2*1.d3 
c      write(*,'(a,3f20.4)') 'xyz ',(xyz(i)*1.d3,i=1,3)
c      stop
c Convert baseline units to meters
      call filsig(sigs,sigxyz,xyz,1)
c
c Check for free coordinates
      iflg=.true.
      do 10 i=1,3
         if(free(ix1+i).eq.1) go to 11
         if(free(ix2+i).eq.1) go to 11
   10 continue
      iflg=.false.
      baseln=baseln*1.0d3
      return
   11 continue
c
c Compute jacobian from geocentric to cartesian coordinates
c   form jacobian for both stations
c
      call zero2d(2,3,2,3,sigcoo)
               
      itemp = 0
      do 50 i=1,2
         ind=3*(i-1)
         if(i.eq.1) itemp=ix1
         if(i.eq.2) itemp=ix2
         kstat=itemp/3+1
         sinlat=dsin(postvl(itemp+1))
         coslat=dcos(postvl(itemp+1))
         sinlon=dsin(postvl(itemp+2))
         coslon=dcos(postvl(itemp+2))
         radius=postvl(itemp+3)

         call live(free,nsite,kstat,3,ilive2)

c        fill jacobian matrix
c        mode = 1: d(x,y,z)/d(lat,lon,rad)
         call fillj(g,gg,sinlat,coslat,sinlon,coslon,radius,ind,
     .                   ilive2,1)

c        mode = 2: d(u,v,w)/d(x,y,z)
         if (i.eq.1)
     .   call fillj(g,gg,sinlat,coslat,sinlon,coslon,radius,ind,
     .                   ilive2,2)

c        no factor 3, no rms scaling
c
         do 20 ico=1,3
            if (ilive2(ico).eq.0) goto 20
            sigcoo(i,ico)=sigma(ilive2(ico))
 20      continue
c
         if (i.eq.1) call copy1i(1,3,0,ilive2,ilive1)
         if (i.eq.2) call copy1i(1,3,3,ilive2,ilive1)
c
 50   continue
c
c     construct geocentric coordinates covariance matrix
c     units in meters
c
c variance and covariances - station 1 & station 2
      do 80 i=1,2
         ind=3*(i-1)
         do 30 k=1,3
            k1=k+ind
            k2=k1*(k1+1)/2
            ilk=ilive1(k1)
            siggeo(k2)=sigcoo(i,k)**2
            if(ilk.eq.0.or.k.eq.3) go to 30
            indk=ind+k
            do 40 j=k+1,3
               ilj=ilive1(j+ind)
               if(ilj.eq.0) go to 40
               indj=ind+j
               indxs=(indj*indj-indj)/2+indk
               indxa=ilj*(ilj-1)/2+ilk
               siggeo(indxs)=a(indxa)*sigcoo(i,k)*sigcoo(i,j)
               siggeo(indxs)=siggeo(indxs)/(sigma(ilk)*sigma(ilj))
 40         continue
 30      continue
 80   continue
c
c cross-covariances : station 1 & station 2
      do 90 ico1=1,3
         il11=ilive1(ico1)
         il12=ilive1(ico1+3)
         if (il11.ne.0.and.il12.ne.0)
     .    call dosig(siggeo,ico1,ico1,il12,il11,sigma,sigcoo)
         do 100 ico2=1,3
            if(ico1.eq.ico2) goto 100
            il21=ilive1(ico2)
            il22=ilive1(ico2+3)
            if (il11.ne.0.and.il22.ne.0)
     .    call dosig(siggeo,ico2,ico1,il22,il11,sigma,sigcoo)
            if (il12.ne.0.and.il21.ne.0)
     .    call dosig(siggeo,ico1,ico2,il12,il21,sigma,sigcoo)
 100     continue
 90   continue
c
c correct sign for first site partials
      do 55 i=1,3
      do 55 j=1,3
 55   g(i,j)=-g(i,j)
c
c     propagate cartesian coordinates covariance matrix
c
      call gpgt(g,siggeo,3,6,sigxyz,temp1,temp2)
c
      call filsig(sigs,sigxyz,xyz,2)
c
c  propagate to local coordinates

      alat = postvl(ix1+1)
      along = postvl(ix1+2)
      aradius = postvl(ix1+3)
      call dxyz_to_neu (xyz, sigxyz,
     .                  alat,along,aradius, 
     .                  1, uvw, siguvw)

c      call gpgt(gg,sigxyz,3,3,siguvw,temp1,temp2)
c
c  compute local coordinates (in case needed)
c      call mmply(gg,xyz,uvw,3,3,1)
c
      call filsig(siglcl,siguvw,uvw,3)
c
c     compute jacobian from cartesian coordinates to distance
      h(1,1)=xyz(1)/baseln
      h(1,2)=xyz(2)/baseln
      h(1,3)=xyz(3)/baseln
c     propagate distance variance
      call gpgt(h,sigxyz,1,3,sigd,temp1,temp2)
c     change units to meters
      baseln=baseln*1000.d0
      dsig=1000.d0*dsqrt(sigd)
c
      return
      end
c--------------------------------------------------------------------
      subroutine filsig(sigs,sig,xyz,mode)
c
      implicit none
                 
      integer mode

      real*8 sigs(3,3),sig(6),xyz(3),s1,s2,s3

      if (mode.eq.2) goto 10
c
c convert to meters
      sigs(1,1)=xyz(1)*1000.d0
      sigs(2,1)=xyz(2)*1000.d0
      sigs(3,1)=xyz(3)*1000.d0
      if (mode.eq.1) goto 20
c
c  variances
 10   s1=dsqrt(sig(1))
      s2=dsqrt(sig(3))
      s3=dsqrt(sig(6))
c  change units to meters
      sigs(1,2)=s1*1000.d0
      sigs(2,2)=s2*1000.d0
      sigs(3,2)=s3*1000.d0
c  covariances
      sigs(1,3)=sig(2)/(s1*s2)
      sigs(2,3)=sig(4)/(s1*s3)
      sigs(3,3)=sig(5)/(s2*s3)

 20   continue
      return
      end
c--------------------------------------------------------------------
      subroutine dosig(siggeo,ico1,ico2,il1,il2,sigma,sigcoo)
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
                 
      integer ico1,ico2,il1,il2,indxa,ind
      real*8 sigma(maxprm),siggeo(21),sigcoo(2,3)

       indxa=il1*(il1-1)/2+il2
       ind=(ico1+3)*(ico1+2)/2+ico2
       siggeo(ind)=a(indxa)*sigcoo(1,ico2)*sigcoo(2,ico1)
       siggeo(ind)=siggeo(ind)/(sigma(il2)*sigma(il1))

      return
      end

