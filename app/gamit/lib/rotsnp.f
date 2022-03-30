Copyright (c) Massachusetts Institute of Technology and The University of
California at San Diego, 1994. All rights reserved.

      Subroutine rotsnp( idir,jd,t,tdtgpst,iut1pol
     .                 , frame,precmod,iut1,ipole,inut
     .                 , rot,rotdot,sidtm,xpole,ypole )
c
C     Y. Bock and R. King 1986-1994; arguments changed RWK 170412 .
C
C     Calculate the rotation matrix between earth-fixed and inertial
c
c     Input:
C       idir = 1   Earth-fixed to inertial
C             -1   Inertial to Earth-fixed
c       jd         Julian day (PEP JD, beginning at midnight, = MJD + 1)
c       t          Seconds-of-day  (GPST)
c       tdtgpst    TDT - GPST  (seconds)
c       iut1pol    Binary coded control for diurnal pole (bit 1) and UT1 (bit 2)
c       iut1, ipole, inut - unit numbers for tables 
c       frame,precmod - Inertial frame and precession model to be used in rotations

c      Output:
c        rot(3,3)     Rotation matix (rad)
c        rotdot(3,3)  Derivative of rotation matrix (rad/s)
c        sidtm        Sidereal time (rad)
c        xpole        X-pole position (rad)
c        ypole        Y-pole position (rad)

      implicit none
C
      character*5 frame, precmod

      integer*4 idir,jd,jds,iut1pol,iut1,ipole,inut,i,j,k

      real*8 t,ts,tdtutc,tdtgpst,rot,rotdot,prec,rnut,srot,sdrot
     .     , snpmat,pnsmat,snpdmt,pnsdmt,eqe,sidtm,xpole,ypole
     .     , gpstutc,taiutc,sidten,prot
     .     , era,erarot

      dimension prec(3,3),rnut(3,3),srot(3,3),sdrot(3,3),
     1          snpmat(3,3),pnsmat(3,3),snpdmt(3,3),pnsdmt(3,3),
     2          rot(3,3),rotdot(3,3),sidten(3,3),prot(3,3),
     3          erarot(3,3)

cd       write(*,'(a,1x,a5,1x,a5)')' ROTSNP frame precmod ',frame,precmod
C
c Check if old ICRF (IAU76/IAU2000) or new SOFA IAU200A or IAU2006A required 
      if ( precmod .ne. 'IAU76' .and. precmod .ne. 'IAU68' ) then        
      
         call rotsnp_sofa( idir,jd,t,tdtgpst,iut1pol
     .                   , frame,precmod,iut1,ipole,inut
     .                   , rot,rotdot,sidtm,xpole,ypole
     .                   , prec,rnut,prot,sidten,era,erarot )
c CIO based IAU SOFA earth attitude model used - Copy era to sidtm and erarot to sidten
         if ( precmod .ne. 'IAU0A' ) then
            sidtm = era
            sidten = erarot
         endif  
      else
c
C Precession and nutation transformations
         call pnrot(inut,jd,t,tdtgpst,eqe,prec,rnut,frame,precmod)

C Sidereal rotation and polar motion
cd         print *,'after pnrot srotat jd,t,tdtutc,eqe,iut1pol,sdrot,sidtm'
cd     .          , jd,t,tdtgpst,eqe,iut1pol,sdrot,sidtm
cd         stop 
c     srotat (unique in GAMIT) expects UTC
         jds = jd
         ts = t
         gpstutc = taiutc(jds) - 19.d0
         call timinc(jds,ts,-gpstutc)
c     catch possible leap second
         if( jds.ne.jd ) then
           gpstutc = taiutc(jds) - 19.d0
           jds = jd
           ts = t
           call timinc(jds,ts,-gpstutc)
         endif
         tdtutc = taiutc(jds) + 32.184d0
cd         print *,'ROTSNP jd t jds ts gpstutc tdtutc '
cd     .                  ,jd,t,jds,ts,gpstutc,tdtutc
         call srotat( jds,ts,tdtutc,eqe,iut1pol,srot,sdrot,sidtm
     .           , xpole,ypole,sidten,prot,iut1,ipole,precmod)

c Combine the matrices
         if(idir.eq.-1 ) then
           call snp(prec,rnut,srot,snpmat)
           call snp(prec,rnut,sdrot,snpdmt)

c          write(6,69) jd,t,tdtgpst,idir,sidtm
c          write(6,70) ((prec(k,j),j=1,3),k=1,3)
c          write(6,71) ((rnut(k,j),j=1,3),k=1,3)
c          write(6,72) ((srot(k,j),j=1,3),k=1,3)
c          write(6,73) ((sdrot(k,j),j=1,3),k=1,3)
c          write(6,74) ((snpmat(k,j),j=1,3),k=1,3)
c          write(6,75) ((snpdmt(k,j),j=1,3),k=1,3)
c          write(13,70) ((prec(k,j),j=1,3),k=1,3)
c          write(13,71) ((rnut(k,j),j=1,3),k=1,3)
c          write(13,72) ((srot(k,j),j=1,3),k=1,3)
c          write(13,73) ((sdrot(k,j),j=1,3),k=1,3)
c          write(13,74) ((snpmat(k,j),j=1,3),k=1,3)
c          write(13,75) ((snpdmt(k,j),j=1,3),k=1,3)
c  69      format(/,' In ROTSNP, JD, T, TDTGPST, idir, sidtm '
c     .          ,i8,f13.6,f10.6,i3,d22.14)
c  70      format(/,' Precession matrix',/,3(1X,3D22.14,/))
c  71      format(/,' Nutation matrix',/,3(1X,3D22.14,/))
c  72      format(/,' S-matrix',/,3(1X,3D22.14,/))
c  73      format(/,' SDOT-matrix',/,3(1X,3D22.14,/))
c  74      format(/,' SNP-matrix',/,3(1X,3D22.14,/))
c  75      format(/,' SNPDOT-matrix',/,3(1X,3D22.14,/))
           do 100 j=1,3
           do 100 i=1,3
           rot(I,J)=snpmat(I,J)
100        rotdot(I,J)= snpdmt(I,J)

         else
           call pns(prec,rnut,srot,pnsmat)
           call pns(prec,rnut,sdrot,pnsdmt)
c          write(6,69) jd,t,tdtgpst,idir,sidtm
c          write(6,70) ((prec(k,j),j=1,3),k=1,3)
c          write(6,71) ((rnut(k,j),j=1,3),k=1,3)
c          write(6,72) ((srot(k,j),j=1,3),k=1,3)
c          write(6,73) ((sdrot(k,j),j=1,3),k=1,3)
c          write(6,76) ((pnsmat(k,j),j=1,3),k=1,3)
c          write(6,77) ((pnsdmt(k,j),j=1,3),k=1,3)
c          write(13,70) ((prec(k,j),j=1,3),k=1,3)
c          write(13,71) ((rnut(k,j),j=1,3),k=1,3)
c          write(13,72) ((srot(k,j),j=1,3),k=1,3)
c          write(13,73) ((sdrot(k,j),j=1,3),k=1,3)
c          write(13,76) ((pnsmat(k,j),j=1,3),k=1,3)
c          write(13,77) ((pnsdmt(k,j),j=1,3),k=1,3)
c  76      format(/,' PNS-matrix',/,3(1X,3D22.14,/))
c  77      format(/,' PNSDOT-matrix',/,3(1X,3D22.14,/))

            do 200 J=1,3
            do 200 I=1,3
            rot(I,J)=pnsmat(I,J)
200         rotdot(I,J)=pnsdmt(I,J)

         end if

      end if   

      return
      end
