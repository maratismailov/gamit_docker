      subroutine cpsoln(mode)
c
c     copy velocity solution and covariance submatrix
c     logic:
c       1. due to the rank deficiency, we fix some parameters to 
c       get regular normal matrix.  This solution therefore 
c       can be transformed to other kind of solution.  We need
c       all velocity parameters including the 'fixed' parameters
c       in the 'solution' array. 
c       2. We only deal with the velocity field.  So, the site
c       coordinate solution can be discarded and we get a compressed
c       'solution' array
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer mode,isum,igd,i,ik,id,id1,ikv,j
      integer jk,j1,i1,i2,jk1,ija
      character*8 stnnam
      real*8 tmp(3),tmp1(6),cov1(6),cjac(9),tmp2(6)
c
      if (mode.eq.0) goto 100
c     
c     velocity solution
      isum = 0
c     print*,'  Site    No. Comp.  solution' 
      do 20 igd = 1,gdsit
         i = itoj(igd)
         stnnam = sname(i)
         do 40 ik = 1,3
            id1 = (i-1)*6+ik
            id = map(id1)
            ikv = id1+3
            isum = isum+1
            solutn(isum) = 0.0d0
            id1 = map(ikv)
            if (id1.gt.0) solutn(isum) = bnorm(id1)
c           print'(1x,a8,1x,i4,1x,i4,1pg12.4)',stnnam,i,ik,solutn(isum)
 40      continue
 20   continue
c
      if (iaux.gt.0) then
         ija = iaux
         do i = 1,ija
            i1 = nlive-ija+i
            print*,' AUX',i,' (second): ',bnorm(i1)
         enddo
      endif
c
      if (jaux.gt.0) then
         ija = iaux+jaux
         do i = nlive-ija+1,nlive-iaux
            i1 = i-nlive+ija
            id1 = i1-(i1-1)/6*6
            if (id1.le.3) print*,' AUX',i1,' (meter): ',bnorm(i)
            if (id1.gt.3) print*,' AUX',i1,' (radius): ',bnorm(i)
         enddo
      endif
c
      if (jeaux.gt.0) then
         ija = iaux+jaux
         do i = 1,jeaux
            i1 = nlive-jeaux-ija+i
            print*,' Episodic parameter',i,' (meter): ',bnorm(i1)
         enddo
         if (iomode(10).le.0) goto 100
         call pline(24,64,'=',1)
         write(24,'(10x,a)') 'coseismic deformation estimate'
         call pline(24,64,'.',1)
         write(24,'(4x,a4,2x,3(5x,a5,6x,a2,2x),
     .         2(6x,a3),7x,a3,5x,a4)')
     .      'time',' east','+-','north','+-',' up  ','+-','rho',
     .      'lon','lat','site'
         i1 = 0
         do 50 i = 1,iq_sit
            if (quake_use(i).le.0) goto 50
            j1 = iq_ind(nquake)
            if (i.ge.j1) then
               tq = quake_time(nquake)
            else
               jk = 1
               do j = 2,nquake
                  ik = iq_ind(j)
                  if (i.ge.jk.and.i.lt.ik) tq = quake_time(j-1)
                  jk = ik
               enddo
            endif
            i1 = i1+1
            do ik = 1,3
               jk = nlive-jeaux-ija+(i1-1)*3+ik
               tmp(ik) = bnorm(jk)
               do i2 = 1,ik
                  jk1 = jk*(jk+1)/2-ik+i2
                  cov1(ik*(ik-1)/2+i2) = anorm(jk1)
               enddo
            enddo
            j1 = quake_sit(i)
            call getjac(slon(j1),slat(j1),a3,cjac,5)
            call atwa(3,3,cjac,cov1,tmp2,tmp1,1)
            cov1(1) = dsqrt(tmp2(1))
            cov1(2) = dsqrt(tmp2(3))
            cov1(3) = dsqrt(tmp2(6))
            rho = tmp2(2)/cov1(1)/cov1(2)
            if (rho.gt.1.0d0) rho = 1.0d0
            if (rho.lt.-1.0d0) rho = -1.0d0
            call sph_ca(slat(j1),slon(j1),tmp(1),tmp(2),tmp(3),
     .         tmp1(1),tmp1(2),tmp1(3),2)
            write (24,60) tq,tmp1(1),cov1(1),tmp1(2),cov1(2),
     .         tmp1(3),cov1(3),rho,slon(j1)*rtod,
     .         slat(j1)*rtod,sname(j1)
 50      continue
      endif

 60   format (2x,f8.3,3(2x,f8.3,2x,f8.4),3(2x,f8.3),2x,a8)
 100  continue
      return
      end

