      subroutine putctr(key,job)
c
c     force specified parameters to have the same adjustments
c
c     control keys: (key)
c       1. free all parameters (no rank deficiency)
c       2. fix some station or velocity or baseline orientation
c       3. fix network center and no rotation (both location and velocity)
c       4. inner coordinates solution
c       5. outer coordinate solution
c       6. fix network center and no rotation (velocity only)
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      dimension indx(12),coef(12)
      dimension a2(maxprm*6),cl2(maxprm),cob(200),tmp2(maxprm)
      dimension tmp3(maxprm),tmp4(maxprm)

      if (key.le.0.or.job.eq.0) goto 200
c
c     temporary weight 
c     position: 10 cm     velocity: 1 mm/year
      wfac = 1.0d2 
      wvfac = 1.0d6
      npar = gdsit*6
      print*,'npar: ',npar,icnum
c  
c     geodetic height constraint
      if (ncht.gt.0) then
c     constrain cos(lat)*cos(lon)*dx+cos(lat)*sin(lon)*dy+sin(lat)*dz=0
c     constrain cos(lat)*cos(lon)*Vx+cos(lat)*sin(lon)*Vy+sin(lat)*Vz=0
      f = 1.0d0/finv
      e2 = 2.0d0*f-f*f
      do 120 is = 1,ncht
         i2 = ihgt(is)
         i1 = jtoi(i2)
         if (i1.le.0) goto 120
         sl = dsin(slon(i1))
         cl = dcos(slon(i1))
         sf = dsin(slat(i1))
         cf = dcos(slat(i1))
         fac = 1.0d0/(1.0d0-e2*sf**2)
         coef(1) = cf*cl*fac
         coef(2) = cf*sl*fac
         coef(3) = sf*(1.0d0-e2)*fac
c        insert position weighting submatrix
         indx(1) = (i1-1)*6+1
         indx(2) = indx(1)+1
         indx(3) = indx(1)+2
         cob(1) = 0.0d0
         cl2(1) = 1.0d0/wfac
         chio = chi
         write (6,'(6f16.5)') (anorm(ik*(ik+1)/2),ik=1,60)
         call upsln1(npar,cob,a2,cl2,chi,3,indx,coef)
         print*,'cght x: ',i1,chi-chio
c        insert velocity weighting submatrix
         indx(1) = i1*6-2
         indx(2) = indx(1)+1
         indx(3) = indx(1)+2
         cob(1) = 0.0d0
         cl2(1) = 1.0d0/wvfac
         chio = chi
         write (6,'(6f16.5)') (anorm(ik*(ik+1)/2),ik=1,60)
         call upsln1(npar,cob,a2,cl2,chi,3,indx,coef)
         print*,'cght v: ',i1,chi-chio
 120  continue
      endif   
c
      i0 = icnum
      if (i0.le.0) goto 130
c
c     constraint on a specified orientation
c     in 3-D case, this constraint only limit 2 degree freedom
c     we need another constraint to remedy rotation deficiency
      do 10 ic = 1,i0
         write (6,'(6f16.5)') (anorm(ik*(ik+1)/2),ik=1,60)
         i1 = icstr(1,ic)
         i2 = icstr(2,ic)
         is1 = jtoi(i1)
         is2 = jtoi(i2)
         if (is1.le.0.or.is2.le.0) goto 10
c        icord = 1: geodetic constraint (horizontal only) (coor.)
c        icord = 2: geodetic constraint (horizontal only) (velo.)
c        icord = 3: geodetic constraint (horizontal only) (both)
c        icord = 4: geocentric constraint (coor. only)
c        icord = 5: geocentric constraint (velo. only)
c        icord = 6: geocentric constraint (both)
         icord = icstr(3,ic)
         if (icord.ge.1.and.icord.le.3) then
c           constrain dVn - k1*dVe = 0
            npa = 3
            sl = dsin(slon(is1))
            cl = dcos(slon(is1))
            sf = dsin(slat(is1))
            cf = dcos(slat(is1))
            fac1 = -(x(is2)-x(is1))*cl*sf-(y(is2)-y(is1))*sl*sf
     .             +(z(is2)-z(is1))*cf
            fac1 = fac1/(-(x(is2)-x(is1))*sl+(y(is2)-y(is1))*cl)
            indx(1) = is1*6-5
            indx(2) = indx(1)+1
            indx(3) = indx(1)+2
            indx(4) = is2*6-5
            indx(5) = indx(4)+1
            indx(6) = indx(4)+2
            coef(1) = sf*cl-fac1*sl
            coef(2) = sf*sl+fac1*cl
            coef(3) = -cf
         endif
         if (icord.ge.4.and.icord.le.6) then
            fac1 = (y(is2)-y(is1))/(x(is2)-x(is1))
            npa = 2
c           constrain dy-k1*dx = 0
c           constrain dVy-k1*dVx = 0
            indx(1) = is1*6-5
            indx(2) = indx(1)+1
            indx(3) = is2*6-5
            indx(4) = indx(3)+1
            coef(1) = fac1
            coef(2) = -1.0d0
         endif
         do i = 1,npa
            coef(i+npa) = -coef(i)
         enddo
         npa2 = 2*npa
c        insert position constraint
         if (icord.ne.2.and.icord.ne.5) then
            cob(1) = 0.0d0
            cl2(1) = 1.0d0/wfac
            chio = chi
            print*,'before upsln1 cl2: ',cl2(1)
            call upsln1(npar,cob,a2,cl2,chi,npa2,indx,coef)
            print*,'cdir x: ',is1,is2,npa2,chi-chio
         endif
         write (6,'(6f16.5)') (anorm(ik*(ik+1)/2),ik=1,60)
c        insert velocity constraint
         if (icord.ne.1.and.icord.ne.4) then
            cob(1) = 0.0d0
            cl2(1) = 1.0d0/wvfac
            chio = chi
            print*,'before upsln1 cl2: ',cl2(1)
            call upsln1(npar,cob,a2,cl2,chi,npa2,indx,coef)
            print*,'cdir v: ',is1,is2,npa2,chi-chio
         endif
         if (icord.le.3) goto 10
c        constrain dz-k2*dx = 0
c        constrain dVz-k2*dVx = 0
         fac2 = (z(is2)-z(is1))/(x(is2)-x(is1))
         indx(1) = is1*6-5
         indx(2) = is1*6-3
         indx(3) = is2*6-5
         indx(4) = is2*6-3
         coef(1) = fac2
         coef(2) = -1.0d0
         coef(3) = -fac2
         coef(4) = 1.0d0
c        insert position weighting submatrix
         if (icord.ne.5) then
            npa2 = 4
            cob(1) = 0.0d0
            cl2(1) = 1.0d0/wfac
            chio = chi
            call upsln1(npar,cob,a2,cl2,chi,npa2,indx,coef)
            print*,'cdir x: ',is1,is2,npa2,chi-chio
         endif
         if (icord.eq.4) goto 10
c        insert velocity weighting submatrix
         do i = 1,4
            indx(i) = indx(i)+3
         enddo
         npa2 = 4
         cob(1) = 0.0d0
         cl2(1) = 1.0d0/wvfac
         chio = chi
         call upsln1(npar,cob,a2,cl2,chi,npa2,indx,coef)
         print*,'cdir v: ',is1,is2,npa2,chi-chio
 10   continue  
c
 130  i0 = ifnum
      if (i0.le.0) goto 150
c
c     force coordinate correction or velocity along a specified
c     direction (geodetic frame only)
c        icord = 1: force geodetic coordinate correction
c        icord = 2: force site velocity
      do 60 ic = 1,i0
         i1 = icstr(1,ic)
         i2 = icstr(2,ic)
         icord = icstr(3,ic)
         is1 = jtoi(i1)
         is2 = jtoi(i2)
         if (is1.le.0.or.is2.le.0) goto 60
         print*,'fdir: ',is1,is2,icord
c        later consider dx = 0
c           constrain dVn - k1*dVe = 0
         npa = 3
         sl = dsin(slon(is1))
         cl = dcos(slon(is1))
         sf = dsin(slat(is1))
         cf = dcos(slat(is1))
         fac1 = -(x(is2)-x(is1))*cl*sf-(y(is2)-y(is1))*sl*sf
     .             +(z(is2)-z(is1))*cf
         fac1 = fac1/(-(x(is2)-x(is1))*sl+(y(is2)-y(is1))*cl)
         if (icord.eq.1) then
            indx(1) = is1*6-5
            indx(2) = indx(1)+1
            indx(3) = indx(1)+2
            wfc = wfac
         endif
         if (icord.eq.2) then
            indx(1) = is1*6-2
            indx(2) = indx(1)+1
            indx(3) = indx(1)+2
            wfc = wvfac
         endif
         coef(1) = sf*cl-fac1*sl
         coef(2) = sf*sl+fac1*cl
         coef(3) = -cf
c        insert weighting submatrix
         call zero1d(1,maxprm,a2)
         cob(1) = 0.0d0
         cl2(1) = 1.0d0/wvc
         chio = chi
         call upsln1(npar,cob,a2,cl2,chi,npa,indx,coef)
         print*,'fdir : ',icord,chi-chio
         if (icord.eq.1) then
            indx(1) = is2*6-5
            indx(2) = indx(1)+1
            indx(3) = indx(1)+2
         endif
         if (icord.eq.2) then
            indx(1) = is2*6-2
            indx(2) = indx(1)+1
            indx(3) = indx(1)+2
         endif
c        insert weighting submatrix
         cob(1) = 0.0d0
         cl2(1) = 1.0d0/wvc
         chio = chi
         call upsln1(npar,cob,a2,cl2,chi,npa,indx,coef)
         print*,'fdir : ',icord,chi-chio
 60   continue
c
c     force coordinate or velocity to be the same
 150  if (lkx.gt.0) then
         do 160 ib = 1,lkx
            i1 = lkxstr(1,ib)
            i2 = lkxstr(2,ib)
            is1 = jtoi(i1)
            is2 = jtoi(i2)
            if (is1.le.0.or.is2.le.0) goto 160
            indx(1) = (is1-1)*6+1
            indx(2) = indx(1)+1
            indx(3) = indx(1)+2
            indx(4) = (is2-1)*6+1
            indx(5) = indx(4)+1
            indx(6) = indx(4)+2
            coef(1) = 1.0d0
            coef(2) = 1.0d0
            coef(3) = 1.0d0
            coef(4) = -1.0d0
            coef(5) = -1.0d0
            coef(6) = -1.0d0
            cob(1) = 0.0d0
            cob(2) = 0.0d0
            cob(3) = 0.0d0
            call zero1d(1,6,cl2)
            cl2(1) = 1.0d0/wvfac
            cl2(3) = 1.0d0/wvfac
            cl2(6) = 1.0d0/wvfac
            call zero1d(1,npar*3,a2)
            a2(indx(1)) = coef(1)
            a2(indx(4)) = coef(4)
            a2(indx(2)+npar) = coef(2)
            a2(indx(5)+npar) = coef(5)
            a2(indx(3)+npar*2) = coef(3)
            a2(indx(6)+npar*2) = coef(6)
            print*,'index: ',(indx(k),k=1,6)
            chio = chi
            call seqsln(npar,3,bnorm,anorm,cob,a2,cl2,chi,gvm,tmp2,
     .                  tmp3,tmp4,1)
            print*,'tie x: ',is1,is2,chi,chi-chio
 160     continue
      endif
      if (lkv.gt.0) then
         do 170 ib = 1,lkv
            i1 = lkvstr(1,ib)
            i2 = lkvstr(2,ib)
            is1 = jtoi(i1)
            is2 = jtoi(i2)
            if (is1.le.0.or.is2.le.0) goto 170
            indx(1) = (is1-1)*6+4
            indx(2) = indx(1)+1
            indx(3) = indx(1)+2
            indx(4) = (is2-1)*6+4
            indx(5) = indx(4)+1
            indx(6) = indx(4)+2
            coef(1) = 1.0d0
            coef(2) = 1.0d0
            coef(3) = 1.0d0
            coef(4) = -1.0d0
            coef(5) = -1.0d0
            coef(6) = -1.0d0
            cob(1) = 0.0d0
            cob(2) = 0.0d0
            cob(3) = 0.0d0
            call zero1d(1,6,cl2)
            cl2(1) = 1.0d0/wvfac
            cl2(3) = 1.0d0/wvfac
            cl2(6) = 1.0d0/wvfac
            call zero1d(1,npar*3,a2)
            a2(indx(1)) = coef(1)
            a2(indx(4)) = coef(4)
            a2(indx(2)+npar) = coef(2)
            a2(indx(5)+npar) = coef(5)
            a2(indx(3)+npar*2) = coef(3)
            a2(indx(6)+npar*2) = coef(6)
            chio = chi
            print*,'index: ',(indx(k),k=1,6)
            call seqsln(npar,3,bnorm,anorm,cob,a2,cl2,chi,gvm,tmp2,
     .                  tmp3,tmp4,1)
            print*,'tie v: ',is1,is2,chi,chi-chio
 170     continue
      endif
c
      i0 = ibnum
      if (i0.le.0) goto 200
c
c     constraint on baseline length
      do 110 ib = 1,i0
         i1 = ibstr(1,ib)
         i2 = ibstr(2,ib)
         is1 = jtoi(i1)
         is2 = jtoi(i2)
         if (is1.le.0.or.is2.le.0) goto 110
         bl = baslen(is1,is2,0.0d0,3,2)
c        coefficients (xyz coordinate)
         coef(1) = -(x(is2)-x(is1))/bl
         coef(2) = -(y(is2)-y(is1))/bl
         coef(3) = -(z(is2)-z(is1))/bl
         coef(4) = -coef(1)
         coef(5) = -coef(2)
         coef(6) = -coef(3)
c        insert position weighting submatrix
         indx(1) = (is1-1)*6+1
         indx(2) = indx(1)+1
         indx(3) = indx(1)+2
         indx(4) = (is2-1)*6+1
         indx(5) = indx(4)+1
         indx(6) = indx(4)+2
         cob(1) = 0.0d0
         cl2(1) = 1.0d0/wfac
         call upsln1(npar,cob,a2,cl2,chi,6,indx,coef)
c        insert velocity weighting submatrix
         do id = 1,6
            indx(id) = indx(id)+3
         enddo
         cob(1) = 0.0d0
         cl2(1) = 1.0d0/wvfac
         call upsln1(npar,cob,a2,cl2,chi,6,indx,coef)
 110  continue     
c
 200  continue
c
      return
      end

c-------------------------------------------------------------
      subroutine upsln1(npar,cob,a2,cl2,chi,np,indx,coef)
c
c     update solution with one condition equation
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      dimension indx(12),coef(12)
      dimension a2(npar),cl2(npar),cob(200),tmp2(maxprm)
      dimension tmp3(maxprm),tmp4(maxprm)

      call zero1d(1,npar,a2)
      do i = 1,np
         a2(indx(i)) = coef(i)
      enddo 
      print*,'before seqsln cl2: ',cl2(1)
      call seqsln(npar,1,bnorm,anorm,cob,a2,cl2,chi,scale,tmp2,
     .            tmp3,tmp4,1)

      return
      end
