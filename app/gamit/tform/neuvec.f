      program neubas

c     Purpose:
c        Convert Cartesian xyz baselines to neu baselines, given the
c        spherical or geodetic latitude and longitude, Cartesian
c        xyz vector and covariance
c
      integer i,j
      real*8 radcon,alat,alon,hght,slat,slon,clat,clon
      real*8 d(3),cd(3,3),nd(3),cnd(3,3)
      real*8 jac(3,3),work(3,3),val(9)

      radcon = atan(1.d0)/45.d0

      write(6,'(a)') 'Input location (dec. deg.) and vector in form'
      write(6,'(a)')
     .   '  lat, long, x y z sigx sigy  sigz corxy corxz coryz'
      write(6,'(a)') '  <CNTL> Z to stop'

 10   continue
      read(5,*,end=999) alat,alon,(val(i),i=1,9)

      d(1) = val(1)
      d(2) = val(2)
      d(3) = val(3)
      cd(1,1)  = val(4)**2
      cd(2,2)  = val(5)**2
      cd(3,3)  = val(6)**2
      cd(1,2)  = val(4)*val(5)*val(7)
      cd(1,3)  = val(4)*val(6)*val(8)
      cd(2,3)  = val(5)*val(6)*val(9)
      cd(2,1)  = cd(1,2)
      cd(3,1)  = cd(1,3)
      cd(3,2)  = cd(2,3)

      slat = dsin(alat*radcon)
      slon = dsin(alon*radcon)
      clat = dcos(alat*radcon)
      clon = dcos(alon*radcon)

c     (n,e,u) = J(x,y,z)
      do j = 1, 3
      do i = 1, 3
         jac(i,j) = 0.d0
      enddo
      enddo

      jac(1,1) = -slat*clon
      jac(1,2) = -slat*slon
      jac(1,3) =  clat
      jac(2,1) = -slon
      jac(2,2) =  clon
      jac(2,3) =  0.d0
      jac(3,1) =  clat*clon
      jac(3,2) =  clat*slon
      jac(3,3) =  slat

c     convert positions and velocities to neu

      call dgemv ( 'N', 3, 3
     .           , 1.d0, jac, 3
     .           , d, 1
     .           , 0.d0, nd, 1)

      call dgemm ( 'N', 'N', 3, 3, 3
     .           , 1.d0, jac, 3
     .           , cd, 3
     .           , 0.d0, work, 3 )

      call dgemm ( 'N', 'T', 3, 3, 3
     .           , 1.d0, work, 3
     .           , jac, 3
     .           , 0.d0, cnd, 3 )

c     remove negative covariances due to round-off problems

      do i = 1, 3
         if (cnd(i,i)+1.d0 .le. 1.d0) cnd(i,i) = 0.d0
      enddo

*     output the neu vector, sigmas, and correlations

      write(6,'(/,3f15.4,6f8.4)')
     .  nd(1), nd(2), nd(3)
     ., dsqrt( cnd(1,1) ), dsqrt( cnd(2,2) ), dsqrt( cnd(3,3) )
     ., cnd(1,2)/dsqrt( cnd(1,1) )/dsqrt( cnd(2,2) )
     ., cnd(1,3)/dsqrt( cnd(1,1) )/dsqrt( cnd(3,3) )
     ., cnd(3,2)/dsqrt( cnd(3,3) )/dsqrt( cnd(2,2) )

      goto 10

999   continue
      end
