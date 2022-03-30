      subroutine rotsln(mode)
c
c     transfer 3-D geocentric solution and covariance matrix into
c     3-D geodetic frame
c     assuming all coordinates and velocities are estimated
c
c     if X(enu) = T X(xyz), then  C(enu) = T C(xyz) T~
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer mode,igd,i,ix,jgd,i1,j1,j
      real*8 cjac1(9),cjac2(9),temp(12),tcov(12)
c
      do 20 igd = 1,gdsit
         i = itoj(igd)
c        geodetic latitude and longitude
c        **** consider velocity correction later
         a1 = slat(i)
         a2 = slon(i)
         a3 = srad(i)
c        get Jacob matrix
         call getjac(a2,a1,a3,cjac1,5)
         do ix = 1,3
            tcov(ix) = bnorm(igd*6-6+ix)
         enddo
         call axb(3,3,1,cjac1,tcov,temp,1,0)
         do ix = 1,3
            bnorm(igd*6-6+ix) = temp(ix)
            tcov(ix) = bnorm(igd*6-3+ix)
         enddo
         call axb(3,3,1,cjac1,tcov,temp,1,0)
         do ix = 1,3
            bnorm(igd*6-3+ix) = temp(ix)
         enddo
         do 10 jgd = 1,igd
c           special case: diagonal submatrix (symmetric)
            if (jgd.eq.igd) then
c              coordinate
               i1 = igd*6-5
               j1 = i1
               call submtx(i1,j1,3,3,anorm,tcov,1,1)
               call atwa(3,3,cjac1,tcov,temp,cjac2,1)
               call submtx(i1,j1,3,3,anorm,temp,1,2)
c              vellocity
               i1 = igd*6-2
               j1 = i1
               call submtx(i1,j1,3,3,anorm,tcov,1,1)
               call atwa(3,3,cjac1,tcov,temp,cjac2,1)
               call submtx(i1,j1,3,3,anorm,temp,1,2)
c              cross correlation between coor. and velocity
               i1 = igd*6-2
               j1 = i1-3
               call submtx(i1,j1,3,3,anorm,tcov,2,1)
               call axb(3,3,3,cjac1,tcov,temp,1,0)
               call axb(3,3,3,temp,cjac1,tcov,1,0)
               call submtx(i1,j1,3,3,anorm,tcov,1,2)
               goto 10
            endif
c           off-diagonal submatrix
            j = itoj(jgd)
            a1 = slat(j)
            a2 = slon(j)
            a3 = srad(j)
c           get Jacob matrix
            call getjac(a2,a1,a3,temp,5)
c           get transform
            call trans(3,3,temp,cjac2)
c           copy submatrix (coordinate)
            i1 = igd*6-5
            j1 = jgd*6-5
            call submtx(i1,j1,3,3,anorm,tcov,2,1)
c           get submatrix: cjac1*tcov*cjac2
            call axb(3,3,3,cjac1,tcov,temp,1,0)
            call axb(3,3,3,temp,cjac2,tcov,1,0)
c           copy submatrix (X unit = meter)
            call submtx(i1,j1,3,3,anorm,tcov,2,2)
c           copy submatrix (velocity)
            i1 = igd*6-2
            j1 = jgd*6-2
            call submtx(i1,j1,3,3,anorm,tcov,2,1)
c           get submatrix: cjac1*tcov*cjac2
            call axb(3,3,3,cjac1,tcov,temp,1,0)
            call axb(3,3,3,temp,cjac2,tcov,1,0)
c           copy submatrix (V unit = meter/year)
            call submtx(i1,j1,3,3,anorm,tcov,2,2)
 10      continue
 20   continue
c
 100  continue
      return
      end
