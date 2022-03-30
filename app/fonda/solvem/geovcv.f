      subroutine geovcv(mode)
c
c     transfer 3-D geocentric velocity covariance matrix to
c     3-D geodetic velocity covariance matrix.  
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      real*8 cjac1(9),cjac2(9),temp(12),tcov(12)
      integer j1,i,i1,i2,id,igd,jgd,j,mode,id1,id2
c
      covf = 1.0d-10
c
      do igd = 1,nsit
         cova(igd*(igd+1)/2) = covf
      enddo
      do 20 igd = 1,gdsit
         i = itoj(igd)
c        get geocentric coordinate
         x1 = x(i)
         y1 = y(i)
         z1 = z(i)
c        get geodetic coordinate
         call geoxyz(radius,finv,tx,ty,tz,a1,a2,a3,x1,y1,z1,2,hght)
c        get Jacob matrix
         call getjac(a2,a1,a3,cjac1,5)
      do 10 jgd = 1,igd
c        special case: diagonal submatrix (symmetric)
         if (jgd.eq.igd) then
            i2 = 0
            do i1 = 1,3
               id1 = map((i-1)*6+i1+3)
               do j1 = 1,i1
                  id2 = map((i-1)*6+j1+3)
                  i2 = i2+1
                  if (id1.eq.0.or.id2.eq.0) then
                     tcov(i2) = 0.0d0
                     if (j1.eq.i1) tcov(i2) = covf
                  else
                     if (id1.gt.id2) id = id1*(id1-1)/2+id2
                     if (id1.le.id2) id = id2*(id2-1)/2+id1
                     tcov(i2) = anorm(id)
                  endif
               enddo
            enddo
            call atwa(3,3,cjac1,tcov,temp,cjac2,1)
            i2 = 0
            do i1 = igd*3-2,igd*3
               do j1 = igd*3-2,i1
                  i2 = i2+1
                  id = i1*(i1-1)/2+j1
                  cova(id) = temp(i2)
               enddo
            enddo
            goto 10
         endif
c        get geocentric coordinate
         j = itoj(jgd)
         x1 = x(j)
         y1 = y(j)
         z1 = z(j)
c        get geodetic coordinate
         call geoxyz(radius,finv,tx,ty,tz,a1,a2,a3,x1,y1,z1,2,hght)
c        get Jacob matrix
         call getjac(a2,a1,a3,temp,5)
c        get transform
         call trans(3,3,temp,cjac2)
c        copy submatrix
            i2 = 0
            do i1 = 1,3
               id1 = map((i-1)*6+i1+3)
               do j1 = 1,3
                  id2 = map((j-1)*6+j1+3)
                  i2 = i2+1
                  if (id1.eq.0.or.id2.eq.0) then
                     tcov(i2) = 0.0d0
                     if (id1.eq.id2) tcov(i2) = covf
                  else
                     if (id1.gt.id2) id = id1*(id1-1)/2+id2
                     if (id1.le.id2) id = id2*(id2-1)/2+id1
                     tcov(i2) = anorm(id)
                  endif
               enddo
            enddo
c        get submatrix: cjac1*tcov*cjac2
         call axb(3,3,3,cjac1,tcov,temp,1,0)
         call axb(3,3,3,temp,cjac2,tcov,1,0)
c        copy submatrix (V unit = m/year)
            i2 = 0
            do i1 = igd*3-2,igd*3
               do j1 = jgd*3-2,jgd*3
                  i2 = i2+1
                  id = i1*(i1-1)/2+j1
                  cova(id) = tcov(i2)
               enddo
            enddo
 10      continue
 20   continue
c
 100  continue
      return
      end

