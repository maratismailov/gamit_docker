      subroutine imf_grad(val1,val2,val3,val4,val5,dlat,dlon,el,az)

c  subroutine to calculate (by least squares) a plane to pass
c  through four corners of a grid. Borrowed from code by 
c  Arthur Niell.
c
c  P. Tregoning
c  12 January 2004
c
c  Input:
c   val[1-4]   :  the height of the z200 Hpa surface at the
c                 four grid nodes
c   val5       : the interpolated z200 height at our site coordinates
c   dlat, dlon :  size of the grid in lat/long degrees.
c
c  Output:
c    el        : the angle between the normal to the gradient surface
c                and the zenith direction
c    az        : the azimuth in the horizontal plane (clockwise from north) 
c                of the normal

      implicit none

      integer*4 nobs,maxsize,nparm
      parameter (nobs=4,maxsize=4,nparm=3)
      real*4 val1,val2,val3,val4,val5,dlat,dlon,el,az,z0
      real*8 zc(5),xc(nobs),yc(nobs),apriori(nparm),A(nobs,nparm)
     .      ,o_min_c(nobs)
     .      ,W(nobs,nobs),zc0(nobs),ATW(nparm,nobs),ATWA(nparm,nparm)
     .      ,parm_adj(nparm),At(nparm,nobs),ATWAinv(nparm,nparm)
     .      ,ATWomc(nparm,1),hxz,hyz,rxy,rxy2,hz,pi,za,sigma
      integer k
      integer i,j

      pi = 3.14159265d0

c debug:
c      print*,'input parameters: ',val1,val2,val3,val4,val5,dlat,dlon

c  AEN says that 20m is a reasonable sigma value for the z200 height estimates
      sigma = 20.d0

c  form up the z-matrix. Following Niell, 
c z-coordinates of indices and reference point = heights
c zc(1)=zSW, zc(2)=zSE, zc(3)=zNW, zc(4)=zNE, zc(5)=zRP
                                                                        
      do i=1,maxsize
        zc(i) = 0.d0
      end do
                                                                       
      zc(1) = val1
      zc(2) = val2
      zc(3) = val4
      zc(4) = val3
      zc(5) = val5    

      xc(1) = 0.d0
      xc(2) = dlon
      xc(3) = 0.d0
      xc(4) = dlon

      yc(1) = 0.d0
      yc(2) = 0.d0
      yc(3) = dlat
      yc(4) = dlat

       apriori(1) = val5    
       apriori(2) = 0.d0
       apriori(3) = 0.d0

c calculate the partials matrix, A, and the initial values for zc.
       j = 1
       do i=1,nobs
         A(i,j) = 1.d0
       end do

c  xc must be the x coordinate   
       j = 2
       do i=1,nobs
         A(i,j) = -xc(i)
       end do
                                                                        
c  yc must be the y coordinate                                    
       j = 3
       do i=1,nobs
         A(i,j) = -yc(i)
       end do

c initialize zc0
       do i=1,nobs
         zc0(i) = 0.d0
       end do
                                                                        
c calculate fit values for zc
       do i=1,nobs
         do j=1,nparm
           zc0(i) = zc0(i) + A(i,j) * apriori(j)
         end do
       end do

c       print*,'zc0 ',zc0

c define the vector o_min_c that is the difference between
c   the observables, zc, and the values zc0 calculated with apriori
c   parameter values
       do i=1,nobs
           o_min_c(i) = zc(i) - zc0(i)
       end do
        
c       print*,'Apriori ',apriori 
c       print*,'o-c',o_min_c

c calculate the weighting matrix
      do j = 1,nobs
           do k = 1,nobs
                 if (j.eq.k) then
                    W(j,k) = 1.d0/(sigma*sigma)
                 else
                    W(j,k) = 0.d0
                 end if
           end do
c      print*,'W ',(w(j,i),i=1,4)
       end do

c  transpose A
      call transp(A,At,nobs,nparm)

c      print*,' A'
c      do i=1,nobs
c          print*,(A(i,j),j=1,nparm)
c      enddo
c
c      print*,' At'
c      do i=1,nparm
c          print*,(At(i,j),j=1,nobs) 
c      enddo
c
c      print*,'Wmat'
c      do i =1,nobs
c        print*,(w(i,j),j=1,nobs)
c      enddo

c  AtW
      call matmult(At,W,ATW,nparm,nobs,nobs)
c      print*,' AtW'
c      do i=1,nparm
c          print*,(Atw(i,j),j=1,nobs) 
c      enddo


c  AtWA
      call matmpy(ATW,A,ATWA,nparm,nobs,nparm)
c      print*,' AtWA'
c      do i=1,nparm
c          print*,(Atwa(i,j),j=1,nparm)
c      enddo



c invert it. Routine inver2 will return the inverted matrix into ATWA
      call invert(ATWA,atwainv,nparm,nparm)
c      print*,' AtWA inverted'
c      do i=1,nparm
c          print*,(Atwainv(i,j),j=1,nparm)
c      enddo

c     now find adjustment to parameters;
c     remember : ATW = Atr*W
c     first calculate ATW*o_min_c and call it ATWomc
      call matmpy(ATW, o_min_c, ATWomc, nparm,nobs, 1)

c     next calculate invATWA*ATWomc to get the adjustment
c     to the apriori values.  Call the output parm_adj
      call matmpy(ATWAinv, ATWomc, parm_adj, nparm, nparm, 1)
c      print*,'solution is',parm_adj

c The solution for the parameters can be found by adding the vectors
c  apriori and parm_adj.  This solution will be stored in the vector
c  that currently holds the apriori values so that iterative improvement
c  of the apriori values can be done.
       do i=1,nparm
          apriori(i) = apriori(i) + parm_adj(i)
       end do

       z0 = apriori(1)
       hxz = apriori(2)
       hyz = apriori(3)

c Convert gradients to azimuth and zenith angles (degrees).
       rxy = sqrt(hxz * hxz + hyz * hyz)
       rxy2 = rxy * rxy
       hz = 1.d0 / sqrt(1.d0 + rxy2)
       el = 180.d0 / pi * atan2(1.d0,rxy)
       az = 180.d0 / pi * atan2(hxz,hyz)
       za = 90.d0 - el

       return
       end

c -----------------------------------------------------------------------
 
        subroutine matmult(M1,M2,M3,l,m,n)
 
c  multiplies lxm matrix M1 by a mxn matrix M2 to give a lxn matrix M3
c  written by Paul Tregoning
c  9th June 1992
 
        implicit none
        integer i,j,k,l,m,n
        real*8 M1(l,m),M2(m,n),M3(l,n),temp
 
        temp = 0.d0
 
        do 10 i=1,l
          do 30 k=1,n
            do 20 j=1,m
              temp = M2(j,k) * M1(i,j) + temp
c              print*,i,k,j,M2(j,k) * M1(i,j),temp
20          continue
 
          M3(i,k) = temp
          temp = 0.d0
30       continue
10      continue
 
 
        return
        end

c -----------------------------------------------------------------------
      subroutine invert(a,b,n)
 
c  subroutine stolen from someone else
 
      implicit none
 
      integer n,i,j,k
      real*8 a(n,n),b(n,n),z
 
      do  i=1,n
        do  j=1,n
          b(i,j) = 0.0
        enddo
        b(i,i) = 1.0
      enddo
 
      do k = 1,n
        do i = 1,n
               if(i.ne.k)then
                 z = a(i,k)/a(k,k)
            do j=1,n
              a(i,j) = a(i,j)-a(k,j)*z
              b(i,j) = b(I,j) - b(k,j)*z
            enddo
               endif
             enddo
             z = a(k,k)
             do j=1,n
               a(k,j) = a(k,j)/z
               b(k,j) = b(k,j)/z
             enddo
           enddo
 
      return
      end
 
 























