      subroutine delaunay (ndim,itri,ntri,idst)
      
c     Find delaunay triangles
c     
c     Based on ACORD, by
c     C.F. Watson
c     Computers and Geosciences, vol. 8, p 97-101, 1982
c     
c     Typed and tested by Kurt Feigl MIT
c     
c     Read data points and form all 3-tuples such that no other point lies
c     within that 3-tuple's circumcircle.

      include 'solvem.fti'

c     INPUT:
c     array dimensions
      integer ndim,ig
c     OUTPUT:
c     index of stations in triangles
      integer*4 itri(3*ndim)
c     number of triangles
      integer ntri,idst(ndim)
       
c     maximum number of points
      integer maxpnt
      parameter (maxpnt = 3*maxsit)

c     maximum number of triangles
      integer maxcon
      parameter (maxcon = maxsit)

c     Mercator X,Y of data points in rubber units
      real*8 pnt(maxpnt+3,2)                       
c     circumcircle center and radius squared for each 3-tuple.
      real*8 tetr(2*maxpnt+1,3)
c     data point indices in input order of each 3-tuple
      integer itetr(2*maxpnt+1,3)
c     last-in, first-out stack of indices of vacant 3-tuples.
      integer istack(2*maxpnt+1)
c     temporary list of edges of deleted 3-tuples
      integer ktetr(maxcon,2)
c     pointer to istack
      integer id
c     pointer to tetr and itetr
      integer jt      

c     other variables
      real*8 xpnt(3,3), det(2,3), xmin,ymin,xmax,ymax
     .,    datax,datay,dsq,dd
     .,    rat,rat0,ron,ron0
     .,    b0,rscale,scalec



      real*8 dlonmin,dlatmin,dlonmax,dlatmax 

      integer itemp(3,2)
      integer isp,i,j,k,l,n,km,kmt,nuc,kt,k1,l1,l2,jz,i2

      data itemp /1,1,2,2,3,3/
      data xpnt  /-1.,5.,-1.,-1.,-1.,5.,2.,2.,18./
      data isp   /1/
      data id    /2/
      data xmin,ymin,xmax,ymax /1.e30,1.e30,-1.e30,-1.e30/       
      data dlonmin,dlatmin,dlonmax,dlatmax /1.e30,1.e30,-1.e30,-1.e30/       

c     total number of points is 3 more than number of stations
      n = ndim + 3

      if (gdsit .gt. maxpnt) then
         print *,' DELAUNAY: gdsit, maxpnt ',gdsit,maxpnt
         stop ' in DELAUNAY'
      endif

c     scan the data for extrema
      do 10 ig = 1,ndim
         i = idst(ig)
         dlatmax = max (dlatmax,slat(i))
         dlatmin = min (dlatmin,slat(i))
         dlonmax = max (dlonmax,slon(i))
         dlonmin = min (dlonmin,slon(i))
  10  continue              
         
c     map into Transverse Mercator coordinates
c     origin is at (rat0,ron0)
      rat0 = (dlatmin+dlatmax)/2.0d0
      ron0 = (dlonmin+dlonmax)/2.0d0
      rscale = 6378.2064d0
      do 15 i = 4,n                                  
         ig = idst(i-3)
         rat = slat(ig)
         ron = slon(ig)
         b0 = cos(rat) * sin(ron-ron0)
         pnt(i,1) = 0.50d0 * rscale * log ((1.0d0+b0)/(1.0d0-b0))
         pnt(i,2) = rscale*(datan(dtan(rat)/dcos(ron-ron0))-rat0)

c        find extrema 
         xmax = max(xmax,pnt(i,1))
         xmin = min(xmin,pnt(i,1))
         ymax = max(ymax,pnt(i,2))
         ymin = min(ymin,pnt(i,2))
  15  continue


c     figure out limits
      datax = xmax - xmin
      datay = ymax - ymin
      
      scalec = max(datax,datay)

c     scale the X,Y coordinates
c     Moving the origin will ensure that the 3 exterior points in xpnt
c     are indeed exterior to the data points.

      do i = 4,n
         pnt(i,1) = (pnt(i,1) - xmin)/scalec
         pnt(i,2) = (pnt(i,2) - ymin)/scalec
      enddo

c     load 3 points exterior to the net       
      do 3 i = 1,3
         itetr(1,i) = i
         tetr(1,i) = xpnt(i,3)
         do 2 j = 1,2
            pnt(i,j) = xpnt(i,j)
   2     continue
   3  continue

      do 6 i = 2, 2*maxpnt+1
         istack(i) = i
   6  continue                                         

      do 50 nuc = 4,n
         km = 0
c        Loop through established 3-tuples
         do 30 jt =1,isp
c           Test if new data point is within the jt circumcircle
            dsq = tetr(jt,3) - (pnt(nuc,1) - tetr(jt,1))**2
            if (dsq .lt. 0.0) goto 30
            dsq = dsq - (pnt(nuc,2) - tetr(jt,2))**2
            if (dsq .lt. 0.0) goto 30
c           delete this 3-tuple, but save its edges
            id = id - 1
            istack(id) = jt
c           Add edges to ktetr but delete if already listed
            do 28 i = 1,3
               l1 = itemp(i,1)
               l2 = itemp(i,2)
               if (km .gt. 0) then
                  kmt = km
                  do 24 j = 1,kmt
                     if (itetr(jt,l1) .ne. ktetr(j,1)) goto 24
                     if (itetr(jt,l2) .ne. ktetr(j,2)) goto 24
                     km = km - 1
                     if (j .gt. km) goto 28
                     do 22 k = j,km
                        k1 = k + 1
                        do 20 l = 1,2
                           ktetr(k,l) = ktetr(k1,l)
  20                    continue
  22                 continue
                     goto 28
  24              continue           
               endif
  26           km = km + 1
               ktetr(km,1) = itetr(jt,l1)
               ktetr(km,2) = itetr(jt,l2)
  28        continue
  30     continue

c        form new 3-tuples
         do 48 i = 1,km
            kt = istack(id)
            id = id + 1
c           calculate the circumcircle and radius squared 
c           of point ktetr(i,*) and place in tetr(kt,*)
            do 44 jz = 1,2
               i2 = ktetr(i,jz)
               det(jz,1) = pnt(i2,1) - pnt(nuc,1)
               det(jz,2) = pnt(i2,2) - pnt(nuc,2)
               det(jz,3) = det(jz,1) * (pnt(i2,1) + pnt(nuc,1))/2.0
     .                   + det(jz,2) * (pnt(i2,2) + pnt(nuc,2))/2.0
  44        continue
            dd = det(1,1) * det(2,2) - det(1,2) * det(2,1)
            tetr(kt,1) = (det(1,3)*det(2,2) - det(2,3)*det(1,2))/dd
            tetr(kt,2) = (det(1,1)*det(2,3) - det(2,1)*det(1,3))/dd
            tetr(kt,3) = (pnt(nuc,1) - tetr(kt,1))**2
     .                 + (pnt(nuc,2) - tetr(kt,2))**2
            itetr(kt,1) = ktetr(i,1)
            itetr(kt,2) = ktetr(i,2)
            itetr(kt,3) = nuc
  48     continue
         isp = isp + 2
  50  continue

      ntri = 0
c     loop over all triangles
      do 90 jt = 1,isp
         if (itetr(jt,1) .ge. 4 .and. tetr(jt,3) .le. 1.) then
            id = ntri*3
            ntri = ntri + 1
            itri(id+1) = jtoi(idst(itetr(jt,1) - 3))
            itri(id+2) = jtoi(idst(itetr(jt,2) - 3))
            itri(id+3) = jtoi(idst(itetr(jt,3) - 3))
         endif
   90 continue

      return

      end





