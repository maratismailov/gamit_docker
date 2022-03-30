      subroutine calstr(ifile,idim,mode)
 
*     root: getstr by Kurt Feigl March, 1991
*     modified to fit FONDA frame by Dong, 1991.07.30
*
*     modified by Gilbert Ferhat May 11, 1994
*     to have correct triangle name with 8-char
*
*     This program will calculate strain rates from a set of velocities.
*     mode = 1:
*     Triplets of stations are chosen to form Delaunay triangles.
*     In each triangle, the horizontal strain rates are uniquely
*     determined and printed out in terms of their eigenvelues.
*     mode = 2:
*     get averaged horizontal strain rate over whole network
*     mode = 3:
*     get averaged horizontal strain rate over specified sub-networks
*
*     ifile: index of strain rate file
*
 
      include 'solvem.fti'
      integer mnsit
      parameter (mnsit=maxsit)
 
*     VARIABLES
 
*     list of indices for triangle vertices
      integer itri (21*mnsit)
      integer ntri, kk ,k2
c     latitude, longitude in degrees
      real*8 dlat3(mnsit),dlon3(mnsit)
c     covariances in (mm/yr)(e1,n1,e2,n2,e3,n3)
      real*8 vcov(maxprm*6),xv(6),cx(21)
c     principal components and their uncertainties
      real*8 eps1,eps1sig,eps2,eps2sig
c     clockwise spin rate and sigma in rad/yr
      real*8 spin,spinsig
c     azimuth of more compressive eigenvector (degrees) and uncertainty
      real*8 theta,thetasig
      real*8 gamma1,gamma2,g1sig,g2sig
      real*8 ch2
 
      real*8 dlatc,dlonc
      character*8 name1,name2,name3,namea,nameb,namec
 
      integer jt,j,idim,mode,ifile,isit
      equivalence (vcov,aprm)
 
*   rcpar       - Function to read runstring
*   iel, indx   - Pointers to positions in strings
*   ns          - Number of sites found
*   nstn        - number of stations in list
 
 
      integer*4 nstn
c     list of pointers to stations
 
      character*96 pname
 
      integer idst(maxsit),isubt(3*mnsit)
 
      integer i,k,i1,i2,i3,k1,j1,j2,ipr,istn,ic,inet
      character*128 line
      integer count_arg,lift_arg,match_name,lbad
      character*8 sitnam
 
*     true if disaster
      logical lerror
 
c     read mode
      read (18,*) mode
      if (mode.le.0) goto 300
      if (mode.le.2) then
         read (18,*) istn 
         if (istn.le.2) goto 300
         isit = 0
         lbad = 0
         do 20 i1 = 1,istn
            read (18,'(a)',end=30) line
            ic = count_arg(line)
            if (ic.le.0) goto 20
            do 10 j = 1,ic
               i2 = lift_arg(line,sitnam,j)
               if (i2.le.0) goto 10
               i3 = match_name(nsit,i2,sname,sitnam)
               if (i3.le.0) then
                  print*,' CALSTR name mismatch:',sitnam
                  lbad = lbad+1
                  goto 10
               endif
               k = map(i3*6)
               if (k.le.0) goto 10
               isit = isit+1
               idst(isit) = i3
               if (mode.eq.2) itri(isit) = jtoi(i3)
               if (isit.ge.istn-lbad) goto 30
 10         continue
 20      continue
      endif
      if (mode.eq.3) then
         read (18,*) inet 
         if (inet.le.0) goto 300
         isit = 0
         j1 = 0
         do 130 k = 1,inet
         read (18,*) istn 
         if (istn.le.2) goto 300
         k1 = 0
         lbad = 0
         do 120 i1 = 1,istn
            read (18,'(a)',end=30) line
            ic = count_arg(line)
            if (ic.le.0) goto 120
            do 110 j = 1,ic
               i2 = lift_arg(line,sitnam,j)
               if (i2.le.0) goto 110
               i3 = match_name(nsit,i2,sname,sitnam)
               if (i3.le.0) then
                  print*,' CALSTR name mismatch:',sitnam
                  lbad = lbad+1
                  goto 110
               endif
               j2 = map(i3*6)
               if (j2.le.0) then
                  lbad = lbad+1
                  goto 110
               endif
               isit = isit+1
               k1 = k1+1
               itri(isit) = jtoi(i3)
               if (k1.ge.istn-lbad) goto 125
 110        continue
 120     continue
 125     if (k1.eq.0) goto 130
         j1 = j1+1
         isubt(j1) = k1
 130     continue
         ntri = j1
      endif
 
 30   close (18)
      if (isit.le.2) goto 300
      inet = isit
 
      if (mode.eq.1) then
c        get delaunay triangles
         write (*,*) ' Calculating strain rate over Delaunay triangles'
         call delaunay (isit,itri,ntri,idst)
         do i = 1,ntri
            isubt(i) = 3
         enddo
         write (ifile,*) 
     .      '* Strain rates estimated in Delaunay triangles'
      endif
      if (mode.eq.2) then
         write(*,*) ' Calculating averaged strain rate over the network'
         isubt(1) = isit 
         ntri = 1
         write (ifile,*) '* Strain rates estimated in whole network'
         pname(1:7) = 'NETWORK'
         ipr = 7
      endif
 
      if (mode.eq.3) then
         write(*,*) ' Calculating averaged strain rate over subnetwork'
         write (ifile,*) '* Strain rates estimated in subnetwork'
         pname(1:7) = 'SUBNET_'
         ipr = 9
      endif
 
      write (ifile,17) 'LAT','LON','EPS1','EPS1SIG','EPS2'
     .,'EPS2SIG','CW-SPIN','SPINSIG','THETA','THETASIG'
     .,'TRIANGLE'
      write (ifile,17) '(deg)','(deg)','(1/yr)','(1/yr)','(1/yr)'
     .,'(1/yr)','(1/yr)','(1/yr)','(deg)','(deg)','corners'
  17  format (12(a9,1x))
      IF (IOMODE(14) .EQ. 1) THEN
         WRITE (27,18) 'LAT','LON','GAMMA1','GAM1SIG','GAMMA2',
     .   'GAM2SIG','dVe/de','dVe/dn','dVn/de','dVn/dn','TRIANGLE'
         write (27,18) '(deg)','(deg)','(1e-6/yr)','(1e-6/yr)',
     .   '(1e-6/yr)','(1e-6/yr)','(1e-6/yr)','(1e-6/yr)',
     .   '(1e-6/yr)','(1e-6/yr)','corners'
      ENDIF
  18  format (11(a9,1x))
 
c     loop over all triangles
      ic = 0
      do 90 jt = 1,ntri
         nstn = isubt(jt)
         do 60 i = 1,nstn
            idst(i) = itri(ic+i)
            isit = itoj(idst(i))
            dlat3(i) = slat(isit)*rtod
            dlon3(i) = slon(isit)*rtod
            j1 = idst(i)
c              covariances in (mm/yr)**2
               do k = 1,i
                  j2 = idst(k)
                  k1 = j2*3-2
                  i1 = (3*j1-2)*(j1-1)*3/2+k1
                  if(j2.gt.j1) i1 = (3*j2-2)*(j2-1)*3/2+j1*3-2
                  i2 = i1+1
                  if(j2.gt.j1) i2 = (3*j2-2)*(3*j2-1)/2+j1*3-2
                  i3 = (2*i-1)*(i-1)+k*2-1
                  vcov(i3) = cova(i1)*1.d6
                  if(k.lt.i) vcov(i3+1) = cova(i2)*1.d6
                  i1 = (3*j1-2)*(3*j1-1)/2+k1
                  if(j2.gt.j1) i1 = (3*j2-2)*(3*j2-3)/2+j1*3-1
                  i2 = i1+1
                  if(j2.gt.j1) i2 = (3*j2-2)*(3*j2-1)/2+j1*3-1
                  i3 = i*(2*i-1)+2*k-1
                  vcov(i3) = cova(i1)*1.d6
                  vcov(i3+1) = cova(i2)*1.d6
               enddo
           
c              name of polygon
               i2 = (i-1)*9
               if (mode.eq.1) 
     .            write (pname(i2+1:i2+9),'(a8,''-'')') sname(isit)(1:8)
 
 60      continue
CKURT    Pad subnet number with zeros to make 2 digits...
         if (mode.eq.3) write(pname(8:9),'(i2.2)') jt
         ic = ic+nstn
c         write (*,*) 'Working on polygon: ',pname(1:5*nstn)
c
         call strain (idst,vcov,gamma1,gamma2,g1sig,g2sig
     .                ,xv,cx,eps1,eps1sig,eps2,eps2sig
     .                ,theta,thetasig
     .                ,spin,spinsig
     .                ,nstn,lerror,ch2)
c
         if (lerror) goto 90
 
c        get centroid
         call centroid (dlat3,dlon3,dlatc,dlonc,nstn)
 
c        print out one line:
         if (.not. lerror) then
            if (mode.ge.2) goto 50
            ipr = 9*nstn
c           alphebetic order for name
            read (pname,'(a8,1x,a8,1x,a8)') name1,name2,name3
            if (llt(name1,name2)) then
               namea = name1
               nameb = name2
            else
               namea = name2
               nameb = name1
            endif
            if (llt(nameb,name3)) then
               namec = name3
            else
               namec = nameb
               if(llt(namea,name3)) then
                  nameb = name3
               else
                  nameb = namea
                  namea = name3
               endif
            endif
            
            write (pname,'(a8,''-'',a8,''-'',a8)') namea,nameb,namec
c           print *,pname,namea,nameb,namec
            if (mode.eq.1) pname(ipr:ipr) = ' '
 50         write (ifile,200) dlatc,dlonc
     .                     ,eps1,eps1sig,eps2,eps2sig
     .                     ,spin,spinsig,theta,thetasig
     .                     ,pname(1:ipr)
c        Someone need this information, and someone does not like it.
c        What should I do ???
c        Write it to a different "gamma file"
            write (27,221) dlonc,dlatc,
     .      gamma1,g1sig,gamma2,g2sig,(xv(k),k=3,6),
     .      pname(1:ipr)
            do kk = 3,6
               k1 = kk*(kk-1)/2+2
               k2 = kk-2
               write (27,'(60x,4(1x,1pe9.2))') (cx(k1+k),k=1,k2)
            enddo
            write (27,223) ch2,2*nstn-6
         endif
 
   90 continue
200   format (2(1x,f9.4),6(1x,1pe9.2),2(1x,0pf9.4),1x,a)
221   format (2(1x,f9.4),8(1x,1pe9.2),1x,a)
223   format (4x,'chi-2:',1pe11.4,5x,'# of freedom:',i5)
 
 
      print *,' Number of stations  = ',inet
      print *,' Number of polygons  = ',ntri
      
 300  continue
      return
      end

