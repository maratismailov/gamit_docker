Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.
                            
      Subroutine COVST1( covkep )
c
c   propagate state-vectors to keplerian elements including covariances
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'     
      include 'parameters.h'

      integer ind2,indx1,indx2,igood
     .      , isat,indxs1,indxs2,igrav2,jj1,jj2,jj3,ind1,ielem,i,j,k

      real*8 x(6),e(6),covx(21),cove(21),dx(6,6),temp1(6),temp2(6)
     .     , covkep(maxcov),xi,yi,zi,xdoti,ydoti,zdoti,gm

c     old code--no longer used:
c     dimension be(6),dxt(6,6),workf(6,6),works(21),check(6,6)
c     real*8 rcond    

      logical debug/.false./
c
c  find first live orbital element in entire parameter list(indx1)
c     and in live parameter list (indx2)
c
      indx1=0
      indx2=0
      do 10 i=1,ntpart
      if(islot1(i).gt.500.and.islot1(i).le.1100) go to 11
      if (free(i).eq.1) indx2=indx2+1
      indx1=indx1+1
   10 continue
   11 continue

c  loop over all satellites
      igood=0
      do 100 isat=1,nsat
       indxs1=indx1+(norb)*(isat-1)
c  skip unadjusted satellite (assuming 6 orbit elements at same status)
       if(free(indxs1+1).eq.0) go to 100
c  count number of live non-gravitational parameters
      igrav2=0          
c** rwk 190519: The upper index in the loop below is designed to
c**      go beyond the live SRP parameters for this SV without leaking
c**      into the SRP parameters of the next SV.  This works when 
c**      maxorb is 18 but not when it's 22.  The logic is flawed.
c**      For now I'm going to kluge a fix by reducing the index by 1.
c      do 25 i=indxs1+7,indxs1+maxorb  
       do 25 i=indxs1+7,indxs1+maxorb -1 

      if(islot1(i).le.1100.or.islot1(i).gt.2400) go to 25
      if(free(i).eq.0) go to 25
      igrav2=igrav2+1
   25 continue
       igood=igood+1
       indxs2=indx2+(6+igrav2)*(igood-1)
c
c  form keplerian covariance matrix (diagonal)
      do 200 j=1,6
c  save state-vector for this satellite
      x(j)=preval(indxs1+j)+adjust(indxs1+j)
  200 continue
c  extract state-vector covariances from covariance matrix
      do 300 j=1,6
      jj1=(j*j-j)/2
      jj2=j+indxs2
      jj3=(jj2*jj2-jj2)/2
      do 300 k=1,j
      ind1=k+jj1
      ind2=k+indxs2+jj3
      covx(ind1)=a(ind2)
      if(debug) print *
     .,'COVST1 isat igrav2 indx1 indx2 indxs1 indxs2 ind1 ind2 covx(1) '
     .       , isat,igrav2,indx1,indx2,indxs1,indxs2,ind2,ind1,covx(1) 
  300 continue
c
c  compute Jacobian of e wrt x
      xi=x(1)
      yi=x(2)
      zi=x(3)
      xdoti=x(4)
      ydoti=x(5)
      zdoti=x(6)
      gm=398603.46D0
      call elem(xi,yi,zi,xdoti,ydoti,zdoti,e,3,dx,gm,0)
c  convert state vector to keplerian elements
c      call xyzkep(x,e)
c
c  compute jacobian for transformation of e to x
c      t=0.d0
c      nv=6
c      call kepxyz(t,e,nv,x,dx)
c
c  compute inverse transformation (find inverse of dx)
c   transpose dx
c       do 250 i=1,6
c       do 250 j=1,6
c       dxt(i,j)=dx(j,i)
c       check(i,j)=dx(i,j)
c  250 continue
c  form symmetric matrix dxt*dx
c       call mmply(dxt,dx,workf,6,6,6)
c  transfer to symmetric storage mode
c       do 260 i=1,6
c       ii=(i*i-i)/2
c       do 260 j=1,i
c       ind=ii+j
c       works(ind)=workf(i,j)
c  260  continue
c
c       call inver2(works,be,1,6,rcond,ier)
c  transfer back to full storage mode
c       do 261 i=1,6
c       do 261 j=1,6
c       if(i-j) 3,4,4
c3      ind=i+(j*j-j)/2
c       go to 5
c4      ind=j+(i*i-i)/2
c5      continue
c       workf(i,j)=works(ind)
c  261  continue
c dx now contains its inverse
c       call mmply(workf,dxt,dx,6,6,6)
c
c   convert to keplerian covariances
      call gpgt(dx,covx,6,6,cove,temp1,temp2)
      if(debug) print *,'COVST1 covx(1-42) ',(covx(i),i=1,42) 
      if(debug) print *,'COVST1 cove(1-42) ',(cove(i),i=1,42)

c
c backup keplerian covariances
      ielem=21*(igood-1)
      call copy1d(1,21,ielem,cove,covkep)
c
  100 continue
c
      return
      end
