Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      Subroutine ABNEW(nb12,ianew,x,nded)
c
c     calculate ANEW, BNEW (see description of SOLVE2)
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h' 
      include 'parameters.h'
c
      integer nb12,ianew,ijnew,indx1,indx2,ipos,nded,idi,idj
     .      , ii,jj,kk,ll,i,j,k
      real*8 x(maxprm),dtemp,dtemp0,sm
c
c statement function to compute positions in a symmetric matrix
c  stored in diagonal form
      ipos(i,j)=min0(i,j) + max0(i,j)*(max0(i,j)-1)/2
c
C     ANEW  = A22 - A21 INV(A11) A12   (A21  =  TRANSPOSE A12)
c
      sm = 1.d-13
      idi = 0
      ijnew = ianew
      if (nlive.eq.0) go to 601
c
      do 400 ii = 1,ntpart
         if(free(ii).ne.0) go to 400
         idi = idi+1
         idj = 0
         do 410 jj = 1,ii
            if(free(jj).ne.0) go to 410
            idj = idj+1
c  the position of anew in the a array
            ijnew = ijnew+1
C
C  IJNEW is the position of ANEW  =  the position of A22
C
C  IDI and IDJ are the indices of ANEW if it were
C  Written as a matrix of just "dead" indices.
C
C    NOW DO THE INNER PRODUCTS.	(A21 inv(A11) A12)
c
            dtemp0 = 0.0d0
            if((dabs(x(ii)).lt.sm).or.(dabs(x(jj)).lt.sm)) go to 390
c
            do 300 kk = 1,nlive
               dtemp = 0.0d0
               do 200 ll = 1,nlive
                  indx1 = nb12+nlive*(idj-1)+ll
                  dtemp = dtemp+a(ipos(kk,ll))*a(indx1)
200            continue
               indx2 = nb12+nlive*(idi-1)+kk
               dtemp0 = dtemp0+dtemp*a(indx2)
300         continue
c
390         continue
c ianew is essentially a new a22
            a(ijnew) = a(ijnew)-dtemp0
c
410      continue
400   continue
C
C-----------------------------------------------------------------
C We could squeeze the Bnew calculation in above .
C ( do ( A21 *inv(a11) )	 |B1>	 as that is what DTEMP gives you
C but seeing as how b  contains inv(A11)|B1>
C we can do  A21 (inv(A11)|B1>) as easily. Which is
C faster depends on the relative numbers of live
C and fixed parameters. Let's do the latter as it's
C more readable)
C----------------------------------------------------------------------
C
C |BNEW> = (  |B2> - A21 [  inv(A11) |B1>  ] )
C
C BNEW is stored in b after b. (which is where B22 is already)
c
      if (nded.eq.0) go to 601
c
      do 600 i = 1,nded
         dtemp = 0.0d0
         do j = 1,nlive
            indx1 = nb12+nlive*(i-1)+j
            dtemp = dtemp+a(indx1)*b(j)
         enddo
         k = nlive+i
         b(k) = b(k)-dtemp
600   continue
c
601   continue
c----------------------------------------------------------------------
      return
      end
