      Subroutine ADDWGT( indx,ndim,n,wgtmat )

c     Add an nxn (a priori) weight matrix to the lower triangular normal matrix
c     Written by R. King 980916 using code extracted from existing weighting
c     routines (wsat, wstat, etc) written by Y. Bock and others
      
c    Input
c    -----
c     indx         :   index in normal equations of last element before the current one 
c     ndim         :   dimension of weight matrix to be added
c     n            :   # rows and columns of the weight matrix to be added
c     wgtmat(ndim) :   weight matrix (lower triangular) to be added 

c     a(maxnrm)    :   Normal matrix, in common /nmat/, in solve.h
c     alc(nlcom)   :   LC normal matrix (N11 part), in common /lcnom1/ in solve.h
c          

c    Output:  Updated a- and alc-matrices
c    ------
c           

      implicit none

      include '../includes/dimpar.h'  
      include 'solve.h'
      
      integer*4 indx,ndim,n,j,jj1,jj2,jj3,k,ind1,ind2

      real*8 wgtmat(ndim)
           
      do j=1,n
        jj1=(j*j-j)/2
        jj2=j+indx
        jj3=(jj2*jj2-jj2)/2 
c        print *,'j jj1 jj2 jj3 ',j,jj1,jj2,jj3
        do k=1,j
          ind1=k+jj1
          ind2=k+indx+jj3                    
          if( ind1.gt.ndim.or.ind2.gt.maxnrm) then
c            print *,'ADDWGT k jj3 indx ndim ind1 ind2 '
c     .          ,k,jj3,indx,ndim,ind1,ind2
c            print *,'ADDWGT: data  ind2 alc ',ind2,alc(ind2) 
          endif
          a(ind2) = a(ind2) + wgtmat(ind1)
          alc(ind2) = alc(ind2)+wgtmat(ind1)
c      print *,'ind1 ind2 covsvant(ind1) ',ind1,ind2,covsvant(ind1) 
c       print*,'ADDWGT: indx ind1 ind2 alc ',indx,ind1,ind2,alc(ind2)
 
         enddo
      enddo 

      return
      end

