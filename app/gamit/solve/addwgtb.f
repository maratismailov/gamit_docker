      Subroutine ADDWGTB( indx,ndim,n,wgtmat )   

c     Add weight to the A and AN22 matrices for biases
c     (different from ADDWGT in filling the AN22 rather than ALC matrix)
c      R. King 070402
      
c    Input
c    -----
c     indx         :   index in normal equations of last element before the current one 
c     ndim         :   dimension of weight matrix to be added
c     n            :   # rows and columns of the weight matrix to be added
c     wgtmat(ndim) :   weight matrix (lower triangular) to be added 

c     a(maxnrm)    :   Normal matrix, in common /nmat/, in solve.h
c     an22(maxwm1) :   Bias part of normal matrix  in common /lcnom4/ in solve.h
c          

c    Output:  Updated a- and an22-matrices
c    ------
c           

      implicit none

      include '../includes/dimpar.h'  
      include 'solve.h'
      
      integer*4 indx,ndim,n,j,jj1,jj2,jj3,k,ind1,ind2

      real*8 wgtmat(ndim)
             
c  Add weights to A matrix (code same as ADDWGT) 
      do j=1,n
        jj1=(j*j-j)/2
        jj2=j+indx
        jj3=(jj2*jj2-jj2)/2 
c        print *,'j jj1 jj2 jj3 ',j,jj1,jj2,jj3
        do k=1,j
          ind1=k+jj1
          ind2=k+indx+jj3     
          if( ind1.gt.ndim.or.ind2.gt.maxnrm) then
            print *,'ADDWGTB k jj3 indx ndim maxnrm ind1 ind2 '
     .          ,k,jj3,indx,ndim,maxnrm,ind1,ind2
c            print *,'ADDWGTB: data  ind2 a ',ind2,a(ind2)   
             stop
          endif
          a(ind2) = a(ind2) + wgtmat(ind1)
c       print *,'ind1 ind2 covsvant(ind1) ',ind1,ind2,covsvant(ind1) 
c       print*,'ADDWGTB: indx ind1 ind2 alc ',indx,ind1,ind2,a(ind2)
         enddo
      enddo       

c  Add weights to AN22 matrix

      do j=1,n
        jj1=(j*j-j)/2
        jj2=j
        jj3=(jj2*jj2-jj2)/2 
c         print *,'j jj1 jj2 jj3 ',j,jj1,jj2,jj3
        do k=1,j
          ind1=k+jj1
          ind2=k+jj3    
          if( ind1.gt.ndim.or.ind2.gt.maxwm1) then
            print *,'ADDWGTB k jj3 indx ndim maxwm1 ind1 ind2 '
     .          ,k,jj3,indx,ndim,maxwm1,ind1,ind2
            stop
c            print *,'ADDWGT: data  ind2 a22 ',ind2,an22(ind2) 
          endif    
          an22(ind2) = an22(ind2) + wgtmat(ind1)
c          print *,'ind1 ind2 covsvant(ind1) ',ind1,ind2,covsvant(ind1) 
c          print*,'ADDWGTB: indx ind1 ind2 an22 '
c     .          ,indx,ind1,ind2,an22(ind2)
        enddo
      enddo 

      return
      end


