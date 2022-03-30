      SUBROUTINE WSAT

c  Written by Yehuda Bock
c  Mod:  SV antenna offsets added and made implicit none-- rwk 980910


c     ADD SATELLITE WEIGHTS TO NORMAL MATRIX

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

c     in parameters.h: parameter values, slots, and uncertainties
c     real*8 preval(maxprm),sigma(maxprm)
c     integer*4 idms(maxprm),islot1(maxprm),free(maxprm)
c     character*20 rlabel(maxprm)

c       local
      integer*4 i,j,k,ii,jj,indx,indxs,ind,ier,icnt,nv
      real*8 cove,covx,x,e,dx,temp1,temp2,be,t
      real*8 covsvant(6),rcond
  
     

C MOD TAH 910111: Changed dimension of work1 and work2 from
C     9 to 45 to allow full storage of covariance matrix
c     rwk 950513: set with variable dimension (maxorb now = 15)
c     rwk 950819: save cartesian covariances for h-file
      DIMENSION COVE(21),COVX(maxorb_cov),X(6),E(6)
      DIMENSION DX(6,6),TEMP1(maxorb_cov),TEMP2(maxorb_cov),BE(maxorb)
c     note: temp1 and temp2 are dimensioned maxorb_cov here only because they are
c           used in the debug printout; in sb LWSAT, they are maxorb


c      print *,'WSAT: weighting satellite parameters'

c  find first orbital element
c  assumption is that orbital elements are adjacent in parameter list
c   with non-gravitational parameters interleaved as well
      INDX=0
      DO 10 I=1,NTPART
       IF(ISLOT1(I).GT.500.AND.ISLOT1(I).LE.1100) GO TO 11
       INDX=INDX+1
   10 CONTINUE
   11 CONTINUE


C  loop over all satellites

      do 200 i=1,nsat
        indxs=indx + (norb)*(i-1)
c        print *,'orbit i indx norb indxs ',i,indx,norb,indxs

        do j=1,norb
C         save the state-vector for this satellite
          IF (J.LE.6) X(J)=PREVAL(INDXS+J)
          JJ=(J*J-J)/2
          do k=1,j
            IND=K+JJ
c           set the covariance matrix of the constraints - nominally Keplerian elements
            IF (J.LE.6) COVE(IND)=0.D0
            IF(J.LE.6.AND.J.EQ.K) COVE(IND)=sat_apr(I,J)**2
c           write(*,*) sat_apr(i,j)**2
c           Non-gravitational force constraints
            IF (J.GT.6) COVX(IND)=0.0D0
            IF(J.GT.6.AND.J.EQ.K) COVX(IND)=sat_apr(I,J)**2
          enddo
        enddo

c       WRITE(6,'(/,a)') 'Covariance matrix for orbits'
c        ICNT=0
c        do II=1,6
c          do jj=1,ii
c            ICNT=ICNT+1
c            TEMP1(JJ)=COVE(ICNT)
c          enddo
c          WRITE(6,'(1X,9D12.5)') (TEMP1(IJ),IJ=1,II)
c        enddo

c       If the input constraints were Keplerian (normal case), need to convert
c       to non-diagonal Cartesian

        if( satwt_type.eq.'KEP' ) then
C         need first to get Keplerian elements of state vector
          call xyzkep(x,e)
c          write(6,'(1x,f16.9)') (E(K),K=1,6)
c         now compute the Jacobian for transformation of covariances from Keplerian to Cartesian
          T=0.D0
          NV=6
          CALL KEPXYZ(T,E,NV,X,DX)
c          WRITE(6,'(a,i3)') 'Orbit transformation -- satellite',i
c          WRITE(6,'(1X,F16.9)') (X(K),K=1,6)
c          WRITE(6,'(6(1X,D12.5))') ((DX(K,J),J=1,6),K=1,6)
C         convert to cartesian covariances
c          write(6,*) (cove(ijk),ijk=1,21)
          CALL GPGT(DX,COVE,6,6,COVX,TEMP1,TEMP2)
        else if( satwt_type.eq.'CAR' ) then
          do ii=1,21
            covx(ii) = cove(ii)
          enddo
        else
          call report_stat('FATAL','SOLVE','wsat',' '
     .       ,'Illegal input covariance type',0)
        endif

c       save the cartesian covariance for the h-file
        do ii=1,maxorb_cov
          covorbx(ii,i) = covx(ii)
        enddo

C
C       Now invert the state-vector covariance matrix for the normal equations
C
c      WRITE(6,'(/,a,i3)') 'Orbit transformation -- satellite',i
c      ICNT=0
c      DO  II=1,norb
c        DO JJ=1,II
c          ICNT=ICNT+1
c          TEMP1(JJ)=COVX(ICNT)
c        enddo
c        WRITE(6,'(1X,15D12.5)') (TEMP1(IJ),IJ=1,II)
c      enddo
         CALL NRMSC2(COVX,BE,TEMP1,6,0,0)
CD       do ii = 1, 6
CD          icnt = ii*(ii-1)/2
CD          write(6,'(1X,15D12.5)') (covx(icnt+jj),jj=1,ii)
CD       enddo
         CALL INVER2(COVX,BE,1,6,rcond,IER)
         if( ier.ne.0 )call report_stat('WARNING','SOLVE','wsat',' '
     .    , 'Bad inverse of a priori state-vector covariance matrix',0)
c        write(*,*) 'I,IER',i,ier
CD       do ii = 1, 6
CD          icnt = ii*(ii-1)/2
CD          write(6,'(1X,9D12.5)') (covx(icnt+jj),jj=1,ii)
CD       enddo
         CALL NRMSC2(COVX,BE,TEMP1,6,1,0)
         do ii = 7,norb
           icnt = ii*(ii+1)/2
           if( covx(icnt).gt.0.d0 ) covx(icnt) = 1.d0/covx(icnt)
         enddo
c         WRITE(6,'(/,a)') 'Covariance matrix for orbits'
c         ICNT=0
c         DO II=1,norb
c           do jj=1,ii
c             ICNT=ICNT+1
c             TEMP1(JJ)=COVX(ICNT)
c           enddo
c           WRITE(6,'(1X,15D12.5)') (TEMP1(IJ),IJ=1,II)
c         enddo

C      Add the weight matrix for orbital parameters to lower triangle normal matrix
            
        call addwgt( indxs,maxorb_cov,norb,covx ) 
c 
c        do j=1,norb
c          JJ1=(J*J-J)/2
c          JJ2=J+INDXS
c          JJ3=(JJ2*JJ2-JJ2)/2
c          DO K=1,J
c            IND1=K+JJ1
c            IND2=K+INDXS+JJ3
c            a(IND2)=a(IND2)+COVX(IND1)
c            ALC(IND2)=ALC(IND2)+COVX(IND1)
c          enddo  
c        enddo
                                                    
  200 continue
c-----end of loop on satellites for orbital parameters     

c      print *,'WSAT last orbit ind2 ',ind2


c     Loop over satellites for SV antenna offsets
     
      if (svantwgt ) then 

        do 320 i=1,nsat 

          indxs=indx + nsat*norb + 3*(i-1)
c          print *,'SV ant i indx norb indxs ',i,indx,norb,indxs

c         create lower triangle of weight matrix 
          do 300 j=1,3
            jj=(j*j-j)/2
            do 300 k=1,j
              ind=k+jj
              covsvant(ind)=0.d0
              if(j.ne.k) goto 300
              covsvant(ind)=sat_apr(i,norb+j)**2
c             print *,'j ind covsvant(ind)',j,ind,covsvant(ind)
c             invert to get weight
              if( covsvant(ind).gt.0.d0 ) then 
                 covsvant(ind) =1.d0/covsvant(ind)
              else
                 covsvant(ind)=0.d0
              endif
  300     continue

c         add weight matrix to lower triangular normal matrix   
          call addwgt( indxs,6,3,covsvant ) 
c*        do j=1,3
c*            jj1=(j*j-j)/2
c*            jj2=j+indxs
c*            jj3=(jj2*jj2-jj2)/2 
c*c            print *,'j jj1 jj2 jj3 ',j,jj1,jj2,jj3
c*            do k=1,j
c*              ind1=k+jj1
c*              ind2=k+indxs+jj3                    
c*c              print *,'k jj3 indxs ind1 ind2 ',k,jj3,indxs,ind1,ind2
c*              a(ind2) = a(ind2) + covsvant(ind1)
c*              alc(ind2) = alc(ind2)+covsvant(ind1)
c*c          print *,'ind1 ind2 covsvant(ind1) ',ind1,ind2,covsvant(ind1)  
c*             enddo
c*          enddo 

  320 continue
 
      endif
            
      RETURN
      END




