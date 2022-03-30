Copyright 1995 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      SUBROUTINE LWSAT

C     ADD SATELLITE WEIGHTS TO NORMAL MATRIX (change to loose constraint)

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer*4 isat(maxsat)
      common/sats/isat

c       local   
      integer*4 i,j,k,ii,jj,indx,indxs,ind,ier,icnt,nv,norb_cov
c     debug
c      integer*4 jj1,jj2,jj3,ind1,ind2
      real*8 cove,covx,x,e,cove1,covx1,dx,temp,temp1,temp2,be,t,covsvant
     .     , rcond
      dimension covsvant(6)
      DIMENSION COVE(21),COVX(maxorb_cov),X(6),E(6)
      DIMENSION COVE1(21),COVX1(maxorb_cov)
      DIMENSION DX(6,6),TEMP1(maxorb),TEMP2(9),BE(maxorb)
      LOGICAL WOLD

c      Intialize the covariance-matrix length for the orbit model
 
      norb_cov = norb*(norb+1)/2


c       Determine the number of orbital parameters

c         find first orbital element
c         assumption is that orbital elements are adjacent in parameter list
c         with non-gravitational parameters interleaved as well
      INDX=0
      DO 10 I=1,NTPART
       IF(ISLOT1(I).GT.500.AND.ISLOT1(I).LE.1100) GO TO 11
       INDX=INDX+1
   10 CONTINUE
   11 CONTINUE


c         Write the a priori values to the Q-file and log file

c        print *,'LWSAT satwt_type ',satwt_type
c        print *,'LWSAT apr2: ',(sat_apr2(1,j),j=1,maxorb)
        if( satwt_type.eq.'KEP' ) then
          if( logprt ) write(6,118)
          if ( iqflag.eq.1) write(10,118)
  118     format(/,
     1' Keplerian a priori orbital errors (dimensionless except',
     2   ' semi-major axis (km))',/,
     3' Sat#             Semiaxis   Eccen.  Inclin.   Asc.node  '
     4,'Perigee   M.anom.     rad1      rad2      rad3      rad4'
     5,'      rad5      rad6      rad7      rad8      rad9',/)
         do i=1,nsat
           if( logprt ) write(6,119) i,isprn(i),(sat_apr2(i,j),j=1,norb)
           if (iqflag.eq.1)
     .      write(10,119) i,isprn(i),(sat_apr2(i,j),j=1,norb)
         enddo
  119     format(i3,2x,'PRN ',i2,2x,f13.4,1p19e10.1)
        else if( satwt_type.eq.'CAR' ) then
          if( logprt ) write(6,120)
          if ( iqflag.eq.1 ) write(10,120)
  120     format(/,
     1 ' Cartesian a priori orbital errors (km km/s) + force paramters'
     a ,' (fraction))',/,
     2   '   X          Y          Z       Xdot      Ydot      Zdot   ',
     3   '  rad1   rad2   rad3   rad4   rad5   rad6   rad7   rad8   rad9
     4 ')
          do i=1,nsat
            if( logprt ) write(6,121) i,(sat_apr2(i,j),j=1,norb)
            if (iqflag.eq.1) write(10,119) i,(sat_apr2(i,j),j=1,norb)
          enddo
  121     format(i3,6f6.3,1p9e10.1)
        else
           call report_stat('FATAL','SOLVE','lwsat',' '
     .                     ,'Illegal input covariance type',0)
        endif    

c       SV antenna offsets

      if( svantwgt ) then
         if( logprt ) write( 6,125)
         write(10,125)
 125     format(/,' Apriori SV antenna offset errors (m)',/,
     .          ' Sat#    dX       dY       dZ ',/)
         do i=1,nsat
           if( logprt ) 
     .      write( 6,'(i4,3x,3f8.3)') i,(sat_apr2(i,norb+j),j=1,3)
           write(10,'(i4,3x,3f8.3)') i,(sat_apr2(i,norb+j),j=1,3)
         enddo
      endif


c   Loop over satellites for integrated orbital parameters

      do 500 i=1,nsat

       INDXS=INDX+(norb)*(I-1)
       WOLD = .true.
C
C  form keplerian covariance matrix (diagonal)
C
       DO 200 J=1,norb
c         see if constraints were applied in the first solution
c         checks to see if a priori variances < 1.d-10 (0.1 ppb)
          IF (DABS( sat_apr(I,J)).LT.1.0D-20) WOLD = .false.
c          print *,'wold ',wold
C  SAVE STATE-VECTOR FOR THIS SATELLITE
          IF (J.LE.6) X(J)=PREVAL(INDXS+J)
          JJ=(J*J-J)/2
          do 200 K=1,J
           IND=K+JJ
           IF (J.LE.6) COVE(IND)=0.D0
           if ( j.le.6 .and. j.eq.k ) cove(ind) = sat_apr2(i,j)**2
           IF (J.GT.6) COVX(IND)=0.0D0
           if ( j.gt.6 .and. j.eq.k ) covx(ind) = sat_apr2(i,j)**2
           IF (J.LE.6) COVE1(IND)=0.D0
           IF (J.GT.6) COVX1(IND)=0.0D0
           IF (.not. WOLD) GOTO 200
           IF(J.LE.6.AND.J.EQ.K) COVE1(IND)=sat_apr(I,J)**2
           IF(J.GT.6.AND.J.EQ.K) COVX1(IND)=sat_apr(I,J)**2
  200  CONTINUE

c      If the input constraints were Keplerian (normal case), need to
c      convert to non-diagonal Cartesian

       if( satwt_type.eq.'KEP' ) then
c        need first to get Keplerian elements of state vector
         call xyzkep( x,e )
c         if( logprt ) write(6,'(1X,F16.9)') (e(k),k=1,6)
C        now compute the Jacobian for transformation of e to x
         t=0.D0
         nv=6
         call kepxyz( t,e,nv,x,dx )
c         WRITE(6,'(1X,F16.9)') (X(K),K=1,6)
c         WRITE(6,'(6(1X,D12.5))') ((DX(K,J),J=1,6),K=1,6)
C        convert to cartesian covariances
         CALL GPGT(DX,COVE,6,6,COVX,TEMP1,TEMP2)
         IF (WOLD) CALL GPGT(DX,COVE1,6,6,COVX1,TEMP1,TEMP2)
       else if( satwt_type.eq.'CAR' ) then
         do ii=1,21
           covx(ii) = cove(ii)
           covx1(ii) = cove1(ii)
         enddo
       endif
c      save cartesian covariances for h-file
       do ii= 1,norb_cov
           covorbx2(ii,i) = covx(ii)
       enddo
c       print *,'covx ',covx
c       print *,'covx1',covx1

C      Invert the state-vector covariance matrix

       CALL NRMSC2(COVX,BE,TEMP1,6,0,0)
       CALL INVER2(COVX,BE,1,6,rcond,IER)
       if( ier.ne.0 )call report_stat('WARNING','SOLVE','lwsat',' '
     .    ,'Bad inverse of loose-soln state-vector covariance matrix',0)
       CALL NRMSC2(COVX,BE,TEMP1,6,1,0)

       IF (WOLD) THEN
          CALL NRMSC2(COVX1,BE,TEMP1,6,0,0)
          CALL INVER2(COVX1,BE,1,6,rcond,IER)
         if( ier.ne.0 )call report_stat('WARNING','SOLVE','lwsat',' '
     .    ,'Bad inverse of tight-soln state-vector covariance matrix',0)
          CALL NRMSC2(COVX1,BE,TEMP1,6,1,0)
       ELSE
          do ii = 1,norb_cov
             covx1(ii) = 0.0d0
          end do
       ENDIF

c      Invert the non-gravitational parameter covariance matrix

       do ii = 7,norb
          icnt = ii*(ii+1)/2
          if( covx(icnt).gt.0.d0 ) covx(icnt) = 1.d0/covx(icnt)
          if(WOLD.and.covx1(icnt).gt.0.d0 )
     .       covx1(icnt) = 1.d0/covx1(icnt)
       end do
c       print *,'inverted covx ',covx
c       print *,'inversted covx1 ',covx1

C      Add the differences of the loose and tight weight matrices to the normal matrix
    
      do  ii=1,norb_cov
        covx(ii) = covx(ii) - covx1(ii)
      enddo
      call addwgt( indxs,maxorb_cov,norb,covx ) 
                
c      DO 300 J=1,norb
c         JJ1=(J*J-J)/2
c         JJ2=J+INDXS
c         JJ3=(JJ2*JJ2-JJ2)/2
c         DO 300 K=1,J
c            IND1=K+JJ1
c            IND2=K+INDXS+JJ3  
c          print *,'LWSAT: data  ind2 alc ',ind2,alc(ind2)
c            if (l2flag.eq.1.or.l2flag.eq.3.or.l2flag.eq.4 ) then
c               ALC(IND2) = ALC(IND2)+COVX(IND1)-COVX1(IND1)
c      print*,'LWSAT: indx ind1 ind2 alc ',indx,ind1,ind2,alc(ind2)
c            else
c               A(IND2) = A(IND2)+COVX(IND1)-COVX1(IND1)
c            endif
c      print*,'LWSAT: j,k,ind1,covx,covx1 ='
c     .               ,j,k,ind1,covx(ind1),covx1(ind1)
c  300 CONTINUE

  500 CONTINUE


c     Loop over satellites for SV antenna offsets

      do 700 i=1,nsat
     
       do j = 1,6
         covsvant(j) = 0.d0
       enddo

       indxs = indx + nsat*norb + 3*(I-1)
       wold = .true.
                       
c      creat lower triangle of svant weight matrix.
       do 650 j=1,3
         JJ=(J*J-J)/2
         DO 640 K=1,J
           IND=K+JJ
           covsvant(IND)=0.D0
           IF(J.NE.K) GO TO 640
           covsvant(IND)=sat_apr2(i,norb+j)**2
c           print*, 'j ind covsvant(ind)',j,ind,covsvant(ind)
c   temp = (new constraint)/(old constraint)
c   (remember that weight = 1.0/(constraint**2) )
c   weight increment = (new weight) - (old weight)
c                    = 1.0/(new variance) - 1.0/(old variance)
c                    = (1.0 - temp*temp)/(new variance)
c   new weight = old weight + weight increment
c
           if(sat_apr(i,norb+j).gt.0.d0) then
             temp = sat_apr2(i,norb+j)/sat_apr(i,norb+j)
           else
C   no svant constraints in constrained solutions
             temp=1.d0
           endif
C   invert to get weight          
c           print *,'LWSAT i j temp ',i,j,temp
           if(covsvant(ind).gt.0.d0) then
             covsvant(IND)=(1.D0-TEMP*TEMP)/covsvant(IND)
           else
             covsvant(ind)=0.0d0
           endif
  640    continue
  650  continue

c   Add weight matrix to lower triangle normal matrix
     
       call addwgt( indxs,6,3,covsvant )

c      do j=1,3
c        jj1=(j*j-j)/2
c        jj2=j+indxs
c        jj3=(jj2*jj2-jj2)/2
c        do k=1,j
c          ind1=k+jj1
c          ind2=k+indxs+jj3     
c          print *,'SVANT: data  ind2 alc ',ind2,alc(ind2)
c          a(ind2)=a(ind2)+covsvant(ind1)
c          alc(ind2)=alc(ind2)+covsvant(ind1)
c          print *,'ind1 ind2 convsvant ',IND1,IND2,covsvant(IND1) 
c        print*,'SVANT: indx ind1 ind2 alc ',indx,ind1,ind2,alc(ind2)
c        enddo
c      enddo

  700 continue        

      RETURN
      END







