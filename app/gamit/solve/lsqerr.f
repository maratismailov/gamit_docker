C
      Subroutine LSQERR( covkep,constraints,free_fix,phase_obs )

C     Compute error correlation (and covariance) matrix and
c     formal standard errors (scaled and unscaled)

c     Input
c       nobs  : # observations
c       ncoef : # of coefficients for the fit
c       A     : inverse of normal equation coefficient matrix
c               (covariance matrix)      

        
c     Flags for program flow (common /flags/)
c       iqflag : =1 write the q-file, =2 do not write
c       ioflag : =1 write the o-file, =2 do not write

      Implicit none

      include '../includes/dimpar.h'
      include 'solve.h'    
      include 'parameters.h' 

      INTEGER IDIAG,IOLD
      integer jelf
   

      logical orbpart

      character*4 free_fix,phase_obs
      character*5 constraints

      CHARACTER HEADER*80,message*80,soltyp*29

      integer*4 iparm(maxprm)

      real*8 temp(maxprm),covkep(maxcov)
                          
c     local    
 
      integer*4 iyear,imonth,iday,ihr,imn,isec,ihnsec,numkey
     .        , nlivestn,ind,kl,ij4,ijump,ij,i,j 
      real*8 oldrms

      logical debug/.false./   

      DATA HEADER(1:50)/
     .'       USER  SOLN  DIFF  PHASE CONSTRAINTS BIASES '/
      DATA HEADER(51:80)/
     .'      PARAMETERS      H-FILE  '/
               

c     Compute the prefit and postfit goodness of fit

      OLDRMS=DSQRT(R2SUM/(DBLE(NOBS-NLIVE)))
      if(debug) print *,'LSQERR nobs nlive r2sum oldrms '
     .                        , nobs,nlive,r2sum,oldrms 
      CHI2=SCLERR**2*DBLE(NOBS-NLIVE)
      if(debug) print *,'LSQERR nobs nlive r2sum sclerr oldrms chi2 '
     .       ,nobs,nlive,r2sum,sclerr,oldrms,chi2
         
c     Write the key words to the screen and the q- and o-files

      call keywrd( constraints,free_fix,phase_obs,numkey )
      if( logprt ) WRITE(6,'()')
      if( logprt ) WRITE(6, '(/,A,/,A,15(1X,A5))')
     .      HEADER,' KEYS:',(KEYWORD(I),I=1,NUMKEY)
      if( logprt ) WRITE(6,'()')
      IF (iqflag.EQ.1 ) THEN
c         WRITE(10,'(A,/,A,15(1X,A5))')
          write(10,'(a1,/,a1,15(1x,a5))')
     .    HEADER,' KEYS:',(KEYWORD(I),I=1,NUMKEY)
         WRITE(10,'()')
      ENDIF
      if (IOFLAG.EQ.1 ) THEN
         WRITE(15,'(/,A,/,A,15(1X,A5))')
     .    HEADER,' KEYS:',(KEYWORD(I),I=1,NUMKEY)
         WRITE(15,'()')
      ENDIF

    
c     Write the version, date, time, and input files to screen, q-file, and o-file

      CALL GETDAT (IYEAR,IMONTH,IDAY)
      CALL GETTIM (IHR,IMN,ISEC,IHNSEC)
      if( logprt ) WRITE(6,100)  qfiln,IYEAR,IMONTH,IDAY,IHR,IMN,ISEC
      IF (iqflag.EQ.1) THEN
         WRITE(10,100) qfiln,IYEAR,IMONTH,IDAY,IHR,IMN,ISEC
      ENDIF
      IF(IOFLAG.EQ.1) THEN
         WRITE(15,100) qfiln,IYEAR,IMONTH,IDAY,IHR,IMN,ISEC
      ENDIF
  100 FORMAT (1X,'Ephemeris and survey data files',
     .10X,'(',A16,I4,'/',I2,'/',I2,2X,I2,':',I2,':',I2,')')  
      if (logprt ) write(6,105) tfiln,obfiln(1),cfiln(1) 
      if (iqflag.eq.1 ) write(10,105) tfiln,obfiln(1),cfiln(1) 
      if (ioflag.eq.1 ) write(15,105) tfiln,obfiln(1),cfiln(1) 
  105 FORMAT((1X,A16,1X,A16,1X,A16))
      DO J = 2, nsite
        if( logprt ) WRITE (6,115) obfiln(J),CFILN(J)
         IF (iqflag.EQ.1) WRITE(10,115) obfiln(J),CFILN(J)
         IF (IOFLAG.EQ.1) WRITE(15,115) obfiln(J),CFILN(J)
  115    FORMAT((18X,A16,1X,A16))
      enddo
      if( logprt ) WRITE(6,125) mfiln
      IF (iqflag.EQ.1) WRITE(10,125) mfiln
      IF (IOFLAG.EQ.1) WRITE(15,125) mfiln 
  125 FORMAT (1X,'MERGE File: ',A16)


c     Write the observation summary and solution type to the screen and files

      if( logprt ) WRITE(6,130) (I,I=1,NSAT)
      if( logprt ) WRITE (6,135) (ISEEN(I),I=1,NSAT)
130   FORMAT(/,' Channels used: ',32(4X,I2))
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
135   FORMAT(17X,50I6)
      if(keyword(12).eq.' GCR ') soltyp= 'Constrained bias-free nrms ='
      if(keyword(12).eq.' GCX ') soltyp='Constrained bias-fixed nrms ='
      if(keyword(12).eq.' GLR ') soltyp= 'Loose bias-free nrms ='
      if(keyword(12).eq.' GLX ') soltyp= 'Loose bias-fixed nrms ='
      write(message,'(a,d10.3)') soltyp,sclerr
      call report_stat('STATUS','SOLVE','lsqerr',' ',message,0)
      if( logprt ) write(6,150) nobs,istart,iend,inter,idecim
      if(iqflag.eq.1 )  write(10,150) nobs,istart,iend,inter,idecim
      if(ioflag.eq.1 )  write(15,150) nobs,istart,iend,inter,idecim
  150 format(/,
     .   ' Double-difference observations: ',i6,/
     .,  ' Epoch numbers ', i4, ' to',i5,2x
     .,   'Interval: ',i5, ' s   decimation: ',i3)
      if( logprt ) write(6,151) itor(3,1),itor(1,1),itor(2,1)
     .           , int(tor(1,1)),int(tor(2,1)),tor(3,1)
      if(iqflag.eq.1) write(10,151) itor(3,1),itor(1,1),itor(2,1)
     .                          , int(tor(1,1)),int(tor(2,1)),tor(3,1)
      if(ioflag.eq.1) write(15,151) itor(3,1),itor(1,1),itor(2,1)
     .                          , int(tor(1,1)),int(tor(2,1)),tor(3,1)
  151 format(' Start time: ',4i4,2x,i3,2x,f7.3)
      if( logprt ) write(6,152)    ntpart,nlive,oldrms,sclerr
      if(iqflag.eq.1.or.gloprt) write(10,152) ntpart,nlive,oldrms,sclerr
      if(ioflag.eq.1) write(15,152)           ntpart,nlive,oldrms,sclerr
  152 format(/,
     .   ' Total parameters: ', i5,3x
     .,   'live parameters: ', i5,/
     .,  ' Prefit nrms: ',e12.5,4x
     .,   'Postfit nrms:',e12.5,/
     .,   ' -- Uncertainties not scaled by nrms' )
      IF(iqflag.EQ.1) WRITE(10,130) (I,I=1,NSAT)
      IF(iqflag.EQ.1) WRITE(10,135) (ISEEN(I),I=1,NSAT)
      IF(IOFLAG.EQ.1) WRITE(15,130) (I,I=1,NSAT)
      IF(IOFLAG.EQ.1) WRITE(15,135) (ISEEN(I),I=1,NSAT)


c     Compute the uncertainties

      IF(NLIVE.LE.0) GO TO 990
      DO I=1,NLIVE
        IDIAG=I*(I+1)/2
        SIGMA(I)=DSIGN(DSQRT(DABS(A(IDIAG))),A(IDIAG))
CD      WRITE(6,'(1x,3(e22.15,1x)') SIGMA(I)
      enddo


c     Compute the Keplerian covariances and save for LSQDO1
c     (necessary here because the covariance matrix is normalized in the
c     correlation matrix computation)
c     Do if any orbit partial present
      orbpart = .false.
      do ind=1,ntpart
       IF((ISLOT1(IND).GT.500).AND.(ISLOT1(IND).LE.2400)) orbpart=.true.
      enddo
      if(debug) print *,'LSQERR covkep(1-42) ',(covkep(i),i=1,42)
      if( orbpart )  call covst1( covkep )
          

c     Compute the correlation matrix

      IJ=0
      NLIVESTN = 0
      DO I=1,NTPART
         KL=ISLOT1(I)
         IF (KL.GT.0 .AND. KL.LE.300 .AND. FREE(I).EQ.1 )
     .     NLIVESTN = NLIVESTN+1 
         If (free(i).eq.1 ) then
          ij = ij + 1
          iparm(ij) = i
         endif
      enddo
c     write to the o-file the correlation matrixs for only the live station parameters
      IF(IOFLAG.EQ.1) WRITE(15,460)
  460 FORMAT (/,1X,'Error correlation matrix (5-decimal accuracy):',/)
      IOLD=1
      DO I=1,NLIVESTN  
c***  temp write out all correlations 
c*      do i=1,nlive
*       MOD TAH 980310: Replaced calc by jelf call
        IDIAG=jelf(I+1)           !  = I*(I+1)/2
        DO  J=1,I
          IJ4=IOLD+J-1
          TEMP(J)=A(IJ4)/(SIGMA(I)*SIGMA(J))
CD        WRITE(6,'(1x,d22.15)')) A(IJ4)
        enddo
CD       WRITE(6,510) IPARM(I),(TEMP(J),J=1,I)
        IF(IOFLAG.EQ.1) WRITE(15,510) IPARM(I),(TEMP(J),J=1,I)
  510   FORMAT (1X,I3,'. ',10(10F9.5,:,/,6x))
        IOLD=IDIAG+1
      enddo           
c     write to the q-file the largest correlation coefficients
      if( constraints.ne.'loose') then 
        if( logprt ) WRITE(6,620) correl_prt
        IF(iqflag.EQ.1) WRITE(10,620) correl_prt
        IF(IOFLAG.EQ.1) WRITE(15,620) correl_prt
 620    FORMAT(/' Correlation coefficients greater than',F9.6,':'/)
        IJ4=0
        IJUMP=0
        DO 660 I=1,NLIVE
          DO 650 J=1,I
            IJ4=IJ4+1
            IF(I.EQ.J) GO TO 650
            TEMP(J)=A(IJ4)/(SIGMA(I)*SIGMA(J))
            IF(DABS(TEMP(J)).LT.correl_prt) GO TO 650
            IJUMP=IJUMP+1
            if( logprt ) WRITE(6,640) RLABEL(IPARM(I)),
     1      RLABEL(IPARM(J)),TEMP(J)
            IF(iqflag.EQ.1) WRITE(10,640) RLABEL(IPARM(I)),
     1                             RLABEL(IPARM(J)),TEMP(J)
            IF(IOFLAG.EQ.1) WRITE(15,640) RLABEL(IPARM(I)),
     1                             RLABEL(IPARM(J)),TEMP(J)
  640       FORMAT(2X,A20,2X,A20,' : ',F8.5)
  650     CONTINUE
  660   CONTINUE
c       or indicate that no correlations are greater than the threshold
        IF(IJUMP.EQ.0) THEN
          if( logprt ) WRITE(6,670)
  670    FORMAT(' None')
          IF(iqflag.EQ.1) WRITE(10,670)
          IF(IOFLAG.EQ.1) WRITE(15,670)
        END IF
c     endif on whether tight or loose solution
      endif  
                            
c     Write the h-file header block

      if (gloprt .and.  
c        standard mode: write h-file only for final loose solutions
     .  ( (ihmode.eq.0.and.constraints.eq.'loose'.and.hfiln(6:6).eq.'a')
c        second mode: write h-file only for (all) constrained solutions
     . .or. (ihmode.eq.1.and.constraints.eq.'tight') 
c        third mode: write h-file for all solutions
     .   .or. (ihmode.ge.2) ) ) then      
       call glohed( numkey,oldrms,iparm )  
      endif

  990 continue
              
      RETURN
      END
