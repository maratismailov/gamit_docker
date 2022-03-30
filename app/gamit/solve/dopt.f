      Subroutine DOPT( mbias0, last_nonbias )
        
C     Construct the bias-D operator by reading the choice of baselines and satellites
c     from those selected by AUTCLN, or ordering by baseline length and selecting a
c     base satellite by number of observations.  Original routine written by Danan Dong  
c     December 1987; revised included multi-session by D. Dong February 1987; revised
c     to include AUTCLN input by R. King December 2003.
       

c     Dong's logic for constructing the bias-D operator:
 
c        +---> (clean dead stations and dead bias perameters)
c        |
c        +---> (rearrange live stations and live biases)
c        |
c        +---> (calculate and sort out all baselines by their lengths)
c        |
c        +---> (determine base-satellite with most observations) (ibasat)
c        |
c        +---> (construct mapping operator with base-satellite)
c           |
c           +---> (loop over other satellites) (isa)
c           |
c           +---> (find first shortest baseline with isa)
c           |                                              <------+
c           +---> (find second shortest baseline with isa)        |
c           |                                              (no)   |
c           +---> (check the independence) ----------------------->
c              |  (yes)
c              +---> (get one effective mapping operator)
c        |
c        +---> (there are stations no observation with ibasat ?) (miss)
c        |   (yes)
c        |   +--> (construct mapping operator with miss) (same logic as above)
c        |---|
c            +--> (no)
c        |
c        +---> DBADST (construct operator with bad stations) (not used)
c        |
c        +---> REMENU (reorder effective dd-biases by baseline lengths)
c        |
c        +---> (fill indicators)
      

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h' 
      include 'parameters.h'  

      integer*4 mbias0,ibnum,ix1,nbad,nlen,i2,i3
     .        , iold,isa,isum0,id2,last_nonbias,ibsat,id1,llbad
     .        , i1,li,ifirst,isum,ngood,lbias1,lbias1dd,ix2,isum1,ibasat
     .        , miss,iswtch,iloop,ilen,i10,j1o,j1,j2,iau,ii,jj,mm,i,j,k 
              

      INTEGER MAXBAS
      real*8 baslin,bln   
      PARAMETER (MAXBAS=MAXSIT*(MAXSIT-1)/2)
      LOGICAL FIND
      character*256 message

      integer ISTAT(MAXSIT)
      integer  NDEX(MAXBAS*2),IBAS(MAXSIT)
      integer  LSTAT(MAXSIT)
      real*8 BASLEN(MAXBAS)
      integer IBAD(MAXSIT),JBSAT(MAXSAT)

      logical debug/.false./
                          

c---Initialize the counters

C---- Number of bad stations
      LLBAD = LBFRE(1)
      I1 = nsite*NSAT 
C---- Find bias parameters and their indices
      IFIRST =  IDXB(MBIAS0+1)
      IBNUM =  IDXB(1)
      IF (IFIRST.LT.0) IFIRST = -IFIRST
      IF (IBNUM.LT.0) IBNUM = -IBNUM
      IBNUM = last_nonbias - IBNUM+1  
c---- Kill extra biases for L1 only and LC mode.
      if( l2flag.le.1 ) then  
        ix1 = ifirst + i1
        do i= ix1, ix1+i1-1
          if( free(i).ne.0 ) then
            free(i) = 0
            nlive = nlive -1
          endif
        enddo
       endif
       if(debug) print *,'DOPT before bias setting 1-way nlive = ',nlive 
       if(debug) print *,'  from earlier call to DOUBLE and OPERA ',
     .      ' irowd irowdt ipntd ipntdt all zero'

C---- Find stations without effective observations (dead stations)
      NBAD = 0  
      DO 10 I = 1,nsite
         ISUM = 0
         DO 20 J = 1,NSAT
            K = IFIRST+NSAT*(I-1)+J-1
C---- for manually fixed biases, iuse equivalent to zero
            IF (FREE(K).EQ.0) iuse(I,J) = 0
            ISUM = ISUM+iuse(I,J)
            IF(iuse(I,J).GT.0) GO TO 20
C---- Kill bias with no observation
            FREE(K) = 0
            IF(IBAND.EQ.2) FREE(K+I1) = 0
 20      CONTINUE
         IF (ISUM.EQ.0) THEN
            NBAD = NBAD+1
            GOTO 10
            ENDIF
C---- Kick out bad stations temporarily.
         IF (LLBAD.GT.0) THEN
           DO LI = 1,LLBAD
             IF (I.EQ.LBAD(LI)) THEN
               NBAD = NBAD+1
               GOTO 10
             ENDIF   
           enddo
         ENDIF
      ISTAT(I-NBAD) = I  
  10  continue
      NGOOD = nsite-NBAD



C---- Calculate the length-ordered sequence of all possible independent baselines
                
c     ndex   = list of all independent baselines, sorted by length 
c              e.g., 1 5  2 6   means BL 1 is site1/site5, BL 2 is site2/site6 
c              This array is filled in DOPT because the full list is needed
c              in forming the normal equatioins for the DD biases from the one-way
c              biases; that is, there seems to be an implicit order of independent
c              baselines that includes all stations.
                                 
      call zero1i(1,maxsit,ibas)
      call zero1i(1,maxbas*2,ndex)
      id1 = 0 
      ISUM = 0        
      if(debug) print *,'DOPT ngood ',ngood
      DO I = 1,NGOOD-1
         II = ISTAT(I)  
         IX1 = 3*(II-1)+1
         DO J = I+1,NGOOD
            JJ = ISTAT(J)
            IX2 = 3*(JJ-1)+1
            I1 = ISUM*2
            ISUM = ISUM+1
            BASLEN(ISUM) = BASLIN(COORDS(IX1),COORDS(IX2))
            NDEX(I1+1) = II
            NDEX(I1+2) = JJ  
           if(debug) print *,'  i j jj i1 x1 ix2 coords isum baslen  '
     .      ,i,j,jj,i1,ix1,ix2,coords(ix1),coords(ix2),isum,baslen(isum)
         enddo
       enddo
      NLEN = ISUM      
      if(debug) then
         print *,'DOPT initial nlen ngood isum ndex '
     .    ,nlen,ngood,isum
         write(*,'(10i7)') (ndex(i),i=1,50)
       endif
C---- Sort out all baselines in sequences of length.
      IF (nlen.GT.1) CALL SORTBL(NLEN,2,BASLEN,NDEX)
      if(debug) then 
         print *,'Sorted nlen ndex ',nlen
         write(*,'(10i7)') (ndex(i),i=1,50)
       endif


c---Get the bias D-operator from the combinations choosen by AUTCLN
               
      if(debug) print *,'DOPT l2flag ',l2flag
      if( l2flag.eq.4 ) then    

        if( nfiln(1:1).eq.' ')  call report_stat('FATAAL','SOLVE','dopt'
     .   ,' ','LC_AUTCLN requested but no N-file (pre solution?)',0)
        call read_biases( last_nonbias,lbias1,lbias1dd,nlen,ndex )  

c        returns the number of NL or WL DD biases (lbias1dd), the baseline
c        indices (ndex), and (in common) the 'live' bias list (ipnt2d),
c        the bias operator pointers (ipntd), and the pointers to
c        baselines for these biases (ipntdt)   
c**      returns also the nubmer of live DD combinations (lbias1 - nsite*nsat - dead)
c**      but this may be from unnecessary code, and in any case is not a logical
c**      calculation for read_biases; revisit this later
c       print *,'At end of new-style nlive= ',nlive   
c      store the number of viable DD biases in common /bbii/ for get_widelane 
       l1bias = lbias1dd

      if(debug) then
        print *,'After read_biases irowdt '
        write(*,'(10i7)') (irowdt(i),i=1,50)     
      endif 


      else


c---Construct the bias D-operator here

          
C---- Initialize
      CALL ZERO1I(1,MAXND,IROWDT)
      CALL ZERO1I(1,NSAT,JBSAT)    

C---- Sort out live bias parameters.
      ISUM = 0
      IX1 = 0
      DO I = 1,nsite
         DO 320 J = 1,NSAT
            ISUM = ISUM+1
            ISIGMA(ISUM) = 0
            IF (iuse(I,J).LT.1)  GOTO 320
            IX1 = IX1+1
            ISIGMA(ISUM) = IX1
            IPNT2D(IX1) = ISUM
            JBSAT(J) = JBSAT(J)+1
 320     CONTINUE
      enddo    

C---- LBIAS1 : Number of live 1-way L1 bias parameters.
      LBIAS1 = IX1    
      if(debug) then 
        print *,'DOPT 1 ways '
        print *,'ngood nbad lbias1 ',ngood,nbad,lbias1
        print *,'isum isigma ',isum
        write(*,'(10i7)') (isigma(i),i=1,isum)
        print *,'ipnt2d ix1 ',ix1
        write(*,'(10i7)') (ipnt2d(i),i=1,50)
      endif

      


C---- Find satellite with most observations
      I2 = 0
      I3 = 0 
      ibsat = 0 
      FIND = .false.
      DO J = 1,NSAT
         IF (JBSAT(J).EQ.nsite) FIND = .true.
         ISUM = 0
         ISUM1 = 0
         DO I = 1,NGOOD
            I1 = ISTAT(I)  
            ISUM = ISUM+iuse(I1,J)
         enddo
         IF (JBSAT(J).EQ.nsite) ISUM1 = ISUM   
         IF (ISUM.GT.I2) THEN
            I2 = ISUM
            IBASAT = J
         ENDIF
         IF (ISUM1.GT.I3) THEN
            I3 = ISUM1
            IBSAT = J
         ENDIF
      enddo

C---- If the base-satellite missed some one-way observations?
      MISS = 0   

      IF (FIND.AND.I3.GT.0) THEN
         IBASAT = IBSAT
         GOTO 525
      ENDIF         
c      print *,'DOPT istat ',(istat(i),i=1,7)
      DO I = 1,NGOOD
         I1 = ISTAT(I)
         IF (iuse(I1,IBASAT).EQ.0) THEN
            MISS = MISS+1
            ISTAT(I) = -I1 
         ENDIF
      enddo    
      if(debug) print *,' ngood ibasat miss iuse '
     .  ,ngood,ibasat,miss,(iuse(i,ibasat),i=1,7)
 
 525  MM = 0  
      IROWD(2) = 0
      ISWTCH = 2
      ILOOP = 0
      IOLD = 0  
      

C---- DO loop with all other satellites---------(long loop )

      DO 590 ISA = 1,NSAT       

      IF (ISA.EQ.IBASAT) GOTO 590
      ILOOP = ILOOP+1
C---- Compare observations to ISA with that to previous satellite
      IF (ILOOP.GT.1) THEN
         DO 530 I = 1,NGOOD
           I1 = ISTAT(I)  
           IF (I1.LE.0) GOTO 530  
           IF (iuse(I1,ISA).GT.0.AND.iuse(I1,IOLD).EQ.0) GOTO 535
           IF (iuse(I1,ISA).EQ.0.AND.iuse(I1,IOLD).GT.0) GOTO 535
 530     continue
         GOTO 585
      ENDIF
C---- Pick out stations with no observations to ISA
 535  NBAD = 0
      DO 540 I = 1,NGOOD 
         I1 = ISTAT(I)      
         IF (I1.LE.0) GOTO 540
         IF (iuse(I1,ISA).EQ.0) THEN
         NBAD = NBAD+1
         IBAD(NBAD) = I1
         ENDIF
 540  CONTINUE
C---- For satellite-selecting option. 
      IF (NBAD.EQ.NGOOD-MISS) THEN
        write(message,'(a,i2)') 'No observation to satellite ',isa
        call report_stat('WARNING','SOLVE','dopt',' ',message,0)
        ILOOP = ILOOP-1
        GOTO 590
      ENDIF
C---- Find first short baseline with observation to ISA
      DO 545 I = 1,NLEN
         JJ = (I-1)*2  
         I1 = NDEX(JJ+1)
         I2 = NDEX(JJ+2)
         IF (iuse(I1,ISA).EQ.0.OR.iuse(I2,ISA).EQ.0)
     *   GOTO 545
C---- Check if there is a station belongs to MISS.
         IF (iuse(I1,IBASAT).EQ.0.OR.iuse(I2,IBASAT).EQ.0)
     *   GOTO 545
         IBAS(1) = I
         ISUM0 = 2*I
         GOTO 550
 545  CONTINUE
      write(message,'(a,i2,a,i2)') 'No double difference for satellite '
     .     ,ibasat,' and satellite ',isa
      call report_stat('WARNING','SOLVE','dopt',' ',message,0)
      ILOOP = ILOOP-1
      GOTO 590
C---- II : Number of selected stations
 550  ISUM = ISUM0                         
c      print *,' isum ibas ',isum
c      write(*,'(10i7)') (ibas(i),i=1,8)
      II = 2
      LSTAT(1) = I1
      LSTAT(2) = I2
C---- Skip baseline searching in single baseline case
      IF (II.GE.NGOOD-MISS-NBAD) GOTO 180
C---- Find nearest station to all selected stations
C---- I1,I2 : index of two stations of next searched baseline
 150  I1 = NDEX(ISUM+1)
      I2 = NDEX(ISUM+2)
C---- Check if there is a station belongs to MISS.
      IF (iuse(I1,IBASAT).EQ.0.OR.iuse(I2,IBASAT).EQ.0) THEN 
         ISUM = ISUM+2
         GOTO 150
      ENDIF
      IF (NBAD.EQ.0) GOTO 560
      DO 570 I = 1,NBAD
         JJ = IBAD(I)
         IF (I1.NE.JJ.AND.I2.NE.JJ) GOTO 570
         ISUM = ISUM+2
         GOTO 150
 570  CONTINUE
 560  DO 130 IX1 = 1,II
         IF (I1.EQ.LSTAT(IX1).OR.I2.EQ.LSTAT(IX1)) GOTO 140
 130  CONTINUE
C---- No overlaped station. Shift to next baseline
      ISUM = ISUM+2
      GOTO 150
C---- Second baseline is selected
 140  IF (I1.EQ.LSTAT(IX1)) ID1 = I2
      IF (I2.EQ.LSTAT(IX1)) ID1 = I1
      DO 160 ID2 = 1,II
C---- This baseline has been used already.
         IF (ID1.EQ.LSTAT(ID2)) THEN  
         ISUM = ISUM+2
         GOTO 150
         ENDIF
 160  CONTINUE
C---- New independent baseline is found.
      IBAS(II) = ISUM/2+1    
      II = II+1
      LSTAT(II) = ID1
      IF (ISUM.EQ.ISUM0) ISUM0 = ISUM0+2
      ISUM = ISUM0  
c      print *,' New isum ibas ',isum
c      write(*,'(10i7)') (ibas(i),i=1,15)

C---- All independent baselines have been found
      IF (II.EQ.NGOOD-MISS-NBAD) GOTO 180
      GOTO 150
 180  CONTINUE    

C---- Sort out all selected baselines according to their lenghth
         ISUM = NGOOD-MISS-NBAD-1  
         CALL SORT1I(ISUM,IBAS) 
c      print *,'nbad miss ngood isum ',nbad,miss,ngood,isum    
c      print *,' sorted ibas ',isum
c      write(*,'(10i7)') (ibas(i),i=1,8)
C---- Construct D-opt operator
 585  DO K = 1,ISUM
         ilen = ibas(k)
         IX1 = (Ilen-1)*2
         I1 = NDEX(IX1+1)
         I2 = NDEX(IX1+2)
C---- J1 :  base-satellite
C---- J2 :  second live satellite
         J1 = IBASAT
         J2 = ISA
c---- fill D-operator  
c         print *,' call FILLD mm mrowd ',mm,iswtch
         CALL FILLD(I1,I2,J1,J2,ISWTCH,MM,Ilen) 
      enddo
      IOLD = ISA
 590  CONTINUE
c      print *,'D-operator for biases iswtch mm ipntd ',iswtch,mm
c      write(*,'(10i7)') (ipntd(i),i=1,mm)  
c      print *,'  irowd '
c      write(*,'(10i7)') (irowd(i),i=1,mm)
c      print *,'  d '
c      write(*,'(10f7.3)') (d(i),i=1,150)    
c      print *,' irowdt '
c      write(*,'(10i7)') (irowdt(i),i=1,50)     

c--------END OF LOOP OVER OTHER SATELLITES


                          
c---- Reconstruct the D-opt operator (add independent bias parameters) 
c    if some base-satellite observations missed

      IF (MISS.EQ.0) GOTO 600
      DO 610 I = 1,NGOOD
         I1 = ISTAT(I)
         IF (I1.GT.0) GOTO 610
         I1 = -I1
         I10 = I1
C----    Find first live satellite for station I1
         DO  J = 1,NSAT
            IF (iuse(I1,J).GT.0) GOTO 630
         enddo
 630     J1  =  J
 635     J1O  =  100
C---- Second live satellite for station I1
         DO 640 J2 = J1+1,NSAT
         I1 = I10
         IF (iuse(I1,J2).EQ.0) GOTO 640
C---- Search shortest baseline with observation to J1 and J2
         FIND  =  .false.
         DO 650 II = 1,NLEN
            ID1 = NDEX(II*2-1)
            ID2 = NDEX(II*2)
C---- The baseline must be connected with station I1.
            IF (ID1.NE.I1.AND.ID2.NE.I1) GOTO 650
            IF (ID1.EQ.I1) I2 = ID2
            IF (ID2.EQ.I1) I2 = ID1
            IF (iuse(I2,J1).EQ.0.OR.iuse(I2,J2).EQ.0) GOTO 650
C---- Check if station I2 is missed by satellite IBASAT
            IF (I2.GT.I1.AND.I.LT.NGOOD) THEN
               DO 660 IX1 = I+1,NGOOD
                  ID1 = ISTAT(IX1)
                  IF(ID1.LT.0.AND.I2.EQ.-ID1) GOTO 650
 660           CONTINUE
            ENDIF
C---- Match the order to NDEX
            IF (I1.GT.I2) THEN
               ID1 = I1
               I1 = I2
               I2 = ID1
            ENDIF
c---- fill D-operator
            CALL FILLD(I1,I2,J1,J2,ISWTCH,MM,II)
            GOTO 640
 650     CONTINUE
c** the following statement is redundant because FIND is always false
c   if you've branched to here from above 660
         IF (.not.FIND) THEN
            IF (J2.LT.J1O) J1O  =  J2
         ENDIF
 640  CONTINUE
      IF (J1O.LT.NSAT) THEN
         J1  =  J1O
         GOTO 635
      ENDIF
 610  CONTINUE  
c      print *
c     . ,'After handling MISS D-operator for biases iswtch mm ipntd '
c     .         ,iswtch,mm
c      write(*,'(10i7)') (ipntd(i),i=1,mm) 
c      print *,'  irowd '
c      write(*,'(10i7)') (irowd(i),i=1,mm)
c      print *,'  d '
c      write(*,'(10f7.3)') (d(i),i=1,mm) 
c      print *,' irowdt '
c      write(*,'(10i7)') (irowdt(i),i=1,50)     

c---  end branch for case miss > 0 

 
C---- Total number of independent double-difference live L1 bias parameters
 600  lbias1dd = MM/4  
c      print *,'Enlarged D-operator for biases mm lbias1dd ipntd '
c     .   ,mm,lbias1dd
c      write(*,'(10i7)') (ipntd(i),i=1,mm)    
c      print *,'   ipntdt '
c      write(*,'(10i7)') (ipntdt(i),i=1,40)  
c      print *,' irowdt '
c      write(*,'(10i7)') (irowdt(i),i=1,50)     


C---- Reorder D-optimal by length of baselines
      CALL REODRD(lbias1dd,4)    

c      print *,'Reordered D-operator for biases mm ipntd ',mm
c      write(*,'(10i7)') (ipntd(i),i=1,mm) 
c      print *,'   ipntdt '
c      write(*,'(10i7)') (ipntdt(i),i=1,40)

C---- CONSTRUCT BAD STATION BASELINE BIAS OPERATOR
      IF (LLBAD.GT.0) THEN
         call dbadst(llbad,lbias1dd,ngood,istat,baslen,coords,ndex
     .              ,iswtch,last_nonbias)   
c      print *,'After DBADSTiswtch mm ipntd ',iswtch,mm
c      write(*,'(10i7)') (ipntd(i),i=1,mm) 
c      print *,'  irowd '
c      write(*,'(10i7)') (irowd(i),i=1,mm)
c      print *,'  d '
c      write(*,'(10f7.3)') (d(i),i=1,mm) 
c      print *,' irowdt '
c      write(*,'(10i7)') (irowdt(i),i=1,50)     

      endif   

      NRD = ISWTCH-2
      IROWD(1) = ISWTCH-2   

c      print *,'At end of old-style lbias1dd nlive ',lbias1dd,nlive      

c*** end IF on whether or not to read biases from external file
      endif
              

c--- Find the index of the longest baseline to be used for narrow-lane resolution

        limitb = 30000
        do i = 1,lbias1dd
           i1 = ipntdt(i)
           bln = baslen(i1)
           if(bln.gt.nldmax) then
              limitb = i
              goto 700
           endif
        enddo   
  700   continue



C---- Recalculate number of live parameters and reset the bias parts of the free and adjust arrays          

c       = original # from all parameters (including 1-way biases)
c         corrected by the difference between original 1-way biases and actual DD biases  

      I1 = nsite*NSAT
      NLIVE = NLIVE-I1+lbias1dd  
      IF (L2FLAG.GE.2) NLIVE = NLIVE-I1+lbias1dd
      II = last_nonbias + lbias1dd*IBAND
      IF (L2FLAG.EQ.1) II = II-lbias1dd
      IX1 = IBNUM
c     put all live bias parameters on the top of B.
      DO I = last_nonbias+1,II
         FREE(I) = 1
         IX1 = IX1+1
         IDXB(IX1) = I
      enddo
c     set remaining bias parameters fixed.
      I2 = IFIRST+I1*IBAND-1
      DO I = II+1,I2
         FREE(I) = 0
         IX1 = IX1+1
         IDXB(IX1) = -I
         ADJUST(I) = 0.0D0
      enddo     
c      print *,'DOPT Bias indices idxb 1-98 '
c      write(*,'(10i7)') (idxb(i),i=1,98)   
c      print *, 'free(211-270) '
c      write(*,'(10i7)') (free(i),i=211,270)
 

C-- Construct the double-difference bias labels for Q-file. 

c      print *,'Before REMENU lbias1dd ndex ',lbias1dd
c      write(*,'(10i7)') (ndex(i),i=1,50) 
c      print *,'  ipntdt '
c      write(*,'(10i7)') (ipntdt(i),i=1,140) 
c      print *,'  ipnt2d '
c      write(*,'(10i7)') (ipnt2d(i),i=1,50)   
c      print *,'  ipntd '
c      write(*,'(10i7)') (ipntd(i),i=1,140)
      if(debug) then 
        print *,'DOPT bef REMENU lbias1dd last_nonbias ndex(1-50) '
     .   ,lbias1dd,last_nonbias
        write(*,'(10i7)') (ndex(i),i=1,50) 
      endif  
      call remenu( lbias1dd,ndex,last_nonbias ) 
      if(debug)  then 
        print *,'DOPT aftr REMENU last_nonbias rlabel',last_nonbias
        do i=last_nonbias+1,last_nonbias + 2*lbias1dd
         write(*,'(a20)') rlabel(i)
        enddo   
      endif 
            

c      print *,'At end of DOPT,  lbias1dd nlive ',lbias1dd,nlive      
      mm = lbias1dd*4
              

C---- Locate index of D-optimal operator by column     
c      print *,' DOPT Before ISIGMA reset mm ',mm     
c      print *,'  isigma '
c      write(*,'(10i7)') (isigma(i),i=1,140) 
c      print *,'  ipntdt '
c      write(*,'(10i7)') (ipntdt(i),i=1,140) 
      NCD = LBIAS1   
      mm = lbias1dd*4
      IROWDT(1) = NCD  
c      print *,'lbias1 ncd mm ',lbias1,ncd,mm 
c      print *,'Before summing irowdt '
c      write(*,'(10i7)') (irowdt(i),i=1,50)     
      DO  I = 1,NCD
         ISIGMA(I) = 1
         IROWDT(I+2) = IROWDT(I+2)+IROWDT(I+1)
      enddo                    
      
      DO I = 1,MM
         IX1 = (I-1)/4+1
         ID1 = IPNTD(I)
         IAU = IROWDT(ID1+1)+ISIGMA(ID1)
         IPNTDT(IAU) = IX1
         DT(IAU) = D(I)  
         ISIGMA(ID1) = ISIGMA(ID1)+1     
c         print *,'i ix1 id1 iau ipntdt dt isigma '
c     .     ,i,ix1,id1,iau,ipntdt(iau),dt(iau),isigma(id1)  
      enddo  

c       print *,'lbias1dd d ',lbias1dd
c       write(*,'(10f7.0)') (d(i),i=1,150)
c       print *,'   dt '
c       write(*,'(10f7.0)') (dt(i),i=1,150)

c       print *,' nsite nsat iuse ',nsite,nsat
c       write(*,'(10i7)') ((iuse(i,j),j=1,nsat),i=1,nsite)                                        
c       print *,' accumulated irowdt '
c       write(*,'(10i7)') (irowdt(i),i=1,50)     
c   *** I think that isigma is a temporary array here, reset again
c   *** Test this by setting it to zero  
c      do i=1,maxprm
c       isigma(i) = 0
c      enddo
c***  end test
c      print *,'After isigma ipntdt reset '
c      print *,'  isigma = 0 '
c      write(*,'(10i7)') (isigma(i),i=1,mm+4)     
c      print *,' ipntdt '
c      write(*,'(10i7)') (ipntdt(i),i=1,160)  
c      print *,'  irowd '
c      write(*,'(10i7)') (irowd(i),i=1,50)

C---- L1BIAS: Number of total independent L1 biases.
      L1BIAS = lbias1dd
c* 300  CONTINUE
       continue
c      write(6,*) 'l1bias',lbias1dd
c      do 301 icnt=1,ntpart
c       write(6,*) rlabel(icnt)
c  301 continue    

      RETURN
      END

