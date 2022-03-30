Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1995.   All rights reserved.
C Update EJP 8 July 09 

            subroutine obsmod( iscrn,iuj,ischan,nchan,ichan,jd0,t0
     .      , jfiln,none,sum2,ilast,itimel,freql1,freql2,vel_light
     .      , jdrcvr,trcvr,delay,atmdel,phctdel,reldel,dipoldel
     .      , iongrpdel,ionphsdel,idcb
     .      , dcbcorr,svcepc,svcrat,svcacc,svcl1,clkrat,clkacc,clkcub
     .      , inter,jdc,tc,phase1,phase2,range1,range2 
     .      , simulation,simdel,simdobs )
c
c PURPOSE: subroutine to compute the modelled carrier beat phase
c
c PARAMETERS:
c         IN: iscrn   : unit number for screen output               I*4
c             iuj     : unit number for jfile                       I*4  
c             ischan  : vector of satellite numbers                 I*4(maxsat)
c             nchan   : number of satellites                        I*4
c             ichan   : satellite array index                       I*4
c             jd0     : julian day of initial epoch, also receiver
c                       polynomial ref time                         I*4
c             t0      : time of day in sec of initial epoch, also
c                       receiver polynomial ref time                R*8
c             jfiln   : file name of jfile                          C*16
c             none    : debug flag                                  C*4
c             sum2    : accumulated phase integeral                 R*8(maxsat)
c             ilast   : start epoch for phase integration           I*4(maxsat)
c             itimel  : epoch number                                I*4
c             freql1,freql2 : carrier frequencies constants         R*8
c             vel_light: velocity of light m/sc                     R*8
c             jdrcvr  : PEP Julian Day of receiver observation epoch
c                       corrected for rec clock error.              R*8
c             trcvr   : time of day in seconds of receiver obs
c                       epoch corrected for rec clock error         R*8
c             delay   : theoretical delays (geometry only)          R*8(maxdat,maxsat)
c             atmdel   : delay caused by the dispersive atmosphere   R*8(maxsat)
c             phctdel : delay due to variable antenna phase centre  R*8(2,maxsat)
c             dipoldel: delay due to orientation of receiver
c                       antenna and transmitter                     R*8(2)
c             iongrpdel: 2d & 3rd order ionspheric group delay (L1,L2 R*8 (2)
c             ionphsdel: 2d & 3rd order ionspheric phase delay (L1,L2 R*8 (2)
c             idcb    : flag indicating type of C1/P2 PR correction I*4   
c             dcbcorr : DCB correction for this SV                  R*8
c             svcrat,svcacc : rate and acceleration of satellite
c                             clock                                 R*8
c             clkrat,clkacc,clkcub : rate, acceleration, and cubic
c                      coefficients of the receiver clock poly      R*8
c             inter   : sampling interval of the receiver.
c             simulation: T/F for simulation mode                   L*4
c             simdel  : differential geometric delay at L1 and L2
c                       due to a simulated station displacement     R*8(2)
c
c        OUT: jdc     : Julian Day day for the ref time of sat clock poly I*4
c             tc      : time of day in secs for the ref time
c                       of satellite clock polynomial               R*8
c             phase1  : modeled carrier beat phase L1 (cycles)      R*8
c             phase2  : modeled carrier beat phase L2 (cycles)      R*8
c             range1  : modelled pseudorange  L1      (cycles)      R*8
c             range2  : modelled pseudorange  L2      (cycles)      R*8 
c             svcepc,svcrat,svcacc : offset, rate and acceleration 
c                       of satellite clock                          R*8
c             svcl1   : L1 phase correction due to satellite clock  R*8
c             simdobs : L1,L2,P1,P2 differences (cycles) due to 
c                       simulated station displacement              R*8(4)

c
c SUBROUTINES CALLED: readj, timinc, timdif
c
c CREATED: 27th DEC 1993               LAST MODIFIED: 27th DEC 1993
c
c AUTHOR: rwk, kf, put in SR by S McClusky.
c
            implicit none
c
            include '../includes/dimpar.h'
c
            character*16 jfiln
            character*4 none
c
            integer*4 i,ichan,itimel,inter,icall,ilast(maxsat)
     .      ,iscrn,iuj,ischan(maxsat),nchan,jd0,jd1,jdc
     .      ,jd2,jdrcvr,idcb
c
            real*8 trm1(2),tt1,tt2,t0,tc,svdt,svcepc,svcrat,svcacc,svcl1
     .      ,valid,trmtc,trmat1,trmat2,tint,sum2(maxsat),freql1,freql2
     .      ,phscor,trm2(2),trm3(2),trm4(2),trm5(2),phctdel(2,maxsat)
     .      ,delay(maxdat,maxsat),atmdel(maxsat),fdev,phase1,phase2
     .      ,vel_light,trcvr,clkrat,clkacc,clkcub,timdif,trmt0,trmt2
     .      ,reldel,dipoldel(2),range1,range2,simdel(2),simdobs(4)
     .      ,dcbcorr,iongrpdel(2),ionphsdel(2)

            logical simulation

c           First term is arbitrary initial phase in satellite
            TRM1(1) = 0.d0
            TRM1(2) = 0.d0

c**debug
c           print *,' '
c           print *,'PHSCOR: ichan,jfiln,none,jd0,t0,jdrcvr,trcvr'
c    .             ,  ichan,jfiln,none,jd0,t0,jdrcvr,trcvr
c           print *,'ischan(1-nchan) : ',(ischan(i),i=1,nchan)

C           Second term models effect of drifting satellite oscillator.
C           We must do the integral because SVCRAT and SVCACC
c           are functions of time.
c           Start at the nearest 5s tick to t0 and step along
c           by INTER seconds to handle any missing data.

c           print *,'integrating on chan, t1-1,t2 :'
c    .               ,ichan,ilast(ichan),itimel

            do 333 i = ilast(ichan)+1,itimel
c              update time 1
c              start at an arbitrary (but consistent) time
               jd1 = jd0
               tt1 = 5.0d0*dnint(t0/5.0d0)
               call timinc (jd1,tt1,dble(inter*(i-2)))

c            Get J-file record for time 1

c              secret back door for debugging:
c              Can avoid use of J-file by specifying NONE for name
c              print *,' At t1 jfiln,none: ',jfiln,none
               if(index(jfiln,none) .eq. 0) then
                 icall= 1   
                 call readj ( iuj,ischan,nchan,ischan(ichan)
     .                      , jd1,tt1, svdt,icall,jdc,tc
     .                      , svcepc,svcrat,svcacc,valid )
c                print *,'read JFILE 1, jd1 tt1,tc,svdt,svcepc,svcrat: '
c     .                      ,          jd1,tt1,tc,svdt,svcepc,svcrat
               endif

c              frequency deviation at time 1
               trmtc = timdif(jd1,tt1,jdc,tc)
               trmat1 = svcrat + 2.0d0*svcacc*trmtc

c              update time 2
c              start at an arbitrary (but consistent) time
               jd2 = jd0
               tt2 = 5.0d0*dnint(t0/5.0d0)
               call timinc (jd2,tt2,dble(inter*(i-1)))

c            Get J-file record for time 2

c              secret back door for debugging:
c              print *,' At t2 jfiln,none: ',jfiln,none
c              Can avoid use of J-file by specifying NONE for name
               if(index(jfiln,none) .eq. 0) then
                 icall= 1  
                 call readj ( iuj,ischan,nchan,ischan(ichan)
     .                      , jd2,tt2, svdt,icall,jdc,tc
     .                      , svcepc,svcrat,svcacc,valid )
c                print *,'read JFILE 2, jd2,tt1,tc,svdt,svcepc,svcrat: '
c     .                      ,          jd2,tt1,tc,svdt,svcepc,svcrat

               endif
                                                     

c            Frequency deviation at time 2

               trmtc = timdif(jd2,tt2,jdc,tc)
               trmat2 = svcrat + 2.0d0*svcacc*trmtc

c              time of integration between 1 and 2
               tint  = timdif(jd2,tt2,jd1,tt1)

c              integrate using trapezoidal rule
               sum2(ichan) = sum2(ichan) + tint*(trmat1+trmat2)/2.d0

c              print *,'trmtc,trmat1,trmat2,tint,sum2 :'
c    .                , trmtc,trmat1,trmat2,tint,sum2(ichan)

c              change the epoch number for debugging
               if (itimel .lt. 0 .and. ichan .eq. 1) then
                  write(*,fmt=('(2i4,7(1pe16.8,1x))')) itimel
     .                 , ischan(ichan),svcrat,tt1,tt2,tint,
     .            trmat1*freql1,trmat2*freql1,sum2(ichan)*freql1
               endif
 333        continue

c           now correct for any difference between time 2
c           and the actual receiver time (jdrcvr,trcvr)
c           This is the critical term for Selective Availability
            trmt2 = timdif(jdrcvr,trcvr,jd2,tt2)
            phscor =  trmt2 * (svcrat + 2.0d0*trmt2*svcacc)
c           print *,'trmt2,phscor ',trmt2,phscor

c           change the epoch number for debugging
            if (itimel .lt. 0 .and. ichan .eq. 1) then
               write(*,fmt=('(2i4,3(1pe16.8,1x))')) itimel,ischan(ichan)
     .              ,  svcrat,trmt2,phscor
            endif

C           update most recent time and value
            ilast(ichan) = itimel

            TRM2(1) = FREQL1*(SUM2(ICHAN) + PHSCOR)
            TRM2(2) = FREQL2*(SUM2(ICHAN) + PHSCOR)

c           Third term models the delay
c           time since reference epoch for sat clock poly
            TRMTC = timdif(JDRCVR ,TRCVR ,JDC,TC)
c 
            TRM3(1) = -FREQL1 * (1.0d0 + SVCRAT + 2.0d0*SVCACC*TRMTC)
     .                        * (DELAY(1,ICHAN) + ATMDEL(ichan) +
     .                          PHCTDEL(1,ICHAN) + reldel + dipoldel(1)
     .                          + ionphsdel(1) )
c  must also turn off dipoldel down where range1 and range2 are computed!!!

            TRM3(2) = -FREQL2 * (1.0d0 + SVCRAT + 2.0d0*SVCACC*TRMTC)
     .                        * (DELAY(2,ICHAN) + ATMDEL(ichan) +
     .                          PHCTDEL(2,ICHAN) + reldel + dipoldel(2)
     .                          + ionphsdel(2) ) 
C Liz Debug - Check values of terms
C      Print*,'MODEL\obsmod TRM3(1): ',TRM3(1) ,'reldel:',reldel
C      Print*, 'MODEL\obsmod ATMDEL(satellite):', ATMDEL(ichan),ichan,
C     .  '  Total atmospheric delay (seconds)','Delay1:', DELAY(1,ichan)   
cd           write(*,'(a,i2,4d12.4,f20.4)')
cd     .       'OBSMOD ichan atmdel phctdel dipoldel ionphsdel trm3 '
cd     .    ,ichan,atmdel(ichan),phctdel(1,ichan),dipoldel(1),ionphsdel(1)
cd     .    ,trm3(1)
c  must also turn off dipoldel down where range1 and range2 are computed!!!

c           Fourth term is proportional to square of delay, and
c           is of order 10e-12 cycles.  We neglect it here.
            TRM4(1) = 0.0d0
            TRM4(2) = 0.0d0

c           Fifth term is phase of receiver oscillator.
c           Note that leading order term has already canceled with TRM2
c           Time of nominal time tag after nominal first epoch.
c           (JD0,T0) is the reference epoch for the receiver clock polynomial
c           estimated in S-file.

            TRMT0 = timdif(JDRCVR,TRCVR,JD0,T0)
            FDEV = CLKRAT*TRMT0
     .           + CLKACC*TRMT0*TRMT0
     .           + CLKCUB*TRMT0*TRMT0*TRMT0
            TRM5(1) = -FREQL1 * FDEV
            TRM5(2) = -FREQL2 * FDEV

c           Add them all up. Sign convention is DOPPLER
            PHASE1 = TRM1(1) + TRM2(1) + TRM3(1) + TRM4(1) + TRM5(1)
            PHASE2 = TRM1(2) + TRM2(2) + TRM3(2) + TRM4(2) + TRM5(2)
c           if( itimel.eq.30) then
c            print *,'trm1,trm2,trm3,trm4,trm5,phase1,phase2:  '
c     .             , trm1,trm2,trm3,trm4,trm5,phase1,phase2
c           endif                    

c           save the SV clock effect on the phase to assist AUTCLN 
C           svcl1 = trm2(1) + trm3(1)
* MOD TAH 991004: Only term 2 needed (trm3 contains the light propagation
*           term as well).  Also re-call readj to get the satellite clock
*           at this time and then remove the peice that model has 
*           alreay included.   
            call readj ( iuj,ischan,nchan,ischan(ichan)
     .                 , jdrcvr,trcvr, svdt,icall,jdc,tc
     .                 , svcepc,svcrat,svcacc,valid ) 
            svcl1 = trm2(1)
* MOD TAH 000718: svcepc does not have the expected meaning.  svepc is
*           satellite clock offset at the nearest j-file entry, not at
*           time of this measurement.  The polynomial evaluation at this
*           epoch is actually given by svdt.  Therefore to maintain the
*           intended meaning of svepc (i.e., the part of clock offset
*           not included in the model phase calculation, we should compute
*           it as svdt - trm2(1)/freqL1
C           svcepc = svcepc - trm2(1)/freql1
            svcepc = svdt - trm2(1)/freql1
                                                     
c           For simulation, compute the difference in phase and pseudorange 
c           observations due to a station displacement - enters only through term 3
           
            if( simulation) then
              simdobs(1) = 
     .         -freql1 * (1.d0 + svcrat + 2.d0*svcacc*trmtc) * simdel(1)
              simdobs(2) = 
     .         -freql2 * (1.d0 + svcrat + 2.d0*svcacc*trmtc) * simdel(2) 
              simdobs(3) = simdobs(1)
              simdobs(4) = simdobs(2)
            endif
                
c           For the range, remove the phase wrap (dipoldel) and antenna  phase 
c           centre (phctdel) effects (computed as negative above so add here to remove),
c           and add the correction for the C1-P1 and P2 bias for correlating receivers 
c           rwk 070515: also swap the phase delay for the group delay for 2d & 3rd order
c           ionospheric effects (sign taken care of in sb iondel)
C EJP 8 July 2009 Changed so iongrpdel will be added in the same way as ionphsdel is added in TRM3
C             above (i.e. subtracted). The sign is then converted below as before, so the correction 
C            ends up positive in the pseudorange convention.       
            range1 = phase1 +
     .                       freql1*(dipoldel(1) + phctdel(1,ichan)
     .                       + ionphsdel(1) )-FREQL1*iongrpdel(1)
C              EJP 8 Jul 09 Switched iongrpdel(1) + ionphsdel(1)  for 2 for range2
            range2 = phase2 +
     .                       freql2*(dipoldel(2) + phctdel(2,ichan)
     .                       + ionphsdel(2) )-FREQL2*iongrpdel(2) 
                  
c           3 cases for pseudorange correction due to differential code bias (dcb):
c             (units are nanoseconds from dcb.dat)  
c             0 : No corrections needed (non-cross-correlating receiver)
c             1 : Correct C1 only (bit 28 of data_flag set)
c             2 : Correct C1 and P2 (bit 29 of data_flag set)  
c             sign is positive to 'c' here because of sign reversal below      
c            print *,'idcb freql1 freql2 dcbcorr '
c     .                 ,idcb,freql1,freql2,dcbcorr
c           stop
            if( idcb.gt.0 ) then
              range1 = range1 + freql1*1.d-9*dcbcorr
              if( idcb.eq.2 ) then 
                range2 = range2 + freql2*1.d-9*dcbcorr
              endif
            endif              
 
c           Convert to PSEUDORANGE sign convention.
            PHASE1 = -PHASE1
            PHASE2 = -PHASE2
            range1 = -range1
            range2 = -range2 
            if( simulation ) then
              do i=1,4
                simdobs(i) = -simdobs(i)
              enddo
            endif
                                     
c           print some things

C             Print*, 'phase1', phase1, 'range1',range1
C             Print*, 'phase2', phase2, 'range2',range2
C             Print*, 'iongrpdel1', iongrpdel(1)
C             Print*, 'iongrpdel2', iongrpdel(2)
C             Print*, 'ionphsdel1', ionphsdel(1)
C             Print*, 'ionphsdel2', ionphsdel(2)

c             print *,'svcrat svcacc trmtc ',svcrat,svcacc,trmtc
c             print *,'phctdel ',phctdel
c             print *,'delay atmdel ',delay(1,ichan),atmdel(ichan)
c             print *,'reldel dipoldel ',reldel,dipoldel  
c             print *,'dcbcorr ',dcbcorr

            if (ischan(ichan) .eq. 99 ) then
               write(iscrn,*) ' PRN ',ischan(ichan),' epoch ',itimel
               write(iscrn,*) ' 0    ',jd0   ,t0
               write(iscrn,*) ' RCVR ',jdrcvr,trcvr
               write(iscrn,*) ' C    ',jdc   ,tc
               write(iscrn,*) ' SVCOF',svcepc,svcrat,svcacc
               write(iscrn,'(a,2(1x,1pe21.14))')
     .            'TIMES: TRMT0, TRMTC        ',TRMT0, TRMTC
               write(iscrn,'(a,3(1x,1pe21.14))')
     .            'TRM1',TRM1(1),TRM1(2),TRM1(1)-TRM1(2)
               write(iscrn,'(a,4(1x,1pe21.14))')
     .            'TRM2',TRM2(1),TRM2(2),TRM2(1)-TRM2(2),TRM2(2)/TRM2(1)
               write(iscrn,'(a,4(1x,1pe21.14))')
     .            'TRM3',TRM3(1),TRM3(2),TRM3(1)-TRM3(2),TRM3(2)/TRM3(1)
               write(iscrn,'(a,3(1x,1pe21.14))')
     .            'TRM4',TRM4(1),TRM4(2),TRM4(1)-TRM4(2)
               write(iscrn,'(a,4(1x,1pe21.14))')
     .            'TRM5',TRM5(1),TRM5(2),TRM5(1)-TRM5(2),TRM5(2)/TRM5(1)
               write(iscrn,'(a,4(1x,1pe21.14))')
     .            'PHAS',PHASE1,PHASE2,PHASE1-PHASE2,PHASE2/PHASE1
               write(iscrn,'(a,3(1x,1pe21.14))') 'DEL ',
     .            DELAY(1,ICHAN),DELAY(2,ICHAN),
     .            DELAY(1,ICHAN)-DELAY(2,ICHAN)
               write(iscrn,'(a,3(1x,1pe21.14))') 'RANG',
     .            PHASE1*vel_light/FREQL1,PHASE2*vel_light/FREQL2,
     .            PHASE1*vel_light/FREQL1-PHASE2*vel_light/FREQL2
            endif
            return
            end
