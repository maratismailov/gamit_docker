Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      Subroutine LSQPRT
C                                                        
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'     

      INTEGER IHR,IMN,ISEC,IHNSEC,IYEAR,IMONTH,IDAY,nblen
      CHARACTER*16 NAMSHT
      real*4 hversn

c       Write the header in the Q-file

      WRITE(10,1000) VERS(1:nblen(vers))
1000  FORMAT(//,10X,'Program SOLVE Version ',A,//)

      CALL GETDAT(IYEAR,IMONTH,IDAY)
      CALL GETTIM(IHR,IMN,ISEC,IHNSEC)
      CALL GETUSR(NAMSHT)
CD     WRITE(6,496) IYEAR,IMONTH,IDAY,IHR,IMN,ISEC,OWNER,NAMSHT
      WRITE(10,496) IYEAR,IMONTH,IDAY,IHR,IMN,ISEC,OWNER,NAMSHT
496   FORMAT(' SOLVE Run on ',I4,'/',I2,'/',I2,2X,I2,':',I2,':',I2,/
     .   1x,'OWNER: ',a3,'  OPERATOR: ',a16,//)

c       Write the header in the H-file (single header for tight and loose solutions)   

c         Version number introduced August 2010; prior versions will be read 
c         as 0. by htoglb were mostly changed in a backward-compatible manner, 
c         but for but for documentation sake, assign them (unused) values here
c            0.0 Original version                                 SOLVE 7.31 900901
c            0.1 Add antenna offset                                     8.16 911217  
c            0.2 Add digit to EOP format                                8.17 920112    
c            0.3 Change EOP time format                                 9.13 920611  
c            0.4 Add solution mode (GCR, GCX, GLR, GLX)                 9.21 921230
c            0.5 Add comments                                           9.42 941118
c            0.6 Add reference frame and precession model               9.44 950512   
c            0.7 Change position of elev cut, # zen delays, ant model   9.46 950525
c            0.8 Add radiation pressure model                           9.46 950613     
c            0.9 Add antenna reference point offset                     9.48 950718         
c            1.0 Add error model information and decimation interval    9.56 950831
c            1.1 Change headers for receiver/antenna information        9.57 951002                                                 
c            1.2 Write out a priori covarinaces for SV antenna offsets  9.90 991029
c            1.3 Increase format size for observations                 10.15 041105
c            1.4 Add average atmospheric loading values                10.18 050228
c            1.5 Add ocean-tide and atmospheric tide models            10.23 051115
c            1.6 Add nutation and gravity model                        10.37 071221
c            2.0 Change structure to get full rcvr and SV ant info     10.42 100827
c            3.0 Replace I*4 block number with C*20 svantbody          10.49 141004 
                            
      if (gloprt .and.  
c        standard mode: write h-file only for final solutions
     .  ( (ihmode.eq.0.and.hfiln(6:6).eq.'a')
c        second mode: write h-file only for all constrained solutions
     .   .or. ihmode.ge.1 ) )  then
* MOD TAH 200127: Increased h-file verion to 3.1 to account for change in format
*       of the satellite antenna offsets.
        hversn = 3.1
        write (21,10) hversn,qfiln,VERS(1:nblen(vers))
        write (21,20) iyear,imonth,iday,ihr,imn,isec,owner,namsht
        write (21,40) mfiln    
        write (21,'(1x,a)') 'Datum: geocentric coordinates' 
      endif
 10   format (/,20x,'GAMIT H-file  Version ',f4.1,/
     .       ,1x,'Parallel Q-file: ',a16,/
     .        ,1x,'SOLVE version: ',a)
 20   format(1x,'Running time: ',i4,'/',i2,'/',i2,2x,i2,':',i2,':',i2,
     .   /,1x,'Owner: ',a3,'  Operator: ',a16)
 40   FORMAT (1X,'M-file name: ',A16)

      RETURN
      END









