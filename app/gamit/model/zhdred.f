Copyright (c) Aero Service Division, 1986. All rights reserved.
      SUBROUTINE ZHDRED( IUZ,iprnt,FJDZ1,FJDZ2 )
C
C.....Interpolate brightness temperatures from table
C.....Y. Bock (Modification of ATMRED) - 2/5/87
C     J. Davis - 2/15/87
C     R. King  6/3/87 - split into ZHDRED and WVRRED
C
      implicit none
C
      CHARACTER*4 SCARR(20),SITECD,UPPERC
      character*80 fn

      INTEGER*4 JULDAY,month,ihr1,ihr2,imn1,imn2,iday,idoy1,idoy2
     .        , iyr1,iyr2,iuz,iprnt,i,ioerr

      real*8 fjdt1,fjdt2,fjdz1,fjdz2
C
      ioerr = 0
CD    ISCRN=6

      inquire ( unit=IUZ, name=fn, iostat=ioerr )
      if (ioerr .ne. 0 ) then
         call report_stat('FATAL','MODEL','zhdred',fn,
     1   'Z-file does not exit or is not opened:',ioerr)
      endif
C
C.....FIND THE WORD "END" TO INDICATE THE END OF THE INFO HEADER S
    1 READ(IUZ,2,iostat=ioerr) (SCARR(I),I=1,20)
    2 FORMAT(20A4)
      if (ioerr .ne. 0 ) then
         call report_stat('FATAL','MODEL','zhdred',fn,
     1   'Error reading Z-file:',ioerr)
      endif
      IF (UPPERC(SCARR(1)).EQ.UPPERC('END ')) GO TO 4
CD     WRITE(6,3) (SCARR(I),I=1,20)
      WRITE(50,3) (SCARR(I),I=1,20)
    3 FORMAT (1X,20A4)
      GO TO 1
C
C  SKIP 2 HEADERS
    4 CONTINUE
      READ(IUZ,2) (SCARR(I),I=1,20)
      READ(IUZ,2) (SCARR(I),I=1,20)
C.....READ START AND STOP TIMES OF TABLE AND CONVERT TO JD
      READ(IUZ,2000) IYR1,IDOY1,IHR1,IMN1,IYR2,IDOY2,IHR2,IMN2,SITECD
      WRITE(iprnt,2001) IYR1,IDOY1,IHR1,IMN1,IYR2,IDOY2,IHR2,IMN2,SITECD
CD     WRITE(ISCRN,2001) IYR1,IDOY1,IHR1,IMN1,IYR2,IDOY2,IHR2,IMN2,SITECD
 2000 FORMAT (I4,1X,I3,2(1X,I2),2X,I4,1X,I3,2(1X,I2),2X,A4)
 2001 FORMAT (1X,I4,1X,I3,2(1X,I2),2X,I4,1X,I3,2(1X,I2),2X,A4)
      IYR1=IYR1-1900
      CALL MONDAY(IDOY1,MONTH,IDAY,IYR1)
      FJDT1 = JULDAY(MONTH,IDAY,IYR1)
      FJDT1 = FJDT1+(DBLE(IHR1)+DBLE(IMN1)/60.D0)/24.D0
      IYR2=IYR2-1900
      CALL MONDAY(IDOY2,MONTH,IDAY,IYR2)
      FJDT2 = JULDAY(MONTH,IDAY,IYR2)
      FJDT2 = FJDT2+(DBLE(IHR2)+DBLE(IMN2)/60.D0)/24.D0
      fjdz1 = fjdt1
      fjdz2 = fjdt2
C
CD     WRITE(ISCRN,864) FJD,FJDT1,FJDT2
  864 FORMAT(' FJD,FJDT1,FJDT2 :',3F15.6)
C
C SKIP TWO RECORDS
      READ(IUZ,2) (SCARR(I),I=1,20)
      READ(IUZ,2) (SCARR(I),I=1,20)
C
C
      RETURN
      END
