      SUBROUTINE WHDRED( IUW,IPRNT,FJDW1,FJDW2,H2OTYP,SEALOC )
C
C     Read the weather (W-) file header record
C
      implicit none
C
      CHARACTER*1 SEALOC,H2OTYP,RELHUMID,DEWPOINT,SEA,LOC
      CHARACTER*3 CHKEND,BUF3
      CHARACTER*8 PARAM(3),LEVEL(3),UNITS(3)
      CHARACTER*80 HEADER,fn
C
      INTEGER*4 JULDAY,imnth,ihr1,ihr2,imn1,imn2,iyr1,iyr2,iday
     .        , idoy1,idoy2,iuw,iprnt,i,ioerr

      real*8 fjdw1,fjdw2

      DATA CHKEND/'END'/
      DATA RELHUMID,DEWPOINT,SEA,LOC/'R','D','S','L'/
      call uppers(chkend)
      call uppers(dewpoint)
      call uppers(relhumid)
      call uppers(sea)
      call uppers(loc)
      ioerr = 0

      inquire ( unit=IUW, name=fn, iostat=ioerr )
      if (ioerr .ne. 0 ) 
     .   call report_stat('FATAL','MODEL','whdred',fn,
     .   'W-file does not exit or is not opened:',ioerr)

C        Read and echo to the P-file the header comments   

      write(iprnt,'(a)') ' W-file header:'
10    READ (IUW,15,iostat=ioerr) HEADER
15    FORMAT (A80)
      if (ioerr .ne. 0 ) 
     .   call report_stat('FATAL','MODEL','whdred',fn,
     .   'Error reading weather file header:',ioerr)
      BUF3 = HEADER(1:3) 
      call uppers(buf3)
      IF (BUF3.EQ.CHKEND) GO TO 30   
      write(iprnt,'(1x,a80)') header
      GO TO 10

C        Read the start and stop times of the table and convert to JD

30    READ(IUW,*,iostat=ioerr) IYR1,IDOY1,IHR1,IMN1,IYR2,IDOY2,IHR2,IMN2 
      if (ioerr .ne. 0 ) 
     .   call report_stat('FATAL','MODEL','whdred',fn,
     .   'Error reading weather file times',ioerr)
      call check_y2k(iyr1)  
      WRITE (IPRNT,45) ' Start/Stop times  '
     .     , IYR1,IDOY1,IHR1,IMN1,IYR2,IDOY2,IHR2,IMN2
45    FORMAT (a,I4,1X,I3,2I3,1X,I4,1X,I3,2I3)
      CALL MONDAY(IDOY1,IMNTH,IDAY,IYR1)
      FJDW1 = JULDAY(IMNTH,IDAY,IYR1)    
      FJDW1 = FJDW1+(DBLE(IHR1)+DBLE(IMN1)/60.D0)/24.D0
      CALL MONDAY(IDOY2,IMNTH,IDAY,IYR2)
      FJDW2 = JULDAY(IMNTH,IDAY,IYR2)
      FJDW2 = FJDW2+(DBLE(IHR2)+DBLE(IMN2)/60.D0)/24.D0
C
C
C        Read the weather parameter type, sea level or local, and units
C
      READ (IUW,50,iostat=ioerr) (PARAM(I),I=1,3)   
      if (ioerr .ne. 0 ) 
     .   call report_stat('FATAL','MODEL','whdred',fn,
     .   'Error reading weather file parameters',ioerr)
      call uppers(param(3))
50    FORMAT (16X,A8,1X,A8,3X,A8)
      IF (PARAM(3).EQ.'RELHUMID') H2OTYP=RELHUMID
      IF (PARAM(3).EQ.'DEWPOINT') H2OTYP=DEWPOINT
      READ (IUW,50,iostat=ioerr ) (LEVEL(I),I=1,3) 
      if (ioerr .ne. 0 ) 
     .    call report_stat('FATAL','MODEL','whdred',fn,
     .   'Error reading weather file level',ioerr)
      call uppers(level(1))
      READ (IUW,60,iostat=ioerr ) (UNITS(I),I=1,3) 
      if (ioerr .ne. 0 ) 
     .   call report_stat('FATAL','MODEL','whdred',fn,
     .   'Error reading weather file units',ioerr)
60    FORMAT (16X,A8,1X,A8,3X,A8)
      IF (LEVEL(1).EQ.'SEA     ') SEALOC=SEA
      IF (LEVEL(1).EQ.'LOCAL   ') SEALOC=LOC  
      write(iprnt,'(a,3(1x,a8))') ' Type of measuerments: '
     .    ,(param(i),i=1,3)
      write(iprnt,'(23x,3(1x,a8))') (level(i),i=1,3)
      write(iprnt,'(a,3(1x,a8))') ' Units               : '
     .    ,(units(i),i=1,3)

      RETURN
      END
