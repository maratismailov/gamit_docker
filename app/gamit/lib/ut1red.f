      SUBROUTINE UT1RED( IUT1,JD,FRACT,TAIUT1,UT1DOT,IUTTYP )
C
C     Read UT1 values from an external data set
C     Output variable is TAI-UT1
C     R. Abbot - Nov 1984 from PEP routines
C     R. King -  Jun 1987
C
      implicit none
C
      character*32   afmt
      INTEGER*4 JD,MJD,JD1,JD2,JDT1,JDT2,J00001,iutval(12),itsav
     .        , iut1,int,nvtot,nr,nrec,nrecs,ifirst,npr,nvr
     .        , ioerr,itype,iuttyp,it,nback,nrbmin,lap,i,j,n

      REAL*8 TAB(4),UTVAL(12),Y1(2),Y2(2),taiut1,ut1dot,s,t,f2,fract
     .     , rmult

      integer*4 len,rcpar

      character*80 prog_name
      character*256 message


      SAVE UTVAL,JDT1,JDT2,JD1,JD2,Y1,Y2,
     $     ITSAV,NVTOT,NVR,NREC,IFIRST,NPR,INT,rmult,AFMT,itype

      DATA J00001/2400001/,jd1/0/
      nrecs=0
      nback=0
      nrbmin=0
      lap=0

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)

c     handle problem when JD skips backwards, as in multiple ARC integrations
      if( jd.lt.jd1 ) ifirst = 0

      IF (IFIRST.NE.0) GO TO 130
C
C     Read the Header Records, assuming format goodies
c     on the second line.
C 
      REWIND IUT1
c     READ (IUT1,50) afmt,JDT1,JDT2,NPR,INT
      READ (IUT1,50) afmt,itype,JDT1,JDT2,NPR,INT,rmult
c  50 FORMAT (/,a35,I7,1X,I7,2X,I1,2X,I1)
   50 FORMAT (/,a32,1x,i1,1x,i7,1X,I7,2X,I1,2X,I1,10x,g7.0)

c     find the end of the format format
      i = index(afmt,' ')
      afmt = afmt(1:i-1)
c     choke if you can't find it
      if (i.lt.1) goto 540

      IFIRST = 1
      NVTOT = 0
      NREC=0
      nrbmin= 12/npr + 1
      JDT2 = JDT2 - INT
      ITSAV=-9999

C     Read data record into storage
C     Keep reading until (12/NPR) records are in storage
c     read them as integers, convert them to reals

  100 IF (NVTOT.GT.12-NPR) GO TO 120
      READ (IUT1,fmt=afmt,ERR=520,iostat=ioerr)
     .MJD,(IUTVAL(I),I=1,NPR),NVR
      do 105 i = 1,npr
         utval(nvtot+i) = dble(iutval(i))
  105 continue
      NREC = NREC + 1

C     JD1 and JD2 are the limits of usable values in storage
C     and depend on the interpolation scheme
      IF (NVTOT.EQ.0) JD1 = MJD + J00001 + INT
      IF (NVR.EQ.0) NVR = NPR
      NVTOT = NVTOT + NVR
      ITSAV=-9999
      JD2 = MJD + J00001 + (NVR-2)*INT

      IF (JD2.LT.JDT2) GO TO 100
C        Save the first date just in case the header JDT1 is wrong
  120 IF (IFIRST.EQ.0) JDT1 = JD1
C
C     Is JD within the range of the table?
C
  130 IF (JD.LT.JDT1.OR.JD.GT.JDT2) GO TO 500
C     Is JD too close to JD1 ?
      IF (JD.GE.JD1) GO TO 170
C        If so, backspace and try again
      if( nback.gt.0 .and. nrec.ge.nrecs ) lap = nrbmin+1+nrec-nrecs
      nrecs = nrec
      nback = nback + 1
      if( nback.gt.100) then
c       get calling program name and m-file name for report_stat
        len = rcpar(0,prog_name)
        write(message,'(a,i8,a,2i8)') 'Infinite loop in UT1RED, jd='
     .       ,jd,' jd1,jd2=',jd1,jd2
        call report_stat('FATAL',prog_name,'lib/ut1red',' ',message,0)
      endif
      N = (JD1-JD)/INT/NPR + lap
      IF (N.GT.NREC) N = NREC
      DO 160 I=1,N
  160 BACKSPACE IUT1
      NREC = NREC - N
      NVTOT = 0
      GO TO 100
C
C        Is JD too close to JD2 ?
  170 IF (JD.LT.JD2) GO TO 200
C        If so, shift storage and read another record
      NVTOT = NVTOT - NVR
      DO 180 I=1,NVTOT
  180 UTVAL(I) = UTVAL(NPR+I)
      JD1 = JD1 + NPR*INT
      GO TO 100
C
C        Calculate interpolation times and value of tabular points
  200 T = JD - JD1
      T = (T + FRACT)/INT
      IT = T
      T = DMOD(T,1.D0)
      S = 1.D0 - T
      IF (IT.EQ.ITSAV) GO TO 230
      DO 210 I=1,4
      J = IT + I
  210 TAB(I) = UTVAL(J)
C
C
C        Calculate interpolation Y-vector
C
      DO 220 I=1,2
      NR = I + 1
      F2 =     0.166666666666667D0 * (TAB(NR+1)+TAB(NR-1))
      Y1(I) =  1.333333333333333D0 * TAB(NR) - F2
  220 Y2(I) = -0.333333333333333D0 * TAB(NR) + F2
      ITSAV = IT
C
C
C        Second difference interpolation
C
c 230 TAIUT1 = (T* (Y1(2)+T*T*Y2(2)) + S * (Y1(1)+S*S*Y2(1)))*1.D-5
  230 TAIUT1 = (T* (Y1(2)+T*T*Y2(2)) + S * (Y1(1)+S*S*Y2(1)))*rmult
      ut1dot = ( y1(2)+3.d0*t*t*y2(2) - y1(1)-3.d0*s*s*y2(1) )*rmult/int
      GO TO 999
C
C
C        Out of range of table
C
  500 WRITE(message,510) JD,JDT1,JDT2
  510 FORMAT('JD= ',I7,' out of range of ut1 table, JDT1= '
     .,I7,' JDT2= ',I7)
      call report_stat('FATAL',prog_name,'lib/ut1red',' ',message,0)
  520 WRITE(message,530) MJD,JD2,JDT2
  530 FORMAT('File error in UT1 table, MJD = '
     .,I5,' JD2 = ',I7,' JDT2 = ',I7)
      call report_stat('FATAL',prog_name,'lib/ut1red',' ',message,ioerr)
  540 call report_stat('FATAL',prog_name,'lib/ut1red',' '
     .,'Cannot find format statement in ut1 file:',0)
  999 iuttyp = itype

      RETURN
      END
