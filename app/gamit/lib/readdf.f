      SUBROUTINE READDF( BATCH,LUIN,SNAMES,DAYNUM,istat )
C
C     Main functions :
C      1. Read D-file, assumed open on LUIN.
C      2. Check if all necessary files really exist.
C      3. restore the informations of D-file in arrays.
C      4. reorder stations according to the alphabat order.
C
C     Input :
C        BATCH: Batch or interactive (not currently used in CTOX)
C        LUIN  : D-file unit number
C
C     Output :
C        SNAMES : 4-chracter site codes for XFILS
C        DAYNUM : 3-character day-of-year for XFILS
C
      implicit none

      include '../includes/dimpar.h'

      CHARACTER*1 XORX,BATCH
      character*3 daynum
      CHARACTER*4 SNAMES(maxsit),SNAME,snames1(maxsit)
      CHARACTER*16 XFILS(maxsit),XFILS2(maxsit),XJUNK(maxsit)
     .             ,FNAME,LOWERC,TFILE,XFILE,IFILE,JFILE,LFILE    
      character*80  prog_name
      character*256 message
 

      integer luin,icol,msessn,isessn,istat1,kstat
     .      , jstat,istat,ioerr,len,rcpar,i


c        Get calling program name and m-file name for report_stat

      len = rcpar(0,prog_name)
                       
c        Dummy statement to avoid compiler warning for unused variable
       if( batch.eq.'B' ) then
         continue
       endif
     
         
c        Read the D-file one line at a time till the end

c     first two lines are the GNSS code (not needed) and the # sessions (obsolete)
      read(luin,'(1x)')
      read(luin,'(1x)')

c        Read L, T, I, and J file names but do not use

      READ(LUIN,30,iostat=ioerr,err=900,end=901) LFILE
         CALL  LOWERS( LFILE )
         IF( LFILE(1:1) .ne. 'l' )  THEN
           call report_stat('FATAL',prog_name,'lib/readdf',lfile,
     .      'First entry in D-file should be L-file. Check order',ioerr)
         endif
      READ(LUIN,30,iostat=ioerr,err=900,end=901) TFILE
      READ(LUIN,30,iostat=ioerr,err=900,end=901) IFILE
 30   FORMAT(A16)
      READ(LUIN,'(A16)',iostat=ioerr,err=900,end=901) JFILE
C     Read the number of stations
      READ(LUIN,*,iostat=ioerr,err=900,end=901) ISTAT
         
c     Site loop 
      DO JSTAT=1,ISTAT
   
C       read the x-file nameE
        READ(LUIN,30,iostat=ioerr,err=900,end=901) XFILE
        fname = lowerc(xfile)
        call uppers(xfile)
c       IS IT AN C-, OR X-FILE ?
C       FIND THE PERIOD IN FILENAME
        ICOL=INDEX(XFILE,'.')
C       FIRST LETTER OF FILENAME
        XORX=XFILE(ICOL-6:ICOL-6)
c       Day number of filename
        daynum=xfile(icol+1:icol+3)
        if(xorx.ne.'X'.and.xorx.ne.'C') then
           call report_stat('FATAL',prog_name,'lib/readdf',xfile,
     .       'Bad X- or C-file name on D-file ',0)
        else
          XFILS2(jstat)=XFILE
        endif  
c     endof station loop
      enddo
C
C     Reorder X(C) file according to the alphabet, removing duplicates
      call sort_string(maxsit,istat,xfils2,0,xjunk,istat1,xjunk)
      istat = istat1
      do i=1,istat
        xfils2(i) = xjunk(i)
      enddo

C     Assign the reordered X(C)-file name to array XFILS
      KSTAT=0
      DO JSTAT=1,istat
        KSTAT=KSTAT+1
        XFILS(JSTAT)=XFILS2(KSTAT)
        SNAME=XFILS2(KSTAT)(ICOL-5:ICOL-2)
        snames(kstat) = sname
c     end of site loop
      enddo         
      GOTO 999

c     come here on file error
  900 continue
      if (ioerr.ne.0) then
        call report_stat('FATAL',prog_name,'lib/readdf',' ',
     .  'Error while reading D-file. ',ioerr)
      endif
      return

c     come here on unexpected end-of-file
  901 call report_stat('FATAL',prog_name,'lib/readdf',' ',
     .'Unexpected EOF on D-file ',ioerr)

  999 return
      end
