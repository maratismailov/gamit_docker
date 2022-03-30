      SUBROUTINE RDSITT( SITE, NIN, INCHR, NOUT, OUTCHR, LUN, reqd,ill )

C     subroutine to read site table file
c     S. Shimada 1990, mods by R. King Sep 93

      implicit none

      logical reqd,found_all

      INTEGER  *  4  NIN, NOUT, LUN, ill, ICHRX, ICHR0, ICHR1, I

      CHARACTER*(*)  INCHR, OUTCHR
      character*4 site,stat
      character*80 message
      CHARACTER*200  LINE, HEAD

c   Input:
c       LUN      I*4      Logical unit number of sittbl. file
c       SITE     C*4      Site id to be found
c       NIN      I*4      Number of characters in the table entry question
c       INCHR    CHAR     Table entry question to be read
c       NOUT     I*4      Number of characters in the table entry answer
c       REQD     LOGICAL  .TRUE. if the calling program requires the entry
c                         .FALSE. if it's ok to have the entry missing
c   Output:
c       OUTCHR   CHAR     Table entry answer
c       ILL      I*4      Return code = 0 if ok;
c                                     = -1 if header entry not found, unless
c                                           REQD = .false.
c                                     = -2 if site not found


           
      ill = 0
      ichr0 = 0
      ichr1 = 0
      nout = 0              
      found_all = .false.
      rewind (lun)
      READ( LUN, '(A200)' )  HEAD
      CALL  UPPERS( HEAD )
      ICHRX = INDEX( HEAD, INCHR(1:NIN) )
c      print *,'RDSITT ichrx ',ichrx
      IF( ICHRX .EQ. 0 )  THEN
         if( reqd ) then
           write(message,'(2a)') 'Table term not found: ',inchr(1:nin)
           call report_stat('WARNING','FIXDRV','rdsitt',' ',message,0)
           ill = -1
           return
         else
c          blank for tests on blank or zero
           nout = nin    
           outchr = ' '
           do i=1,nout
              outchr(i:i)= ' '
           enddo
           ill = 0
           return
         endif
      ENDIF
C
      DO  10  I = ICHRX-1, 1, -1
         IF( HEAD(I:I) .EQ. ' ' )  THEN
            ICHR0 = I+1
            GO TO 20
         ENDIF
   10 CONTINUE
   20 CONTINUE
      DO  30  I = ICHRX+NIN, 200
         IF( HEAD(I:I) .EQ. ' ' )  THEN
            ICHR1 = I-1
            GO TO 40
         ENDIF
   30 CONTINUE
   40 CONTINUE
      NOUT = ICHR1 - ICHR0 + 1
C
      STAT = SITE
      CALL  UPPERS( STAT )
   50 CONTINUE
      READ( LUN, '(A200)', END=70 )  LINE
      CALL  UPPERS( LINE )  
      if( line(1:4) .eq. 'ALL ' ) then
         outchr(1:nout) = line(ichr0:ichr1)
         found_all = .true.
      endif
      IF( LINE(1:4) .NE. STAT )  GO TO 50
      OUTCHR(1:NOUT) = LINE(ICHR0:ICHR1)
      ill = 0
      RETURN
           
   70 if( .not.found_all ) then
        write(message,'(2a)') 'Site name not found: ',stat
        call report_stat('WARNING','FIXDRV ','rdsitt',' ',message,0)
        do i=1,nout
          outchr(i:i)= ' '
        enddo
        ill = -2
      endif
      RETURN
      END
