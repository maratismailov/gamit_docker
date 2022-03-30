      SUBROUTINE RDSEST( NIN, INCHR, NOUT, OUTCHR, LUN, REQD, ILL )

C     subroutine to read SESTBL. file
c     S. Shimada 1990, mods by R. King Sep 92, Jun 95, Oct 95

      implicit none

      logical       REQD

      INTEGER*4  NIN, LUN, NOUT, ill, ioerr, ieq, isc, end_of_line
     .         , len, rcpar

      CHARACTER    LINE*120, INCHR*(*), OUTCHR*(*)
     .           , in_string*120, table_string*120
     .           , message*256, prog_name*80

c   Input:
c       LUN      I*4      Logical unit number of SESTBL. file
c       NIN      I*4      Number of characters in the table entry question
c       INCHR    CHAR     Table entry question to be read
c       NOUT     I*4      Number of characters in the table entry answer
c       REQD     LOGICAL  .TRUE. if the calling program requires the entry
c                         .FALSE. if it's ok to have the entry missing
c   Output:
c       OUTCHR   CHAR     Table entry answer
c       ILL      I*4      Return code (=0 if ok; = -2 if format error;
c                            = -1 if entry not found, unless REQD = .false.

c   Rules:  The input string is matched with the table keyword without
c             case sensitivity.
c           The table entry must follow an '=' and one blank.

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)

      REWIND( LUN )
   10 CONTINUE
      READ(LUN,'(A80)',END=30,IOSTAT=IOERR)  LINE
      in_string = inchr(1:nin)
c      print *,'in_string ',in_string(1:nin)
      call uppers(in_string)
c      print *,'in_string ',in_string(1:nin)
      table_string = line(1:nin)
c      print *,'table_string ',table_string(1:nin)
      call uppers(table_string)
c      print *,'table_string ',table_string(1:nin)
c      IF( LINE(1:NIN) .NE. INCHR(1:NIN) )  GO TO 10
      if( in_string(1:nin) .ne. table_string(1:nin) ) goto 10
      ieq = INDEX( LINE(NIN+1:80), '=' )
      isc = index( line(nin+1:80), ';' )
      if( isc.ne.0 ) then
        end_of_line = nin + isc  -1
      else
        end_of_line = nin + ieq + nout + 1
      endif
      outchr = ' '
      if (ieq .NE. 0 )  THEN
         if( (nin+ieq+nout+1).gt.120 )
     .     call report_stat('FATAL',prog_name,'lib/rdsest',' '
     .                      ,'Output line too large',0)
         OUTCHR(1:NOUT) = LINE(NIN+ieq+2:end_of_line)   
c        convert to uppercase except for Unix file name for local I/O
         if( in_string(1:5).ne.'LOCAL' )  CALL  UPPERS( OUTCHR ) 
         ill = 0
      ELSE
         write(message,'(2a)') 'sestbl. format improper: ',inchr(1:nin)
         call report_stat( 'WARNING',prog_name,'lib/rdsest',' '
     .                   , message,0)
         ill = -2
      ENDIF

      RETURN
C
  30  if ( reqd ) then
        write(message,'(2a)') 'Table term not found: ',inchr(1:nin)
c        print *,'message: ',message
        call report_stat('WARNING',prog_name,'lib/rdsest',' ',message,0)
        ill = -1
      else
c       blank for tests on blank or zero
        outchr = ' '
        ill = 0
      endif
      RETURN
      END

