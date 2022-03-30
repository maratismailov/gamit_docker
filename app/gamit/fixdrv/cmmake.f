      SUBROUTINE CMMAKE( istat, CFILE, MFILE, SNAMES, BFIL2
     .                 , EXPRMT, ISAT, TOTSAT, LSITE, LSESS 
     .                 , nzen, gradest, ngrad, ierr_sestbl, ierr_sittbl)
C     Subroutine to make CFMRG batch file
C
C     S.Shimada              original at NRCDP
C     S.Shimada   01/26/90   modified at IGPP
C     S.Shimada   02/24/90   reviced to 4-character site codes
C
C     input arguments           
c        istat : number of sites in the session
C        CFILE  : C-file name list
C        MFILE  : M-file name
C        SNAMES : 4-chracter codes for XFILS
C        BFIL2  : secondary batch file name
C        EXPRMT : type of experiment        
C        ISAT   : number of satellites
c        TOTSAT : list of satellites
C
C     sample CFMRG secondary batch file
C     BATCH       : Batch file
C     1ABT        : Station code
C     1B75        : Station code
C                   (blank line)
C     C1ABTD.278  : C-file name
C     C1B75D.278  : C-file name
C     END         : End of stations
C     EEEEEEE     : Explicit biases
C     MTESTA.278  : M-files
C     Y           : Coordinate partials?  ('Y', Unless EXPRMT = 'ORBIT' )
C     Y           : Atmospheric partials? This should be hard wired
c      n n n n .. : Number of zenith delay parameters per station (1st session; others follow)
C     Y           : Orbital partials?  ('Y', Unless EXPRMT = 'BASEL' or 'KINEM')
C     Y           : SV antenna offset partials?  ('Y' if orbit partials included) 
c     Y           : Atmospheric gradient parameters estimated? (Y/N)
c      n n n ..   : Number of gradient parameters per station (1st session; others follows)

      implicit none

      include '../includes/dimpar.h'

      logical reqd,gradest

      integer*4    totsat(maxsat),istat,icnt,isat
     .           , ill,ierr_sestbl,ierr_sittbl,n,i,j,nout
     .           , lsess,lsite,nzen,ngrad,ioerr
                                                                  
      CHARACTER * 1   gradient
      CHARACTER * 4   SITNAM,SNAMES(maxsit),ansnum,anshr
      character * 5   exprmt
      CHARACTER * 16  BFIL2,CFILE(maxsit),MFILE
      CHARACTER * 42  COMMENT,AFMT
      CHARACTER * 120 line
      character * 256 message
   


c  Check for maximum dimensions

      if( istat.gt.maxsit ) then 
        write(message,'(a,i3,a,i3)')
     .    'Network size (ISTAT',istat,' > MAXSIT (=',maxsit      
         call report_stat('FATAL','FIXDRV','cmmake',' ',message,0)
      endif


C     Write the execution line in the main batch file

      WRITE( 17, '(A,A16)' )  'cfmrg  < ', BFIL2


C     Open the secondary batch file and indicate batch (only mode now allowed)

      OPEN( 21, FILE=BFIL2, STATUS='UNKNOWN' )
      WRITE( 21, '(A)' )  'BATCH'

C     Write the site codes, ending with a blank line

      do N = 1, istat
         WRITE( 21, '(A4,16X,A)' )
     .        SNAMES(N), '4 letter site code'
      enddo
      WRITE( 21, '(1X)' )


C     Write the complete satellite list

      WRITE(AFMT,30) ISAT
   30 FORMAT ('(',I2,'I3,5X,A18)')
      comment='Total PRN Numbers '
      WRITE( 21, AFMT) (TOTSAT(ICNT),ICNT=1,ISAT),COMMENT


C     Write C-file names, separating session with blanks and ending with END

      do i = 1, istat
        write( 21, '(a16,4x,a)' ) cfile(i), ' C-file'
      enddo
      write( 21, '(a)' ) 'END'


C     Write whether explicit or implicit biases
C
      DO  N = 1, istat
         LINE(N:N) = 'E'
      enddo
      WRITE( 21, '(A)' )  LINE(1:istat)


C     Write the M-file name

      WRITE( 21, '(A16,4X,A)' )  MFILE, 'M-file'


C     Indicate whether there are station coordinate partials
C
      IF( EXPRMT .EQ. 'ORBIT' )  THEN
         WRITE( 21, '(A1,19X,A)' ) 'N', 'coordinate partials?'
      ELSE
         WRITE( 21, '(A1,19X,A)' ) 'Y', 'coordinate partials?'
      ENDIF


C     Indicate atmosphere partials (This should be hard wired)

      WRITE( 21, '(A1,19X,A)' )
     .     'Y', 'atmospheric partials? Now hard-wired'


c     Write out the number of zenith delay parameters

c     Hard wire the partial flag and number of parameters, even if not
c     planning to estimate zenith delay in SOLVE.
      write(afmt,40) istat
   40 format ('(',i3,'i4,5x,a42,i2)')
      comment='Number zenith delay parameters  '
      write( 21, afmt) (nzen,j=1,istat),comment

C     Indicate orbit partials and SV antenna offset partials

      IF( EXPRMT .EQ. 'BASEL' .OR. EXPRMT .EQ. 'KINEM') THEN
         WRITE( 21, '(A1,19X,A)' ) 'N', 'orbital partials?' 
         WRITE( 21, '(A1,19X,A)' ) 'N', 'SV antenna offset partials?' 
      ELSE
         WRITE( 21, '(A1,19X,A)' ) 'Y', 'orbital partials?'  
         WRITE( 21, '(A1,19X,A)' ) 'Y', 'SV antenna offset partials?' 
      ENDIF

c     Write out the number of gradient parameters
                          
      if( gradest ) then  
        write( 21, '(a1,19X,a)' )
     .     'Y', 'gradient parameters estimated? (Y/N)'   
        write(afmt,40) istat
        comment='Number gradient parameters  '
        write( 21, afmt) (ngrad,j=1,istat),comment
      else
        write( 21, '(a1,19X,a)' )
     .     'N', 'gradient parameters estimated? (Y/N)'
      endif

      CLOSE( 21 )
      RETURN
      END

