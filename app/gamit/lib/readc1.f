      Subroutine readc1 (lunit,ntext,text)

      implicit none

      include '../includes/dimpar.h'

      integer*4          nversn,iflag,ntext
      character*4        buf4
      character*16       cfname
      character*80       text(maxtxt),prog_name
      character*256      message
      integer*4          lunit,ioerr,inqerr,len,rcpar,i

c     Read the first part of a C-file which is record #, version #, and textual history
c     K. Feigl/R. King  July 89; last modified to version 9.30ff by R. King May 95
c     Modfied for release 9.80ff (norb added to record 2, preval assignments changed)
c         R. King Sep 98 
c     Modified for release 10.2 (model parameters added to record 2). R. King Jan 2005 
c      (no non-comment changes to this routine) 
c     Modified for release 10.4 (antenna model descriptors lengthened on record 2 and
c       daily-average atmospheric loading added to record 4)     
c     Modified for release 10.41 (dryzen, wetzen added to record 2)  R. King jan 2013
c     Modified to verson 1042 for addition of eradmod and antradmod on record 2. R. King 27 Mar 2014
c     Modfied to version 1060 for addition of GNSS character in record 2.  R. King 29 August 2014
c     T. Herring Jan 2020: Added L1/L2 satellite PCO. version 10.71.

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)

c     first read only the first two integers, to make sure we have the right version
                     
      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag,nversn                    
c rwk 060905: check for binary compatibility
      if( iflag.ne.1 ) 
     .   call report_stat('FATAL',prog_name,'lib/readc1',' ',
     .'Error reading first record of c-file; binary incompatibility?',0)
     
c rwk 050131: no longer check for ancient (<1995) versions, but rather for 9.80 vs 10.2.
* MOD TAH 200126: Check to see if 1071 version.
      if( nversn.ne.1071 ) then    
        write(message,'(a,i4,a)') 'Old version of C-file (',nversn,')'
        call report_stat('FATAL',prog_name,'lib/readc1',' ',message,0)
      endif
      rewind ( lunit )


c     now read the first three integers, checking that ntext isn't too big for dimensions

      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag,nversn,ntext
      if( ntext.gt.maxtxt ) then
c     Get the calling module name for report_stat
        len =rcpar(0,prog_name)
        write(message,'(a,i4,a,i4)')
     .    'READC1: ntext (',ntext,') > maxtxt (',maxtxt
        call report_stat('FATAL',prog_name,'lib/readc1',' ',message,0)
      endif
      rewind ( lunit )

c     now read the whole record for real

      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag,nversn,ntext,(text(i),i=1,ntext)

 1000 if (ioerr .ne. 0) then
           inquire ( unit=lunit,name=cfname,iostat=inqerr )
           call report_stat('FATAL',prog_name,'lib/readc1',cfname
     .                     ,'Error reading first record of C-file--is it
     . empty from a previous step?',ioerr)
      endif

      if (iflag .ne. 1) then
         write(buf4,'(i4)') iflag
         call report_stat('FATAL',prog_name,'lib/readc1',buf4
     .                   , 'Wrong iflag: ',0)
      endif
c      print *,'READC1: ntext = ',ntext,' text follows: '
c      do 2000 i=1,ntext
c         print *,text(i)
c 2000 continue
      return
      end
