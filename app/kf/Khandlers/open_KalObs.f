CTITLE OPEN_KALOBS
 
      subroutine open_KalObs(KO_name,ierr)

      implicit none
 
c     Routine to open KalObs
C MOD TAH 870128 Modified to use FmpOpen routine rather than OPEN call
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
 
 
      character*(*) KO_name
 

*   ierr            - FMP error flag

      integer*4 ierr 
 
***** Open file, first append appeend file type if not already there
      call FmpOpen(ko_dcb, ierr, KO_name, 'RWO', 1)
 
      call report_error('FmpOpen', ierr,'open',KO_name,0,'OPEN_KALOBS')
 
 
***** Now use FmpShortName to get the name of the file.  In particular
*     this will return the cartridge on which this is located
C MOD TAH 891117:  Following code is not needed on a UNIX system...
C     although we could use something like this to get path name
C     and the function GETCWD(dirname) to see where we are.
 
C     if( ierr.ge.0 ) then
C         ierr = FmpShortName( ko_dcb, ierr, KO_name )
C         call report_error('FmpShortName',ierr,'nam',KO_name,0,
C    .                      'OPEN_KALOBS')
C     end if
 
***** Set value_read false
 
*                             ! Force the first record of file to be
      values_read = .false.
*                             ! read, so that we will get the record
*                             ! lengths and start record numbers in
*                             ! RW_KalObs_Block
      end
 
