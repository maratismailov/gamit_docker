CTITLE CORRUPT_KALOBS
 
      subroutine corrupt_KalObs(type,error,record,action,abort)

      implicit none
 
*     J.L. Davis                   3:27 PM  TUE.,  2  JUNE, 1987
*
*     This routine will tell the user (LU=1) that something bad
*     has happened while reading or writing the KalObs file.
*     The following actions occur:
*
*     0. Checks is error=0.  If so, return.
*
*     1. An error message is written saying that an error has
*        occurred.
*
*     2. The KalObs file is closed.
*
*     3. REPORT_ERROR is called with the terminate flag set
 
 
*       abort               - 1=abort, 0=don't
*   ,   close_error         - Error flag on close of KalObs
*   ,   error               - Error flag passed by user
*   ,   record              - The record number in KalObs
 
      integer*4 abort, close_error, error, record
 
*       action              - What happened (e.g. 'read')
*   ,   type                - Type of error (e.g. 'IOSTAT')
 
      character*(*) action, type
 
*       prog                - Routine name
 
      character*14 prog
 
      data prog / 'CORRUPT_KALOBS' /
 
***** Return if no error occurred
      if (error .eq. 0) return
 
***** Write a message
      write(*,100) record
  100 format(/,' ERROR: occurred at record #',I3,' of KalObs',
     .       /,' **** WARNING **** KalObs may be corrupt')
 
***** Try to close KalObs
      call close_KalObs(close_error)
 
***** Check for error on close
      call report_error('FMP',close_error,'clos','KalObs',0,prog)
 
***** Report origninal error
      call report_error(type,error,action,'KalObs',abort,prog)
 
      end
 
