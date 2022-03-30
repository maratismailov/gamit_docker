CTITLE WWRIT
 
      subroutine WWRIT(dcb,ierr,data,words, record)

      implicit none
 
 
*     Routine which acts like VWRIT but the number of words to write
*     (WORDS) is integer*4.  This routine mulitply calls VWRIT with
*     and I*2 length.  NOTE:  The record to be written is passed
*     optionally as an I*2.  If larger record numbers are required then
*     EPOSN should be used as with VWRIT
*     UNIX: Just use the writd routine
 
C    .    data(1j)    ! Data to be write (EMA)
*   data(131072j)    - Data to be write (EMA)
*   dcb(16)     - Dcb buffer to be used
*   ierr        - VWRIT error return
*   record      - I*2 record number to be write (optional, next
*               - records are write if not given)
*   words       - Number of words to be write
 
      integer*4 data(*), dcb(16),  ierr, record, words
 
****  Simple use the writd routine
 
      Call writd(dcb,ierr,data,words, record)
 
***** Thats all
      return
      end
