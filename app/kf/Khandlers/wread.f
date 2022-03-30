CTITLE WREAD
 
      subroutine WREAD(dcb,ierr,data,words, length, record)
 

      implicit none
 
*     Routine which acts like VREAD but the number of words to read
*     (WORDS) is integer*4.  This routine mulitply calls VREAD with
*     and I*2 length.  NOTE:  The record to be written is passed
*     optionally as an I*2.  If larger record numbers are required then
*     EPOSN should be used as with VREAD
*     UNIX: Now we can just use readd as before.
 
C    .    data(1j)    ! Data to be read (EMA)
*   data(131072j)    - Data to be read (EMA) (Use large arbirtary
*                   - number
*   dcb(16)     - Dcb buffer to be used
*   ierr        - VREAD error return
*   length      - I*2 length of last segment read (optional)
*   record      - I*2 record number to be read (optional, next
*               - records are read if not given)
*   words       - Number of words to be read
 
      integer*4 data(*), dcb(16),  ierr, length, words, record
 
****  Simply call readd emulation routine
 
 
      call readd(dcb,ierr,data,words, length, record)
 
 
***** Thats all
      return
      end
 
