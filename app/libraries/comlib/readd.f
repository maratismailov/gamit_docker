ctitle readd
 
      subroutine readd(dcb,ierr,ibuf, numr, lenread, recnum)

      implicit none
 
*     Routine to emulate the HP 1000 readd subroutine.
*     The input are:
*     DCB - DCB buffer as defined in FmpOpen
*     ierr - error return:
*         0 -- all OK
*         -1 -- EOF reached, lenread returns 0
*         -2 -- File not open
*     ibuf- integer*4 array for return
*     numr - number of I*4 words to read
*     lenread - actual number of words read
*     recnum  - number of record to read. 0 for next record
 
*   dcb(16)     - DCB buffer
*   ierr        - Fmp error
*   ibuf(*)     - Buffer for read
*   numr        - Number of words to read
*   lenread     - number read
*   recnum      - Record number to read
 
      integer*4 dcb(16), ierr, ibuf(*), numr, lenread, recnum
 
*     LOCAL VARIABLE
 
*   actrecn     - actual record number to read
*   i,j         - counter for reading records
*   numw        - actual number of records to read
 
      integer*4 actrecn, i, j , numw
 
****  MAke sure file is open
      if( dcb(1).lt.700 ) then
          ierr = -2
          lenread = 0
          return
      end if
 
****  Get record to read
      if( recnum.eq.0 ) then
          actrecn = dcb(3) + 1
      else
C         actrecn = (recnum-1)*dcb(2) + 1
          actrecn = recnum
      end if
 
****  Now do the read
      if( dcb(5).le.2 ) then
          numw = (numr-1)/dcb(2) + 1
          do j = 1, numw
             read(dcb(1), iostat=ierr, rec=actrecn+j-1 ) 
     .           (ibuf((j-1)*dcb(2)+i),i=1,dcb(2))
          end do
          actrecn = actrecn+numw-1
      else
          read(dcb(1), iostat=ierr ) (ibuf(i),i=1,numr)
      end if
 
*                             ! Update DCB
      if ( ierr.eq.0 .or. ierr.eq.110 ) then
          dcb(3) = actrecn
          lenread = numr
      else if ( ierr.gt.0 ) then
          ierr = -ierr
      end if
 
      if ( ierr.eq.-1 ) then
          lenread = 0
      end if
 
****  Thats all
      return
      end
