CTITLE    ...................................................................
 
      subroutine trim_lead(buffer,ierr)
c
c     Routine to remove leading blanks from string.
c     See description of ierr returns
c
c Variables
c ---------
c buffer -- buffer to have its leading blanks removed
c ierr   -- a returned error number:
c           ierr = 0 if OK
c           ierr = -2 if string is all blanks
c           ierr = -3 if string is too long to be copied
c
      character*(*) buffer
 
c
      integer*4 ierr
 
c
c new_buffer -- dummy buffer used to copy buffer needed becuase of bug
c     in fortran complier
c ilen  -- length of buffer
c ib    -- counter used to find first non-blank character
c
      character*80 new_buffer
 
c
      integer*4 ilen, ib
 
c
c scratch common
c --------------
c assign local variables to scratch common
c
      common new_buffer
 
c
c Functions used
c --------------
c Trimlen -- HP utility
c
      integer*4 trimlen
 
c
c.... Firstly make sure buffer is not too long to be copied
c     through 'new_buffer'
      ierr = 0
      ilen = trimlen(buffer)
c
      if( ilen.eq.0 )then
*                              ! set string empty error
         ierr = -2
         return
*                              ! see if too long
      else
*                              ! string too long
         if( ilen.gt.80 ) then
            ierr = -3
            return
         end if
      end if
c
c.... String is fine now find leading blanks
      ib = 1
*                                          ! check for blank
      do while ( buffer(ib:ib).eq.' ' )
         ib = ib + 1
      end do
c
c.... Now copy string back into itself. Due to compiler bug we must
c     use the intermediate buffer new_buffer
*                           ! kill off leading blanks
      if( ierr.ge.0 ) then
         new_buffer = buffer(ib:ilen)
         buffer     = new_buffer
      end if
c
      return
      end
 
