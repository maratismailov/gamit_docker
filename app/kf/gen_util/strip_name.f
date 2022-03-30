ctitle
 
      subroutine strip_name(buffer,names,isize,iel,decode)

      implicit none 
 
c
c     routine to find the occurence of the first 8 non-blank
c     characters in buffer in the list names. If no match
c     is found then iel is returned equal to -1. If decode is non-zero
c     then an error message will be printed if no
c     match is found.
c     If bit 16 is set in decode the string 'all' will decoded
c     If 'all' is passed it will be interpretted
c     then as being applied to all vals associated with this quantity.
c     In this case iel is returned with value 999999.
c
c     Buffer is returned with the matching characters stripped
c     from it so that it can be further decoded if necessary.
c     If the buffer is empty (i.e., all blanks), iel is returned
c     as -2.
c
c Variables
c ---------
c buffer  -- an <=256 character buffer which will have its first 8
c        characters matched to names
c names   -- the list of names which will be searched
c isize   -- the number of name strings in names
c iel     -- the number of the string in names which matched the
c        characters in buffer.  Returned -1 if no match found.
c        Returned -2 is string is all blanks
c decode  -- bit 1 set will cause an error message to be printed when
c        no match is found between buffer and names.
c            bit 2 will allow the string 'ALL' to be decoded
c Local variables
c ---------------
c new_buffer -- 256 character buffer used to remove the matching
c        string from buffer (needed because of bug in ftn77)
c ib      -- integer variable used to count blanks in buffer
c ibmax   -- length of the buffer array.
c actual_len -- the actual length of the string (trailing blanks not
c        not counted)
c kbit    -- SOLVE logical function to tell if bit is on
 
      character*(*) buffer
 
      character*256 new_buffer
c
      character*8 names(*)
 
c
      integer*4 isize, iel, decode, ib, ibmax, actual_len
 
*   i,j,k   - Loop counters
      integer*4 i, cand
 
c
C     logical kbit
 
c
c.... Remove any leading blanks in buffer
 
      ibmax = len(buffer)
      ib = 1
      do while ( buffer(ib:ib).eq.' ' )
         ib = ib + 1
c
c....    see if we have reached the end of the string
*                               ! string is empty
         if( ib.gt.ibmax ) then
            iel = -2
            return
         end if
c
      end do
c
c.... now try to match the first 8 characters of buffer with strings
c     in name
      iel = -1
      call casefold(buffer(ib:ib+7))
      do i = 1,isize
         if( index(names(i),buffer(ib:ib+7)).ne.0 ) iel = i
      end do
c
c.... see if we can match 'all' to the string if bit 2 of decode is on
c mod simon 2/22/2002 (layhey compiler compat problem) convert oct 100000 to integer 32768.
c      if( cand(decode,o'100000').ne.0 ) then
      if( cand(decode,32768).ne.0 ) then
         if( buffer(ib:ib+2).eq.'ALL' ) iel = 999999
      end if
c
c.... see if we found a match
*                           ! no match found
      if( iel.eq.-1 ) then
c
c....    see if user wants to find out that there was no match
c mod simon 2/22/2002 (layhey compiler compat problem) convert oct 77777 to integer 32767.
c         if( cand(decode,o'77777').ne.0 ) then
         if( cand(decode,32767).ne.0 ) then
c
c....       get the length of 'actual' length of the string
            actual_len = ibmax+1
            do  i = ibmax, ib+7, -1
               if( buffer(i:i).eq.' ' .and. actual_len.eq.i+1  )
     .            actual_len = i
            end do
c            
            write(6,'(3a)') 'Unable to match'
     .            , buffer(ib:actual_len),' : Ignoring this line'
         end if
c
*            ! a match found so strip name from string
      else
c
c....    now strip out the characters which have been matched.
         new_buffer = buffer(ib+8:)
         buffer = new_buffer
c
      end if
c
      return
      end
 
c.........................................................................
