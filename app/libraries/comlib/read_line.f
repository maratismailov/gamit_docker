CTITLE 'READ_LINE'
 
      subroutine read_line(line,indx,type,err,value,cvalue)

      implicit none
 
c
c     Routine to read the next 'thing' in string LINE, starting
c     with character indx.  The 'thing' is returned in VALUE,
c     and the type of VALUE is given by the user in TYPE:
c
c              TYPE             type
c              ----             ----
c              'CH'             character
c              'I2'             integer*2
c              'I4'             integer*4
c              'R4'             real*4
c              'R8'             real*8
c              'R16'            real*16
c
c     The IOSTAT error is returned in ERR.  indx is updated to point
c     at the next blank character.  If type is 'CH', the result is
c     returned in CVALUE, otherwise in VALUE.
c
c
c     J.L. Davis 870729 Transported to Vax from HP.  Added REAL*16
c     feature.
c     T.A. Herring 890301 Removed I*2 references (except for value)
c     so that routine is more generic Fortran.
c     R. W. King 970111 Dimensioned value to 8 rather than 1 to allow
c     bounds checking
 
      character*(*) line, type, cvalue
 
c
*         value(1)   - Passed as shortest length so that all
*                    - will match.
 
      integer*2 value(8)
 
      integer*4 indx, err, length, indx_end, num_words
 
      integer*4 i
 
c
      integer*2 i2_value(8)
 
      integer*4 i4_value
 
      real*4 r4_value
 
      real*8 r8_value
 
c     real*16   r16_value ! Real*16 not supported on the convex.  Subs
c                         ! Real*8 as fix.
 
      real*8 r16_value
 
*       option_found  - Indicates that the type passed was valid
*                     - Input value is left unchanged.
      logical option_found
 
c
      equivalence (i2_value, i4_value)
      equivalence (i2_value, r4_value)
      equivalence (i2_value, r16_value)
      equivalence (i2_value, r8_value)
c
c.... Initialize results
      num_words = 0
*                      ! Default no error
      err       = 0
c
c.... Determine length of string
      length = len(line)
      if( indx.le.0 ) indx = 1
* MOD TAH 951026: Check to make sure that indx does not get too long
      if( indx.gt.length ) indx = length
c
c.... Find next nonblank character
      do while (line(indx:indx) .eq. ' ' .and. indx .lt. length)
c
c....   Increment indx
        indx = indx + 1
c
      end do
c
c.... Did we reach the end of the line?  At this point we do not know if
c     we exited the above line because we were are the end of the line
c     or because the next non-blank character was the last character in
c     the line.  Check now by seeing if line(indx:indx) is blank.  If it
c     we have reached end of line so set ERR=-1 (EOF)
 
*                                         ! We did reach end of line
      if ( line(indx:indx).eq.' ' .and. indx.eq.length ) then
          err = -1
* MOD TAH 921102: We are not going to do anything so say we found the
*         option, so that an error message is not generated.
* MOD TAH 951026: Set the end index value and added check on indx itself.
          option_found = .true.
          indx_end = indx
*                                         ! We still have characters to
      else
*                                         ! read
 
c
c....   Find the next blank character
        indx_end = indx
        do while (line(indx_end:indx_end) .ne. ' '
     .            .and. indx_end .lt. length)
c
c....     Increment end indx
          indx_end = indx_end + 1
c
        end do
 
        option_found = .false.
        call casefold( type )
c
c....   Do we want character data?
        if (type .eq. 'CH') then
c
c....     Simply assign character substrings
          cvalue     = line(indx:indx_end)
          option_found = .true.
c
        end if
c
c....   Integer*2?
        if (type .eq. 'I2') then
c
c....     Read into I2 variable
          i2_value(1) = 0
          read (line(indx:indx_end),*,iostat=err) i2_value(1)
          if (err .eq. 0) num_words = 1
          option_found = .true.
c
        end if
c
c....   Integer*4?
        if (type .eq. 'I4') then
c
c....     Read into I4 variable
          i4_value = 0
          read(line(indx:indx_end),*,iostat=err) i4_value
          if (err .eq. 0) num_words = 2
          option_found = .true.
c
        end if
c
c....   Real*4?
        if (type .eq. 'R4') then
c
c....     Read into R4 variable
          r4_value = 0
          read(line(indx:indx_end),*,iostat=err) r4_value
          if (err .eq. 0) num_words = 2
          option_found = .true.
c
       end if
c
c....   Real*16?
        if (type .eq. 'R6') then
c
c....     Read into R16 variable
          r16_value = 0
          read(line(indx:indx_end),*,iostat=err) r16_value
          if (err .eq. 0) num_words = 8
          option_found = .true.
c
        end if
c
c....   Real*8?
        if (type .eq. 'R8') then
c
c....     Read into R8 variable
          r8_value = 0
          read(line(indx:indx_end),*,iostat=err) r8_value
          if (err .eq. 0) num_words = 4
          option_found = .true.
c
        end if
 
      end if
 
c.... Transfer to VALUE, if necessary
      if (num_words .gt. 0) then
 
          do i = 1, num_words
              value(i) = i2_value(i)
          end do
 
      end if
c
c.... Update indx
      indx = indx_end
 
c.... Make sure we founf option
      if( .not.option_found ) then
          call bad_option( type, 'READ_LINE')
          err = 1
      end if
c
      END
