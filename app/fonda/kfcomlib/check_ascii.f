CTITLE CHECK_ASCII

      subroutine check_ascii(string)

*     Routine to check a string for non-ascii characters.  Any 
*     non-ascii characters are replaced with the @ character

* PASSED Variables

* string  - The string to checked.  (ichar should be between 32-127
*      for each character in the string.  If non-ascii characters
*      found, all charcaters after 32 are set to blanks.

      character*(*) string

* LOCAL variables

*   i  - loop counter
*   ascii_ichar   -  Value of ascii character (nulls are set to blank)

       integer*4 i, ascii_ichar

*   non_ascii - Set true if non-ascii characters are found

       logical non_ascii

*****  Set non_ascii false to start
       non_ascii = .false.
       do i = 1, len(string)
           ascii_ichar = ichar(string(i:i))
          if( ascii_ichar.eq.0 )  then
              string(i:i) = ' '
              ascii_ichar = 32
           else if ( ascii_ichar.lt.32 .or.
     .               ascii_ichar.ge.128 ) then
              string(i:i) = '@'
              non_ascii = .true.
           end if
       end do

*****  If we had non-ascii characters then clear the end of the
*      string
       if( non_ascii .and. len(string).gt.32 ) then
           string(32:)  = ' '
       end if

*****  Thats all
       return
       end

 

