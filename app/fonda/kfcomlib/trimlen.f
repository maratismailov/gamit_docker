CTITLE 'TRIMLEN'
 
 
      integer*4 function trimlen(string)
 
 
*     Routine to return the length of the used portion of string.
*     Length is with trailing blanks removed.
 
*         len_string    - declared length of string
*         i             - Loop counter
 
      integer*4 len_string, i
 
*       blanks          - Indicates blanks are still being
*                       - found
 
      logical blanks
 
*             string    - the string whose length we want
 
      character*(*) string
 
***** Get full length of string, and work backwards
 
      len_string = LEN(string)
      i = len_string
 
      blanks = .true.
 
*     Scan string from end to front
      do while ( blanks )
          if( i.gt.0 ) then
              if( string(i:i).eq.' ' ) then
                  i = i - 1
              else
                  blanks = .false.
              end if
          else
*                                     ! Force exit at end of string
              blanks = .false.
          end if
      end do
 
***** Save length of string
      trimlen = i
 
***** Thats all
      return
      end
 
 
