CTITLE OPEN_LU
 
      subroutine open_lu( unit, name, ierr, status)
 
      implicit none

*     This routine will either open the file NAME to unit UNIT
*     or, if name is a numerical value, it will return the
*     the number in UNIT.
 
*   ierr    - IOSTAT error on open
*   unit    - the unit number (either for opening or to be
*           - set to the number in name
 
      integer*4 ierr, unit
 
*   name    - Name of the file or the unit number
*   status  - Menthod of opening file.
 
      character*(*) name, status
 
* LOCAL VARIABLES
 
*   lu_num  - NUmber of unit read from name
*   jerr    - error in decoding lu_num from name
 
      integer*4 lu_num, jerr

*   Final_status - Status acutally used to open file
*   final_access - Final access for casefolding

      character*20 final_status, final_access
 
****  See if name contains an lu number
      if( name(1:1).ne.'/') then
          call check_num(name, jerr)
          if( jerr.eq.0 ) then
              read(name,*, iostat=jerr) lu_num
          end if 
      else
          jerr = 1
      end if
*                             ! Looks like an lu number
      if( jerr.eq.0 ) then
          unit = lu_num
          ierr = 0
*                             ! Must be a file
      else

****      See if we a appening to file
          final_access = status
          call casefold( final_access )
          if( final_access(1:6).eq.'APPEND' ) then
              final_status = 'unknown'
          else
              final_status = status
              final_access = 'sequential'
          end if
              
          open(unit, file=name, iostat=ierr, status=final_status,
     .               access=final_access)
      end if
 
****  Thats all
      return
      end
 
 
 
CTITLE CHECK_NUM

      subroutine check_num( string, err )

*     Routine to check if a string containes only numeric values.
*     Test each character for [0-9]; lead +- and one decimical point.

* PASSED VARIABLES
* string  - String to be checked
* err     - Error return -- 0 if numberic;
*                           non-zero for non-numeric.

      integer*4 err
      character*(*) string

* LOCAL VARIABLES
* num_dp     - Number of decimal points (only one allowed)
* num_pm     - Number of + - signs
* num_ch     - Number of non-blank lead characters
* num_dp     - Number of deciminal points 
* trimlen    - Function to return Lengh of string.
* len_str    - Length of string 
* i          - Loop counter

* asc_0, asc_9  - Ascii codes for 0 and 9

      integer*4 num_dp, num_pm, num_ch, trimlen, len_str, i
      integer*4 asc_0, asc_9

****  Start; initialize counters
      err = 0

      len_str = trimlen(string)
      if( len_str.eq.0 ) err = -1
      num_dp = 0
      num_pm = 0
      num_ch = 0
      asc_0  = ichar('0')
      asc_9  = ichar('9')

****  Loop over string see if numeric
      i = 0
      do while ( err.eq.0 .and. i.lt.len_str )
         i = i + 1

*        Check the next character
         err = 1
         if( string(i:i).eq.' ' .and. num_ch.eq.0 ) then
             err = 0
         end if 

*        Check for first +- sign
         if( num_pm.eq.0 .and. num_ch.eq.0 ) then
             if( string(i:i).eq.'+' .or. 
     .           string(i:i).eq.'-'     ) then
                 num_pm = num_pm + 1
                 num_ch = num_ch + 1
                 err = 0
             end if
         end if 

*        check to see if digit
         if( ichar(string(i:i)).ge.asc_0 .and.
     .       ichar(string(i:i)).le.asc_9  ) then
             err = 0
             num_ch = num_ch + 1
         end if 

*        check for decimimal point (one allowed)
         if( string(i:i).eq.'.' .and. num_dp.eq.0 ) then
             num_dp = num_dp + 1
             num_ch = num_ch + 1
             err = 0
         end if

      end do

****  Thats all
      return
      end

