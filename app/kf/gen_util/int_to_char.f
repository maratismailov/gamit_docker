ctitle
 
      subroutine int_to_char(ch_string, ival, num_char)

      implicit none 
c
c     routine to convert the integer 'ival' into a character
c     string 'ch_string'.
c
c Variables
c ---------
c ch_string  -- character*2 string which will have 'ival' put into
c     it
c ival   -- integer value which has two characters stored in it.
c num_char -- the number of charaters to converted
c
c Local variables
c ---------------
c jval   -- the high order byte of ival
c kval   -- the low order byte of ival
c jchar  -- starting character number in word
c num_wrd -- number of words to be converted
c
      character*(*) ch_string
 
c
      integer*4 ival(1), num_char
 
c
*   i,j,k       - Loop counters
*   iel         - Set to 1 or 2 for first or second byte
*   kval        - Low order byte
	
      integer*4 jval, jchar, num_wrd, i, iel, kval
 
c
c.... find out number of words to be converted
      num_wrd = (num_char-1)/2 + 1
c
c.... check to see if too many words
      if( num_wrd.gt. len(ch_string) ) then
         write( *  ,100) num_wrd, len(ch_string), (ival(i),i=1,num_wrd)
  100    format(/" ** warning ** ",i3," words attempted to be ",
     .      ' converted to a string with ',i3,' characters',/,
     .      ' Offending string was ',20a2)
         num_wrd = (len(ch_string)-1)/2 + 1
      end if
c
c.... split each word of ival into two characters
      do i = 1,num_wrd
c
c....    compute beginning character number
         jchar = 2*(i-1) + 1
c
c....    see if only one character in last word
         iel = 2
         if( i.eq.num_wrd ) then
            if( num_wrd*2.ne.num_char ) then
               iel = 1
*                 ! even number of characters
            else
               iel = 2
            end if
         end if
c
c....    compute the ascii codes of the two characters in ival(i)
         jval  = ival(i)/256
         kval  = ival(i) - jval*256
c
c....    concatinate the two characters (char is a FTN77 function, //
c        concatinates strings).
         if( iel.eq. 2 ) then
            ch_string(jchar:) = char( jval ) // char( kval )
 
*             ! only one character needed
         else
            ch_string(jchar:) = char( jval)
         end if
c
      end do
c
c.... thats all
      return
      end
 
c.........................................................................
