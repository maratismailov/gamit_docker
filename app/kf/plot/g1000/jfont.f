CTITLE JFONT
 
      subroutine jfont(type_char)
 
*     This routine saves the font type in the g1000 common.
*     and sent a PRIWTX call to actually set the size of font.
*     The character size as sent is doubled for cartographic
*     because of its small size. Only works with High quality
*     PWRITX text.
 
      include 'g1000.h'

*   type_char - Character string containing the font type or
*               number. 1 - Simplex, 3 - Times, 5 - symbol
*   type_loc  - Local copy of type_char

      character*(*) type_char
      character*20  type_loc
 
*   type    - Font to be set 1 - Cartographics, 2 - Principle
*             size
 
 
      integer*4 type, ierr
 
c     write(*,100)
c 100 format(' Fonts do not need to be set anymore.  High quality'
c    .    ' text is output using',/,
c    .    'PWRITX routine, medium quality text using PWRITY and',
c    .    ' low quality text',/,
c    .    'using PWRIT.',/,
c    .    'To invoke greek letters in a string include ',
c    .    '''PGL''',' to produce ',/,
c    .    '   lower case greek letters at the principal size',/,
c    .    'See NCAR Graphics manual for other cases.')

*     See if numerical type passed.

      read(type_char,'(I2)', iostat=ierr ) type
      if( ierr.ne.0 ) then    ! Try to get string
          type_loc = type_char
          call casefold( type_loc)
          type = 1
          if( type_loc(1:2).eq.'TI' ) type = 2
      end if
 
      jfnt = type
      jfntset = .false.
 
 
      return
      end
 
