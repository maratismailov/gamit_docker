Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      subroutine NEWCHR(chr)
C
C update the year/iteration code
c
      CHARACTER*1 chr

      if (ichar(chr).ge.ichar('0') .and.
     +    ichar(chr).le.ichar('9')) then
         chr = 'A'
      else if (ichar(chr).ge.ichar('A') .and.
     +        ichar(chr).lt.ichar('Z')) then
         chr = char(ichar(chr)+1)
      else if (ichar(chr).ge.ichar('a') .and.
     +         ichar(chr).lt.ichar('z')) then
         chr = char(ichar(chr)+1)
      else
         chr = 'A'
      endif
C
      return
      END
