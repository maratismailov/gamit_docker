      character*80 function unpad (filnam, suffix)
c
c     This function joins strings together for file name use so that
c     the blank characters in the string are elliminated.
c     This allows suffix cancatennation.
c
c     by peter morgan for the Apollo, January 1987.
c
c     The routine assumes that blanks are NOT part of file name strings

      character*(*) suffix
      character*(*) filnam
      character*2 two_blanks
      integer*4   no,no2

      data two_blanks/'  '/

      no=index(filnam,two_blanks)
      no=no-1
      no2 = len(filnam)
      no = min(no,no2,80)
      unpad=filnam(1:no)//suffix

      return
      end
