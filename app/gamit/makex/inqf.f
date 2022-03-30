      subroutine inqf(nam,lexist)
c     inquire to see if a file exists

      include '../includes/makex.h'
      character*(*)  nam
      logical        lexist

      inquire (file   =  nam,
     .         exist  =  lexist)

      return
      end

