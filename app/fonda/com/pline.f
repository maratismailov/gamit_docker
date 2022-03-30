      subroutine pline(ifile,length,char,mode)
c
c     print a line of symbol
c        mode = 1: line of symbol only
c        mode = 2: add an empty line before symbol line
c        mode = 3: add an empty line after symbol line
c        mode = 4: add empty lines both before and after symbol line
c
      integer ifile,length,mode,i
      character*1 char
      character*256 line
c
      if (mode.eq.2.or.mode.eq.4) write (ifile,'(3x)')
      do 10 i = 1,length
         line(i:i) = char
 10   continue
      write (ifile,'(a)') line(1:length)
      if (mode.ge.3) write (ifile,'(3x)')
c
      return
      end
