      program change_ps
 
*     Program to change the ps.tah file.
 
*               line    - Longest line expected.
 
      character*14400 line
      integer*4 iline(3600)

      character*1 ans
 
*   trimlen - Length of string
*   ierr    - IOSTAT error.
 
      integer*4 trimlen, ierr, lenl, nrec, first

      write(*,'(a)') char(49)
 
      open(100, file ='ps.tah', iostat=ierr , access='direct',
     .      recl=14400)
      open(200, file ='ps.th2', access='direct', recl=14400)
 
*     Start reading file
      nrec = 0
      read(100,iostat=ierr, rec=1) iline

      write(*,'('' Words 228,229 '',2i5, 2a)') iline(229), iline(230),
     .        char(iline(229)), char(iline(230))
      write(*,'('' Words 316,317 '',2i5, 2a)') iline(317), iline(318),
     .        char(iline(317)), char(iline(318))

      iline(229) = ichar('5')
      iline(317) = ichar('5')
      iline(230) = ichar('0')
      iline(318) = ichar('0')

      write(*,'('' Words 228,229 '',2i5, 2a)') iline(229), iline(230),
     .        char(iline(229)), char(iline(230))
      write(*,'('' Words 316,317 '',2i5, 2a)') iline(317), iline(318),
     .        char(iline(317)), char(iline(318))

      write(*,'('' Words 2909 is '',i6)') iline(2909) 
      iline(2909) = iline(2909)/2

      write(200, rec=1 ) iline
 
*     Close
      close(100)
      close(200)
      end
 
 
