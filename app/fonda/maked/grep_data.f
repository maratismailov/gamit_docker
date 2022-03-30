      subroutine grep_data(fil1,in_opt,mode)
c
c     grep multiple input data files to get 
c     compressed data file
c
      implicit real*8(a-h,o-z)

      character*80 subfile
      character*80 in_opt
      character*80 command
      character*3 fmt2
      integer fil1,mode,ifile,i,len,ierr 
      integer nblen

      print*,' Begining to extract baseline vectors from o-files ....'
      mode = 25
c     new ofile format
c      if (in_opt(4:4).eq.'N'.or.in_opt(4:4).eq.'n') mode = 27
      if (in_opt(3:3).eq.'N'.or.in_opt(3:3).eq.'n') mode = 27
c     remove underscore from in_opt
      if (in_opt(2:2).eq.'_') then
          in_opt(2:2)=' '
      end if

      ifile = 0
      read (in_opt,'(a3)') fmt2
c     fmt2(6:6) = ' '
      if (nblen(fmt2).eq.0) then
         print*,'missing in_opt file'
         goto 1000
      end if
 
      rewind (fil1)
c     get input file name list
      do 20 i = 1,200
         read (fil1,'(a)',err=100,end=1000) subfile
         ifile = ifile+1
         call blank(command)
c        temporary file collect all necessary lines from o-file
         command(1:8) = 'grep -e '
         command(9:13) = ''''//fmt2//''''
         len = nblen(subfile)
         command(14:14+len+1) = ' '//subfile(1:len)
         if (ifile.eq.1) command(15+len+1:26+len+1) = ' > tmpofile'
         if (ifile.gt.1) command(15+len+1:26+len+1) = ' >> tmpofile'
         ierr = system(command)
 20   continue
      goto 1000

 100  print *,'GREP_DATA: error in reading data file: ',subfile
      stop 'MAKED: GREP_DATA ...'
         
 1000 continue
      return
      end
