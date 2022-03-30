      subroutine merge_list(fil1,fil2,mode)
c
c     merge several maked.out file to create a list file
c                                Danan Dong 07/26/95
c
      implicit real*8(a-h,o-z)

      character*80 subfile
      character*8  st1,sitnam(500)
      character*80 line,note
      integer fil1,mode,fil2,ifile,isit,nsit
      integer i,i1,i2,j
      integer match_name,nblen

      print*,' Begining to merge file list ......'
      mode = 25
      ifile = 0
      nsit = 0
      rewind (fil1)
c     get input file name list
      do 20 i = 1,200
         read (fil1,'(a)',err=100,end=70) subfile
         if (subfile(1:1).eq.'*') goto 20
         i1 = nblen(subfile)
         if (i1.lt.1) goto 20
         ifile = ifile+1
         open (25,file=subfile(1:i1),status='old',err=110)
c        get site number
         read (25,'(a)') line
         read (line,*) isit
         if (isit.le.0) goto 10
         do 30 j = 1,isit
            read (25,'(a)') st1
            if (ifile.eq.1) then
               nsit = nsit+1
               sitnam(nsit) = st1
            else
               i2 = match_name(nsit,8,sitnam,st1)
c              new site, not match current name list
               if (i2.le.0) then
                  nsit = nsit+1
                  sitnam(nsit) = st1
               endif
            endif
 30      continue
 10      close (25)
 20   continue
c
 70   print*,' total sites =',nsit
      if (nsit.le.0) goto 1000
c
c     create a list file
      rewind (fil1)
      note(1:16) = '{number of sites'
      write (fil2,'(i5,50x,a16)') nsit,note(1:16)
      do i = 1,nsit
         write (fil2,'(a8)') sitnam(i)
      enddo
      note(1:19) = '{number of subfiles'
      write (fil2,'(i5,50x,a19)') ifile,note(1:19)
      do 50 i = 1,200
         read (fil1,'(a)',err=50,end=80) subfile
         if (subfile(1:1).eq.'*') goto 50
         i1 = nblen(subfile)
         if (i1.lt.1) goto 50
         write (fil2,'(a)') subfile(1:i1)
 50   continue
c      
 80   continue
      close (fil1)
      goto 1000

 100  print *,'MERGE_LIST: error in reading data file: ',subfile
      stop 'MAKED: MERGE_LIST ...'
 110  print *,'MERGE_LIST: error in opening file: ',subfile
      stop 'MAKED: MERGE_LIST ...'
         
 1000 continue
      return
      end
