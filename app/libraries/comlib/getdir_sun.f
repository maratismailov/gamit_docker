      subroutine getdir (list,maxn,files,n)

      implicit none
c
c     return a listing of the directory

c     input:   list   = UNIX wild card
c     output   files  = files meeting wild card spec
c              n      = number of files

      character*(*) list
      character*256 message
      character*80 command,prog_name
      character*80 files(*)
      integer i, n, maxn, luin, lunit
      integer istat,istat2,system,ioerr,unlink,ftell,len,rcpar
      integer ioff6
      logical okay
      integer trimlen

*     Character string with a quote in to make it easier to add
*     quotes to written strings. 
      character*(1) quote

      EXTERNAL lunit,unlink,system,fcheck

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)
              
*      Set the quote string (Ascii 39 is a quote)
      quote = char(39)

c     delete the temporary file if it exists.
      istat = 0   
      inquire (file='tmp.getdir',exist=okay)
      if( okay ) istat = unlink('tmp.getdir')
c     don't use fcheck since it's in gamit/lib but getdir is in /libraries/comlib
c      if (fcheck('tmp.getdir'))  istat = unlink('tmp.getdir')
      if (istat .ne. 0) then
         call report_stat('WARNING',prog_name,'lib/getdir','tmp.getdir '
     .          , 'Problem deleting temporary file',istat)
         call ferror(istat,6)
         n = 0
         return
      endif

* MOD TAH 980204: force the ls to be done under csh so that we can
*     use the re-direction features of csh with out knowing what sh
*     is being used for the run.  
      write (command,10) quote, list(1:trimlen(list)), quote
 10   format('/bin/csh -c ',a1,
     .       '\ls -1 ',a,' |& egrep -v match > tmp.getdir',a1)

c     write the shell command
C     write (command,10) list(1:trimlen(list))
* CRON DEBUG removed & since this does not work with sh-posix
C10   format ('\ls -1 ',a,' | egrep -v match > tmp.getdir') 
C10   format ('\ls -1 ',a,' |& egrep -v match > tmp.getdir') 
c rwk 940620: This form of 'ls' is designed to avoid any user-defined aliases and
c     system-dependent locations for 'ls'.  It will work as advertised with the
c     HP compiler, pre-1997 Gnu compilers, and on newer Sun compilers.  With Sun 
c     compilers 1.4 and older and with the newer Gnu compiler, the '\' is interpreted  
c     as an escape and ignored, leaving just 'ls'.  If 'ls' is aliased to a longer 
c     form (e.g. ls -a), getdir will fail.  Users confronting this problem should 
c     either remove the alias or use the proper absolute for their system (e.g. 
c     '/usr/bin/ls')
c
      n = trimlen(command)      
      command = command(1:n)

c     remember where we are in unit 5
c     ioff5 = ftell(5)
c     if (ioff5.lt.0) then
c        write (6,*) 'GETDIR: error on ftell. IOFF5 = ',ioff5
c        call ferror(-ioff5,6)
c     endif

c     remember where we are in unit 6
      ioff6 = ftell(6)
      if (ioff6.lt.0) then
        call report_stat('WARNING',prog_name,'lib/getdir',' '
     .          , 'Problem with unit 6',ioff6)
        call ferror(-ioff6,6)
      endif

c     execute the shell command
      istat = system(command)

c     restore file pointer on unit 5, because system screws it up.
c     istat2 = fseek(5,ioff5,0)
c     if (istat2.ne.0) then
c        write (6,*) 'GETDIR: error on seek on unit 5.'
c        call ferror(istat2,6)
c     endif

c     restore file pointer on unit 6, because system screws it up.
c** function       istat2 = fseek(6,ioff6,0)
c** intrinsic      call fseek(6,ioff6,0)
      call fseekg(6, ioff6, 0, istat2)
      if (istat2.ne.0) then
        call report_stat('WARNING',prog_name,'lib/getdir',' '
     .          , 'Error on seek on unit 6',istat2)
         call ferror(istat2,6)
      endif

      if (istat .eq. 0) then
         luin = lunit()
         open (unit = luin,file = 'tmp.getdir',iostat = ioerr
     .        , status='unknown' )
         if (ioerr .eq. 0) then
            n = 0
            do i =1,maxn
               read (luin,'(a)',err=30,end = 30) files(i)
               n = n + 1
            enddo
 30         continue
            close(luin)
         else
            call report_stat('WARNING',prog_name,'lib/getdir'
     .          ,'tmp.getdir','Error opening temporary file',ioerr)
            call ferror(ioerr,6)
            n = 0
         endif
      else if (istat .eq. 256) then   
c        no longer print this message here--put warning in calling routine -rwk 970212
c        write(message,'(2a)') 'No files matching list: ',list
c        call report_stat('WARNING',prog_name,'lib/getdir',' ',message,0)
         n = 0
      else
         write(message,'(2a)') 'Error executing shell command ',command
         call report_stat('WARNING',prog_name,'lib/getdir',' '
     .                   , message,istat)
         call ferror(istat,6)
         n = 0
      endif

      return
      end
