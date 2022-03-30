      character*(*) function pickfn (list,flen)

      implicit none

c     pick a file name from the available choices
c     rwk 970831:  use '0' from list to set file to 'none' (allow skipping entry in some programs)
c     rwk 980721: Change from a subroutine to a function to match the new C-routine written by P. Fang
c                 The C function has been substituted for the Fortran for most systems, but it uses
c                 the function fnmatch which is not available under Sun OS/4 CC and older versions of GCC
c                 The subroutine getdir is in the same catagory (.f on os4, .c elsewhere)

c     UNIX wild card             
      character*(*) list
      character*80 files(1000),buff80,prog_name,file
      character*1 letter
      integer nch,plen,rcpar,i
      integer ipick,nfiles,len,flen,trimlen
      integer ioerr,nerr
c      logical fcheck
      logical okay

c     get calling program name and m-file name for report_stat
      plen = rcpar(0,prog_name)

      nerr = 0  
      call getdir (list,1000,files,nfiles)

      if (nfiles .gt. 0) then
         write (6,'(/,a)') 'Available files:'
         do 10 i = 1,nfiles
            write (*,'(i3,1x,a)') i,files(i)
  10     continue

  20     continue
         if (nerr .gt. 5) then
           call report_stat('FATAL',prog_name,'lib/pickfn',' '
     .                      ,'Tried 5 times to find file',0)
         endif

         write (6,'(/,1x,a,1x,$)') 
     .           'Enter a file name or pick a number (0 for none):'
         nch = len(file)
         read (5,'(a)') buff80(1:nch)

         call ljust (80,buff80)
         letter = buff80(1:1)
         if (lge(letter,'1') .and. lle(letter,'9')) then
            read (buff80,'(i3)',iostat=ioerr) ipick
            if (ioerr .eq. 0) then
               if (ipick .lt. 1 .or. ipick .gt. 1000) then
                  nerr = nerr + 1
                  goto 20
               else
                  file = files(ipick)(1:nch)
               endif
            else
               nerr = nerr + 1
               goto 20
            endif
         elseif ( letter.eq.'0' ) then  
c               Curiously, setting the filename to blank causes the result of 
c               inquire ( file=filename,exist=lexist ) to return lexist as T,
c               at least with Sun OS/4.  Hence I'm using the filename 'none'.--rwk 970902
            file = 'none' 
            goto 99
         else
            file = buff80(1:nch)
         endif 
         flen = trimlen(file)
         pickfn = file(1:flen)
c        double check existence of file 
c        don't use fcheck since it's in gamit/lib but pickfn is in /libraries/comlib
c        if (.not.fcheck(file)) then
         inquire (file=file,exist=okay)
         if( .not.okay ) then
            call report_stat('WARNING',prog_name,'lib/pickfn',file
     .                    , 'Cannot find file: ',0)
            nerr = nerr + 1
            goto 20
         endif
      else
         call report_stat('WARNING',prog_name,'lib/pickfn',list
     .                   ,'Cannot find list: ',0)
         file = ' '
      endif

   99 return
      end       



