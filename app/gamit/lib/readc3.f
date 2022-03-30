      subroutine readc3 (lunit
     .,          nlabel,nparam
     .,          islot,idms,rlabel,preval)

c     read the third part of a C-file which tells about the parameters
c     K. Feigl/R. King July 89

      implicit none

      include '../includes/dimpar.h'

      integer            lunit,nlabel,nparam
      integer            iflag,i
      integer*4          ioerr,inqerr,len,rcpar


      integer            islot(maxlab),idms(maxlab)
      real*8             preval(maxprm)
      character*20       rlabel (maxlab)

      character*4        buf4
      character*16       cfname
      character*80       prog_name

c     get the calling module name for report_stat
      len =rcpar(0,prog_name)

      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag,nlabel,nparam
     .,     (islot (i),i=1,nlabel)
     .,     (idms  (i),i=1,nlabel)
     .,     (rlabel(i),i=1,nlabel)
     .,     (preval(i),i=1,nparam)

 1000 if (ioerr .ne. 0) then
         inquire ( unit=lunit,name=cfname,iostat=inqerr )
         call report_stat('FATAL',prog_name,'lib/readc3',cfname
     .                   ,'Error reading C-file',ioerr)
      endif

      if (iflag .ne. 3) then
         write(buf4,'(i4)') iflag
         call report_stat('FATAL',prog_name,'lib/readc3',buf4
     .                    ,'Wrong iflag: ',0)
      endif

CD     print *,'READC3: nparam ',nparam
CD     print *,'READC3: nlabel ',nlabel
CD     print *,'READC3: islot  ',(islot (i),i=1,nlabel)
CD     print *,'READC3: idms   ',(idms  (i),i=1,nlabel)
CD     print *,'READC3: rlabel ',(rlabel(i),i=1,nlabel)
CD     print *,'READC3: preval ',(preval(i),i=1,nparam)

      return
      end
