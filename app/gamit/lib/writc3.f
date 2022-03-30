      subroutine writc3 (lunit
     .,          nlabel,nparam
     .,          islot,idms,rlabel,preval)

c     write the third part of a C-file which tells about the parameters
c     K. Feigl/R. King July 89

      implicit none

      include '../includes/dimpar.h'

      integer            lunit,nlabel,nparam
      integer            iflag,i
      integer*4          ioerr
      integer*4          len,rcpar


      integer            islot(maxlab),idms(maxlab)
      real*8             preval(maxprm)
      character*20       rlabel (maxlab)
      character*16       cfname
      character*80       prog_name

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)
      inquire( unit=lunit, name=cfname, iostat=ioerr )
      if( ioerr.ne.0 ) goto 1000

      iflag = 3

      write (unit    =  lunit
     .,     iostat  =  ioerr
     .,     err     =  1000)
     .      iflag,nlabel,nparam
     .,     (islot (i),i=1,nlabel)
     .,     (idms  (i),i=1,nlabel)
     .,     (rlabel(i),i=1,nlabel)
     .,     (preval(i),i=1,nparam)

 1000 if (ioerr .ne. 0) then
         call ferror (ioerr,6)
         call report_stat('FATAL',prog_name,'lib/writc3',cfname
     .                      ,'Error writing C-file header',ioerr)
      endif

CD     print *,'WRITC3: nparam ',nparam
CD     print *,'WRITC3: nlabel ',nlabel
CD     print *,'WRITC3: islot  ',(islot (i),i=1,nlabel)
CD     print *,'WRITC3: idms   ',(idms  (i),i=1,nlabel)
CD     print *,'WRITC3: rlabel ',(rlabel(i),i=1,nlabel)
CD     print *,'WRITC3: preval ',(preval(i),i=1,nparam)

      return
      end
