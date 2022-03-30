      subroutine oscrat (isunit,jdb)
c     open scratch space
      implicit none
      integer*4 isunit,ioerr,jdb
      character*14   ftmp
                
c      Assign the names so that multiple arcs can run on a single machine
      
      ftmp = 'tmp_'
      write(ftmp(5:14),'(i2,a1,i7)') isunit,'_',jdb
      open( unit=isunit,file=ftmp,status='unknown',form='unformatted'
     .    , iostat=ioerr )
c      OPEN(UNIT=ISUNIT,STATUS='SCRATCH',FORM='UNFORMATTED',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','ARC','oscrat','scratch',
     .  'Error could not open scratch file necessary for ARC run',ioerr)
      endif

      return
      end
