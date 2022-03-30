      Subroutine readm2 (lunit,
     $                   it0,t00,
     $                   nepch,inter,nskip,
     $                   ncfile,cfname,stawgt,
     $                   nsat,satwgt )
c
c Purpose:
c     Read the second m-file header, corresponding to
c     session subheaders
c
c M. Murray 890703;  last modified for version 9.21ff by R. King 950616
c
c-----------------------------------------------------------------------

      implicit none

      include '../includes/dimpar.h'

c input (all unchanged on output):
c                                     lunit = m-file logical unit
      integer       lunit
c                                     it0   = first epoch time tag: year, month, day
      integer       it0(3)
c                                     t00   = first epoch time tag: hour, minute, second
      real*8        t00(3)
c                                     nepch = total number of epochs (le.maxepc)
      integer       nepch
c                                     inter = time interval between epochs (seconds)
      integer       inter
c                                     nskip = decimation factor from SOLVE (initially=1)
      integer       nskip
c                                     ncfile= number of cfiles 
      integer       ncfile
c                                     cfname= c-file names in session
      character*16  cfname(maxsit)
c                                     stawgt= station weights
      real*8        stawgt(maxsit)
c                                     nsat  = number of satellites (le.maxsat)
      integer       nsat
c                                     satwgt= satellite weights
      real*8        satwgt(maxsat)
c
c internal:
c                                     iflag = m-file header flag
      integer       iflag
c                                     nwds  = number of words in header following iflag and nwds
      integer       nwds
c                                     i     = index
      integer       i             
c                                     idy   = dummy session number (always = 1)
      integer       idy

c                                     ioerr = system iostat variable (must be i*4)
      integer*4     ioerr

      integer*4     len,rcpar
      character*16  mfile
      character*80  prog_name
      character*256 message

      logical debug/.false./

c subroutines and functions:  report_stat
c
c-----------------------------------------------------------------------

       inquire( unit=lunit, name=mfile, iostat=ioerr )
       if( ioerr.ne.0 ) goto 1000

       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,
     $       idy,it0,t00,
     $       nepch,inter,nskip,
     $       ncfile,(cfname(i),i=1,ncfile),(stawgt(i),i=1,ncfile),
     $       nsat,(satwgt(i),i=1,nsat)
c                   

 1000 if (ioerr .ne. 0) then
c     Get calling program name and m-file name for report_stat
         len = rcpar(0,prog_name)
         call report_stat('FATAL',prog_name,'lib/readm2',mfile
     .                      ,'Error reading M-file',ioerr)
      endif
c
      if (iflag .ne. 2) then
        write(message,'(a,i4)') 'Wrong iflag: ',iflag
        call report_stat('FATAL',prog_name,'lib/readm2',mfile,message,0)
      endif

      if(debug) then 
        print '(1x,14i5)', iflag,nwds,nsat,idy
        print '(1x,14i5)', it0
        print '(1x,3(d22.15,1x))', t00
        print '(1x,14i5)', nepch,inter,nskip,ncfile
        print '(/,1x,40a16)', (cfname(i),i=1,ncfile)
        print '(/,1x,40d7.1)',(stawgt(i),i=1,ncfile)
      endif 

      return
      end
