      Subroutine readm1 (lunit,
     $                   nepch,mtpart,
     $                   alabel,idms,islot1,aprval,adjust,
     $                   nsat,isat,
     $                   nsite,sitet,
     $                   nrfile,rfname,
     $                   ntfile,tfname,
     .                   norb )
c
c Purpose:
c     read the main m-file header, corresponding to
c     all sessions
c
c M. Murray 890703; last modified for CFMRG version 9.40  by R. King 980914
c
c---------------------------------------------------------------------------
c
      implicit none

      include '../includes/dimpar.h'
c
      integer       lunit
c                                     nepch = total number of epochs (le.maxepc)
      integer*4     nepch
c                                     ndy   = total number of sessions (now always = 1, read but not passed)
      integer*4     ndy
c                                     mtpart= total number of partial parameters (le.maxprm)
      integer*4     mtpart
c                                     alabel= parameter labels displayed in solve menu
      character*20  alabel(maxprm)
c                                     idms  = degree, minute, sec flag for parameters
      integer*4     idms(maxprm)
c                                     islot1= parameter codes (see fills1)
      integer*4     islot1(maxprm)
c                                     aprval= a priori parameter values
      real*8        aprval(maxprm)
c                                     adjust= adjustments to a priori parameter values
      real*8        adjust(maxprm)
c                                     nsat  = total number of satellites (le.maxsat)
      integer*4     nsat
c                                     isat  = satellite prn numbers
      integer*4     isat(maxsat)
c                                     nsite = total number of sites 
      integer*4     nsite
c                                     sitet = monument labels
      character*12  sitet(maxsit)
c                                     nrfile= total number of c-files (le.maxsit)
      integer*4     nrfile
c                                     rfname= all x-file names
      character*16  rfname(maxsit)
c                                     ntfile= total number of t-files (le.maxtfl)
c                                             (currently = ndy)
      integer*4     ntfile
c                                     tfname= t-file names
      character*16  tfname(maxtfl)
c                                     norb= number of orbit parameters (ICs + non-grav)
      integer*4     norb

c Internal:
c                                     iflag = m-file header flag
      integer       iflag
c                                     nwds  = number of words in header following iflag and nwds
      integer*4     nwds
c                                     i     = index
      integer       i
c                                     nversn = version number
      integer*4     nversn
c                                     ioerr = system iostat variable (must be i*4)
      integer*4     ioerr
c                                     len = dummy argument for integer function rcpar, which returns the calling program
      integer*4     len,rcpar


      character*16  mfile
      character*80  prog_name
      character*256 message

c subroutines and functions:  report_stat, rcpar

c-----------------------------------------------------------------------
          
c  Get the calling program name for report_stat
       len = rcpar(0,prog_name)

  
       inquire( unit=lunit, name=mfile, iostat=ioerr )
       if( ioerr.ne.0 ) goto 1000

       do i=1,maxprm
         alabel(i) = ' '
         idms(i) = 0
         islot1(i) = 0
         aprval(i) = 0.d0
         adjust(i) = 0.d0
       enddo

c       first read only the first three integers, to make sure we have the right version
                                  
       read( unit = lunit ) iflag
       rewind(lunit)
       read( unit = lunit ) iflag,nwds
       rewind (lunit)
       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,nversn

c      if this is an earlier version, nversn will be the old ndy and be < 100
      if( nversn.lt.940 ) then
         call report_stat('FATAL',prog_name,'lib/readm1',' '
     .              ,'Old version of M-file ( <9.31 )--rerun CFMRG',0)
      endif
      rewind ( lunit )


c       next read the number of parameters and compare with the dimensions

       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,nversn,
     $       ndy,nepch,mtpart
      if( mtpart.gt.maxprm ) then 
         write(message,'(a,i6,a,i6)') '# parameters on m-file=',mtpart
     .       ,' > maxprm=',maxprm
         call report_stat('FATAL',prog_name,'lib/readm1',' ',message,0)
      endif           
      if( ndy.ne.1 ) call report_stat('FATAL',prog_name,'lib/readm1',' '
     .    ,'# sesssions ne 1, multisession no longer supported ',0 ) 
      rewind ( lunit )
     
c       next read the number of satellites and compare with the dimensions
                             
       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,nversn,
     $       ndy,nepch,mtpart,
     $       (alabel(i),i=1,mtpart),
     $       (idms(i),i=1,mtpart),
     $       (islot1(i),i=1,mtpart),
     $       (aprval(i),i=1,mtpart),
     $       (adjust(i),i=1,mtpart),
     $       nsat
      if( nsat.gt.maxsat ) then 
         write(message,'(a,i2,a,i2)') '# sats on m-file=',nsat
     .       ,' > maxsat=',maxsat
         call report_stat('FATAL',prog_name,'lib/readm1',' ',message,0)
      endif 
      rewind ( lunit )

c       next read the number of sites and compare with the dimensions

       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,nversn,
     $       ndy,nepch,mtpart,
     $       (alabel(i),i=1,mtpart),
     $       (idms(i),i=1,mtpart),
     $       (islot1(i),i=1,mtpart),
     $       (aprval(i),i=1,mtpart),
     $       (adjust(i),i=1,mtpart),
     $       nsat,(isat(i),i=1,nsat),
     $       nsite
      if( nsite.gt.maxsit ) then 
         write(message,'(a,i3,a,i3)') '# sites on m-file=',nsite
     .       ,' > maxsit=',maxsit
         call report_stat('FATAL',prog_name,'lib/readm1',' ',message,0)
      endif 
      rewind ( lunit )


c       Now read the record for real

       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,nversn,
     $       ndy,nepch,mtpart,
     $       (alabel(i),i=1,mtpart),
     $       (idms(i),i=1,mtpart),
     $       (islot1(i),i=1,mtpart),
     $       (aprval(i),i=1,mtpart),
     $       (adjust(i),i=1,mtpart),
     $       nsat,(isat(i),i=1,nsat),
     $       nsite,(sitet(i),i=1,nsite),
     $       nrfile,(rfname(i),i=1,nrfile),
     $       ntfile,(tfname(i),i=1,ntfile),
     $       norb
c
 1000 if (ioerr .ne. 0) then
         call report_stat('FATAL',prog_name,'lib/readm1',mfile
     .                      ,'Error reading M-file',ioerr)
      endif

      if (iflag .ne. 1) then
        write(message,'(a,i4)') 'Wrong iflag: ',iflag
        call report_stat('FATAL',prog_name,'lib/readm1',mfile,message,0)
      endif
c
c     print '(14i5)', iflag,nwds,nversn,mtpart,nsat,nsite,ndy,nepch
c     print '(a20)', (alabel(i),i=1,mtpart)
c     print '(14i5)', (idms(i),i=1,mtpart)
c     print '(14i5)', (islot1(i),i=1,mtpart)
c     print '(3(d22.15,1x))', (aprval(i),i=1,mtpart)
c     print '(3(d22.15,1x))', (adjust(i),i=1,mtpart)
c     print '(14i5)', (isat(i),i=1,nsat)
c     print '(a)', (sitet(i),i=1,nsite)
c     print '(14i5)', nrfile  
      
c     print '(a16)', (rfname(i),i=1,nrfile)
c     print '(14i5)', ntfile
c     print '(a16)', (tfname(i),i=1,ntfile)
c     print '(i3)',norb
c
      return
      end
