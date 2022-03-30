      Subroutine writm1 (lunit,
     $                   nepch,mtpart,
     $                   alabel,idms,islot1,aprval,adjust,
     $                   nsat,isat,
     $                   nsite,sitet,
     $                   nrfile,rfname,
     $                   ntfile,tfname,
     .                   norb )
c
c purpose:
c     write the main m-file header, corresponding to
c     all sessions
c
c M. Murray 890703; last modified for version 9.40ff by R. King 980914
c                   
c
c-----------------------------------------------------------------------
c
      include '../includes/dimpar.h'
c
c input (all unchanged on output):
c                                     lunit = m-file logical unit
      integer       lunit
c                                     nepch = total number of epochs (le.maxepc)
      integer       nepch
c                                     ndy   = total number of sessions (le.maxtfl)
      integer       ndy
c                                     mtpart= total number of partial parameters (now always = 1, 
c                                        written but not passed) 
      integer       mtpart
c                                     alabel= parameter labels displayed in solve menu
      character*20  alabel(maxprm)
c                                     idms  = degree, minute, sec flag for parameters
      integer       idms(maxprm)
c                                     islot1= parameter codes (see fills1)
      integer       islot1(maxprm)
c                                     aprval= a priori parameter values
      real*8        aprval(maxprm)
c                                     adjust= adjustments to a priori parameter values
      real*8        adjust(maxprm)
c                                     nsat  = total number of satellites (le.maxsat)
      integer       nsat
c                                     isat  = satellite prn numbers
      integer       isat(maxsat)
c                                     nsite = total number of sites 
      integer       nsite
c                                     sitet = monument labels
      character*12  sitet(maxsit)
c                                     nrfile= total number of c-files (le.maxsit)
      integer       nrfile
c                                     rfname= all x-file names
      character*16  rfname(maxsit)
c                                     ntfile= total number of t-files (le.maxtfl)
c                                             (currently = ndy)
      integer       ntfile
c                                     tfname= t-file
      character*16  tfname
c                                     norb= nubmer of orbit parameters (ICs + non-grav)
      integer*4     norb

c
c internal:
c                                     iflag = m-file header flag
      integer       iflag
c                                     nwds  = number of words in header following iflag and nwds
      integer       nwds
c                                     i     = index
      integer       i
c                                     nversn = version number
      integer*4     nversn
c                                     ioerr = system iostat variable (must be i*4)
      integer*4     ioerr
c                                     len = (dummy) length of program name for rcpar call
      integer*4     len,rcpar
c                                     m-file name
      character*16  mfile
c                                     prog_name = name of calling program
      character*80  prog_name

      logical debug/.false./                              

c
c subroutines and functions:  ferror, report_status
c
c-----------------------------------------------------------------------
c

      iflag=1
      nwds=20+20*mtpart+2*nsat+6*nsite+8*nrfile+8*ntfile
* MOD TAH 200126: Increased version to 1071 from 1061
C     nversn = 1071
* MOD TAH 210701: Increased version to 1072 from 1071 to account
*     for changes in slot numbers to accomodate up to 45 Beidou satellites.
*     The new assigments limit the max number of satellites to 45 with
*     99 stations.  
      nversn = 1072

c            
c     multi-session no longer supported
         ndy = 1 
      write (unit   = lunit,
     $       iostat = ioerr,
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
     $       ntfile,tfname,
     $       norb
c
 1000 if (ioerr .ne. 0) then
c       get calling program name and m-file name for report_stat
        len = rcpar(0,prog_name)
        inquire( unit=lunit, name=mfile, iostat=ioerr )
        call ferror (ioerr,6)
        call report_stat('FATAL',prog_name,'lib/writm1',mfile
     .                  ,'Error writing 1st M-file record',ioerr)
      endif       

      if(debug) then 
   
        print *,'WRITM1: iflag,nwds,nversn ',iflag,nwds,nversn
        print *,'WRITM1: mtpart,nsat,nsite,ndy,nep   h '
     .              , mtpart,nsat,nsite,ndy,nep   h
        print *,'WRITM1: alabel ',(alabel(i),i=1,mtpart)
        print *,'WRITM1: idms ', (idms(i),i=1,mtpart)
        print *,'WRITM1: islot1 ', (islot1(i),i=1,mtpart)
        print *,'WRITM1: aprval ', (aprval(i),i=1,mtpart)
        print *,'WRITM1: adjust ', (adjust(i),i=1,mtpart)
        print *,'WRITM1: isat ',(isat(i),i=1,nsat)
        print *,'WRITM1: sitet ', (sitet(i),i=1,nsite)
        print *,'WRITM1: nrfile rfname ',nrfile,(rfname(i),i=1,nrfile)
        print *,'WRITM1: ntfile tfname ',ntfile,tfname
        print *,'WRITM1: norb ',norb
      endif 

      return
      end
