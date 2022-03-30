      Subroutine writm2 (lunit,
     $                   it0,t00,
     $                   nepch,inter,nskip,
     $                   ncfile,cfname,stawgt,
     $                   nsat,satwgt )
c
c Purpose:
c     Write the second m-file header, corresponding to
c     session subheaders
c
c M. Murray 890703; last modified for version 9.21ff by R. King 950616
c
c-----------------------------------------------------------------------
c
      include '../includes/dimpar.h'
c
c Input (all unchanged on output):
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
c                                     ncfile= number of cfiles in session (le.maxnet)
      integer       ncfile
c                                     cfname= c-file names in session
      character*16  cfname(maxsit)
c                                     stawgt= station weights
      real*8        stawgt(maxsit)
c                                     nsat  = number of satellites (le.maxsat)
      integer       nsat
c                                     satwgt= satellite weights
      real*8        satwgt(maxsat)


c Internal:
c                                     iflag = m-file header flag
      integer       iflag
c                                     idy dummy session 3 (always = 1)
      integer       idy 
c                                     nwds  = number of words in header following iflag and nwds
      integer       nwds
c                                     i     = index
      integer       i
c                                     ioerr = system iostat variable (must be i*4)
      integer*4     ioerr
c                                     len = (dummy) length of program name returned by rcpar
      integer*4     len,rcpar
c                                      prog_name = name of calling program
      character*80   prog_name
c                                      mfile = M-file name
      character*16   mfile
            
      logical debug/.false./
c
c subroutines and functions:  ferror, report_status
c
c-----------------------------------------------------------------------
c
      iflag=2
      nwds = 17 + 6*ncfile + 2*nsat
c            
      idy = 1    

      write (unit   = lunit,
     $       iostat = ioerr,
     $       err    = 1000)
     $       iflag,nwds,
     $       idy,it0,t00,
     $       nepch,inter,nskip,
     $       ncfile,(cfname(i),i=1,ncfile),(stawgt(i),i=1,ncfile),
     $       nsat,(satwgt(i),i=1,nsat)

c
 1000 if (ioerr .ne. 0) then
c       get calling program name and m-file name for report_stat
        len = rcpar(0,prog_name)
        inquire( unit=lunit, name=mfile, iostat=ioerr )
         call ferror (ioerr,6)
        call report_stat('FATAL',prog_name,'lib/writm2',mfile
     .                  ,'Error writing 2d M-file header record',ioerr)
      endif
             
      if(debug) then 
        print *,'WRITM2: iflag,nwds,nsat,idy ', iflag,nwds,nsat,idy
        print *,'WRITM2: it0 t00 ',it0,t00
        print *,'WRITM2: nepch,inter,nskip ',nepch,inter,nskip
        print *,'WRITM2: ncfile ',ncfile
        print *,'WRITM2: cfname ',(cfname(i),i=1,ncfile)
        print *,'WRITM2: stawgt ',(stawgt(i),i=1,ncfile)
      endif 

      return
      end
