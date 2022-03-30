      Subroutine writm3 (lunit,
     $                   it0,t00,
     $                   kpart,islot2,
     $                   elvcut,zenmod,numzen,idtzen,
     $                   gradmod,numgrad,idtgrad )
c
c Purpose:
c     Write the third M-file header, corresponding to a
c     single C-file subheader
c
c M. Murray 890703; last modified for version 9.20ff by R. King 950524
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
c                                     kpart = number of c-file partials in islot2
      integer       kpart
c                                     islot2= c-file parameter codes (see fills2)
      integer       islot2(maxprm)
c                                     elvcut= elevation cutoff in SOLVE (degrees)
      real*8        elvcut
c                                     zenmod= zenith delay model ('CON' or 'PWL')
      character*3   zenmod
c                                     numzen= number of zenith delay parameters
      integer*4     numzen
c                                     idtzen= tabular points of piecewise linear zenith model
      integer*4     idtzen(maxatm)        
c                                     gradmod= gradith delay model ('CON' or 'PWL')
      character*3   gradmod
c                                     numgrad= number of gradient parameters
      integer*4     numgrad
c                                     idtgrad= tabular points of piecewise linear gradient model
      integer*4     idtgrad(maxatm)

c 

c  Note:  Current SOLVE requires that the elevation cutoff and number of zenith-delay and gradient
c        parameters be the same for all sites, but we made want to make these site-specific in the future.

c Internal:
c                                     iflag = m-file header flag
      integer       iflag
c                                     nwds  = number of words in header following iflag and nwds
      integer       nwds
c                                     i     = index
      integer       i
c                                     ioerr = system iostat variable (must be i*4)
      integer*4     ioerr
c                                     len = dummy for length of prog_name returned from recpar
      integer*4 len,rcpar
c                                     prog_name = name of calling program
      character*80 prog_name
c                                     mfile = name of M-file
      character*16 mfile

      logical debug/.false./
c
c subroutines and functions:  ferror, report_stat
c
c-----------------------------------------------------------------------
c
      iflag = 3
      nwds = 20 + 2*kpart
c
      write (unit   = lunit,
     $       iostat = ioerr,
     $       err    = 1000)
     $       iflag,nwds,
     $       (it0(i),i=1,3),(t00(i),i=1,3),
     $       kpart,(islot2(i),i=1,kpart),
     $       elvcut,zenmod,numzen,(idtzen(i),i=1,numzen),
     $       gradmod,numgrad,(idtgrad(i),i=1,numgrad)

 1000 if (ioerr .ne. 0) then
c        get calling program name and m-file name for report_stat
         len = rcpar(0,prog_name)
         inquire( unit=lunit, name=mfile, iostat=ioerr )
         call ferror (ioerr,6)
         call report_stat('FATAL',prog_name,'lib/writm3',mfile
     .                      ,'Error writing 3rd M-file record',ioerr)
      endif
        
      if(debug) then 
        print *,'WRITM3: iflag,nwds = ',iflag,nwds
        print *,'WRITM3: it0,t00 = ',it0,t00
        print *,'WRITM3: kpart islot2 = ', kpart, (islot2(i),i=1,kpart)
        print *,'WRITM3: elvcut,zenmod = ',elvcut,zenmod
        print *,'WRITM3: numzen,idtzen = ',numzen,(idtzen(i),i=1,numzen)
        print *,'WRITM3: numgrad,idtgrad = ',numgrad
     .                  ,(idtgrad(i),i=1,numgrad)                       
      endif

      return
      end
