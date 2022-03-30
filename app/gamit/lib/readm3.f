      Subroutine readm3 (lunit,
     $                   it0,t00,
     $                   kpart,islot2,
     $                   elvcut,zenmod,numzen,idtzen,
     $                   gradmod,numgrad,idtgrad )

c
c Purpose:
c     Read the third M-file header, corresponding to a
c     single C-file subheader
c
c M. Murray 890703; last modified for version 9.40ff by R. King 980914
c
c-----------------------------------------------------------------------
      implicit none

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
      integer*4      idtzen(maxatm)  
c                                     gradmod= gradient model ('CON' or 'PWL')
      character*3   gradmod
c                                     numgrad= number of gradient parameters
      integer*4     numgrad
c                                     idtgrad= tabular points of piecewise linear gradient model
      integer*4      idtgrad(maxatm)



c  Note:  Current SOLVE requires that the elevation cutoff and number of zenith-delay parameters
c         be the same for all sites, but we made want to make these site-specific in the future.
c
c Internal:
c                                     iflag = m-file header flag
      integer       iflag
c                                     nwds  = number of words in header following iflag and nwds
      integer       nwds
c                                     i     = index
      integer       i
c                                     ioerr = system iostat variable (must be i*4)
      integer*4     ioerr
c                                     len = dummy argument for integer function rcpar which returns to the calling module
      integer*4     len,rcpar

      character*16  mfile
      character*80  prog_name
      character*256 message

      logical debug/.false./

c subroutines and functions:  report_stat

c-----------------------------------------------------------------------

c      get calling program name and m-file name for report_stat
       len = rcpar(0,prog_name)
       inquire( unit=lunit, name=mfile, iostat=ioerr )
       if( ioerr.ne.0 ) goto 1000

c      first read the array sizes to make sure they're not too large

      read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,
     $       (it0(i),i=1,3),(t00(i),i=1,3),
     $       kpart
       if( kpart.gt.maxprm ) then
        write(message,'(a,i6,a,i6,a)')
     $    'kpart (=',kpart,') for islot2 on m-file > maxprm (='
     $    ,maxprm,') in dimpar.h'
        call report_stat('FATAL',prog_name,'lib/readm3',mfile,message,0)
       endif
       backspace lunit
       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,
     $       (it0(i),i=1,3),(t00(i),i=1,3),
     $       kpart,(islot2(i),i=1,kpart),
     $       elvcut,zenmod,numzen
       if( numzen.gt.maxatm ) then
         write(message,'(a,i2,a,i2,a)')
     $    'numzen (=',numzen,') on m-file > maxatm (='
     $   ,maxatm,') in dimpar.h'
        call report_stat('FATAL',prog_name,'lib/readm3',mfile,message,0)
       endif 

       backspace lunit
       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,
     $       (it0(i),i=1,3),(t00(i),i=1,3),
     $       kpart,(islot2(i),i=1,kpart),
     $       elvcut,zenmod,numzen,(idtzen(i),i=1,numzen),
     $       gradmod,numgrad
       if( 2*numgrad.gt.maxgrad ) then
         write(message,'(a,i2,a,i2,a)')
     $    '2*numgrad (=',2*numgrad,') on m-file > maxgrad (='
     $   ,maxgrad,') in dimpar.h'
        call report_stat('FATAL',prog_name,'lib/readm3',mfile,message,0)
       endif
 
       backspace lunit

c      now read the entire record for real

       read (unit   = lunit,
     $       iostat = ioerr,
     $       end    = 1000,
     $       err    = 1000)
     $       iflag,nwds,
     $       (it0(i),i=1,3),(t00(i),i=1,3),
     $       kpart,(islot2(i),i=1,kpart),
     $       elvcut,zenmod,numzen,(idtzen(i),i=1,numzen),
     $       gradmod,numgrad,(idtgrad(i),i=1,numgrad)


 1000 if (ioerr .ne. 0) then
         call report_stat('FATAL',prog_name,'lib/readm3',mfile
     .                      ,'Error reading M-file',ioerr)
      endif

      if (iflag .ne. 3) then
        write(message,'(a,i4)') 'Wrong iflag: ',iflag
        call report_stat('FATAL',prog_name,'lib/readm3',mfile,message,0)
      endif
      if(debug) then 
        print *,'READM3: iflag,nwds = ',iflag,nwds
        print *,'READM3: it0,t00 = ',it0,t00
        print *,'READM3: kpart islot2 = ', kpart, (islot2(i),i=1,kpart)
        print *,'READM3: elvcut,zenmod = ',elvcut,zenmod
        print *,'READM3: numzen,idtzen = ',numzen,(idtzen(i),i=1,numzen)
        print *,'READM3: numgrad,idtgrad = ',numgrad
     .      ,(idtgrad(i),i=1,numgrad)                                   
      endif 

      return
      end
