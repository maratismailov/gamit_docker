      subroutine writc1 (lunit,ntext,text)

      implicit none

      include '../includes/dimpar.h'

      integer*4          nversn,iflag,ntext
      character*80       text(maxtxt)
      integer*4          lunit,ioerr,i
      integer*4 len,rcpar
      character*80 prog_name
      character*16 cfname

c     write the first part of a C-file which is a textual history      
c     K. Feigl/R. King  July 1989    
c     R. King Sep 1998: Modfied to add 'norb' parameter   
c     R. King Feb 2005  Modified for version 10.2 to add more model names, avg loading 
c                       ocean/atm loading values, and more extra slots    
c     R. King Aug 2010  Modified for releaase 10.4 to length the names of the
c                       antenna models   
c     R. King Jan 2013  Modified for release 10.41 (dryzen, wetzen added to record 2)  
c     R. King Mar 2014  Added E-radiation and antenna-radiation models to record 2
c     R. King Aug 2014  Added gnss character to record 2, version 10.60
c     T. Herring Jan 2020: Added L1/L2 satellite PCO. version 10.71.

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)
      inquire( unit=lunit, name=cfname, iostat=ioerr )
      if( ioerr.ne.0 ) goto 1000

      iflag = 1
* MOD TAH 200126: Increased version to 1071 from 1060.
      nversn = 1071

      write (unit    =  lunit
     .,     iostat  =  ioerr
     .,     err     =  1000)
     .      iflag,nversn,ntext,(text(i),i=1,ntext)

 1000 if (ioerr .ne. 0) then
         call ferror (ioerr,6)
         call report_stat('FATAL',prog_name,'lib/writc1',cfname
     .                   ,'Error writing first C-file record',ioerr)
      endif

c     print *,'WRITC1: ntext = ',ntext,' text follows: '
c     do 2000 i=1,ntext
c        print *,text(i)
c2000 continue
      return
      end
