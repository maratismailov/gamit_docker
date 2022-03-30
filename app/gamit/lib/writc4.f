      subroutine writc4 (lunit
     .,                  msat,mprn
     .,                  iepoch,iyr,idoy,sod,rclock
     .,                  okmet,zendel,pres,temp,relhum,atmlod  
     .,                  sitecd,latr_sph,lonr,radius
     .,                  l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,                  nsave,save )
  

c     Write the fourth part of a C-file which has the time tag for the epoch
c     K. Feigl/R. King July 89  
c     R. King May 1995: Modified to pass kinematic variables (formerly in common)
c         and to add antaz for C-file version 9.30ff.  
c     R. King August 2010: Add atmospheric loading values. 

c     Notes:  The met data are written as real*4 on the file but passed
c             around in GAMIT software as real*8.

      implicit none

      include '../includes/dimpar.h'

c     ***passed***
      integer*4          lunit,msat,iepoch,iyr,idoy,nsave,mprn(maxsat)
     .,                  okmet
      real*4             atmlod(3)
      real*8             sod,rclock,zendel,pres,temp,relhum,save(maxsav)
     .,                  latr_sph,lonr,radius
     .,                  l1z,l1n,l1e,l2z,l2n,l2e,antaz

c     ***C-file and internal to this routine***
      integer            iflag
      integer*4          ioerr,i
      integer*4 len,rcpar
      real*4             pres4,temp4,relhum4  
      character*4        sitecd
      character*16       cfname
      character*80       prog_name

      iflag = 4

c     convert r*8 variables to r*4
      pres4   = pres
      temp4   = temp
      relhum4 = relhum

      write (unit    =  lunit
     .,     iostat  =  ioerr
     .,     err     =  1000)
     .      iflag
     .,     msat,(mprn(i),i=1,msat)
     .,     iepoch,iyr,idoy,sod,rclock
     .,     zendel,okmet,pres4,temp4,relhum4,atmlod
     .,     sitecd,latr_sph,lonr,radius
     .,     l1z,l1n,l1e,l2z,l2n,l2e,antaz  
     .,     nsave,(save(i),i=1,nsave) 

 1000 if (ioerr .ne. 0) then
         call ferror (ioerr,6)
c        get calling program name and m-file name for report_stat
         len = rcpar(0,prog_name)
         inquire( unit=lunit, name=cfname, iostat=ioerr )
         call report_stat('FATAL',prog_name,'lib/writc4',cfname
     .                   ,'Error writing C-file data record',ioerr)
      endif
    
         
      if( iepoch.lt.0 ) then
       print *,'WRITC4: msat  = ',          msat
       print *,'WRITC4: mprn  = ',          (mprn(i),i=1,msat)
       print *,'WRITC4: iepoch= ',          iepoch
       print *,'WRITC4: iyr,idoy,sod= ',    iyr,idoy,sod
       print *,'WRITC4: rclock= ',          rclock
       print *,'WRITC4: okmet=',            okmet
       print *,'WRITC4: pres,temp,relhum=', pres4,temp4,relhum4 
       print *,'WRITC4: atmlod=',           atmlod
       print *,'WRITC4: nsave = ',          nsave
       print *,'WRITC4: save  = ',          (save(i),i=1,nsave)
       print *,'WRITC4: l1z,l1n,l1e = ',      l1z,l1n,l1e
       print *,'WRITC4: l2z,l2n,l2e = ',      l2z,l2n,l2e
       print *,'WRITC4: antaz = ',            antaz
      endif

      return
      end
