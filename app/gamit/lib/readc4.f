      subroutine readc4 (lunit
     .,                  msat,mprn
     .,                  iepoch,iyr,idoy,sod,rclock
     .,                  okmet,zendel,pres,temp,relhum,atmlod   
     .,                  sitecd,latr_sph,lonr,radius
     .,                  l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,                  nsave,save )

c     Read the fourth part of a C-file which has the time tag for the epoch
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
     .,                  okmet,kflag,iep
      real*4             atmlod(3)
      real*8             sod,rclock,zendel,pres,temp,relhum,save(maxsav)
     .,                  latr_sph,lonr,radius 
     .,                  l1z,l1n,l1e,l2z,l2n,l2e,antaz

c     ***C-file and internal to this routine***
      integer            iflag
      integer*4          ioerr,inqerr,len,rcpar,i
      real*4             pres4,temp4,relhum4
      character*4        buf4,sitecd
      character*16       cfname
      character*80       prog_name

c     get the calling module name for report_stat
      len = rcpar(0,prog_name)
                                        
cd      if(iep.eq.953) then 
cd         print *,'READC4 iep lunit ',iep,lunit
cd         read(lunit,iostat=ioerr) iflag,msat
cd         if(ioerr.eq.0) then
cd           print *,'READC4 read iflag msat ',iflag,msat
cd          backspace(lunit)
cd         else 
cd           print *,'READC4 error reading iflag msat '
cd         endif 
cd      endif 

      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag  
     .,     msat,(mprn(i),i=1,msat)
     .,     iepoch,iyr,idoy,sod,rclock
     .,     zendel,okmet,pres4,temp4,relhum4,atmlod
     .,     sitecd,latr_sph,lonr,radius
     .,     l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,     nsave,(save(i),i=1,nsave)  
c     convert r*4 variables to r*8 (except for atmlod)
      pres   = pres4
      temp   = temp4
      relhum = relhum4                    
          
 1000 if (ioerr .ne. 0) then
          inquire ( unit=lunit,name=cfname,iostat=inqerr )
          call report_stat('FATAL',prog_name,'lib/readc4',cfname
     .                    ,'Error reading C-file',ioerr)
      endif

      if (iflag .ne. 4) then
         write(buf4,'(i4)') iflag
         call report_stat('FATAL',prog_name,'lib/readc4',buf4
     .                   ,'Wrong iflag: ',0)
      endif
 
      if( iepoch.lt.0 ) then
       print *,'READC4: msat  = ',            msat
       print *,'READC4: mprn  = ',            (mprn(i),i=1,msat)
       print *,'READC4: iepoch= ',            iepoch
       print *,'READC4: iyr,idoy,sod= ',      iyr,idoy,sod
       print *,'READC4: rclock= ',            rclock
       print *,'READC4: okmet=',              okmet
       print *,'READC4: pres4,temp4,relhum4=',pres4,temp4,relhum4
       print *,'READC4: atmlod=',             atmlod
       print *,'READC4: nsave = ',            nsave
       print *,'READC4: save  = ',            (save(i),i=1,nsave)
       print *,'READC4: l1z,l1n,l1e = ',      l1z,l1n,l1e
       print *,'READC4: l2z,l2n,l2e = ',      l2z,l2n,l2e
       print *,'READC4: antaz = ',            antaz
      endif
      return
      end
