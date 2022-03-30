Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1994.   All rights reserved.

      Subroutine CFOUT 


C     Write the observation records on the output C-File
C       R. W. King   June 1987
C
C   IEPOCH: epoch number for printout           
c   JDOBS : PEP JD of observation epoch
c   TOBS   : seconds-of-day of observation epoch 
c   ZENDEL: zenith delay at the observation epoch   
c   ATMLOD: Atmospheric loading (N E U, mm) 
C   NCHAN : number of channels (satellites) (model.h)
C   NDAT  : number of observation types (2=L1,L2 phase, 4=pseudorange) (model.h)
C   IER=0   if station has good data and non-zero o-c s
c   data_flag: AUTCLN data flag--passed from input C-file
C   ISNR  : signal to noise ration (RINEX convention (0-9))
C   OBS   : observed oneway phases and pseudoranges
c   ampl1, ampl2 : amplitude of L1, L2 signals
C   OMC   : observed minus theoretical oneway phases
C   ELEV  : elevation angles of observed satellites
C   ELVDOT: elevation angle rates
C   AZIM  : azimuth angles of observed satellites
C   AZMDOT: azmimuth angle rates
c   NADANG: nadir angle of SV antennea 
C   ATMDEL: atmospheric delay for each satellite
c   SVCLOCK: SV clock epoch (sec) and L1 phase correction (cyc) for each satellite
C   DELAY : L1 delay for each satellite
C   DRATE : L1 delay rate for each satellite
C   NSAVE : number of extra R*8 values save
C   SAVE  : extra R*8 values for epoch record
C   NPART : number of output partial derivatives  (model.h)
C   TMPART: partial derivatives wrt site and satellite parameters
c   From model.h:   
C   IUC : unit number of the output C-File (in common /units/ of model.h)
C   ISCRN : unit number of the terminal (screen) output
C   IPRNT : unit number of the archive (print) output
c   OKMET : validity mask (binary) for met data
C   PRES  : pressure in mb
C   TEMP  : temperature in deg C
C   WETVAR: relative humidity (fraction of 1.0)    
c   SITECD: 4-character site code
c   L1Z,L1N,L1E,L2Z,L2N,L2E : L1/L2 phase center offsets
c   

      implicit none

      include '../includes/dimpar.h'
      include '../includes/model.h'


      integer*4 iyr,idoy,i,j,k

      real*8 antaz,site_epoch,pos(3)
     .     , decyrs,l1z,l1n,l1e,l2z,l2n,l2e
                          
      character*256 message

      l1z = offl1(1)
      l1n = offl1(2)
      l1e = offl1(3)
      l2z = offl2(1)
      l2n = offl2(2)
      l2e = offl2(3)

      nspare = 1
      spare(1) = 0.d0

      call dayjul( jdobs,iyr, idoy)
      site_epoch = decyrs(iyr,idoy,tobs)
      do i=1,3
        pos(i) = kpos(i) + kvel(i)*(site_epoch-kepoch0)
      enddo 

c     WVR data no longer supported and removed from c-file
C       OKWVR : validity flag for WVR delay (-1=missing, 0=OK)
C       WVRDEL: WVR delay for each satellite                 
c      do i=1,msat  
c       okwvr(i) = -1
c       wvrdel(i) = 0.d0
c      enddo

      call xyz2sph(pos,latr_sph,lonr,radius)
c     C-file units are km
      radius = radius/1000.d0   
cd      print *,'CFOUT msat nsave ndat npart nspare '
cd     .       , msat,nsave,ndat,npart,nspare      
cd      print *,'  rclock,zendel,atmlod '
cd     .       ,   rclock,zendel,atmlod 
      call writc4 (iuc
     .,            msat,mprn
     .,            iepoch,iyr,idoy,tobs,rclock
     .,            okmet,zendel,pres,temp,wetvar,atmlod       
     .,            sitecd,latr_sph,lonr,radius
     .,            l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,            nsave,save )   
                      
      do 100 i=1,msat
         j=1    
 110     if (mprn(i).eq.ischan(j)) then   
cd               print *,'CFOUT prn ndat obs omc '
cd     .        ,ischan(j),ndat,(obs(k,j),k=1,ndat),(omc(k,j),k=1,ndat)
cd              print *,'CFOUT iepoch prn elev azim '
cd     .             ,iepoch,ischan(j),elev(j),azim(j)
            call writc5 (iuc
     .,                  ischan(j)
     .,                  elev(j),elvdot(j)
     .,                  azim(j),azmdot(j),nadang(j)
     .,                  atmdel(j),svclock(1,j),svclock(2,j)
     .,                  delay(1,j),drate(j)
     .,                  ier(j),data_flag(j)
     .,                  ndat
     .,                  obs(1,j)
     .,                  omc(1,j)
     .,                  isnr(1,j)
     .,                  ampl1(j),ampl2(j)
     .,                  nspare,spare
     .,                  npart,tmpart(1,j) )
         else
            j=j+1
            if(j.gt.nchan) then       
               write(message,'(2a,i3,a)') 'Error writing c-file, '
     .           ,'too many channels (',j,')'
               call report_stat('FATAL','MODEL','cfout',' '
     .           ,message,0)
            endif
            goto 110
         endif
 100  continue
      return
      end

