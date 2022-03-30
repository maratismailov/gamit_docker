Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego, 1995. All rights reserved.

      Subroutine CDTRED( iidoy,lpart,allchan,mchan
     .                 , mstart,mstop,idec,ipic )    

C     Read a data record from a C-file
C
C     MOD: MHM 870731  Error flagging from x-files
C          RWK 880607  Add DIMPAR (to allow 7 sats and pseudo-range)
C          MHM 890610  Add option to output data when called from PLOTC
c          RWK 890712  Use library routines to read new C-file format
c          SW  901218  Fix output format for O's and OMC's (for ASCII dump)
c          RWK 141205  Use commons for all C-file variables
c
      implicit none
c
      include '../includes/dimpar.h'  
      include '../includes/units.h'
      include '../includes/model.h'

      CHARACTER*256 MESSAGE


      INTEGER*4  iidoy,mchan,mstart,mstop,idec,ipic

* MOD TAH 101111: Added PRN number to output
      integer*4 prnout(maxsat),iyr,im,id,ihr,min,iprn,ndats,nparts
     .        , is,i,j

      real*8 sod,cl1,cl2,sec,atmtau,antaz,relhum,tau(maxsat)
     .     , radkm, l1z,l1n,l1e,l2z,l2n,l2e

      LOGICAL LGOOD,lpart,allchan

c     Function
      integer*4 julday  

cd DEBUG 
cd      integer*4 iflag,ioerr

C----------------------------------------------------------------------


C        Initialize variables

      do j=1,maxsat
         ier(j) = 1  
         data_flag(j) = 0 
         azim(j) = 0.d0
         elev(j) = 0.d0
         elvdot(j) = 0.d0
         azmdot(j) = 0.d0 
         atmdel(j) = 0.d0
         tau(j) = 0.d0
         drate(j) = 0.d0
         ampl1(j) = 0.d0
         ampl2(j) = 0.d0
         do i=1,maxdat
            isnr(i,j) = 0
            obs(i,j) = 0.D0
         enddo
         do i=1,maxlab
            tmpart(i,j) = 0.d0
         enddo   
      enddo
      do i=1,maxsav
        save(i) = 0.d0
      enddo
           
C     Read the Data from a C-File

      call readc4 (iuc
     .,            msat,mprn
     .,            iepoch,iyr,iidoy,tobs
     .,            rclock,okmet,zendel,pres,temp,relhum,atmlod
     .,            sitecd,latr_sph,lonr,radkm
     .,            l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,            nsave,save  )        
cd      print *,'CDTRED nepoch iepoch ',nepoch,iepoch 

c     Convert the observation time to PEP JD model.h used for writing the x-file 
c     and get HMS from seconds-of-day for dumping the record

      call monday(iidoy,im,id,iyr)
      jdobs = julday(im,id,iyr)      
      call ds2hms( iyr,iidoy,tobs,ihr,min,sec )
      
c     This code no longer needed since kinematic and dynamic not supporte
c     model.h expects radius in meters
c     radius = radkm*1.d3
c     new variable names in common model.h
c      offl1(1)  =  l1z 
c      offl1(2)  =  l1n 
c      offl1(3)  =  l1e 
c      offl2(1)  =  l2z 
c      offl2(2)  =  l2n 
c      offl2(3)  =  l2e 
    
c      read the data, one sat at a time, put into ischan order
c      note that only ier,obs,isnr, and ampl1/l2 are passed back

       do 50 is = 1, msat
         i = 1
   40    if (mprn(is).eq.ischan(i)) then
           call readc5 (iuc
     .,                 iprn
     .,                 elev(i),elvdot(i),azim(i),azmdot(i),nadang(i)
     .,                 atmdel(i),svclock(1,i),svclock(2,i) 
     .,                 tau(i),drate(i)
     .,                 ier(i),data_flag(i)
     .,                 ndats,obs(1,i),omc(1,i),isnr(1,i)
     .,                   ampl1(i),ampl2(i)
     .,                 nspare,spare
     .,                 nparts,tmpart(1,i) )   
*          MOD TAH 101111; Save PRN number for output
           prnout(i) = iprn

cd         write(6,41) is,i,iprn,ier(i),data_flag(i),obs(1,i)
cd  41     format(' CDTRED is,i,iprn,ier,data_flag,nparts,obs:'
cd    .                   ,5i5,e15.7)
cd         if( iepoch.ge.14 ) stop
           if (ndats.ne.ndat .or. iprn.ne.ischan(i)) then
             write(message,'(a,4i5)')
     .          'Error, ndat or iprn missmatch '
     .            ,ndats,ndat,iprn,ischan(i)
             call report_stat('FATAL','CTOX','cdtred',' ',message,0)
           endif
         else
          i=i+1
          if (i.gt.nchan) then
            write(message,'(a,5i5)')
     .          'Error, msat missmatch ',i,nchan,mprn(is),is,msat
                call report_stat('FATAL','CTOX','cdtred',' ',message,0)
          endif
          goto 40
        endif
   50 continue


      IF( iprnt.ne.6 .and.
     .    iepoch.ge.mstart.and.iepoch.le.mstop .and.
     .    idec.eq.ipic ) then

         call ds2hms( iyr,jdobs,tobs,ihr,min,sec )
         write(iprnt,60) iepoch,iyr,jdobs,tobs,ihr,min,sec
   60    FORMAT(/,1X,'Epoch ',I4,3x,i4,i9,f15.7,2x,2i3,f7.2)
         rclock= rclock*1.d6
         write(iprnt,62) rclock
   62    format(3x,'Receiver clock (microsec):',f14.3)
         write(iprnt,64) nsave,(save(i),i=1,nsave)
   64    format(3x,'nsave, save: ',i3,20d16.7)
         write(iprnt,67) (atmlod(i),i=1,3) 
   67    format(3x,'Atm loading (N E U, mm) ',3f10.1)
         write(iprnt,66) okmet,pres,temp,relhum
   66    format(3x,'Met avail:',i3,'  P,T,H: ',3f8.2,/
     .         ,3x,'Elev         Elvdot      Azm          Azmdot
     .  Total atm delay (ns)     PRN')       
         if( nchan.gt.maxsat ) then
            write(message,'(a,i3,a,i3,a)') 'Number of channels = ('
     .                     ,nchan,') > maxsat (=',maxsat,')'
            call report_stat('FATAL','CTOX','cdtred',' ',message,0)
         endif
         do  i=1,nchan
*          MOD TAH 101111: Only output if channel has data in it.
*          if( (allchan.or.i.eq.mchan) .and. data_flag(i).ne.0 ) then
*          MOD TAH 140325: Changed test on data existance to ier=1
*          (data flag is not set until autcln run).
           if( (allchan.or.i.eq.mchan) .and. ier(i).ne. 1 ) then
             atmtau= atmdel(i)*1.d9
             write(iprnt,68) elev(i),elvdot(i),azim(i),azmdot(i)
     .                     , atmtau, prnout(i)
   68        format(3x,4d12.5,4x,f12.2,15x,I3)
c             write(iprnt,69) nspare,(spare(j),j=1,nspare)
c   69        format(3x,'nspare, spare:',i3,8d16.6)
           endif
         enddo
         DO 90 I=1,NCHAN
*         MOD TAH 101111: Only output if channel has data in it.
*         if( (allchan.or.i.eq.mchan) .and. data_flag(i).ne.0 ) then
*         MOD TAH 140325: Changed test on data existance to ier=1
*         (data flag is not set until autcln run).
          if( (allchan.or.i.eq.mchan) .and. ier(i).ne. 1 ) then
            CL1=  OBS(1,I) - OMC(1,I)
            CL2=  OBS(2,I) - OMC(2,I)
            WRITE(IPRNT,70) I,IER(I),data_flag(i), prnout(i)
     .                                        ,OBS(1,I),OMC(1,I),CL1
     1                                        ,OBS(3,I),OMC(3,I)
     2                                        ,OBS(2,I),OMC(2,I),CL2
     3                                        ,OBS(4,I),OMC(4,I)
   70       FORMAT(2X,I2,'  IER=',I2,' dflg=',o12,1x,'PN= ',i3,/,
     .          '   Data     obs              o-c          calc',/
     1            ,3x,'L1    ',3(1X,F14.3),/
     2            ,3X,'P1    ',2(1X,F14.3),/
     3            ,3X,'L2    ',3(1X,F14.3),/
     4            ,3X,'P2    ',2(1X,F14.3) )  
            write(iprnt,74) ampl1(i),ampl2(i)
   74    format(1x,'  Ampl1: ',e12.2,'  Ampl2: ',e12.2)
            write(iprnt,75) tau(i),drate(i),svclock(1,i),svclock(2,i)
   75    format(1x,'  Delay: ',f16.14,'    Delay rate: ',d12.5
     .          ,'    SV offset: ',d12.5,'    SV L1 corr: ',f14.2   )
            if( lpart) WRITE(IPRNT,80) nparts,(TMPART(J,I),J=1,NPARTS)
   80       FORMAT(1X,'  Partials: ',i2,(5D12.5))
           endif
   90       CONTINUE                            
      endif

      RETURN
      END
