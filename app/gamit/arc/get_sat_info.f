Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995/2012.  All rights reserved.

      Subroutine get_sat_info(isat)

c     Read svnav.dat to get the satellite name, mass, and yaw rate and bias
c     R. King - May 2012 - based on R. I. Abbot routine redsat (now obsolete0
 
c     Input parameter: index in satellite array of the satellite being integrated
c     Output values stored in common /raditation/ in arc.h

      implicit none

      include '../includes/dimpar.h' 
      include '../includes/global.h'
      include '../includes/arc.h'

      character*256 message

      integer*4 idms,isat,iyr,imo,iday,ihr,imin
     .        , idir,isvn,idoy,frqchn,svnstart(5),svnstop(5),i,ioerr

      real*8 fomeg0,asat,tics,anomin,wper0,eccen,fincl,temp1,temp2,sec
* MOD TAH 200126: Removed line below.  Not needed
!     real*8 svantdx(3)

      dimension idms(4,4)
      dimension tics(6)

c Get the GNSS and PRN number from the satellite name         

cd      print *,'GET_SAT_INFO isat satnam ',isat,satnam(isat)
      if( satnam(isat)(1:3).eq.'PRN' ) then
         read( satnam(isat)(5:6),'(i2)') iprn                    
         gnss = 'G'
cd         print *,'  gnss iprn ',gnss,iprn     
      else   
         read( satnam(isat)(1:1),'(a1)') gnss
         read( satnam(isat)(2:3),'(i2)') iprn
      endif


c Get the satellite information from the svnav.dat file
c  --call with the IC epoch to avoid PRN changes over day boundaries
c  --year needed to distinguish among reused prn #s  
      call dayjul(jde,iyr,idoy)   
      call monday(idoy,imo,iday,iyr)
      call ds2hms(iyr,idoy,te,ihr,imin,sec)
      idir = -1  
* MOD TAH 190702: Added antpwr to snav_read call
      call svnav_read( idir,iyr,idoy,ihr,imin,gnss,iprn
     .               , isvn,frqchn,antbody,sbmass,bias,yawrate, antpwr
     .               , svnstart,svnstop )
c     convert mass from g to kg
      sbmass = sbmass*1.d-3
      if( antbody(1:1).eq.' ') then
c        blank antenna/body-type means no valid entry in svnav.dat
c        we've made this a fatal error now to avoid subtle errors
c       in the radiation pressure model and (later) the yaw
        write(message,'(a,i3,a)') 'No valid entry for PRN',iprn
     .       ,' in svnav.dat--check #s and dates' 
        write(iarh,'(/,1x,2a)') '**',message
        call report_stat('FATAL','ARC','get_sat_info',' ',message,0)
      endif

c  Write out some information to the output file

      write (iarh,50) iprn
   50 format (/,'******************************* PRN ',i2,' ************
     .**********************')  
      write(iarh,51) iprn,isvn,antbody,sbmass  
   51 format(/,1x,'From svnav.dat:  PRN ',i3,'  SV ',i3
     .  ,'  Ant/Body-type ',a16,' Mass (kg) ',f10.4)
      write (iarh,54) satnam(isat),(satics(i),i=1,nics)
   54 format(/,1x,'Satellite and ICs: ',a16,7(/,3(1x,d22.15)))
      write (iarh,'(a80)')
* MOD TAH 190722: Save the satellite SVN number for use in forces models.
      satsvn(isat) = isvn
      ssvn = isvn   ! Saved value for satellite being integrated.

c Compute keplerian elements

      do  i=1,6
        tics(i)=satics(i)
      enddo
c     check for zero values
      temp1 = tics(1)**2 + tics(2)**2 + tics(3)**2
      temp2 = tics(4)**2 + tics(3)**3 + tics(4)**2
      if( temp1.eq.0.d0 .or. temp2.eq.0.d0 ) then
        call report_stat('FATAL','ARC','get_sat_info',' ',
     .  'Error, satellite position or velocity zero from g-file',0)
      endif     
      call keplr(tics,asat,eccen,wper0,fomeg0,fincl,anomin)    
c     print the elements
      call hofrad(wper0 ,idms(1,1))
      call hofrad(fincl ,idms(1,2))
      call hofrad(fomeg0,idms(1,3))
      call hofrad(anomin,idms(1,4))
      write(iarh,'(2a)') " Sat     Semimjr Ax  Eccen'ty   Perigee      "
     . ,"Inclin'n    Ascen.Node   Mn Anomaly"
      write (iarh,75) satnam(isat)(1:6),asat,eccen,idms
   75 format(1x,a6,1x,f10.3,f11.6,4(1x,i4,i3,i3,'.',i1))
c     bounds check for a, e, and i
c     if (asat.lt.26563.d0.and.asat.gt.26557.d0) go to 80
      if (asat.lt.26585.d0.and.asat.gt.26557.d0) go to 80
      go to 90
   80 if (eccen.lt.0.015) go to 85
      go to 90
   85 if (fincl.lt.(64.d0*cdr).and.fincl.gt.(62.8d0*cdr)) go to 95
c                                  fincl.gt.(63.d0*cdr)) go to 95
   90 continue
   95 continue

      return
      end

