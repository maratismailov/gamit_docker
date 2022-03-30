C Calculate forces from earth radiation using TUM models 
C See http://www.iapg.bv.tum.de/albedo/ and http://acc.igs.org/reprocess2.html (7th Sept 2012)
C Called by ertorb.f
C Calls ERPFBOXW
C EJP Sept 2012 
C Coded to keep TUM code as unchanged as possible in case of further releases.

* MOD TAH 200603: Accounted for M+, K1A and K1B satelites (treated same as M and K1)

      Subroutine earthradTUM(erm,grd,fjd,accel)

C Inputs: 
C		erm: select earth radiation model
C       grd: size of grid: 1 = 2.5 degree grid, 2 = 5.0 degree grid
C                          3 = 10  degree grid (affects computation speed)

C
C Outputs:
C		TUMEvec(3):  accelerations in kms-2 


      implicit none

      include '../includes/dimpar.h'  
      include '../includes/units.h'  
      include '../includes/global.h'
      include '../includes/arc.h'
      
      integer*4 erm,grd,ant,reff,iyr,idoy,imonth,iday,ihr,imin,blknum
     .           , svn,i,iut1pol
      real*8 fjd,mjd,ysat(6),sun(3), rot(3,3),rotdot(3,3),sidtm
      real*8 accel(3)

c     Numerical GPS block code now temporary for UCL routines
      integer*4 iblock
c     Function (also temporary)
      integer*4 gpsblock

      integer*4 cnt, saved_svn 

      character*128 message  ! Message about missing block

      save cnt, saved_svn 

      data saved_svn / 0 / , cnt / 1 / 

      if( ssvn .ne. saved_svn ) then
          cnt = 1
          saved_svn = ssvn 
      endif

C Setup  
C     Set ant to 0 : don't calculate antenna thrust
* COM TAH 200526: Do not change ant=0 value.  Antenna thrust
*     is computed in antpwr from igs_metadata.snx file (if linked)
*     or antpwr_base.f in ertorb.f 
      ant=0    
C     Acceleration in reference frame  (0: inertial, 1: body fixed, 
C      2: sunfixed, 3: orbital-radial,along,crosstrack)
      reff = 0    
C     Satellite posn and vel in m (relative to earth
      do i=1,6
         ysat(i)=sbcor(i)*1000.d0             
      enddo
C     Sun posn vector in m
      do i=1,3
         sun(i)=-ccor(i)*1000.d0
      enddo
      
      imonth=0
      if (erm.eq.2) then
C        calculate the inertial to earthfixed rotation matrix
C*        iut1pol should bet set by the sestbl. according to rwk, but is not in arc.h
* MOD TAH 200505: Read the sestbl. to get the correct value   
C        iut1pol = 11    ! original code
         call get_iut1pol( iut1pol )

         call rotsnp(-1,jde,te,tdtoff,iut1pol,rot,rotdot,sidtm
     .                 , iut1,ipole,inut,frame,precmod)
C        calculate the calendar month for the CERES data  
         call pjdhms(fJD,iyr,idoy,imonth,iday,ihr,imin)
      endif  
c     Retain block numbers locally for GPS temporarily
C     iblock =  gpsblock(antbody)
* MOD TAH 190702: Replaced gpsblock with GNSS name_to_blk  
      
      call name_to_blk( +1, antbody, iblock )

C     Satellite block number (GAMIT has only IIR-A and IIR-B so the 
C     block numbers don't match perfectly.
C MOD TAH 200603 Added new mapping for M+, K1A and K1B
C          TUM code uses:                            MITBlk
CC                             1 = GPS-I          ==  1     
CC                             2 = GPS-II         ==  2
CC                             3 = GPS-IIA        ==  3
CC                             4 = GPS-IIR        ==  Does not exist
CC                             5 = GPS-IIR-A      ==  4
CC                             6 = GPS-IIR-B      ==  5
CC                             7 = GPS-IIR-M      ==  6
CC                             8 = GPS-IIF        ==  7
CC                             9 = GPS-IIIA       ==  8
CC                           101 = GLONASS        == 11
CC                           102 = GLONASS-M      == 12
CC                           102 = GLONASS-M+     == 13  ! Added 200603
CC                           103 = GLONASS-K1     == 14 
CC                           103 = GLONASS-K1A    == 15  ! Added 200603
CC                           103 = GLONASS-K1B    == 16  ! Added 200603
CC                           201 = GALILEO (IOV)  == 23 (GALILEO-1)
CC                           203 = GALILEO (JOV)  == 24 (GALILEO-2)
CC    GPS-IIIA code added based on GFZ acc_albedo_propboxw.f
CC    MITBlk 21, 22 are GALILEO GIOVE-A/B 
      if (iblock .le.3 ) then
          blknum=iblock
      elseif (iblock.ge.4 .and.iblock.le.8) then
          blknum = iblock+1
      elseif( iblock.ge.11 .and. iblock.le.12 ) then
          blknum=(iblock-10) + 100
      elseif( iblock.eq.13 ) then
          blknum = 102
      elseif( iblock.ge.14 .and. iblock.le.16 ) then
          blknum = 103
      elseif( iblock.ge.23 .and. iblock.le.24 ) then
          blknum=(iblock-22) + 200
      else
          if( cnt.eq.1 ) then   ! Tell user there is issue
             cnt = cnt + 1
             write(message,150) ssvn, trim(antbody)
 150         format('No Albedo model for SVN ',I4,' Body ',a)
             call report_stat('WARNING','ARC','earthradTUM',' ',
     .                         message,0)
          endif
          accel = 0.d0  ! Make sure no force
          RETURN
      endif
      
C     SVN: only used for block specific antenna thrust for SVN62/PRN25, and 
C     are not using this TUM code for antenna thrust
C MOD TAH 19092: SVN needed for Galileo mass so a saved versions added to arc.h
C     and set in get_sat_info.
      SVN=ssvn

* MOD TAH 200526: JDE is JD+1 for PEP JD and we need MJD here; Remove 2400001.d0
*     Not needed because only used for ANTPOW 
      mjd=jde-2400001.d0
cd      print *,'earthradTUM jde mjd ',jde,mjd 
      
Cd      print*,erm,ant,grd,reff,imonth,blknum,svn,mjd,'iut1pol',iut1pol
Cd      print*,ysat
Cd      print*,sun
Cd      print*,'rot',rot
C  Call the TUM routine
cd      print *,'Calling ERPFBOXW erm ant grd reff ysat ',erm,ant,grd,reff
cd      print *,'ysat sun rot ',ysat,sun,rot
cd      print *,'month blknum svn mjd ',imonth,blknum,svn,mjd
      call ERPFBOXW(erm,ANT,GRD,REFF,YSAT,SUN,rot,imonth,
     .                    BLKNUM,SVN,MJD,ACCEL)
      if ( cnt.eq.1 ) then 
C         print *,'ERPFBOXW ',BLKNUM,SVN,MJD,ACCEL
          cnt = cnt + 1
      end if
C     convert the accelerations from ms-2 to kms-2     
      do i=1,3
         accel(i)=accel(i)/1000.d0
      enddo
cd      print *,'TUME1 Albedo ',accel
      return
      end
