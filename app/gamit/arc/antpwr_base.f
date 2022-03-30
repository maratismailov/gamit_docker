CTITLE ANTPWR_BASE

      subroutine antpwr_base(antbody, newsv, antpwr )

      implicit none

*     Base routine to return transmitter when the value has not 
*     been read from the igs_metadata.snx file.
*     Routine is called when ANTBK antradmod used and svnav_read
*     read svnav.dat (svnav_read sets power to zero when not found).
*
*     Herring July 2, 2019.

* PASSED VARIABLES
      character*(*) antbody  ! Type of satellite body (IN)
      logical newsv          ! True is this a new satellite. For output (IN)
      real*8 antpwr          ! Transmit power in Watts.
                             ! 1 watt = 1 kg m^2/s^3

* LOCAL
      character*128 message  ! String for report_stat call

****  Look at antbody to see what power to use.
c     Block specific antenna thrust (see http://acc.igs.org/orbits/thrust-power.txt)
c     Compares to Carlos Rodriguez Solano implementation
c     Use Block IIA value for Block I and II
      if( antbody(1:9).eq.'BLOCK I  ' ) then
         antpwr = 76.d0
      elseif( antbody(1:9).eq.'BLOCK II '  .or.
     .        antbody(1:9).eq.'BLOCK IIA' ) then 
         antpwr = 76.d0
      elseif ( antbody(1:11).eq.'BLOCK IIR-A' .or.
     .         antbody(1:11).eq.'BLOCK IIR-B' ) then
         antpwr = 85.d0
      elseif( antbody(1:11).eq.'BLOCK IIR-M') then
c        IIR-M: this includes the M-code power which is nominally switched on 
c        (would be 108W without M code)           
         antpwr = 198.d0 
      elseif( antbody(1:9).eq.'BLOCK IIF' ) then   
c        M-code was switched off for SVN62/PRN25 on 05Apr2011 JD 2455656.5 PJD 2455657) 
         antpwr = 249.d0       
c        (would be 154W without M code)
*        Special code for specific satellite (should be handled in igs_metadata.snx). 
* G062 2010:148:00000 0000:000:00000  240  GPS-IIF; [TP01]; block mean
*        This case is not in igs_metadata.snx.  Not in igs_metadata.snx so don't 
*        implement.
C        if( sname(7:8).eq.'62' .and.fjd.ge.2455657.d0) then
C           antpwr = 154.d0
C        endif     
      elseif( antbody(1:10).eq.'BLOCK IIIA' ) then
**** This is a preliminary estimate of the power based on an IGS AC email
**** from Peter Steigenberg 190104. ---rwk
* G074 2018:357:00000 0000:000:00000  300  GPS-III; guess
         antpwr = 300.d0 
      elseif( antbody(1:9).eq.'GLONASS-M' ) then
         antpwr = 100.d0 
      elseif( antbody(1:9).eq.'GLONASS-K' ) then
* Use R802 as the later satellite (R801 has 135W but seems special).
         antpwr = 105.d0 
      else                 
         if(newsv) then 
            antpwr = 0.d0                       
            write(message,'(a,a20)') 
     .          'Antenna thrust model not available for ',antbody
            call report_stat('WARNING','ARC','ertorb',' ',message,0)
         endif 
      endif

****  Thats all
      return
      end 
