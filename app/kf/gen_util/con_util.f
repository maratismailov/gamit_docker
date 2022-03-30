CTITLE CONOFF

      integer*4 function conoff( gnss )

      implicit none

*     Functiom to assign constellation offset depending on GNSS
*     system.  See svsp3.h
* PRNs are mapped with the additions of constants
* GPS     G    +0
* GLONASS R  +100
* GALIELO E  +200
* BEIDOU  C  +300
* QZSS    J  +400
* SBAS    S  +500
* IRNSS   I  +600

* PASSED VARIBALE 
      character*(*) gnss

*     Assign the offset value
      conoff = -100
      if( gnss(1:1).eq.'G' .or. gnss(1:1).eq.'g' ) then
          conoff = 0
      elseif( gnss(1:1).eq.'R' .or. gnss(1:1).eq.'r' ) then
          conoff = 100
      elseif( gnss(1:1).eq.'E' .or. gnss(1:1).eq.'e' ) then
          conoff = 200
      elseif( gnss(1:1).eq.'C' .or. gnss(1:1).eq.'c' ) then
          conoff = 300
      elseif( gnss(1:1).eq.'J' .or. gnss(1:1).eq.'j' ) then
          conoff = 400
      elseif( gnss(1:1).eq.'S' .or. gnss(1:1).eq.'s' ) then
          conoff = 500
      elseif( gnss(1:1).eq.'I' .or. gnss(1:1).eq.'i' ) then
          conoff = 600
      else
          write(*,220) trim(gnss)
 220      format('**WARNING** Unknown GNSS type ',a,' passed to CONOFF')
      endif

****  Thats all
      return
      end

CTITLE OFFCON

      character*1 function offcon( prn )

      implicit none

*     Functiom to assign constellation offset depending on GNSS
*     system.  See svsp3.h
* PRNs are mapped with the additions of constants
* GPS     G    +0
* GLONASS R  +100
* GALIELO E  +200
* BEIDOU  C  +300
* QZSS    J  +400
* SBAS    S  +500
* IRNSS   I  +600

* PASSED VARIBALE 
      integer*4 prn   ! PRN number with constellatin offset

*     Recover original PRN wiht mod(prn,100)
      offcon = '?'
      if( prn.gt.0 .and. prn.lt. 100 ) then
         offcon = 'G'
      elseif( prn.gt.100 .and. prn.lt. 200 ) then
         offcon = 'R'
      elseif( prn.gt.200 .and. prn.lt. 300 ) then
         offcon = 'E'
      elseif( prn.gt.300 .and. prn.lt. 400 ) then
         offcon = 'C'
      elseif( prn.gt.400 .and. prn.lt. 500 ) then
         offcon = 'J'
      elseif( prn.gt.500 .and. prn.lt. 600 ) then
         offcon = 'S'
      elseif( prn.gt.600 .and. prn.lt. 700 ) then
         offcon = 'I'
      endif

****  Thats all
      return
      end

