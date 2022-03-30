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
      if( gnss(1:1).eq.'G' .or. nss(1:1).eq.'g' ) then
          conoff = 0
      elseif( gnss(1:1).eq.'R' .or. nss(1:1).eq.'r' ) then
          conoff = 100
      elseif( gnss(1:1).eq.'E' .or. nss(1:1).eq.'e' ) then
          conoff = 200
      elseif( gnss(1:1).eq.'C' .or. nss(1:1).eq.'c' ) then
          conoff = 300
      elseif( gnss(1:1).eq.'J' .or. nss(1:1).eq.'j' ) then
          conoff = 400
      elseif( gnss(1:1).eq.'S' .or. nss(1:1).eq.'s' ) then
          conoff = 500
      elseif( gnss(1:1).eq.'I' .or. nss(1:1).eq.'i' ) then
          conoff = 600
      else
          write(*,220) trim(gnss)
 220      format('**WARNING** Unknown GNSS type ',a,' passed to CONOFF')
      endif

****  Thats all
      return
      end


