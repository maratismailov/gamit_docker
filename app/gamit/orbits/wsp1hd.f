c
      subroutine wsp1hd( iungs,iscrn,mjd,fmjd,delt,nepoch,numsv,issat
     .                 , org,orbtyp,icrd )
c
c      R.W. King   November 1991

c      Write the header records for NGS Standard Product #1 orbit format
c      Reference:  B.W. Remondi, "Extending the National Geodetic Survey
c                  Standard GPS Orbit Formats", NOAA Tech. Rep. NOS 133 NGS 46
c                  Rockville, MD, November 1989.

c        iungs   Unit number for output NGS-format file
c        mjd     Modified Julian Day of start
c        fmjd    Fractional day (GPST, from midnight)
c        delt    Epoch interval
c        nepoch  Number of epochs
c        numsv   Number of PRNs
c        idsv    PRN numbers written to NGS-format file (dimension 34)
c        issat   PRN numbers for T-file (dimension maxsat)
c        orbtyp  Orbit Type (F=fitted, E=extrapolated or predicted
c                            B=Broadcast)
c        org     Agency source of ephemeris
c        icrd    Coordinate System of ephemeris
c                    72 = WGS72
c                    84 = WGS84
c                    85 = Earth-fixed 1985 (IERS)
c
c
      implicit none
c
      include '../includes/dimpar.h'
c
      character*1 orbtyp
      character*3 org
      integer*4 iungs, iscrn, mjd, iyr, iday, imon, ihr,imin
     .        , nepoch, numsv, idsv, issat, igpswk, icrd, itflag
     .        , jd, idoy, i
      real*8 fmjd, delt, sec, sod, sow, utcoff

      dimension issat(maxsat),idsv(34)

c  Get Gregorian date and GPS week from MJD, fract

      jd = mjd + 2400001
      call dayjul ( jd,iyr,idoy )
      call monday ( idoy,imon,iday,iyr )
      sod = fmjd*86400.d0
c     avoid roundoff by assuming even second
      sod = anint(sod)
      call ds2hms( iyr,idoy,sod,ihr,imin,sec )
c     determine gps week (input and output is GPST)
      itflag = -4
      call timcon( itflag,igpswk,sow,iyr,idoy,ihr,imin,sec,utcoff )
c
c  Write the first record of the SP1 file header
c
      write(iungs,10) iyr,imon,iday,ihr,imin,sec,delt,mjd,fmjd,nepoch
     .          , orbtyp,org
   10 format(1x,'#  ',i4,4i3,1x,f10.7,1x,f14.7,1x,i5,1x,f15.14
     .      ,1x,i6,a1,1x,a3)
c
c  Fill the IDSV array

      do 15 i=1,34
   15 idsv(i)= 0
      do 16 i=1,numsv
   16 idsv(i) = issat(i)

c  Write the second record of the SP1 file header
c
      write(iungs,20) numsv,(idsv(i),i=1,34),icrd,igpswk
   20 format(1x,'+ ',i2,1x,35i2,1x,i3)
c
c  Echo the SP-1 file headers
c
      write(iscrn,30) orbtyp,icrd,org,igpswk
     1          , iyr,imon,iday,ihr,imin,sec,mjd,fmjd,nepoch,delt
     1          , numsv,(idsv(i),i=1,numsv)
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
   30 format(//
     1      ,1x,'Header written on NGS Standard Product #1 file:',//
     2      ,1x,'Orbit Type: ',a1,'   Coordinate System : ',i2
     3      , '  Organization : ',a3,'   GPS Week : ',i3,/
     4      ,1x,'Start epoch GPST :',i4,4i3,1x,f10.7,/
     5      ,1x,'            MJD  : ',i5,1x,f15.14,/
     6      ,1x,'Number epochs    : ',i6,'  Interval :',f8.2,' sec',//
     7      ,1x,i2,' satellites, PRN #s = ',50i3)
c
      return
      end

