      Subroutine rsp1hd( iungs,iyear,imon,iday,ihr,imin,sec,delt
     .                 , mjd,fmjd,nepoch,numsv,issat )
c
c      R.W. King   March 1988

c      Modified 11 November 1991

c      Read the header records for NGS Standard Product #1 orbit format
c      Reference:  B.W. Remondi, "Extending the National Geodetic Survey
c                  Standard GPS Orbit Formats", NOAA Tech. Rep. NOS 133 NGS 46
c                  Rockville, MD, November 1989.

c        mjd     Modified Julian Day of start
c        fmjd    Fractional day (from midnight)
c        delt    Epoch interval
c        nepoch  Number of epochs
c        numsv   Number of PRNs
c        idsv    PRN numbers from header (dimension 34)
c        issat   PRN numbers for T-file (dimension maxsat)
c        orbtyp  Orbit Type (F=fitted, E=extrapolated or predicted
c                            B=Broacast)
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

      character*1 orbtyp
      character*4 org
      character*256 message,line
      integer*4 iungs,mjd,jd,julday,iyear,iday,imon,ihr,imin,nepoch
     .       , numsv,idsv,issat,igpswk,icrd,i
      real*8 fmjd,fmjdt,delt,ts,sec

      logical debug/.false./

      dimension idsv(34),issat(maxsp3sv)
c

c  Read the first record of the SP1 file header
c                                
      org = ' ' 
      read(15,10) iyear,imon,iday,ihr,imin,sec,delt,mjd,fmjd,nepoch
     .          , orbtyp,org(1:3)  
      call check_y2k(iyear)
   10 format(4x,i4,4i3,1x,f10.7,1x,f14.7,1x,i5,1x,f15.14,1x,i6,a1,1x,a3)
      
c  these removed, along with function FJDTIM, by rwk 990728
c      T(1)=DBLE(IHR)
c      T(2)=DBLE(IMIN)
c      T(3)=SEC
c      IT(1)=IMON
c      IT(2)=IDAY
c      IT(3)=IYEAR
c      CALL FJDTIM(T,IT,FJD,JDN,TS)
c and replaced by the following
      ts =(dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec)/86400.d0
      jd = julday(imon,iday,iyear)

      FMJDT=DBLE(jd)+TS-2400001.D0
C Replace value off header with more precise computed value
      FMJD=TS
c  Read the second record of the SP1 file header
c           
      read(iungs,'(a)') line
      read(line,'(3x,i2)') numsv
      if( numsv.gt.maxsat ) then
          write(message,'(a,i2,a,i2,a)') 
     .        'Number of satellites on SP1 file (',numsv,') > maxsat ('
     .        ,maxsat,')'
        call report_stat('FATAL','NGSTOT','orbits/rsp1hd',' ',message,0)
      endif                                                         
      backspace(iungs)  
      read(15,20) numsv,(idsv(i),i=1,34),icrd,igpswk
   20 format(3x,i2,1x,35i2,1x,i3)
c
c  Echo the SP-1 file headers
       
      if( debug) then
        write(6,30) orbtyp,icrd,org,igpswk
     1          , iyear,imon,iday,ihr,imin,sec,mjd,fmjd,nepoch,delt
     1          , numsv,(idsv(i),i=1,numsv)
   30 format(//,1x,'Header records from NGS Standard Product #1 file:'
     1      , //
     2      ,1x,'Orbit Type: ',a1,'   Coordinate System : ',i2
     3      , '  Organization : ',a3,'   GPS Week : ',i4,/
     4      ,1x,'Start epoch (GPST): ',i4,4i3,1x,f10.7,/
     5      ,1x,'            MJD   : ',i5,1x,f15.14,/
     6      ,1x,'Number epochs     : ',i6,'  Interval :',f7.2,' sec',//
     7      ,1x,i2,' satellites, PRN #s = ',12i3)
       endif

c  Fill the T-file SV array

      do  40 i=1,numsv
  40  issat(i)=idsv(i)
c
      return
      end

