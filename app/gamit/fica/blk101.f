      subroutine blk101
     &   ( debug,iflag,apress,atemp,ahumid,asampl,swveru
     &   , auser,asite,arcvr,antht
     &   , ff,fi,fc,nf,ni,nc )

c     Read a FICA BLK 101
c
c                  GESAR/ASDAP/CORE-ARL:UT FIC Block Definition   -   02-03-88
c
c    FIC Block Name         : GESAR Versions 1.0+ Input Data
c    FIC Block Number       : 101
c    Original Block Number  : 1  (0516 Hexadecimal)
c    Original Block Source  : GESAR
c
c    Number of Floating Point Items : 5
c    Number of Integer Items        : 230
c    Number of Character Items      : 14
c
c    Item Numbers           Item Description                        Units
c    _______________        _________________________________       ______________
c
c    Floating Point Items
c
c    1                            Pressure                          mbar
c    2                            Temperature                       deg C
c    3                            Humidity                          percent
c    4                            Collection rate                   sec
c    5                            Solution rate                     sec
c
c    Integer Items
c    1                            Software version *10
c    2                            First satellite PRN number
c    3                            Doppler aiding value              Hz
c    4-8                          Current date and time             yr,mo,day,hr,min
c    9-80                         Date of scenarios (24 scenarios)  yr,mo,day
c    81-128                       Time of scenarios (24 scenarios)  hr,min
c    129-133                      Scenario end date and time        yr,mo,day,hr,min
c    134-229                      Satellite PRN numbers of
c                                  Scenarios to track
c                                  (SV PRNs,scenarios)
c    230                          External oscillator flag          1=Y/0=N
c
c    Character Items
c    1                            User code
c    2-3                          Receiver serial number
c    4-5                          Site name
c    6-7                          Latitude of initial position      deg
c    8-9                          Longitude of initial position     deg
c    10-11                        Height of initial position        meters
c    12-13                        Antenna height above position     meters
c    14                           Search type                       A=Aided, B=Blind, C=Almanac

c     passed values
C     .true. to print things
      logical     debug
C     1 for data, 2 for ephemeris
      integer*4   iflag

c     scenario arrays
C     24 scenarios + end
C     yr,mo,dy,hr,mn,prn1,prn2,prn3,prn4
      integer*4   iscen(25,9)
c

c     values from receiver operator
      character*8 auser
C     user-input site code
      character*16 asite
C     reciever serial number
      character*16 arcvr
C     operator input antenna ht.
      character*16 antht
C     atmospheric pressure
      real*8       apress
C     atmospheric temperature
      real*8       atemp
C     atmospheric pressure
      real*8       ahumid
c     collection interval
      real*8       asampl
C     software version
      real*8       swveru

      include '../includes/makex.h'

c     FICA ARRAYS
c     the FICA arrays
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
C     number of elements in FICA arrays
      integer*4 nf,ni,nc, mni

      integer*4 i,j,icount
      
c     What is icount for ?  --rwk 970313
      save icount

      data icount/0/
             
      if( debug ) then
        print *,'BLK101 nf,ff : ',nf,(ff(i),i=1,nf)
        print *,'BLK101 ni,fi : ',Mni,(fi(i),i=1,ni)
        print *,'BLK101 nc,fc : ',nc,(fc(i),i=1,nc)
      endif
      icount = icount + 1

      apress = ff(1)
      atemp  = ff(2)
      ahumid = ff(3)
      asampl = ff(4)

      swveru = fi(1)
      swveru = swveru/10.

      auser  = fc(1)
      arcvr  = fc(2)//fc(3)
      asite  = fc(4)//fc(5)
      antht  = fc(12)//fc(13)

c     read in the scenarios

c     yr,mo,dy
      j =  9
      do 10 i=1,24
         iscen(i,1) = fi(j)
         iscen(i,2) = fi(j+1)
         iscen(i,3) = fi(j+2)
         j = j + 3
 10   continue

c     hr,mn
      j = 81
      do 20 i=1,24
         iscen(i,4) = fi(j)
         iscen(i,5) = fi(j+1)
         j = j + 2
 20   continue

c     PRNS
      j = 134
      do 30 i=1,24
         iscen(i,6) = fi(j)
         iscen(i,7) = fi(j+1)
         iscen(i,8) = fi(j+2)
         iscen(i,9) = fi(j+3)
         j = j + 4
  30  continue

c     end yr,mo,dy,hr,mn
      j = 1
      do 40 i=129,133
         iscen(25,j) = fi(i)
         j = j + 1
  40  continue
           
c      if (icount .gt. 2) then
       if( debug ) then
         write (uinfor,113)
         write (uinfor,114) asite
         write (uinfor,115) auser
         write (uinfor,116) arcvr
         write (uinfor,120) antht
         write (uinfor,117) ahumid
         write (uinfor,118) atemp
         write (uinfor,119) apress
         write (uinfor,121) swveru

         write (uscren,113)
         write (uscren,114) asite
         write (uscren,115) auser
         write (uscren,116) arcvr
         write (uscren,120) antht
         write (uscren,117) ahumid
         write (uscren,118) atemp
         write (uscren,119) apress
         write (uscren,121) swveru

 113     format (/,1x,'INFORMATION FROM FICA BLK101:')
 114     format (/,1x,'Operator input site code            ',a16)
 115     format (1x,'Operator initials                   ',a16)
 116     format (1x,'Operator input receiver serial num  ',a16)
 120     format (1x,'Operator input antenna height (m)   ',a16)
 117     format (1x,'Operator input humidity (%)         ',f16.4)
 118     format (1x,'Operator input temperature (deg C)  ',f16.4)
 119     format (1x,'Operator input pressure (Mb)        ',f16.4)
 121     format (1x,'Software version number             ',f16.4)

         write (uscren,148)
         write (uinfor,148)
 148     format (/,1x,'Scenarios: PRN')
         do 200 i=1,24
            if (iscen(i,1) .ne. 0) then
              write (uscren,150) i,(iscen(i,j),j=1,9)
              write (uinfor,150) i,(iscen(i,j),j=1,9)
 150          format (1x,'#',i2,1x,i2.2,'/',i2.2,'/',i2.2,1x,i2.2,':',
     .        i2.2,4(1x,i2))
            endif
 200     continue
         write (uscren,154) (iscen(25,j),j=1,5)
         write (uinfor,154) (iscen(25,j),j=1,5)
 154     format (1x,'END',1x,i2.2,'/',i2.2,'/',i2.2,1x,i2.2,':',i2.2)

      endif

c     write (uscren,158)
c     write (uinfor,158)
c158  format (/,/)

      iflag = 0

      return
      end


