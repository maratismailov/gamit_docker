      subroutine wsp3hd( iungs,iprnt,mjd,fmjd,delt,nepoch
     .                 , numsv,gnss,issat
     .                 , org,orbtyp,crdsys,pcvmod,otlmod, sample_rate )
c
c	P. Fang	March 1993, based on King's wsp1hd
c
c      Write the header records for NGS Standard Product #3 orbit format
c      Reference:  B.W. Remondi, "Extending the National Geodetic Survey
c                  Standard GPS Orbit Formats", NOAA Tech. Rep. NOS 133 NGS 46
c                  Rockville, MD, November 1989.
c  T.Herring modified to allow sample_rate to be passed. 210122

c        iungs   Unit number for output NGS-format file
c        mjd     Modified Julian Day of start
c        fmjd    Fractional day (GPST, from midnight)
c        delt    Epoch interval
c        nepoch  Number /of epochs
c        numsv   Number of PRNs
c        idsv    PRN numbers written to NGS-format file (dimension 85)
c        issat   PRN numbers for T-file (dimension maxsat)
c        orbtyp  Orbit Type (FIT=fitted, EXT=extrapolated or predicted
c                            BRD=Broadcast)
c        org     Agency source of ephemeris
c        crdsys  Coordinate System of ephemeris
c                    WGS72
c                    WGS84
c                    IER85 = Earth-fixed 1985 (IERS)
c                    ITR91
c
c
      implicit none
c
      include '../includes/dimpar.h'

      character*1 otlflg,atlflg,gnss
      character*3 orbtyp,asv,orbmod,clkmod
      character*4 org
      character*5 dused,crdsys
      character*8 otlmod,atlmod      
      character*10 pcvmod
      integer*4 iungs, iprnt, mjd, iyr, iday, imon, ihr,imin
     .        , nepoch, numsv, idsv, issat, igpswk, itflag, isatsgm
     .        , jd, idoy, i
* MOD TAH 210122: Introduced higher sampling rate factor (e.g. for 900-sec sampled
*     t-file to output at 300-seconds, sample_rate = 3
      integer*4 sample_rate   ! Sampling rate multiplier on t-file rate;
                              ! Value must >= 1; default is 1.
      real*8 fmjd, delt, sec, sod, sow, utcoff

      dimension issat(maxsat),idsv(85),isatsgm(85),asv(85)

c isatsgm() are the orbit accuracies, set to zeroes for the time being
      data isatsgm/85*0/

c variable dused has been hardwired as 2-rcvr/2=sat carrier phase
      dused='d    '

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

c  Fill the IDSV array

      do 15 i=1,85
   15 idsv(i)= 0
      do 16 i=1,numsv
   16 idsv(i) = issat(i)

c	Line 1
* MOD TAH 210122: Update the number of epochs for the sample_rate
      write(iungs,910) iyr,imon,iday,ihr,imin,sec,nepoch*sample_rate,
     *                 dused,crdsys,orbtyp,org
c* 910   format('#  ',i4,4i3,1x,f11.8,1x,i7,1x,a5,1x,a5,1x,a3,1x,a4)
c* Fang 010214 change for SP3-B; King 090908 change for SP3-C
910   format('#cP',i4,4i3,1x,f11.8,1x,i7,1x,a5,1x,a5,1x,a3,1x,a4)

c	Line line 2
* MOD TAH 210122: Update the number of epochs for the sample_rate
      write(iungs,920) igpswk,sow,delt/sample_rate,mjd,fmjd
920   format('## ',i4,1x,f15.8,1x,f14.8,1x,i5,1x,f15.13)

c	Lines 3-7
c** mods rwk 090908 to accommodate SP3-C
      do i=1,85     
        asv(i) = '   '
        write(asv(i),'(1x,i2)') idsv(i)
        if( idsv(i).gt.0 ) then
           asv(i)(1:1) = gnss
           if( asv(i)(2:2).eq.' ') asv(i)(2:2) = '0'
        endif
      enddo
c* rwk 090908 change for SP3-C 
      write(iungs,930) numsv,(asv(i),i=1,85)  
  930 format("+   ",i2,3x,17a3,/,4("+        ",17a3,/),$)
c*      write(iungs,930) numsv,(idsv(i),i=1,85)
c* 930   format("+   ",i2,3x,17i3,/,4("+        ",17i3,/),$)  

c	Lines 8-12
      write(iungs,940) (isatsgm(i),i=1,85)
940   format(5("++       ",17i3,/),$)

c	Lines 13-14   
c* rwk 090908 change for SP3-C
      write(iungs,950) 
 950  format("%c  G cc GPS ccc",4(" cccc"),4(" ccccc"),/,$)
      write(iungs,951)
 951  format("%c cc cc ccc ccc",4(" cccc"),4(" ccccc"),/,$)
c* 950   format(2("%c",2(" cc"),2(" ccc"),4(" cccc"),4(" ccccc"),/),$)

c	Lines 15-16
c* rwk 090908 change for SP3-C
      write(iungs,960)
 960  format(2("%f","  1.2500000","  1.025000000","  0.00000000000",
     .         "  0.000000000000000",/),$)
c* 960  format(2("%f","  0.0000000","  0.000000000","  0.00000000000",
c*     .         "  0.000000000000000",/),$)

c	Lines line 17-18
      write(iungs,970)
970   format(2("%i",4("    0"),4("      0"),"         0",/),$)
c	SP3 lines 19-21
      write(iungs,980)
980   format(3("/*                                                ",
     *         "                    ",/),$)      

c  Line  22 
c   The following code was used for earlier versions.   
c*      if( otlmod(1:4).eq.'NONE') then
c*        write(iungs,990) pcvmod,otlmod  
c* temporary proposed format August - 21 November 2006
c* 990     format("/* PCV:",a5," TL:",a8,"   N NONE         ",   
c*      *      "CLK:N/A ORB:CoM") 
c   erroneous code used 22 November 2006 - 14 April 2011
c*990     format("/* PCV:",a10," OL/AL:",a8,"   N NONE         ",   
c*     *      "YN ORB:CoN CLK:CoM") 
c*      else
c*        write(iungs,991) pcvmod,otlmod
c* 991     format("/* PCV:",a5," TL:",a8,"   Y NONE         ",
c*(     *      "CLK:N/A ORB:CoN")   
c*991     format("/* PCV:",a10," OL/ALL:",a8,"   Y NONE         ",
c*     *        "YN ORB:CoN CLK:CoN")  
c*      endif
c   New code correct according to IGS MAIL 5490 22 Nov 2006, 
c   applied 15 April 2011 ---rwk 
      atlmod = 'NONE    '
      otlflg = 'Y'
      atlflg = 'N'
      orbmod = 'CoN'
      clkmod = 'CoN'
      write(iungs,990) pcvmod,otlmod,atlmod,otlflg,atlflg,orbmod,clkmod
  990 format('/* PCV:',a10,' OL/AL:',2(a8,1x),2a1,' ORB:',a3,' CLK:',a3)

c  Echo the SP-3 file headers

      write(iprnt,'(//,a)') '------------------------------------------'
      write(iprnt,'(/,a)') '** SP3 HEADER INFORMATION **'
      write(iprnt,30) orbtyp,crdsys,org,igpswk
     1          , iyr,imon,iday,ihr,imin,sec,mjd,fmjd,nepoch,delt
     1          , numsv  
   30 format(/
     2      ,1x,'Orbit Type: ',a3,'   Coordinate System : ',a5
     3      , '  Organization : ',a4,'   GPS Week : ',i4,//
     4      ,1x,'Start epoch GPST : ',i4,4i3,1x,f10.7,/
     5      ,1x,'            MJD  : ',i5,1x,f15.14,//
     6      ,1x,'Number epochs    : ',i6,'  Interval :',f8.2,' sec',//
     7      ,1x,i2,' satellites ')
      write(iprnt,'(16(1x,a3),/)') (asv(i),i=1,numsv)
      write(iprnt,31) pcvmod,otlmod,atlmod,otlflg,atlflg,orbmod,clkmod
   31 format(/,'Corrections on PCV comment line:',/
     .      ,'  PCV model: ',a10,'  OTL model : ',a8,'  ATL model : ',a8
     .      ,'  OTL CMC included : ',a1,'  ATL CMC included : ',a1
     .      ,'  ORB CMC applied : ',a3,'  CLK CMC applied : ',a3) 

      return
      end


