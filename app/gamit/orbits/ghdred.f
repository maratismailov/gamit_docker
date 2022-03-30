Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.
      subroutine ghdred(gfname,nics,gframe,gfprec,ngics,iuin,program
     .                 ,trun0,trunf,te,tb,tf,jde,jdb,jdf,delt,tdtoff
     .                 ,time_type,srpmod,gfnut,gfgrav,icsnam)
c
c       Reads ic epoch and model parameters from the G-file, and
c       the start,stop times of the integration from the ARC batch file
c       R. King  June 1995 from old timred of March 1987
c
c Input:  gfname :   Name of ics (g-) file
c         nics   :   Number of parameters (ICs + radiation pressure) to be integrated
c         iuin   :   Unit number of input gfile
c         program:   Name of calling program
c
c Output:
c         gfprec: precession model to which IC's refer in gfile
c         gframe: Inertial frame to which IC's refer in gfile
c         gfsrp : Radiation pressure model for gfile coefficients
c         ngics : Number of parameters on G-file
c
c    Output stored in common - PT950727: moved to argument list!
c
c       common/stint/
c            trun0     time of start of ephemeris, in seconds past epoch
c            trunf     time of end of ephemeris, in seconds past epoch
c       common/timxtr/
c            jde       julian day number of ics epoch
c            jdb       julian day number of start time
c            jdf       julian day number of stop time
c            te        seconds of day of ics epoch
c            tb        seconds of day of start time
c            tf        seconds of day of stop time
c            delt      ephemeris tabular interval in seconds
c            tdtoff    time offset between TDT and integration time
c
c       common/srprok/
c            icsnam    parameter names, from G-file and/or input models

      implicit none
c
      include '../includes/dimpar.h'
c
      character*3 upperc
      character*4 sphrc_par(9),srdyz_par(9),srxyz_par(9)
     .          , ecom1_par(15),ecom2_par(19),uclr1_par(15)
     .          , uclr2_par(15),icsgfl(maxorb),g_time_flag,program
      character*5  gfprec,gframe,gfsrp,gfnut,gfgrav
      character*8  scarr(12)
      character*16 gfname,xfname
      character*23 buf23
      character*69 buf69
      character*80 prog_name
      character*256 message

      integer*4 julday,id,ih,im,idoy,min,iy,ilen,nics,ngics,i,j,iuin
      integer*4 len,rcpar,jds

      real*8 taiutc,sec,utcoff,ts

c      common/stint/trun0,trunf
      real*8 trun0,trunf

c      common/timxtr/te,tb,tf,delt,tdtoff,jde,jdb,jdf
      real*8 te,tb,tf,delt,tdtoff
      integer*4 jde,jdb,jdf

c time and frame definition - arc, ghdred, wrthed
c       common/timfrm/time_type
       character*4 time_type

c      common/srprok/srpmod,icsnam
      character*4 icsnam(maxorb)
      character*5 srpmod

      data ecom1_par/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     $ ,'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN','BCOS','BSIN'/  
      data ecom2_par/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     $ ,'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN','BCOS','BSIN'
     . ,'D2CS','D2SN',D4CS','D4SN'/                                      
      data ecomc_par/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     $ ,'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN','BCOS','BSIN'
     . ,'D2CS','D2SN',D4CS','D4SN'/                                      
      data uclr1_par/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     $ ,'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN','BCOS','BSIN'/  
       data uclr2_par/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     $ ,'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS','YSIN','BCOS','BSIN'/  

c        Get calling program name for report_stat
      len =  rcpar(0,prog_name)

c        Read the first line of the G-file for IC epoch and models

      buf69 = ' '
      read(iuin,'(a69)',err=5) buf69
      goto 10
    5 buf69(34:39) = '  UTC '
   10 if( buf69(3:3).eq.' ') then
c       old format, 2-digit year
        read(buf69,11,err=90) iy,idoy,ih,min,sec,g_time_flag
     .                    , gframe,gfprec,gfsrp,gfnut,gfgrav
   11   format(i2,1x,i3,1x,i2,1x,i2,f3.0,20x,a4,5(1x,a5))
        call fix_y2k(iy)
      else
c       new format, 4-digit year
        read(buf69,12,err=90) iy,idoy,ih,min,sec,g_time_flag
     .                    , gframe,gfprec,gfsrp,gfnut,gfgrav
   12   format(i4,1x,i3,1x,i2,1x,i2,f3.0,18x,a4,5(1x,a5))
      endif

c        Check the reference frame, precession, nutation, and gravity models and assign if necessary

      if(gframe.eq.'     ' .or. gframe.eq.'1950.0') gframe = 'B1950'
      if(gframe.ne.'B1950'.and.gframe.ne.'J2000')
     .   call report_stat('FATAL',prog_name,'orbits/ghdred',gframe
     .                      ,'Unknown reference frame:',0)
c     if no precession model specified, assign IAU68 to B1950, IAU76 to J2000
      if(gfprec.eq.'     '.and.gframe.eq.'B1950')gfprec = 'IAU68'
      if(gfprec.eq.'     '.and.gframe.eq.'J2000')gfprec = 'IAU76'
      if( gfnut.eq.'     ' ) gfnut = 'IAU80'
      if( gfgrav.eq.'     ' ) gfgrav = 'IGS92'

c        Read the parameter number and names from the G-file

c  Clear IC parameter arrays (icsgfl and icsnam).
        do j = 1,maxorb
          icsgfl(j) = ' '
          icsnam(j) = ' '
        enddo
        read(iuin,'(i2)') ngics
        if( ngics.gt.maxorb ) then
           write(message,'(a,i2,a,i2)') 'Parameter number on G-file ('
     .         ,ngics,') exceeds MAXORB (',maxorb,')'
          call report_stat('FATAL',prog_name,'orbits/ghdred',' ',message,0)
         endif
         backspace iuin
         read(iuin,'(i2,23(1x,a4))') ngics,(icsgfl(i),i=1,ngics)

c  Check to see that the radiation pressure model selected is consistent
c  with the model specified on the g-file--both with the name on the
c  first line and with the coefficient names on the second line
c  (No longer allow a blank (missing) radiation pressure model since the
c  original defaults (SPHRC) is no longer supported. )
       
      if(gfsrp.ne.'ECOM1'.and.gfsrp.ne.'BERNE'.and.gfsrp.ne.'ECOM2'.and.
     . gfsrp.ne.'ECOMC'.and.gfsrp.ne.'UCLR1'.and.gfsrp.ne.'UCLR2' ) then
        write(message,'(a,a5,a)') 'G-file radiation model (',gfsrp
     .    ,') not supported'
        call report_stat('FATAL',prog_name,'orbits/ghdred',' '
     .   ,message,0)
      endif
      if( gfsrp.eq.srpmod ) then
          do j=1,ngics
             icsnam(j) = icsgfl(j)
          enddo
      else
          if(program.eq.'ARC ')write(message,13) gfsrp,srpmod
          if(program.eq.'ARC ')write(8,13) gfsrp,srpmod
   13     format('G-file solar radiation pressure parameters ('
     .      ,a5,') are not consistant with input file model (',a5,')' )
          call report_stat('WARNING',prog_name,'orbits/ghdred',' '
     .                     ,message,0)
c*          if(program.eq.'ARC ')write(6,14)icsgfl
          if(program.eq.'ARC ')write(8,14)icsgfl
   14     format(1x,'G-file: ',15(a4,1x))
          call uppers(srpmod)
          if( srpmod.eq.'BERNE'.or.srpmod.eq.'ECOM1'.or.
     .            srpmod.eq.'UCLR1' ) then
              do j=1,15
                icsnam(j) = ecom1_par(j)
              enddo
          elseif( srpmod.eq.'ECOM2' ) then
              do j=1,19
                 icsnam(j) = ecom2_par(j)
              enddo
          else
            if(program.eq.'ARC ')
     .        call report_stat('FATAL',prog_name,'orbits/ghdred',srpmod
     .                        , 'Unrecognized SRP model: ',0)
          endif
          if(program.eq.'ARC ')then
c*            write(6,15) (icsnam(j),j=1,nics)
            write(8,15) (icsnam(j),j=1,nics)
   15       format(1x,'Model : ',15(a4,1x))
c*            write(6,'(1x)')
          endif
      endif

c
   19 read (iuin,'(12a8)',end=95) scarr
      if (upperc(scarr(1)(1:3)).eq.'END') go to 20
      go to 19


c        Check the time-type (GPST or UTC) and convert to JD, seconds-of-day

  20  continue

c PT950725: put in a kluge to get around the fact that there is no time_type
c           definition when running GTOG. Assume it is GPST, since this
c           program will only be run using GPST (hopefully!)
      if(program.eq.'GTOG') time_type = 'GPST'

      if( g_time_flag.eq.'    ' ) g_time_flag = 'UTC '
      call monday( idoy,im,id,iy )
      jde= julday( im,id,iy )
      te= ih*3600.d0 + min*60.d0 + sec
      write(8,21) gfname,g_time_flag
  21  format(/,' Epoch of ICs, read from: ',a16,/,1X,
     1           ' Yr  Mn Dy    Hr  Min Sec  JD    Sec of Day (',a4,')')
      write(8,22) iy,im,id,ih,min,sec,jde,te
  22  format(1x,i4,2(1x,i2),3x,2(1x,i2),1x,f3.0,2x,i7,f15.8 )
      utcoff = taiutc(jde) - 19.d0
      if( time_type.eq.'UTC ') then
         jds = jde
         ts = te
         if( g_time_flag.eq.'GPST') call timinc(jds,ts,-utcoff)
c        utcoff will be wrong if GPST JD is a day later than UTC JD at a leap-sec
         if( jds.ne.jde ) then
            utcoff = taiutc(jds) - 19.d0
            call timinc(jds,ts,-utcoff)
         endif
         jde = jds
         te = ts
      elseif( time_type.eq.'GPST') then
         if( g_time_flag.eq.'UTC ') call timinc(jde,te,utcoff)
      else
        call report_stat('FATAL',prog_name,'orbits/ghdred',time_type
     .                      ,'Bad time type',0)
      endif

c
c        Compute the start and stop times for the integration
      if(program.eq.'ARC ')then
c       see if the observation times are to read from an x-file
        read(5,30) xfname
   30   format(a16)
        if (xfname.eq.'                ') goto 40
        call report_stat('FATAL',prog_name,'orbits/ghdred',' ' 
     .      ,'Option to read from X-file not supported',0)
c       epochs can be yy mm dd (old) or ..yy doy (new)
c       read the observation times from the input (batch) file
c       start time:
   40   read (5,'(a23)') buf23
        if( buf23(1:2).eq.'  ' .or. buf23(1:2).eq.'19' .or.
     .      buf23(1:2).eq.'20' ) then
c         new style
          read(buf23,45) iy,idoy,ih,min,sec
   45     format(i4,1x,i3,1x,2(i2,1x),f8.5)
          call monday(idoy,im,id,iy)
        else
c         old style
          read (buf23,50) iy,im,id,ih,min,sec 
   50     format (5(i2,1x),f8.5)
        endif
        ilen = 15 
        call check_y2k(iy)
        jdb= julday( im,id,iy )
        tb= ih*3600.d0 + min*60.d0 + sec
        write(8,60) time_type
   60   format(/,' Input start time for observations: ',/,1X,
     1           ' Yr  Mn Dy    Hr  Min Sec  JD    Sec of Day (',a4,')')
        write(8,61) iy,im,id,ih,min,sec,jdb,tb
   61   format(1x,i4,2(1x,i2),3x,2(1x,i2),1x,f3.0,2x,i7,f15.8 )
c       stop time:
        read (5,'(a23)') buf23
        if( buf23(1:2).eq.'  ' .or. buf23(1:2).eq.'19' .or.
     .      buf23(1:2).eq.'20' ) then
c         new style
          read(buf23,45) iy,idoy,ih,min,sec
        else
c         old style
          read (buf23,50) iy,im,id,ih,min,sec  
          call check_y2k(iy)    
        endif
        call monday(idoy,im,id,iy)
        jdf= julday( im,id,iy )
        tf= ih*3600.d0 + min*60.d0 + sec
        write(8,70) time_type
   70   format(/,' Input stop time for observations: ',/,1X,
     1           ' Yr  Mn Dy    Hr  Min Sec  JD    Sec of Day (',a4,')')
        write(8,61) iy,im,id,ih,min,sec,jdf,tf
c
c        Compute start and stop times for integration by expanding
c        the observation times by 11 tabular intervals and rounding
c        to the nearest tabular interval. Originally 7 tabular intervals
c        was enough, but now we need to be able to interpolate the tfile
c        up to 1 hour before the start of a session to handle satellites
c        that begin the session in eclipse. McClusky 950405   
c        Occasional peculiar eclipse sequence means that we need 13 tabular 
c        intervals. Matt King/Bob King 010615/010717  
c        Another case found where 13 is not enough; increase to 14. R King 030721
c
        trun0= (jdb-jde)*86400.d0 + (tb-te)
        trunf= (jdf-jde)*86400.d0 + (tf-te)
        trun0= trun0 - dmod(trun0,delt) - 14.d0*delt
        trunf= trunf - dmod(trunf,delt) + 14.d0*delt
        jdb= jde
        tb = te
        call timinc( jdb, tb, trun0 )
        jdf= jde
        tf = te
        call timinc( jdf, tf, trunf )
      endif
c        calculate terrestrial dynamical time (TDT) - integration time (UTC or GPST)

      if( time_type.eq.'UTC ' ) then
        tdtoff= 32.184d0 + taiutc(jde)
      elseif( time_type.eq.'GPST') then
        tdtoff= 32.184d0 + 19.0d0
      endif
      write(8,80) time_type,tdtoff
   80 format(/,1x,'TDT-',a4,' at IC epoch   = ',f7.4,' sec')
      return
   90 call report_stat('FATAL',prog_name,'orbits/ghdred',gfname
     .                      ,'Error reading G-file: ',0)
   95 call report_stat('FATAL',prog_name,'orbits/ghdred',gfname
     .                     ,'No END found of header of G-file',0)
      end
