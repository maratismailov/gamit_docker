      Program TESTORBEX

*     Test converting the y-file yaw angles into quaternions and writting an ORBEX file.
*     This can be done entirely within the terrestrial (sp3) system, so separate out the
*     terrestrial-to-inertial transformations (t-file, n-body) that will be needed later
*     to use ORBEX quaternions in MODEL's inertial computations.
*     R. King September 2019

      implicit none

      include '../include/dimpar.h'
                             
* ORBEX variables 
      integer*4 iuobx,nosat 
      real*8  quatern(4)
      character*20 obxfiln 
      character*3 prn

* SP3 variables            
      integer*4 lusp3,nssat
      real*8 
      character*20 sp3filn

* Y-file variables        
      integer*4 luy,nysat,ievent,nyversn
      real*8 yatt
      character*20 yfiln 

* T-file variables
      integer*4 lut,ntsat,iprn,nepcht,ji0,jil,iy1,iy2,jlast,jnow,iendf 
      real*8 satcrd(6)
      character*20 tfiln,svantbody

* Other variables
      integer*4 jdobs,iepoch,ioerr,i  
      real*8 tobs,sun(6),xhat(3),yhat(3),zhat(3),fjd,tdtoff
                 
* Flag on whether or not to do the terrestrial-to-inertial computations
      logical inertial/.false./

     
* Set the input files and test epoch
                                
      iuy = 1
      iuo = 2              
      iut = 3
      iut1 = 4
      ipole = 5
      yfiln = '    '
      orbexfiln = ' '
      tfiln = ' ' 
*     The ORBEX file is for Week 2024 Day 0  2018 October 21 PJD 2458413 
      jdobs = 2458413
      tobs = 1500.d0 
*     In 2018 this is SV66, a Block IIF 
      prn = 'G27'
      iprn = 27
      svantbody = 'BLOCK IIF' 

* Open the files  
                   
      if(inertial) then 

*       T-file
        call topens (tfiln,'old',iut,ioerr)
                                                      
*       EOP files
        open (unit=iut1,file='ut1.',status='old',iostat=ioerr)
        if (ioerr .ne. 0) call report_stat('FATAL','TESTORBEX'
     .   ,'testorbex','ut1','Error opening UT1 table:',ioerr)
        open (unit=ipole,file='pole.',status='old',iostat=ioerr)
        if (ioerr .ne. 0) call report_stat('FATAL','TESTORBEX'
     .     ,'testorbex','pole.','Error opening pole table:',ioerr)
      endif
                        
*       Y-file
      open(unit=iuy,file=yfiln,form='unformatted' 
     .    ,access='sequential',status='old',iostat=ioerr)
      if(ioerr.ne.0) call report_stat('FATAL','TESTORBEX'
     .   ,'testorbex',orbixfiln,'Error opening y-file ',ioerr)

*       SP3 file
      open(unit=iusp3,file='sp3filn,status='old',iostat=ioerr)
      if (ioerr .ne. 0) call report_stat('FATAL','TESTORBEX'
     .    ,'testorbex',sp3filn,'Error opening SP3 file:',ioerr)

*       ORBEX file
      open(unit=iuobx,file='orbexfiln,status='old',iostat=ioerr)
      if (ioerr .ne. 0) call report_stat('FATAL','TESTORBEX'
     .    ,'testorbex',orbexfiln,'Error opening ORBEX file:',ioerr)

cc* Get the Sun coordinates
cc                               
cc      tdtoff =  (32.184d0 + 19.0d0)/86400.d0
cc      ivel = 0 
cc      fjd = jdobs + tobs/86400.d0 + tdtoff
cc      call ephtp(fjd,3,ivel,sun)  

* Get the SV coordinates
  
      if(inertial) then 
        call thdred *  *****  
        isat = iarray( ischan(ichan),itsat,ntsat)                                                 
        satarg = timdif( jdobs,tobs,jdbt,tbt ) 
        call gsatel( 2,satarg,iut,isat,svec,ytrp,yy
     .           , ntsat,sdelt,nintrs,ji0.ji1,iy1,iy2,jlast
     .           , jnow,iendf,nepcht ) 
       else
         call read_sp3



* Get the yaw angle from the y-file

      read(iuy,iostat=ioerr) nyversn
      if( ioerr.ne.0 ) then
           call report_stat('FATAL','MODEL','setup',yfiln
     .      ,'Error reading version from y-file header--old-style file?'
     .      ,ioerr)  
      else
       if(nyversn.lt.1051.or.nyversn.gt.1061  ) then
         write(message,'(a,i4,a)') 'Incompatible y-file version ('
     .      ,nyversn,')  File: '
        call report_stat('FATAL','MODEL','setup',yfiln,message,0)     
      endif
      rewind(iuy)
      if( nyversn.eq.1051 ) then
        read(iuy) nyversn,ytfile,yaw_s,yaw_e,nyinter,nyepoch
     .          , nysat,(ysat(isat),ysvn(isat),yblk(isat),isat=1,nysat)       
      elseif( nyversn.eq.1061 ) then
        read(iuy) nyversn,ytfile,yaw_s,yaw_e,nyinter,nyepoch,nysat
     .           , (ysat(isat),ysvn(isat),svantbody(isat),isat=1,nysat)
      endif
      call satatt( iuy,jd*1.d0+t/86400.d0,ntsat
     .           , yatt,ievent,iepoch)


* Get the SV-fixed coordinates in the SV system from the y-file

      call svbody_coords( iuy,satcrd,sun,xhat,yhat,zhat
     .                  , iprn,svantbody,yatt,ievent,iepoch)
      write(*,'(a,9f8.4)') 'SV coords in SV sys from y-file: '
     .             , xhat,yhat,zhat                                              

* Read the quaternions from the ORBEX file

      read_orbex_head(iuo,nosat)
      read_orbex_att(jdobs,tobs,nosat,prn,quatern)
      write(*,'(a,4f8.5)') 'Quaternion from ORBEX: ',quatern      
      
* Compute the inverse of the M-matrix from the quaternions to get the 
* terrestrial system phase center

  
      write(*,'(    )') '
* 





