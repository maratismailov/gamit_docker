      program erptab

c     Make PEP-type UT1 and pole table from GLOBK pmu values
c     R. King    13 December 1991
c     R. King    22 April 2015 - add another decimal place to match IERS and pmu tables

      character*1 reg,lowerc
      character*80 infn,ut1fn,polefn
      character*80 buf80,title

      integer*4 ndim,iterm,iscrn,in,iut1,ipole,idelt,julday
     .        , itype,ihr,imin,iday,imon,iyr,i,n

      parameter (ndim=20000)

      real*8 xmjd(ndim),x,xp(ndim),sigx,y,yp(ndim),sigy,ut1,sigut1
     .     , step,steps,units,taiut1(ndim)

      data iterm/5/,iscrn/6/,in/10/,iut1/11/,ipole/12/,n/0/
     .   , units/1.d0/,itype/2/


c       Open the input and output files

      write(iscrn,1)
    1 format(//,1x,' Enter the input earth rotation filename:',/)
      read(iterm,2) infn
    2 format(a)
      open (unit=in,file=infn,status='old',err=91)
      write(iscrn,3)
    3 format(//,1x,' Enter the output UT1. filename:',/)
      read(iterm,2) ut1fn
      open (unit=iut1,file=ut1fn,status='unknown',err=92)
      write(iscrn,4)
    4 format(//,1x,' Enter the output pole. filename:',/)
      read(iterm,2) polefn
      open (unit=ipole,file=polefn,status='unknown',err=93)


c       Read in title for PEP-format file

      write(iscrn,5)
    5 format(//,1x,'Enter 80-character title for output file:')
      read(iterm,6) title
    6 format(a80)


c       Is UT1 regularized ?   Default is yes (itype=2)

      write(iscrn,7)
    7 format(//,1x,'Is UT1 regularized (UT1R) ? (Y/N)')
      read(iterm,8) reg
    8 format(a)
      if( lowerc(reg).eq.'n') itype = 4


c       Read and echo comment lines

      write(iscrn,10)
   10 format(//,1x,'Title records for input earth-rotation file:',/)

   20 read(in,21,end=50) buf80
   21 format(a80)
      if( buf80(1:1).eq.'*' ) then
        write(iscrn,22) buf80
   22   format(1x,a80)
        goto 20
      endif

c     Go decode the data lines

c*    read(buf80,31) iyr,imon,iday,ihr,imin,x,sigx,y,sigy,ut1,sigut1
c*   31 format(1x,i4,4i3,1x,f8.5,1x,f6.4,1x,f8.5,1x,f6.4,1x,f13.6,1x,f8.5)
      read(buf80,*) iyr,imon,iday,ihr,imin,x,sigx,y,sigy,ut1,sigut1 
c     call check_y2k(iyr)
* MOD TAH 000218: Do direct check to aviod writing warning file
      if( iyr.le.1900) then 
                  
c        earliest GPS launch date is 1978; earliest space-geodetic data about 1960 
         if( iyr.gt.60 ) then
            iyr = iyr + 1900
         else
            iyr = iyr + 2000
         endif  

      endif
* END no check_y2k replacement
        
      n= n+1
      if( n.gt.ndim ) then
         write(iscrn,32) n
   32    format(1x,'Number of values exceeds dimension',i5,' Stop')
         call report_stat('FATAL','ERPTAB','erptab',' ',
     .                    'Number of values exceeds maximum allowed',
     .                    ndim)

      endif
c     write(iscrn,31) iyr,imon,iday,ihr,imin,x,sigx,y,sigy,ut1,sigut1
      if( ihr.ne.0 .or. imin.ne.0 ) then
         write(iscrn,41) n
   41    format(1x,'Non-zero hr min at line',i5,' Stop')
         stop
         endif
      xmjd(n) = julday(imon,iday,iyr) - 2400000 - 1
      taiut1(n) = -ut1
      xp(n) = x
      yp(n) = y
c     units are seconds of time and arc
      units = 1.d0

c     write(iscrn,8001) xmjd(n)
c8001  format(1x,'xmjd=',f20.3)
      goto 20


c       Check for evenly-spaced data

   50 step = xmjd(2)-xmjd(1)
      do 55 i=3,n
        steps = step
        step  = xmjd(i) - xmjd(i-1)
        if( step.ne.steps ) then
           write(iscrn,51) i
   51      format(1x,'Irregular spacing at line',i10)
           stop
           endif
   55  continue
       idelt = step


c        Write the output PEP-style tables

      call peput1( iscrn,iut1,title,idelt,units,itype,n,xmjd,taiut1 )
      call pepwob( iscrn,ipole,title,idelt,units,n,xmjd,xp,yp )
      write(iscrn,61)
   61 format(1x,'PEP UT1 and POLE tables written')

      stop

   91 stop 'Error opening input earth-rotation file, stop'
   92 stop 'Error opening output PEP UT1 file, stop'
   93 stop 'Error opening output PEP pole file, stop'
      end

      subroutine peput1( iscrn,iut1,title,idelt,units,itype,n,t,u )
c
c  Write a PEP-format output ut1 table from an array of values
c
c   iut1   output unit number
c   title  title for output file
c   idelt  interval (days) of values in output table
c   units  input units (sec) output is 1.e-7 sec)
c   n      number of valuies in input arrays
c   t      array of input dates (mjd)
c   u      array of input values of TAI-UT1

      integer*4 iut1,iscrn,iu(6),mjd(6),npr,idelt,itype,jd1,jd2,i,j,k,n
     .        , k1

      real*8  t(*),u(*),units,outunt,uu,taiutc

      character*80 title
      character*24 varfmt

      data outunt/1.e-7/, npr/6/
     .   , varfmt/'(5x,i5,6(i10,1x),1x,i1) '/
c       itype= 2 for TAI-UT1 (only type allowed at present)

      jd1= int(t(1)) + 2400000
      jd2= int(t(n)) + 2400000

      write(iut1,10) title,varfmt,itype,jd1,jd2,npr,idelt,outunt
   10 format(a80,/,a24,8x,i2,1x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0)

      i= 0
   20 k = 0
   30 i = i + 1
      if( i.gt.n ) goto 50
      k = k + 1
      uu = u(i)*units/outunt
      iu(k) = int(sign(abs(uu)+0.5,uu))
      mjd(k) = t(i)

c  PT 950525: check whether the ut value is tai-utc or utc-at
      if(iu(k).lt.20.d6)then
        iu(k) = iu(k) + taiutc(mjd(k)+2400001)/outunt
      endif

      if ( k.ne.6 ) goto 30
      write(iut1,varfmt) mjd(1),(iu(j),j=1,6)
      goto 20
   50 continue
      if (k.eq.0 ) goto 70             
      k1=(6-k)*11
      write(varfmt(8:24),'(i1,a,i2,a)') k,'(i10,1x),',k1,'x,i1)'
cd    print *,'varfmt ',varfmt
      write(iut1,varfmt) mjd(1),(iu(j),j=1,k),k

   70 return
      end

      subroutine pepwob( iscrn,iwob,title,idelt,units,n,t,x,y )
c
c  Write a PEP-format output pole table from an array of values
c
c   iwob   output unit number
c   title  title for output file
c   idelt  interval (days) of values in output table
c   units  input units (arcsec) output is 1.e-6 arcsec)
c   n      number of valuies in input arrays
c   t      array of input dates (mjd)
c   x      array of input values of x component of pole
c   y      array of input values of y component of pole


      integer*4 iwob,iscrn,ix(6),iy(6),mjd(6),npr,idelt
     .        , jd1,jd2,i,j,k,n

      real*8  t(*),x(*),y(*),units,outunt,xx,yy

      character*80 title
      character*24 varfmt

      data outunt/1.e-6/,  npr/6/
     .   , varfmt/'(5x,i5,12i8,1x,i2)      '/

      jd1= int(t(1)) + 2400000
      jd2= int(t(n)) + 2400000

      write(iwob,10) title,varfmt,jd1,jd2,npr,idelt,outunt
   10 format(a80,/,a24,11x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0)

      i= 0
   20 k = 0
   30 i = i + 1
      if( i.gt.n ) goto 50
      k = k + 1
      xx = x(i)*units/outunt
      yy = y(i)*units/outunt
      ix(k) = int(sign(abs(xx)+0.5,xx))
      iy(k) = int(sign(abs(yy)+0.5,yy))
      mjd(k) = t(i)
      if ( k.ne.6 ) goto 30
      write(iwob,varfmt) mjd(1),(ix(j),iy(j),j=1,6)
      goto 20
   50 continue
      if (k.eq.0 ) goto 70      
      k1=(12-2*k)*8                     
      write(varfmt(8:24),'(i2,a,i2,a)') 2*k,'i8,',k1,'x,i3)'
cd    print *,'varfmt ',varfmt                        
      write(iwob,varfmt) mjd(1),(ix(j),iy(j),j=1,k),2*k

   70 return
      end
                   
