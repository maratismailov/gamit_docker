Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995. All rights reserved.
      subroutine wrthed
c
c     write out the ephemeris (t-) file header
c     R.I. Abbot - november 1984
c     modified for new version of subroutine timred by R. King March 87

      implicit none

      include '../includes/dimpar.h'    
      include '../includes/global.h'
      include '../includes/arc.h'


      integer*4 ih0,ihnsec,iyear,imonth,iday,imin0,ihe,ihf,ihr,imn
     .        , imine,iminf,isec,nrad,i,k

      real*8 sec0,xmine,xminf,he,hf,sece,secf,t,xmin0,h0
                             
      character*1 upperc
      character*2 buf2
      character*4 buf4
      character*80 head,comt(3)
      character*256 message

      integer*4 ime,ide,iye,im0,id0,iy0,imf,idf,iyf,nepoch,nsat1,nintrp
     .        , idoy,nclen

c
      head = ' '
      do i=1,3
       comt(i)= ' '
      enddo


c        set up T-file header
c
      head(1:39)=' Tabular ephemeris file generated from '
      head(40:55)=gfname
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      write (buf2,'(i2)') iday
      read  (buf2,'(a2)') head(57:58)
      head(59:59)='-'
      write (buf2,'(i2)') imonth
      read  (buf2,'(a2)') head(60:61)
      head(62:62)='-'   
      write (buf4,'(i4)') iyear
      read  (buf4,'(a4)') head(63:66)
      write (buf2,'(i2)') ihr
      read  (buf2,'(a2)') head(68:69)
      head(70:70)=':'
      write (buf2,'(i2)') imn
      read  (buf2,'(a2)') head(71:72)
      head(73:73)=':'
      write (buf2,'(i2)') isec
      read  (buf2,'(a2)') head(74:75)
      write (buf2,'(i2)') nics
      read  (buf2,'(a2)') head(79:80)   

c        convert times of ic epoch, start and stop times of
c        the integration to calender day, hr, min, sec
c        prior to GAMIT release 9.9 (ARC 9.50), the ic epoch year was
c        measured from 1900 (e.g. 87 instead of 1987), but the start and 
c        stop years were the full 4 digits; now all are 4 digits
c          
          
c*      call mdyjul( ime,ide,iye,jde )  --routine removed from ARC, use /lib routines
      call dayjul( jde,iye,idoy )
      call monday( idoy,ime,ide,iye )   
c*    changed to avoid non-roundoff  rwk 031111
      call ds2hms( iye,idoy,te,ihe,imine,sece )
cd      print *,'iye idoy te ihe imine sece '
cd     .       , iye,idoy,te,ihe,imine,sece
c*      call hms( te/86400.d0,ihe,imine,sece )
      he= ihe
      xmine= imine
c*      call mdyjul( im0,id0,iy0,jdb ) 
      call dayjul( jdb,iy0,idoy )
cd      print *,'jdb,iy0,idoy ',jdb,iy0,idoy
      call monday( idoy,im0,id0,iy0 )   
c*    changed to avoid non-roundoff  rwk 031111
      call ds2hms( iy0,idoy,tb,ih0,imin0,sec0 )
c*      call hms( tb/86400.d0,ih0,imin0,sec0 ) 
cd      print *,'iy0 idoy tb ih0 imin0 sec0 '
cd     .       , iy0,idoy,tb,ih0,imin0,sec0
      h0= ih0
      xmin0= imin0
c*      call mdyjul( imf,idf,iyf,jdf )   
      call dayjul( jdf,iyf,idoy )
      call monday( idoy,imf,idf,iyf )   
c*    changed to avoid non-roundoff  rwk 031111
      call ds2hms( iyf,idoy,tf,ihf,iminf,secf )
c*      call hms( tf/86400.d0,ihf,iminf,secf ) 
cd      print *,'iyf idoy tf ihf iminf secf '
cd     .       , iyf,idoy,tf,ihf,iminf,secf
      hf= ihf
      xminf= iminf

                                        
c        Compute the number of epochs for the T-file header

      nepoch= ( (jdf-jdb)*86400.d0 + (tf-tb) ) / delt + .000001d0 + 1


c        Write the computed start and stop times to the archive file

      write(iarh,20) time_type,iye,ime,ide,ihe,imine,sece
     1          , iy0,im0,id0,ih0,imin0,sec0
     2          , iyf,imf,idf,ihf,iminf,secf,nepoch,delt,diint
   20 format(/,' Times written to T-file header  (',a4,')',/
     1        ,' IC epoch:  ',i4,4i3,f10.6,/
     2        ,' Start   :  ',i4,4i3,f10.6,/
     3        ,' Stop    :  ',i4,4i3,f10.6,/
     4        ,' No. epochs: ',i4,/
     5        ,' Tabular interval: ',f7.1,'  Integration interval:'
     6        ,f9.4,/)

c       write the T-file header

cd    print *,'WRTHED writing the scratch file header isunit ',isunit
      write(isunit) head,ime,ide,iye,he,xmine,sece
     1                  ,im0,id0,iy0,h0,xmin0,sec0
     2                  ,imf,idf,iyf,hf,xminf,secf
     3                  ,delt,nepoch,ut1utc,xpole,ypole

cd      write(iarh,*) 'header follows:'
cd      write(iarh,*) head
cd      write(iarh,*) ime,ide,iye,he,xmine,sece
cd      write(iarh,*) im0,id0,iy0,h0,xmin0,sec0
cd      write(iarh,*) imf,idf,iyf,hf,xminf,secf
cd      write(iarh,*) delt,nepoch,ut1utc,xpole,ypole

      nsat1=1
      if( upperc(apar).eq."V" ) then
        nintrp = (kount+1)*6
      else
        nintrp=(kount+1)*3
      endif
      nrad=nics-6

      comt(1)(11:14) = time_type
      comt(1)(21:30) = frame_name(1:10)
      comt(1)(36:40) = frame
      comt(1)(42:46) = precmod
      comt(1)(48:52) = nutmod
      comt(1)(54:58) = gravmod                           
      comt(1)(60:64) = srpmod
c** rwk 110114: temporary to get through rest of GAMIT. EJP added UCLR2 4April2012
      if((srpmod.eq.'UCLR1').or.(srpmod.eq.'UCLR2')) 
     .                                    comt(1)(60:64) = 'BERNE'
      comt(1)(66:70) = eradmod
      comt(1)(72:76) = antradmod 
         
      nclen = nics
      if(nics.gt.15) nclen = 15
      k = 1
      do i = 1,nclen
        comt(2)(k:k+3) = icsnam(i)
        k = k+5
      enddo         
      if(nics.gt.15) then
        k = 1
        do i=16,nics
          comt(3)(k:k+3) = icsnam(i)
          k = k+5
        enddo
      endif
      write (isunit) comt,nsat1,nintrp,sname
     1             , (satics(i),i=1,6),(radcon(i),i=1,nrad)          
cd      print *,'Record 2 sname satics ',sname,(satics(i),i=1,6)
cd    write(*,40) sname,head,nintrp,nics,(satics(i),i=1,6)
cd   1          , (radcon(i),i=1,nrad)
cd 40 format(//,1x,'T-file headers for ',a16,/,1x,a80,/
cd 1      ,   1x,'nintrp=',i3,'  nics=',i3,'  satics=',/
cd 2      , 6(1x,3d22.15/) )

      return
      end
