      program countx

c     read an X-file and count the number of good observations

c     Written by Kurt Feigl, the hook.
c     Modified by M.Burc Oral Fri Dec 18 11:11:44 EST 1992
*     changes to xhdred call.
c     Modified by R. King 12 Feb 93 for changes to xhdred call
*     adjusted bin for 24 data set
*     added command line call

      implicit none

      include '../includes/dimpar.h'

      INTEGER*4  DATTYP(MAXDAT),LAMBDA(MAXSAT,MAXDAT),NTEXT
     .         , kount(maxsat),ioerr,luin, nepoch, nblen, ndat
     .         , ierx, idoy,jdoy,jm,jy,jh
     .         , inter, iy0,im0,id0,ihr0,min0,latd,latm,lond,lonm
     .         , kbin,nbin,isat, ischan(maxsat), iepoch
     .         , iamp1,iamp2,icount,jdoy0,nchan,nwide,mtime
     .         , iclarg,ii,i,j,k,m
     .        , ircint,isessn,len

      REAL*8  offarp(3), sec, sec0
     .      , data1,data2,data3,data4,seclon,seclat,height
      real*4 swver

      character*1 latflag,lonflag,letter
      character*3 rcvrsw,rxobtyp(maxdat)
      character*6 antcod
      CHARACTER*16 SITNAM,satnam(maxsat)
      character*20  rctype,rcvnum,anttyp,antnum
      character*60 filnam,wildcard,pickfn
      CHARACTER*80 TEXT(MAXTXT),buff80
      character*100 bargph(maxsat)
      CHARACTER*1  gnss

      logical lgood,lclmf,lmfile
           
***** Remove old versions of the status, warning, and error files

      call report_stat('CLEAR','COUNTX',' ',' ', ' ',0)
      call report_stat('CLEAR','LIB',' ',' ', ' ',0)
 
      ii = iclarg(1,filnam)
      if (ii .gt. 0) then
         lmfile = .true.
         lclmf  = .true.
      else    
         wildcard = 'x*.???' 
         len = 6
         filnam = pickfn ( wildcard,len )
         filnam = filnam(1:len)
      endif    
      write (6,5) filnam
  5   format (1x,'Opening: ',a16)

c     open current X-file
      luin = 7
      open (unit   = luin,
     .      file   = filnam,
     .      form   = 'formatted',
     .      status = 'old',
     .      iostat = ioerr )
      if( ioerr.ne.0 ) call report_stat('FATAL','COUNTX','countx'
     .                ,filnam,'Error opening file',ioerr)

      open (unit = 10,status='scratch')   


c     Read and display the X-file header.
      call XHDRED ( luin,10,10
     1,             NEPOCH,INTER,ircint,isessn
     2,             MTIME,IY0,IM0,ID0,IHR0,MIN0,SEC0
     3,             NCHAN,ISCHAN,satnam
     4,             NDAT,DATTYP,rxobtyp,LAMBDA
     4,             offarp,SITNAM,rcvrsw,swver,antcod
     5,             rctype,rcvnum,anttyp,antnum
     6,             LATFLAG,LATD,LATM,SECLAT
     7,             LONFLAG,LOND,LONM,SECLON,HEIGHT
     8,             NTEXT,TEXT,gnss)
C

c     figure out the width of the bins
      nbin = 100
      nwide = nepoch/nbin + 1 + 1
      do i = 1,maxsat
         write (bargph(i),13)
  13     format (100('|'))
      enddo

      jdoy0 = idoy(iy0,im0,id0)

      do 200 j=1,nepoch
         read (luin,'(/,a34)',iostat=ioerr) buff80 
         if( ioerr.ne.0 ) call report_stat('FATAL','COUNTX','countx'
     .                    ,filnam,'Error reading line on X-file',ioerr)
         if (buff80(8:8) .eq. '0') then 
            read (buff80,150,iostat=ioerr) iepoch,icount    
            if( ioerr.ne.0 ) call report_stat('FATAL','COUNTX','countx'
     .                    ,filnam,'Error reading epoch on X-file',ioerr)
  150       format (2i4)
            jdoy = 0
            jy = 0
            jh = 0
            jm = 0
            sec = 0.0d0
         else  
            read(buff80,151,iostat=ioerr)
     .      iepoch,icount,jy,jdoy,jh,jm,sec
  151       format (2I4,I5,I4,2I3,F11.7)    
            if( ioerr.ne.0 ) call report_stat('FATAL','COUNTX','countx'
     .                    ,filnam,'Error reading epoch on X-file',ioerr)

         endif

         if (iepoch .ne. j) then
            write (6,*) 'READX: epoch mismatch',iepoch,' ',j
         endif

         if (mod(iepoch,nwide) .eq. 1) then
            do k = 1,maxsat
               kount(k)= 0
            enddo
         endif

         do 170 m=1,icount
  160       format (10x,2i2,2(1x,d22.15,1x,i3),2(2x,d22.15))
            if (ndat .eq. 4) then   
               read (luin,160,err=1000,end=1020,iostat=ioerr)
     .         ierx,isat,data1,iamp1,data2,iamp2,data3,data4  
            else if (ndat .eq. 3) then
               read (luin,160,err=1000,end=1020,iostat=ioerr)
     .         ierx,isat,data1,iamp1,data2,iamp2,data3
            else if (ndat .eq. 2) then
               read (luin,160,err=1000,end=1020,iostat=ioerr)
     .         ierx,isat,data1,iamp1,data2,iamp2
            else
               read (luin,160,err=1000,end=1020,iostat=ioerr)
     .         ierx,isat,data1,iamp1
            endif

            if (lgood(ierx)) then
               kount(isat) = kount(isat) + 1
            endif
  170    continue

         if (mod(iepoch,nwide) .eq. 0) then
            do k = 1,maxsat
               if (kount(k) .eq. 0) then
                  write (letter,'(a)') '.'
               else if (kount(k) .le. 9) then
                  write (letter,'(i1)') kount(k)
               else
                  write (letter,'(a)') char(kount(k)-10+ichar('A'))
               endif
               kbin = iepoch/nwide + 1
               bargph(k)(kbin:kbin) = letter
            enddo
         endif
  200 continue

      write (6,600) nwide
 600  format(' Epochs/bin = ',i5,10x,
     &     "A=10",10x,"F=15",10x,"K=20",10x,"P=25",10x,"U=30")

      do i = 1,nchan
         write (6,225) ischan(i),filnam(1:nblen(filnam)),bargph(i)
  225    format (1x,'PRN',1x,i2,1x,a,1x,a100)
      enddo

      close (luin)
 1000 continue  
      if( ioerr.ne.0 ) call report_stat('FATAL','COUNTX','countx'
     .                  ,filnam,'Error reading data on X-file',ioerr)
 1020 continue
      stop
      end


