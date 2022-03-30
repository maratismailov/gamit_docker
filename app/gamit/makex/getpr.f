      subroutine getpr ( debug,iwkn0,sow0,inter,nchan,ischan,ndat,dattyp
     .                  , iwkn,sow,icount,iprn,ierr,pr1,ioerr)

c
c     input:
c            iwkn0          GPS week number at initial epoch
c            sow0           GPS seconds of week at initial epoch
c            inter          interval between epochs
c            nchan          number of channels
c            ischan(maxsat) prn #s of satellites in channels
c            ndat           number of data types
c            dattyp         data types available

c     output
c            iwkn           GPS week number at current epoch
c            sow            GPS seconds of week at current epoch
c            icount         number of satellites at current epoch
c            iprn(maxsat)   PRN numbers of satellites at current epoch
c            ierr(maxsat)   error flags for satellites at current epoch
c            pr1(maxsat)    L1 pseudoranges in km for clock corrections
c            ioerr          return code for X-file read (0= ok; -1=eof)
c
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'


      logical debug  
      character*256 message
      integer dattyp(maxdat),iprn(maxsat),ierr(maxsat)
     1      , icount,isat,iamp1,iamp2,k,m,ioerr,ndat
     2      , ierx,iepoch,inter,ischan(maxsat),iwkn,nchan,iwkn0
      real*8 yl1(maxsat),yl2(maxsat),sow,sow0,dt
     1        , pr1(maxsat),pr2(maxsat),data1,data2,data3,data4


      do 105 k=1,nchan
         ierr(k) = 1
         yl1 (k) = 0.0d0
         yl2 (k) = 0.0d0
         pr1 (k) = 0.0d0
         pr2 (k) = 0.0d0
  105 continue

c       Read the epoch number and number of sats at this epoch

      read (uxfile,150,err=1000,end=1020,iostat=ioerr) iepoch,icount
  150 format (/,2i4)

c     Calculate the current GPS time

      dt = dble(inter)*(iepoch-1)
      call secsum( iwkn0,sow0,dt,iwkn,sow )


      if (debug) write(uscren,155) iepoch,icount,ndat,iwkn,sow
  155 format ('In GETPR: iepoch,icount,ndat,iwkn,sow:'
     .       , 1x,4i5,f12.2)

      do 170 m=1,icount
  160    format (10x,2i2,2(1x,d22.15,1x,i3),2(2x,d22.15))
         if (ndat .eq. 4) then
            read (uxfile,160,err=1000,end=1020,iostat=ioerr)
     .      ierx,isat,data1,iamp1,data2,iamp2,data3,data4
         else if (ndat .eq. 3) then
            read (uxfile,160,err=1000,end=1020,iostat=ioerr)
     .      ierx,isat,data1,iamp1,data2,iamp2,data3
         else if (ndat .eq. 2) then
            read (uxfile,160,err=1000,end=1020,iostat=ioerr)
     .      ierx,isat,data1,iamp1,data2,iamp2
         else
            read (uxfile,160,err=1000,end=1020,iostat=ioerr)
     .      ierx,isat,data1,iamp1
         endif

         iprn(m)=ischan(isat)
         ierr(m)=ierx

         if (dattyp(1).eq.1) yl1(m) = -data1
         if (dattyp(2).eq.2) yl2(m) = -data2
         if( ndat.lt.3 ) then
           write(message,'(a,i3,a)') 
     .       'No pseudoranges on X-file (ndat= ',ndat,')'
           call report_stat('FATAL','MAKEK','getpr',' ',message,0)
         endif
         if (dattyp(3).eq.3 .or. dattyp(3).eq.5) pr1(m) = data3
         if (dattyp(4).eq.4) pr2(m) = data4
         if(debug) write(uscren,161) iprn(m),pr1(m)
  161    format(12x,'iprn,pr1:',i3,d15.7)

  170 continue
      return

c     error reading x-files
 1000 call report_stat('FATAL','MAKEK','getpr',' ',
     .'Error, on X-file ',ioerr)

c     come here on end of file
 1020 continue

      return
      end


