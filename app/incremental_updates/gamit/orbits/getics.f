       subroutine getics( iungs,spver,jde,te,gnss_sel
     .                  , numsp3sv,numsat,itsat,satics )

c       Read an NGS SP file to get ICS near the center
c       R.W. King  March 1988
c      	mod. to allow SP3 format, P. Fang March 1993

c    Input:
c       iungs     : unit number of SP3 file  
c       spver     : version of SP file ('1' 'a' 'b' or 'c'')
c       jde       : PEP JD of IC epoch
c       te        : seconds of day of IC epoch
c       gnss_sel  : GNSS code for requested SVs (G R C E J I)
c       numsp3sv  : number of satellites on the SP# file
c       numsat    : number of satellites for the t-file
c       itsat     : array containing PRN #s of satellites for T-file 

c    Output
c       satics    : initial conditions for requested epoch

      implicit none

      include '../includes/dimpar.h'

      character*1 gnss_sel,spver,asvid
      character*31 string
      character*256 message

      integer*4 iungs,julday,jd,jde
     .        , numsp3sv,numsat,itsat(maxsat)
     .        , iyear,iday,imon,ihr,id,imin
     .        , isvcnt,i,j,k,ioerr

      real*8 te,satics(maxorb,maxsat),t,x(maxorb),sec

c Rewind the SP  file and dummy-read the header records

c      write(6,5)
c    5 format(//,1x,'Reading the NGS file to get T-file ICs',
c     1          1x,' near the center of the span',//)
c     rewind the file and dummy-read the header records
      rewind iungs
      if( spver.eq.'1') then
        read(iungs,'(a31)') string
        read(iungs,'(a31)') string
      else
        do i = 1,100
          read(iungs,'(a31)',end=70,err=70,iostat=ioerr) string
          if( string(1:1).eq.'*' ) goto 8
        enddo
    8   backspace iungs
        goto 10
      endif

c Read through all the epochs until the desired epoch is found

c     come back here for each new epoch until the ic epoch found
   10 if(spver.eq.'1') then
        read(iungs,11,end=60) iyear,imon,iday,ihr,imin,sec
   11   format(4x,i4,4(1x,i2),1x,f10.7)
      else
        read(iungs,12,end=60) iyear,imon,iday,ihr,imin,sec
   12   format(3x,i4,4(1x,i2),1x,f10.7)
      endif
* MOD TAH 210101: Updated 2020 to 2100.
      if( iyear.lt.1982.or.iyear.gt.2100 ) then    
       write(message,'(a,i5)') 'Error on SP file; bad year = ',iyear
       call report_stat('FATAL','NGSTOT','orbits/getics',' ', message,0)
      endif
      jd= julday( imon,iday,iyear )
      t= dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec
c*     all new T-files now written in GPST
c*     call timinc( jd,t,-gpsutc )
cd      print *,'numsv,jd,t,jde,te ',numsv,jd,t,jde,te
      if ( jd.lt.jde .or. jd.eq.jde.and.dabs(t-te).gt.1.d-6 ) then
c        We haven't reached the desired epoch yet,spool thru the rest
c        of the data at this epoch, and go read another
         do i=1,numsp3sv
           read(iungs,'(1x)')
         enddo
         goto 10
      else
c       We have found the desired epoch (within 1 microseconds)
        if (spver.eq."1") then
          do j=1,numsp3sv
            read(iungs,'(3x,i2,3(1x,f12.5),3(1x,f11.8))')
     .                  id,(x(i),i=1,6)
cd          write(6,29) id,(x(i),i=1,6)
            do i=1,6
               satics(i,j) = x(i)
             enddo
          enddo
        else               
          isvcnt = 0 
          do j=1,numsp3sv
            read(iungs,'(1x,a1,i2,3(1x,f13.6))',iostat=ioerr) 
     .                  asvid,id,(x(i),i=1,3)
            if(ioerr.ne.0) call report_stat('FATAL','NGSTOT'
     .         ,'orbits/getics',' '
     .        ,'Error reading decoding coordinates line sp3 file',ioerr)
            if( asvid.eq.' ' ) asvid = 'G'       
            do k=1,numsat
              if( asvid.eq.gnss_sel.and.id.eq.itsat(k) ) then
                do i=1,3
                  satics(i,k) = x(i)
                enddo 
              endif
            enddo
          enddo
        endif
      endif

cd      write(6,'(a,6(/,1x,3f11.3,3f9.4))') 
cd     .      'At end of GETICS, SATICS=',((satics(i,j),i=1,6),j=1,numsv)

c Rewind the SP file and position it to be read for tabular data

      rewind iungs
      if (spver.eq."1") then
        do i=1,2
          read(iungs,'(1x)')
        enddo
      else
c*      do i=1,22
c*      read(iungs,'(1x)')
c*      enddo
        do i = 1,100
          read(iungs,'(a31)') string
          if( string(1:1).eq.'*' ) goto 40
        enddo
   40   backspace iungs
      endif

      return
c
   60 write(message,65)
   65 format('End of NGS file without finding calculated midpoint')
      call report_stat('FATAL','NGSTOT','orbits/getics',' ',message,0)
   70 write(message,75)
   75 format('Error or end on NGS file while looking for data records')
      call report_stat('FATAL','NGSTOT','orbits/getics',' ',message
     .                ,ioerr)
      end

