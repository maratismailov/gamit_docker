	program emerge
c
c	merge a set of broadcast ephemeride files and eliminate duplicates
c
c	by pfang@ucsd.edu 93, modified by S. McClusky Jan 04 2000
c       
c  Usage:  Enter file names in response to the questions or
c          ls ????DDD?.YYn | emerge  
               
      Implicit none

      include '../includes/dimpar.h'
   
      integer maxeph
      parameter (maxeph=1000) 

      logical same 

      character*80 header(50),eph(maxsat*maxeph,8),inpfil,outfil 

      integer inplu,outlu,nlast,nhead,i,j,julday
      integer iyr(maxsat*maxeph),imth(maxsat*maxeph)
      integer iday(maxsat*maxeph), ihr(maxsat*maxeph) 
      integer imin(maxsat*maxeph), iprn(maxsat*maxeph)
      integer indx(maxsat*maxeph)

      real*8 sec(maxsat*maxeph), jday(maxsat*maxeph)   

      data inplu/10/,outlu/20/,outfil/"eph.merged"/     

      open(unit=outlu,file=outfil,status="unknown")

c initialize some variables
      nlast=1
      nhead = 1  

c start loop over the input nav files
      do while (nlast.ge.1)
        print *,'Enter input RINEX navigation file name (Ctrl-D to end)'
        read(*,'(a80)',end=80) inpfil
        open(unit=inplu,file=inpfil,status='old')
c save the header text
        nhead=1
        read(inplu,'(a)') header(nhead)
        do while(header(nhead)(61:63).ne."END")
          nhead=nhead+1
      	  read(inplu,'(a)') header(nhead)
        enddo
	     nhead=nhead-1
c loop over entries within one file
        do while (nlast.le.maxsat*maxeph)

c read in one entry
          read(inplu,'(a)',end=60) eph(nlast,1)  

c extract prn number yr, month day hour min and sec
          read(eph(nlast,1),10) iprn(nlast),iyr(nlast),imth(nlast)
     .    ,iday(nlast),ihr(nlast),imin(nlast),sec(nlast) 
10        format(i2,5i3,f5.1)  

c make 4char year and calculate decimal jd for each record
          call fix_y2k(iyr(nlast))
          jday(nlast) = julday(imth(nlast),iday(nlast),iyr(nlast)) + 
     .    ihr(nlast)/24.d0 + imin(nlast)/1440.d0 + sec(nlast)/86400.d0 

c finish reading the rinex nav record for this epoch and satellite
          do i=2,8
            read(inplu,'(a)',end=60) eph(nlast,i)
          enddo

c check if this record is a duplicated entry 
          same = .FALSE.
          i=1
          do while ( .not.same .and. i.lt.nlast )
            if (jday(i).eq.jday(nlast).and.iprn(i).eq.iprn(nlast)) then
               same = .TRUE.
            endif
            i=i+1
          enddo

c if not same, advance one entry
          if( .not.same) nlast=nlast+1

c end of the loop of all entries within one file
        enddo
60	     close(inplu)

c go to next input nav file
      enddo
c dump out non-duplicate nav entries
80    do i=1,nhead+1
        write(outlu,'(a)') header(i)
      enddo
      nlast = nlast - 1  

c sort entries by time order
      call indexx ( nlast, jday, indx )
      do i=1,nlast
        do j=1,8
          write(outlu,'(a)') eph(indx(i),j)
        enddo
      enddo
      stop
      end
