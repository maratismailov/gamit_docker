Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1997.   All rights reserved.

       Program MERGE_RINEX

c Simon McClusky              Feb 1995  
c Modified by P. Tregoning    Oct 1995  
c Modified by S. McClusky     Jan 1996  
c Modified by P. Tregoning    Oct 1996
c Modfied by R. King          Dec 1997   
c Modified by D. Lau & R King Aug 2009
c
c ######## Program to reformat rinex files into daily sessions ########  

c  Normally run using script sh_merge_rinex (see)

c  Runstring: merge_rinex [list_file] [doy] [output file last character] [sampling] [max_obs]   
c   
c            list_file - file containing a list of rinex files to be reformat.
c            doy - day of year of the output rinex file.
c            last character - last character of the output rinex file. 
c            sampling - output sampling interval in seconds  (optional but omission will
c                 leave the header sampling = 0)
c            max_obs - maximum number of observable to include (L1, L2, P1, P2, C1, D1)
c                 optional; if omitted, take all
c
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'

      logical debug 

      character*80 list_file,rxfile(150),rxfilex(150),outfile
     .           , time_line,time_line2,data_line
      character*3 out_day,out_int
      character*1 out_ext,cmax_obs,gnss 
      character*3 rxobtp(maxobt)
      character*3 rxtime
      character*20 rxpgm,rxusr,rxdat,rxobs,rctype,rcvers,anttyp
     .           , rcvnum,antnum
      character*40 rxagy
      character*60 rxcom(maxlin),rxmrk
      real*4 rxver
      real*8 sec,time,last_time,start_time(150),end_time,dum6,dum12,ttb
      real*8 rxint,apx,apy,apz,anth,ante,antn,rxsec0,rxsec1,ant_ht
     .     ,obs_int
      integer idoy,iext,irnxout,jdoy,odoy,nfiles,i,j,k,l,ioerr,idot
     .        ,yy,yy2,mm,dd,hh,minu,id,num,idum1,idum2,idum3
     .        ,idum4,idum5,idum7,idum8,idum9,idum10,idum11,istart
     .        ,jtb,julday,max_obs,nobtyp_old,lurinex(150)
     .        ,indx(150),isec10,kl,nl
      integer irxcom,nwave1,nwave2,nobtyp
     .       ,irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
     .       ,irxyr1,irxmo1,irxdy1,irxhr1,irxmn1  
c     number of satellites on header 
      integer*4 numsv
c     --satellite numbers for header
      integer*4 totsv(maxsat)
c     --number of obs of each type for each satellite for header
      integer*4 numobs(maxobt,maxsat)

      integer days_long     

      character*256 message
 
      logical decimate,first
      irnxout = 9

c
c maximum rinex file length is 31 days at a sampling rate of 30 seconds
      days_long = 31
c     
***** Remove old versions of the status, warning, and error files

      call report_stat('CLEAR','MERGE_RINEX',' ',' ', ' ',0)
 
c write some junk to the screen!!
      print*, '-----------------------------------'
      print*, ' Program "merge_rinex" Version 1.0 '
      print*, '-----------------------------------'
      print*, ' '
c
c read  command line input
      call rcpar(1,list_file)
      call rcpar(2,out_day)
      call rcpar(3,out_ext)
      call rcpar(4,out_int)
      call rcpar(5,cmax_obs)
c
c check if user knows how to use program.
      if (list_file.eq." ".and.out_day.eq." ".and.out_ext.eq." "
     .     .and.out_int.eq." ") then
        call help_rinex_merge
        stop
      endif
c
c convert output day and output ext to integers
      read (out_day,'(i3)',iostat=ioerr) odoy
      read (out_ext,'(i1)',iostat=ioerr) iext
      print*, 'Looking for RINEX data on DAY OF YEAR (doy): ',odoy
      print*, ' '
      read(out_int,'(f10.0)',iostat=ioerr) obs_int  
      if( obs_int.eq.0.d0 ) then
        print*,'No sampling interval input, output all data'
      else
        print*,'Output rinex file sampling at ',obs_int,' seconds. '
      endif
      print*,' '
      read(cmax_obs,'(i1)')max_obs
      if( max_obs.eq.0 ) then
         print*,'Outputting all observables'
      else
        print*,'Outputting ',max_obs,' observables.'
      endif
      print*,' '
c
c open list file
      open(unit=10, file=list_file,status='unknown',err=5)
      goto 6
5     print*, ' Cannot open the list_file ',list_file
      stop
6     continue
c
c read rinex files names and open rinex files given in the list_file
      nfiles = 0
      do 20 j = 1,150
        read(10,10,end=25) rxfile(j)
        nfiles = nfiles+1
10      format(a80)
        lurinex(j) = 10 + j
        open (unit=10+j, file=rxfile(j),status='old',err=15)
        goto 20
15      print*, ' cannot open the rinex file ',rxfile(j)
20    continue
25    continue   
      if( nfiles.eq.0 ) then
        print *,' No files in input list -- stop '
        stop
      elseif( nfiles.eq.1 ) then
        print *,' **Warning: only 1 file in input list'
      endif
c
c make output rinex file name
      idot = index(rxfile(1),'.')
      outfile = rxfile(1)(1:idot-5) // out_day // out_ext //
     .          rxfile(1)(idot:idot+3)
c
c open output rinex file
      open(unit=irnxout, file=outfile,status='unknown',err=30)
      goto 35
30    print*, ' cannot open the output rinexfile ',outfile
      stop
35     continue
c
c loop over rinex files to determine start and end epochs for the specified day
      istart = 0
      last_time = 0.d0
      ant_ht = 0.0000
      do 200 i = 1,nfiles
        urinex = lurinex(i)
        print*, 'Scanning RINEX file:',rxfile(i)
        call rrxhed( debug,gnss,
     .   rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .   rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz,
     .   anth,ante,antn,nwave1,nwave2,nobtyp,rxobtp,rxint,rxtime,
     .   idum1,idum2,idum3,idum4,idum5,dum6,
     .   idum7,idum8,idum9,idum10,idum11,dum12 )

c This version will not yet support RINEX 3

      if( rxver.gt.2.99 ) 
     .   call report_stat('FATAL','MERGE_RINEX','merge_rinex',' '
     .     ,'Cannot yet support RINEX 3: email rwk@chandler.mit.edu',0)
  
c
c check to make sure antenna heights are the same from file to file....
        if (ant_ht.ne.anth .and. i.gt.1) then
          print*,'Antenna heights different check rinex file ',rxfile(i)
          stop
        endif
        ant_ht = anth
c
c read rinex file data records
        first = .true.           
        do 100 k = 1,maxepc*days_long
          read(urinex,40,end=200) time_line
40        format(a80)   
cd          print *,'time_line ',time_line
          read(time_line,45,iostat=ioerr) yy2,mm,dd,hh,minu,sec,id,num
          if( ioerr.ne.0 ) then
            write(*,'(a80)') 'Bad line reading time line',time_line  
            stop 
          endif 
          yy= yy2    
c          print *,'1st read time yy num ',yy,num
          call fix_y2k(yy)
45        format (1x,5(i2,1x),f10.7,2(1x,i2)) 
c         read continuation of time line if #SVs > 12
          if( num.gt.12 ) then
            read(urinex,40,iostat=ioerr) time_line
cd            print *,'read time_line 2 ',time_line
          endif
c
c compute time of observation for checking if duplicate rinex records
          jtb = julday (mm,dd,yy)
          ttb = hh/24.d0 + minu/1440.d0 + sec/86400.d0
          time = jtb + ttb
c         save the start time for reordering files
          if( first ) then
            start_time(i) = time
            first = .false.
          endif
          if ( time .le. last_time) then
            print*, 'Warning rinex file has duplicate lines: ',rxfile(i)
c           compute the number of data lines for this epoch
            nl = num*(int((nobtyp-1)/5)+1)
cd            print *,'1 nobtyp num nl ',nobtyp,num,nl 
            do kl=1,nl
              read(urinex,340,end=200) data_line
cd              print *,'1 kl data_line kl ',kl,data_line
            enddo
            goto 100
          endif
          last_time = time
c
c compute DOY for current rinex record
          jdoy = idoy(yy,mm,dd)
c          print*, 'time_line ',time_line
c          print*, 'odoy, jdoy, yy, mm, dd ',odoy,jdoy,yy,mm,dd
c
          if ( odoy .eq. jdoy .and. istart .eq. 0 ) then
c            print*,'here first epoch'
            irxyr0 = yy
            irxmo0 = mm
            irxdy0 = dd
            irxhr0 = hh
            irxmn0 = minu
            rxsec0 = sec
            istart = 1
          elseif ( odoy .eq. jdoy .and. istart .eq. 1 ) then
c            print*,'last epoch '
            irxyr1 = yy
            irxmo1 = mm
            irxdy1 = dd
            irxhr1 = hh
            irxmn1 = minu
            rxsec1 = sec
          endif      
c         compute the number of data lines for this epoch
          nl = num*(int((nobtyp-1)/5)+1)
cd          print *,'2 nobtyp num nl ',nobtyp,num,nl 
          do kl = 1,nl
            read(urinex,40,end=200) data_line
cd            print *,'1 kl data_line',kl,data_line
          enddo

100     continue
200   continue
          print*, ' '
          print*, 'start ep: ',irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0
          print*, 'stop  ep: ',irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1

c  save the end time
          end_time = last_time
       print*,' end_time = ',end_time

c sort and rewind RINEX files before starting merge   
      call indexx(nfiles,start_time,indx) 
      do i = 1,nfiles
        rewind(lurinex(i))
        rxfilex(i) = rxfile(indx(i))
      enddo
      write(*,*) " "
      write(*,*) "RINEX files ordered by start_time: "
      write(*,*) " "
      do i=1,nfiles
        rxfile(i) = rxfilex(i)
        write(*,*) rxfile(i)
      enddo

c loop over rinex files merging data as we go!!!
      istart = 0
      last_time = 0.d0
      ant_ht = 0.0000
      do 500 i = 1,nfiles
        print*, ' '
        print*, ' Reading RINEX file: ',rxfile(i)
        print*, ' '
        urinex = lurinex(indx(i))
c        print *, 'debug urinex ',urinex
c       urinex stored in makex.h (--change this eventually)

c read rinex file headers
        call rrxhed( debug,gnss,
     .     rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .     rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz,
     .     anth,ante,antn,nwave1,nwave2,nobtyp,rxobtp,rxint,rxtime,
     .     idum1,idum2,idum3,idum4,idum5,dum6,
     .     idum7,idum8,idum9,idum10,idum11,dum12 )

c  set number of SVs to zero to indicate this (optional) entry not valid
         numsv = 0

c check that the requested output sampling rate is not smaller than the input
c sampling of the rinex file (if rxint=0 [=undefined] test will pass ok)
        if(obs_int.ne.0.d0 .and. obs_int.lt.rxint ) then
       write(message,'(a,f5.1,a,f5.1)') 'Requested output sampling rate'
     .       ,obs_int,' is less than input RINEX file sampling ',rxint  
          call report_stat('FATAL','MERGE_RINEX','merge_rinex',' '
     .       ,message,0)
        endif
c
c check to make sure antenna heights are the same from file to file....
        if (ant_ht.ne.anth .and. i.gt.1) then
        print*, ' Antenna heights different check rinex file ',rxfile(i)
        stop
        endif
        ant_ht = anth

c PT951005: allow option to limit rinex observables to only 5 (ie 1 line per
c           satellite. This is specifically for the SIO archive so that the
c           CORS rinex files are only 600Kb not 1.2Mb
         nobtyp_old = nobtyp
        if(max_obs.gt.0.)then
          nobtyp_old = nobtyp
          nobtyp = max_obs
        endif 
               
c       compute the number lines for each data entry
        nl = int((nobtyp-1)/5)+1  ! Number of lines for data
cd        print *,'2 nobtyp nl ',nobtyp,nl 

c write output rinex file header
        if ( i .eq. 1) then   
          call wrxhed (irnxout,gnss,
     .      rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .      rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz,
     .      anth,ante,antn,nwave1,nwave2,nobtyp,rxobtp,obs_int,rxtime,
     .      irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .      irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1,
     .      numsv,totsv,numobs )
        endif

        nobtyp = nobtyp_old
c                           
c read rinex file data records     
        do 400 k = 1,maxepc*days_long     
          read(urinex,340,end=500) time_line
340       format(a80)
          read(time_line,345) yy2, mm, dd, hh, minu, sec, id, num
          yy = yy2 
          call fix_y2k(yy)
345       format (1x,5(i2,1x),f10.7,2(1x,i2))    
c         if #SVs > 12, read the continuation of the time line
          if( num.gt.12 ) read(urinex,340,end=500) time_line2
c
c compute time of observation for checking if duplicate rinex records
          jtb = julday (mm,dd,yy)
          ttb = hh/24.d0 + minu/1440.d0 + sec/86400.d0
          time = jtb + ttb

c PT951023: To avoid dying on incomplete rinex files, stop merging when
c           the "last_time" is exceeded    
c RWK090820: Comment this out temporarily to avoid losing last record
c          if(time.ge.end_time)then
c            print *,'test last time end_time ',time,end_time
c            goto 500
c          else
            if ( time .le. last_time) then
              print*, 'Warning rinex file has duplicate lines: '
     .                ,rxfile(i)
              do kl = 1,nl
                read(urinex,340,end=500) data_line
              enddo
              goto 400
            endif
            last_time = time
c
c compute DOY for current rinex record
            jdoy = idoy(yy,mm,dd)
c           
c determine whether the epoch is aligned with the required epoch for the 
c           input decimation scheme.  For this purpose, round the time to the 
c           nearest 10 seconds so that RINEX sampling at 59.08 (TI4100) or UTC 
c           times will not  cause the decimation algorithm to miss the point 
c
c PT990608: I want to decimate a 5 second file to 30 seconds. Therefore, the
c           10 second rounding is allowing the 25 sec values through. Tighten
c           this down to 1 second - it will still allow TI data.
            if( obs_int.eq.0.d0 ) then
               decimate = .false.
            else
c 10 sec               isec10 = 10*(nint(sec/10.d0))
               isec10 = 1*(nint(sec/1.d0))
c --old        if (mod(minu*60.d0+nint(sec),obs_int).eq.0)then 
               if (mod(minu*60.d0+isec10,obs_int).eq.0)then
                 decimate = .false.
               else
                 decimate = .true.
               endif  
             endif
c            compute the number of data lines for this epoch
             nl = num*(int((nobtyp-1)/5)+1)
cd             print *,'3 nobtyp num nl ',nobtyp,num,nl 
            if ( odoy .eq. jdoy .and..not.decimate) then
              if (mod(minu,60).eq.0 .and. nint(sec).eq.0 ) then
                write(*,350) yy, mm, dd, hh, minu, sec
350             format(1x,'Writing data at: ',i4,4(i2,1x),f10.7)
              endif
              write(irnxout,340) time_line
              if( num.gt.12) write(irnxout,340) time_line2                    
              do kl=1,nl
                 data_line = ' '
                  read(urinex,340,end=500) data_line
cd                  print *,'2 kl data_line ',kl,data_line
                  write(irnxout,340) data_line
              enddo
            else 
              do kl=1,nl
                data_line = ' '
                read(urinex,340,end=500) data_line  
cd                print *,'2 kl skip data_line ',kl,data_line
              enddo
            endif
c RWK this finishes the commented-out if above
c          endif
400     continue
500   continue
c
      print*, 'Normal end in merge_rinex'
c
      stop
      end

c ## Subroutines ##
      subroutine help_rinex_merge
      write(*,*) "HELP for rinex_merge"
      write(*,*) " "
      write(*,*) "Given: "
      write(*,*) " (1) A file containing a list of rinex files. "
      write(*,*) " (2) A day of year--must be 3 digits"
      write(*,*) " (3) Output rinex file session number (char*1)."
      write(*,*) " (4) Output sampling interval (char*3) <optional> "
      write(*,*) " (5) Max # of observables (char*1) <optional> "
      write(*,*) " Create a rinex file for the day of year specified. "
      write(*,*) " "
      write(*,*) "returns rinex file containing all data from 1 day. "
      write(*,*) " "
      write(*,*) " requirements: "
      write(*,*) "               command line input "
      write(*,*) " "
      write(*,*) "usage: merge_rinex [list_file] [DOY] [ext] [int] "
     .          ,"[max_obs]"
      write(*,*) " "
      write(*,*) "list_file - file containing rinex files to be merged "
      write(*,*) "DOY - day of year for output rinex file "
      write(*,*) "ext - session no. (8th char) of the output rinex file"
      write(*,*) "int - output rinex file interval (in seconds) "
      write(*,*) "      optional to take all data, but header will have"
     .             ," interval=0"
      write(*,*) "max_obs - max # observables (L1, L2, P1, P2, C1, D1)"
     .           ,"        (omit to include all, 5 to exclude Doppler)"
      write(*,*) " "
      write(*,*) "example 1: merge_rinex rxfiles.input 016 a 30"
      write(*,*) " "
      write(*,*) "example 2: merge_rinex rxfiles.input 016 a 30 5"
      write(*,*) " "
      write(*,*) "Warning of current limitation:  RINEX files in list"
      write(*,*) "  must be time-ordered by start time or else some"
      write(*,*) "  data may be missed. "
      return
      end
