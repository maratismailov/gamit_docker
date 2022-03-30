      subroutine read_rcvant( dirn,type,antcod_in,antnam,radome
     .                      , rcvcod,rcvnam,pcncod )
c
c PURPOSE: This subroutine will exchange the antenna code for a full antenna name
c          and vice versa. The definitive place of Gamit where all antennas and receivers
c           are listed is in a data file rcvant.dat in gamit/tables. The
c          intention is that all antennas/receivers ever used will be uniquely identified by
c          their 6 character ant/rec code. Thus, passing a 6 char  code
c          through all gamit routines will suffice, as it can be turned into a
c          complete name at any time by calling this routine.
c
c          This routine will open the file rcvant.dat, read it for the required information
c          then close it again.
c
c VARIABLES:
c     dirn     - direction for information transfer: 1: code name  -> full name
c                                                    2: full name -> code cod   I*4
c     type     - flag for receiver or antenna request 1: antenna; 2: receiver   I*4
c     antcod   - antenna code used in gamit                                     C*6
c     antnam   - full antenna name                                              C*20  
c     radome   - radome name (same for gamit and IGS)                           C*5
c     rcvcod   - antenna code used in gamit                                     C*6
c     rcvnam   - full antenna name                                              C*20 
c     pcncod   - character indicating need for C1 or P2' bias correction        C*1
c                ( see gamit/makex/set_dcb_flag.f )


c P. Tregoning   12th September, 1995
c S. McClusky    20th October, 1995  
c R. King         5th February 2003
c T. Herring     10th November 2003: Fixed up opening rcvant.dat

      implicit none

      integer*4 dirn,type,i,ioerr,irec,indx,len,rcpar
                                    
      character*1 pcncod    
      character*5 radome
      character*6 antcod_in,antcod,rcvcod
      character*20 antnam,rcvnam,antenna(2),receiver(2)
      character*80 line,answer,prog_name
      character*256 message
* NOD TAH 031110: Addeded home_dir so that <home_dir>/gg/tables/rcvant.dat can be opened
*     if their is not a local copy.  Also need trimlen for length of string
      integer*4 trimlen
      character*256 home_dir, home_rcvant

      logical found,stopflag
      parameter (irec = 90)

cd      print *,'READ_RCVANT dirn,type,antcod_in,antnam '
cd     .      , 'radome rcvcod,rcvnam,pcncod '
cd     .      , dirn,type,antcod_in,antnam,radome,rcvcod,rcvnam,pcncod

c     get the module name for report_stat calls
      len = rcpar(0,prog_name)   
c     if mstinf2, just issue warning
      if( prog_name(1:6).eq.'mstinf'.or.
     .    prog_name(1:5).eq.'hfupd' .or.
     .    prog_name(1:6).eq.'conver' ) then
        stopflag = .false.
      else
        stopflag = .true.
      endif                                     
               
c   initialize controls and arrays
      ioerr = 0   
      answer = ' '
      found = .false.
      i = 1

cd      print *,'READ_RCVANT antcod antnam radome '
cd     .       ,antcod_in,antnam,radome
cd       write(*,10) dirn,type,antcod_in,antnam,radome,rcvcod,rcvnam
cd10     format(i2,i2,1x,a6,2x,a20,1x,a5,1x,a6,1x,a20)

c  open the file rcvant.dat
* MOD TAH 031110: Close unit 90 first.  It is opened someplace else and
*     this cause error on HPUX system
c      close(unit=90,iostat=ioerr)
      open(unit = irec,file='rcvant.dat',status='old',iostat=ioerr)     
      if(ioerr.ne.0) then
* MOD TAH 031110: See if we can open a copy through the <home_dir>/gg/tables
*       Get the user's home directory:
        call getenv('HOME',home_dir)
        home_rcvant = home_dir(1:max(1,trimlen(home_dir))) //
     .                '/gg/tables/rcvant.dat'
        open(unit=irec,file=home_rcvant,status='old',iostat=ioerr)
      end if
      if(ioerr.ne.0 ) then
        write(message,'(2a)')
     .     'Cannot open file rcvant.dat'
     .     ,'-- may need link to gamit/tables/rcvant.dat'
        call report_stat('FATAL',prog_name,'lib/read_rcvant',' '
     .                  ,message,ioerr)
      endif
     
c now match up the receiver/antenna information as requested by the arguments
c passed to the subroutine

      if ( type.eq.1 ) then
c        it is an antenna request

c get the full RINEX name from the GAMIT code
        if( dirn.eq.1 ) call ant_alias(antcod_in,antcod)

c move through the file until we find the 'ANTENNA' block
15      do while(answer(1:7).ne.'ANTENNA')
          read(irec,'(a80)',iostat=ioerr) line
          if(ioerr.ne.0)  call report_stat('FATAL',prog_name
     .      ,'lib/read_rcvant',' '
     .      ,'Error looking for ANTENNA in rcvant.dat',ioerr)
          if (line(1:1) .ne. ' ' ) goto 15
          indx = index(line,'ANTENNA')
          read(line(indx:indx+6),'(a7)',iostat=ioerr) answer   
          if(ioerr.ne.0)  call report_stat('FATAL',prog_name
     .      ,'lib/read_rcvant',' '
     .      ,'Error reading line for ANTENNA in rcvant.dat',ioerr)
        enddo

20      do while (.not.found)
          read(irec,'(a)',iostat=ioerr) line   
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/read_rcvant',' '
     .       ,'Error reading antenna lines in rcvant.dat',ioerr)
          if (line(1:1) .ne. ' ' ) goto 20
          if (line(2:4) .eq. 'END' ) then  
            if( dirn.eq.1 )  write(message,'(a,a6,a,a6,a)') 
     .         'Input antenna type ',antcod_in,' with alias '
     .          ,antcod,' not in rcvant.dat' 
            if( dirn.eq.2 )  write(message,'(a,a20,a)')
     .         'Antenna name ',antnam,' not found in rcvant.dat'  
            if( stopflag ) then 
               call report_stat('FATAL',prog_name,'lib/read_rcvant'
     .                       ,' ',message,0)     
             else
               call report_stat('WARNING',prog_name,'lib/read_rcvant'
     .                       ,' ',message,0)    
               if( dirn.eq.1 ) antnam = ' '
               if( dirn.eq.2 ) antcod_in = ' '     
               close(irec)
               return
             endif 
          endif

          read(line,'(1x,a6,8x,a20)',iostat=ioerr) antenna(1),antenna(2)
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .     ,'lib/read_rcvant',' ','Error in antenna table format',ioerr)


          if (dirn.eq.1)then

c  we want to convert codnam into full name
            if(antenna(1)(1:6).eq.antcod) then
              found = .true.
              antnam = antenna(2)               
              if( radome(1:4).eq.'UNKN' .or. radome(1:4).eq.'DOME') then
c               check for radome types built into GAMIT codes
                if( antcod.eq.'ASHDMR' ) radome(1:4)='DOME'
                if( antcod.eq.'ASHDMD' ) radome(1:4)='DOME'
                if( antcod.eq.'LC303R' ) radome(1:4)='DOME'
                if( antcod.eq.'NO503R' ) radome(1:4)='DOME'
                if( antcod.eq.'TPSC3R' ) radome(1:4)='SNOW'
                if( antcod.eq.'TPSCC4' ) radome(1:4)='SNOW'
                if( antcod.eq.'AERCHR' ) radome(1:4)='DOME'
              endif
              antnam(17:20) = radome(1:4) 
c              print *,'found antcod antnam ',antcod,antnam
            endif

          elseif(dirn.eq.2)then

c  want to convert full name into codnam
c  RWK 160419: Including column 16 is a problem if there is a '+' there
c           if(antenna(2)(1:16).eq.antnam(1:16)) then                  
            if(antenna(2)(1:15).eq.antnam(1:15)) then                  
              found = .true.
              antcod_in = antenna(1)(1:6)       
c             for a few cases, assign the radome from the name
              if( antnam.eq.'TPSCR3_GGD SNOW ') radome='SNOW '
              if( antnam.eq.'TPSCR4 SNOW     ') radome='SNOW '
c put the radome name in into the end of the antenna name slot
              antnam(17:20) = radome(1:4)            
            endif
          endif
        enddo

        close(irec)  
        return

      elseif( type .eq. 2 ) then

c  it is a receiver request

c move through the file until we find the 'RECEIVER' block
120     do while(answer(1:8).ne.'RECEIVER')
          read(irec,'(a)',iostat=ioerr) line   
          if( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .        ,'lib/read_rcvant',' '
     .         ,'Error looking for RECEIVER in rcvant.dat file ',ioerr)
          if (line(1:1) .ne. ' ' ) goto 120
          indx = index(line,'RECEIVER')
          read(line(indx:indx+7),'(a8)',iostat=ioerr) answer  
          if( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .       ,'lib/read_rcvant',' '
     .       ,'Error reading line for RECEIVER in rcvant.dat ',ioerr)
        enddo

130     do while (.not.found)
          read(irec,'(a)',iostat=ioerr) line     
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/read_rcvant',' '
     .         ,'Error reading receiver lines in rcvant.dat',ioerr)
          if (line(1:1) .ne. ' ' ) goto 130
          if (line(2:4) .eq. 'END' ) then             
            if(dirn.eq.1) write(message,'(a,a6,a)') 'Receiver code '
     .               ,rcvcod,' not found in rcvant.dat'
            if(dirn.eq.2) write(message,'(a,a20,a)') 'Receiver name '
     .               ,rcvnam,' not found in rcvant.dat'   
            if( stopflag ) then
              call report_stat('FATAL',prog_name,'lib/read_rcvant',' '
     .                      ,message,0)
            else  
              if( dirn.eq.2.and.rcvnam(1:6).ne.'------') 
     .        call report_stat('WARNING',prog_name,'lib/read_rcvant',' '
     .                      ,message,0) 
             if( dirn.eq.1 ) rcvnam = ' '
             if( dirn.eq.2 ) rcvcod = ' ' 
             close(irec) 
             return  
            endif
          endif                              

          read(line,'(1x,a6,8x,a20,1x,a1)',iostat=ioerr) 
     .             receiver(1), receiver(2),pcncod 
          if(ioerr.ne.0) call report_stat('FATAL',prog_name
     .       ,'lib/read_rcvant',' '
     .         ,'Error in receiver tables format',ioerr)

          if (dirn.eq.1)then

c  we want to convert codnam into full name
            if(receiver(1)(1:6).eq.rcvcod)then
              found = .true.
              rcvnam = receiver(2)
            endif

          elseif(dirn.eq.2)then

c  want to convert full name into codnam
            if(receiver(2).eq.rcvnam)then
              found = .true.
              rcvcod = receiver(1)(1:6)
            endif
          endif
        enddo

        close(irec) 
        return

      else

        call report_stat('FATAL',prog_name,'lib/read_rcvant',' '
     .                  ,'Unknown request type',0)

      endif
      end
















