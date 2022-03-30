      subroutine rsesfo (lu,debug,check_flg,iyr,iday,isessn
     .   , hour,min,interval,nepoch,isatcount,satarray )

c     Read the session.info (scenario) file to get the start, stop times
c     and the satellite list.  Prior to August, 1999, the file was rigidly
c     formatted, as shown below.  It now is read free-format, with optional
c     comment lines indicated by a non-blank first-column; in the new-style
c     files, the format is ignored, and the first two tokens are checked to see
c     if they are reasonable values for year and day-of-year

c     Last modified by R. King 060815 to return start time in hr min, not GPS week.
c     This latter used by makex but nothing else, so do conversion there. 
     
c     Input
c     -----    
c     lu         I*4      unit number                                
c     debug      Logical  print debug     
c     check_flg  I*4      = 0  return info for a single-line table, no checking (MAKEXP, FIXDRV/seschk)
c                         = 1  check only the day, not the yr or session (AUTECL, MODEL for simred)
c                         = 2  check yr, day-of-yr, and session (MAKEX, FIXDRV/somake, MODEL, BCTOT)
c     iyr        I*4      4-digit year of session  (output if check_flg = 0 )
c     iday       I*4      day-of-year of session   (output if check_flg = 0 )
c     isessn     I*4      session number           (output if check_flg = 0 )
c
c     Output
c     ------                                           
c     hour             I*4      hr of start time
c     min              I*4      minutes of start time  
c     interval         R*8      Sampling interval (seconds)
c     nepoch           I*4      # epochs
c     isatcount        I*4      number of satellites
c     satarray(maxsat) I*4      PRN #s of satellites

    
c   Original format, used prior to 23 Feb 1992:

c  Yr Day SN Int  Max       HH MM     Satellites
c  (i2,1x,i3,1x,a,1x,i4,x,i4,7x,i2,1x,i2,3x,24i3)
c  91 365  1  10  540        3  8     2  6 11 16 18 19

c   Second format, from 23 Feb 1992:

c  Yr Day SN Int  Max       HH MM     Satellites
c  (i2,1x,i3,1x,i2,1x,i3,x,i4,7x,i2,1x,i2,3x,24i3)
c  91 365  1  10  540        3  8     2  6 11 16 18 19
           
c   This format differs from the original one in reading the session
c   as an integer rather than an alphameric.  As originally designed,
c   both the pre- and post-1992 formats can have different size fields 
c   for any of the variables.  In practice, a typical old format has the 
c   session number as 'a1' (vs 'i2' for the new format), and the sampling 
c   interval 'i4' (vs 'i3' for the later format), so that the positions of 
c   the rest of the fields ('Max' [epochs] thru 'satellites') are the same 
c   as the later format.  With this version of the subroutine, we read both
c   of these formats, as well as the current one with a free format read.
    
c    Current file is free format, e.g. 

c # Session.info : free format, non-blank first column is comment
c # Year Day  Session   Interval  #Epochs  Start hr/min   Satellites
c   1991  365    1        120        225        3   8       2  6 11 16 18 19
         
      implicit none

      include '../includes/dimpar.h'  
      include '../includes/makex.h'
                   
      character*80   ,prog_name
      character*256  line
      character*256  message
      integer*4      lu,iscrn,hour,min,
     .               satarray(maxsat),
     .               isatcount,i,
     .               iyr,iyr1,iday,iday1,isessn,isessn1,ioerr,
     .               nepoch,interval,count_arg,narg,
     .               len,rcpar,check_flg

      logical debug,found

      parameter(iscrn=6)

c   Original comments (1992?): 
c           We input iy,iday,hh,mm in UTC, which the receiver assumes to
c           be in GPS time.  The call to TIMCON converts
c           these to week number and seconds of week format, but because
c           of the receiver's naive assumption, we are now off by the
c           number of leap seconds between UTC and GPS time.
c           To correct, I subtract the UTC offset, and the resulting time will
c           be in the standard form for this program, namely GPS week number
c           and second of week.     
c**  RWK 060815: I don't think that users assume that the input start time is UTC;
c           rather usually start on even minutes, assumed to be the same as the
c           RINEX time-tags (GPST).  So remove the conversion below and simply
c           return the read quantities.

                   
c     get calling program name for report_stat
      len = rcpar(0,prog_name)

c ------ Initialization
          
      found = .false.
      do  i=1,maxsat
         satarray(i)=0
      enddo 
c     With current code, the input session number from MAKEX or FIXDRV
c     should never be zero, but check for safety
      if( isessn.eq.0 ) then
        call report_stat('WARNING',prog_name,'lib/rsesfo',' '
     .          , 'Session number is zero; setting to 1',0)
        isessn = 1
      endif     
      rewind (lu)    
                 

c ----- Begin loop reading lines of the file    

      do while (.not.found ) 
        read(lu,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1) then  
          write(message,'(a,i4,1x,i3,1x,i1,a)') 
     .      'Requested yr day session ',iyr,iday,isessn
     .      ,' not found on session.info file'  
          call report_stat('FATAL',prog_name,'lib/rsesfo'
     .         ,'session.info',message,ioerr)  
        elseif( ioerr.ne.0 ) then
           write (message,'(a,i4,1x,i3,1x,i1)') iyr,iday,isessn
     .       ,'Error reading session.info file for yr day session'
     .         , iyr,iday,isessn 
          call report_stat('FATAL',prog_name,'lib/rsesfo'
     .         ,'session.info',message,ioerr)
          return 
        else     
c         read the line to see if valid or comment  
          iyr1 = 0
          iday1 = 0  
          read(line,*,iostat=ioerr) iyr1,iday1    
c         a non-valid data line is any that doesn't have a reasonable year 
          if( debug ) print *,'RSESFO 1 ',iyr1,iday1
          if( iyr1.lt.80 .or. iyr1.gt.2100   .or.
     .      iday1.le.0 .and. iday1.gt.366 ) then
            continue
          else
c           if valid year, day-of-year, process the line
c           determine the number of satellites by counting tokens on the line
            narg = count_arg(line)     
c           there are 7 non-SV entries (yr day sn int nepoch hh mm)
            if( narg.lt.8 ) call report_stat('FATAL',prog_name
     .         ,'lib/rsesfo',' ', 'Too few session.info entries',0)
            isatcount = narg - 7 
            if(isatcount.gt.maxsat)then
              write(message,'(a,i3,2a,i3,a)')
     .            'Too many satellites:',isatcount,'. Increase the'
     .            ,' dimension of maxsat (',maxsat,') in dimpar.h'
              call report_stat('FATAL',prog_name,'lib/rsesfo',' '
     .                      , message,0)
            endif
            read(line,*,iostat=ioerr) iyr1,iday1,isessn1
     .        ,  interval,nepoch,hour,min,(satarray(i),i=1,isatcount)
            if( debug)  print *,'RSESFO 2'
     .        ,iyr1,iday1,interval,nepoch,hour,min,isatcount
c           always change session 0 to session 1 internally
            if(isessn1.eq.0) isessn1=1  
c           convert to 4-digit years
            call fix_y2k(iyr1)
c           ----  Match the input with the session.info entries
            if (debug) print *,'RSESFO 3 check_flg iyr iday isessn '
     .                    ,check_flg,iyr,iday,isessn
            if (debug) print *,'RSESFO 3 iyr1 iday1 isessn1 '
     .                      ,iyr1,iday1,isessn1
            if( check_flg.eq.1 .and. iday1.eq.iday ) then
              found = .true. 
              iyr = iyr1 
              if( debug ) then
                print *,'Found iyr iday ',iyr,iday
                print *,'check_flg nsat ',check_flg,isatcount
                print *,'sats ',(satarray(i),i=1,isatcount)
              endif
            elseif( check_flg.eq.2 .and. 
     .      (iyr1.eq.iyr.and.iday1.eq.iday.and.isessn1.eq.isessn) ) then
              found = .true.   
c             check number of epochs
              if (nepoch .gt. maxepc) then
                write(message,'(a,i5)')
     .           'Too many epochs; setting nepoch =,',maxepc
                call report_stat('WARNING',prog_name,'lib/rsesfo',' '
     .                        , message,0)
                nepoch = maxepc
              endif  
            elseif( check_flg.eq.0 ) then 
c             no checking but set found = .true. to exit gracefully
              found = .true.
              iyr = iyr1
              iday = iday1
              return
            else
c             wrong day/session, go read another record
            endif  
c         endif for valid data line
          endif
c       endif for error reading line
        endif                       
c     end of loop over records
      enddo
      return    
      end
