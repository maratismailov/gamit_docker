      Program read_sestbl_model

c     Short utility to extract sestbl entries for use in 
c     shell scripts, initially to be used by sh_sp3_fit
c     in getting models for the arc batch file. Code based
c     on fixdrv/rdsest.f   R. King 12 May 2012
           
c     Runstring 
c        read_sestbl_model <'sestbl entry'> <nout> 
c
c       where the 'sestbl entry' must be in quotes if multiple tokens,
c       and nout is the number of character expected in the output
c       string.  The program will write the matching entry to stdout, 
c       which can be redirected within the script. If the entry is not
c       found, the program returns a blank (not null) string.

      implicit none
     
      integer*4 iarg,lus,nin,nout,ieq,isc,end_of_line,ioerr
c      function
     .  , iclarg,nblen,index
                        
      character*7 sestbl
      character*120 arg,line,in_string,table_string,result,message
                
      logical eof,found
c      function
     .   , fcheck

      sestbl = 'sestbl.' 
      lus = 1
      if (fcheck(sestbl)) then
        open( lus, file=sestbl, status='old' )
      else
        call report_stat('FATAL','FIXDRV','read_sestbl',' '
     .               ,'sestbl. not available',0)
      endif         

c  get the entry requested
        
      iarg = iclarg(1,arg)
      if( iarg.le.0 ) then  
        call report_stat('FATAL','FIXDRV','read_sestbl',' '
     .               ,'Requested sestbl string missing',0)
      endif
      nin = nblen(arg)
      in_string = arg(1:nin)
cd      print *,'in_string ',in_string(1:nin)
      call uppers(in_string)   
cd      print *,'in string ',in_string(1:nin)    
      iarg = iclarg(2,arg) 
      if( iarg.le.0 ) then  
        call report_stat('FATAL','FIXDRV','read_sestbl',' '
     .               ,'Runstring missing nout ',0)
       else
         read(arg,*) nout
       endif

c  search the sestbl for the requested entry
      
      eof = .false.                   
      found = .false.                  
      do while( .not.found .and. .not.eof )
        read(lus,'(a)',iostat=ioerr) line 
        table_string = line(1:nin)
cd        print *,'table_string ',table_string(1:nin)
        call uppers(table_string)
cd        print *,'table_string ',table_string(1:nin)
        if( in_string(1:nin) .eq. table_string(1:nin) ) then
cd          print *,'found entry, nin ',nin
          found = .true.
          ieq = index( line(nin+1:80), '=' )
          isc = index( line(nin+1:80), ';' )
          if( isc.ne.0 ) then
            end_of_line = nin + isc  -1
          else
            end_of_line = nin + ieq + nout + 1
          endif
          result = ' '   
cd          print *,'ieq isc end_of_line ',ieq,isc,end_of_line 
          if (ieq .ne. 0 ) then
            if( (nin+ieq+nout+1).gt.120 )
     .      call report_stat('FATAL','READ_SESTBL','read_sestbl',' '
     .                      ,'Output line too large',0)
            result(1:nout) = line(nin+ieq+2:end_of_line) 
            call uppers(result)
          else
            write(message,'(2a)') 'sestbl. format improper: '
     .         ,in_string(1:nin)
            call report_stat( 'WARNING','READ_SESTBL','read_sestbl',' '
     .                      , message,0)   
          endif 
cd          print *,'Writing result, nout ',nout
cd          print *,' ' 
          write(*,'(a)') result(1:nout)
       
        elseif( ioerr.eq.-1 ) then
          eof = .true.
          result = ' ' 
          write(*,'(a)') result 
          stop
        elseif( ioerr.lt.0 ) then
          call report_stat('FATAL','READ_SESTBL','read_sestbl',' '
     .                      ,'Error reading sestbl.',0)
        else
          continue
        endif
      enddo
      
      stop
      end


