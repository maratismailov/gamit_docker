CTITLE    ................................................................
 
      subroutine report_error(type,ierr,operation,file,terminate,prog)
c
c     Routine to report a file manipulation error and possibly stop
c     program from running.  The error meassage is printed to the
c     log lu.
c
c Variables
c ---------
c type -- the type of error (usually FMGR or IOSTAT)
c ierr -- the error which ocurred
c operation -- what we were doing at the time
c file -- the file being manipulated
c terminate -- if non-zero causes the program to stop.
c prog -- the program which was running.
c
      character*(*) type, operation, file, prog
 
c
      integer*4 ierr, terminate
 
c
c Local variables
c ---------------
c trimlen -- HP utility for length of string
c ilog -- the log lu
c error -- logical which indicates if there is an error
c len   -- length of string to be output
 
      integer*4 trimlen, ilog, len
 
c
      logical error
 
c
c.... Get the user terminal
      ilog = 6
c
c.... see if error  (for FMGR less than 0 is error, for anything
c     else a non zero value is error
      error = .false.
      if( type(1:2).eq.'FM' ) then
         if( ierr.lt.0 ) error = .true.
      else
         if( ierr.ne.0 ) error = .true.
      end if
c
*                               ! no error
      if( .not.error ) return
c
c.... Report error
      len = max(1,trimlen(file))
      write(ilog,100) type, ierr, operation, file(1:len),
     .   prog
 100  format(/,1x,a," error ",i4," occurred ",a,"ing ",a,
     .   " in routine ",a)
c
c.... see if we should stop
      if( terminate.gt.0 ) then
         write(ilog,110) prog
 110     format(1x,a," terminating")
         stop ' File error reported by REPORT_ERROR'
      end if
c
      return
      end
