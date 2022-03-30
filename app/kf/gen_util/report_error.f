ctitle
 
      subroutine report_error(type,ierr,operation,file,terminate,prog)
 
 
      implicit none 

c
c     Routine to report a file manipulation error and possibly stop
c     program from running.  The error meassage is printed to the
c     log lu.
c
c MOD JLD 880421 To print a more concise message
c
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
      character*16 error_name
 
c
c Local variables
c ---------------
c trimlen -- HP utility for length of string
c ilog -- the log lu
c error -- logical which indicates if there is an error
c
 
      integer*4 TrimLen, ilog, len_file
 
*   id  - Dummy argument for loglu
*   loglu   - HP function to return user's LU
      integer*4 id, loglu
 
c
      logical error
 
c
c.... Get the user terminal
      ilog = loglu(id)
c
c.... Get the error type in upper case
      error_name = type
      call CaseFold(error_name)
c
c.... See if an FMGR error (<0)
      error = .false.
      len_file = max(1,trimlen(file))
 
      if (error_name(1:2) .eq. 'FM' .and. ierr .lt. 0) then
        error = .true.
        write(ilog,100) type,ierr,operation,file(1:len_file),prog
  100   format(/,X,A,I6.5,' occurred ',A,'ing ',A,' in routine ',A)
      end if
c
c.... Other error?
      if (ierr .ne. 0 .and. error_name(1:2) .ne. 'FM') then
        error = .true.
        write(ilog,150) type,ierr,operation,file(1:len_file),prog
  150   format(/,X,A,X,I6,' occurred ',A,'ing ',A,' in routine ',A)
      end if
c
*                               ! no error
      if( .not.error ) return
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
 
c.......................................................................
