CTITLE 'MULTIREAD'
 
      subroutine multiread( line, indx, type, ierr, ivalue,
     .                      cvalue, num )
 
      implicit none

* MOD TAH 090923: Allow num=0 and read to end of line and 
*     return all values.  num contains number of values returned.
* MOD TAH 091030: Added feature if num<0 then read to that number

* 
****  NOTE: code assumes max string length of 256 characters

*     Routine to multiply call read_line and get 'num' values
*     from line.  Errors are always reported by this routine.
 
*         ivalue(*)     - Values to be read from line.  Dummy
*                       - declarer which is used for position
*                       - only.
*          dum          - Dummy value for read_line
* MOD TAH 200709: Added dimension to dum to make the same as tval
*     (Issue with latest gfortran compiler).
      integer*2 ivalue(*), dum(4)
      integer*2 tval(4)   ! Used to test next value read.
 
*         ierr         - IOSTAT error on read
*          num          - number values to be read.
*          indx         - Current pointer in line
*          inc          - increment in ivalue for each position
*                       - read.
*          i,j          - Loop counter
 
      integer*4 ierr, num, indx, inc, i, j
 
*       error_out       - Indicates error reported
 
      logical error_out
 
*             line      - Line to be read
*          cvalue(1)    - Character values to be read
*          type         - Type of argument.  Values can be
*                       - I2 -- integer*2
*                       - I4 -- integer*4
*                       - R4 -- real*4
*                       - R8 -- real*8
*                       - CH -- Character
      character*(*) line, cvalue(*), type
      character*256 tchar
 
*           cdum        - Dummy argument for readline
 
      character*1 cdum
 
****  Get increment we need in ivalue for type of data
      inc = 0
      call casefold( type )
      if( type(1:2).eq.'I2' ) inc = 1
      if( type(1:2).eq.'I4' ) inc = 2
      if( type(1:2).eq.'R4' ) inc = 2
      if( type(1:2).eq.'R8' ) inc = 4
 
      error_out = .false.
 
***** Do characters and other values separately
*                                     ! Characters
      if( type(1:2).eq.'CH' ) then
* MOD TAH 090923: If num is <=0 then read until error or abs(num) entries
          if( num.gt.0 ) then 
              do i = 1, num
                  call read_line(line,indx, type, ierr, dum, cvalue(i))
                  if( .not.error_out .and. ierr.ne.0 ) then
                      call report_error('IOSTAT',ierr,'decod', line, 0,
     .                                  'MULTIREAD')
                      error_out = .true.
                  end if
              end do
          else  ! loop until error
              ierr = 0
              i = 0
              do while ( ierr.eq.0 ) 
                  call read_line(line,indx, type, ierr, dum,tchar)
                  if( ierr.eq.0 ) then
                      i = i + 1 
                      cvalue(i) = tchar
                  end if  
              end do
              num = i
          end if 
*                                     ! Do numeric data
      else
* MOD TAH 090923: If num is <=0 then read until error
          if( num.gt.0 ) then 
              do i = 1, num
                  call read_line(line,indx, type, ierr,
     .                           ivalue(1+(i-1)*inc), cdum )
                  if( .not.error_out .and. ierr.ne.0 ) then
                      call report_error('IOSTAT',ierr,'decod', line, 0,
     .                                  'MULTIREAD')
                      error_out = .true.
                  end if
              end do
          else
*             Read until we hit an error or i reaches abs(num) if nonzero
              ierr = 0
              i = 0
              do while ( ierr.eq.0 ) 
                  call read_line(line,indx, type, ierr, tval, cdum )
                  if( ierr.eq.0 ) then
                      i = i + 1
                      do j = 1,inc
                         ivalue(j+(i-1)*inc) = tval(j)
                      end do
                      if( i.eq. abs(num) ) ierr = -1
                  endif
              end do
              num = i
          endif
 
      end if
 
***** Thats all
      return
      end
 
