CTITLE READ_DFILE
 
      subroutine read_dfile( runstring, nc, cfiles, max_cfiles )

      implicit none
 
*     This routine will read the names of the cfiles from the
*     d-file.  In this case no path is expected.
 
* PASSED VARIABLES
 
*   nc      - Current number of c-files (normally would be 0)
*   max_cfiles  - Maximum number of cfiles allowed.
 
      integer*4 nc, max_cfiles
 
*   runstring       - Contains the name of d-file
*   cfiles(max_cfiles)  - List of names of cfiles
 
      character*(*) runstring, cfiles(max_cfiles)
 
*     LOCAL VARIABLES
 
*   i       - Loop counter
*   ierr    - IOSTAT error for file reading
*   ncd     - Number of cfiles from the d-file
 
      integer*4 i, ierr, ncd
 
*   dummy       - Dummy lines for reading the d-file
 
      character*20 dummy
 
****  Start open the dfile
 
      open(100, file=runstring, iostat=ierr, status = 'old')
      call report_error('IOSTAT',ierr,'open',runstring,1,'READ_DFILE')
      if( ierr.ne.0 ) RETURN
 
****  Loop over the beginning records of the dfile down to where
*     cfile number is
      do i = 1, 6
          read(100,'(a)', iostat=ierr) dummy
      end do
      call report_error('IOSTAT',ierr,'read',runstring,1,'READ_DFILE')
 
*     If no error get nuber of cfiles
      if( ierr.ne.0 ) RETURN
 
      read(100,*, iostat=ierr) ncd
      call report_error('IOSTAT',ierr,'read',runstring,1,'READ_DFILE')
      if( ierr.ne.0 ) RETURN
 
      do i = 1, ncd
          read(100,'(a)', iostat=ierr) cfiles(nc+i)
          call report_error('IOSTAT',ierr,'read',cfiles(nc+i),
     .                1,'READ_DFILE')
*         Force name to cfile name
          if( ierr.eq.0 ) then
              cfiles(nc+i)(1:1) = 'c'
          else
*             Skip this file
              nc = nc - 1
          end if
      end do
 
      close(100)

* MOD TAH 970502: If the cfile names come only from the dfile
*     i.e., nc =0 on entry, then sort thw list
      if ( nc.le.0 ) call csort(nc+ncd, cfiles) 

      nc = nc + ncd
 
****  Thats all
      return
      end
 
 
CTITLE csort
 
      subroutine csort( num, clist)

      implicit none
 
*     This routine uses an exchande sort algormithm to sort
*     the list clist into ascending order.  There are num values
*     in clist.
 
*   num     - Number of values to be sorted
*   clist(num)  - List to be sorted in to ascending order.
 
      integer*4 num
      character*(*) clist(num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   smallest_one    - Smallest integer in current pass.
*   cswap   - Value used to swap charcter strings 
 
 
      integer*4 i,j, smallest_one
      character*256  cswap
 
****  Start loop using exchange sort
 
      do i = 1, num
          smallest_one = i
          do j = i+1, num
              if( clist(j).lt. clist(smallest_one) ) then
                  smallest_one = j
              end if
          end do 
 
*****     See if we should swap
          if( smallest_one.gt. i ) then
              cswap = clist(smallest_one)
              clist(smallest_one) = clist(i)
              clist(i) = cswap
          end if
      end do
 
***** Thats all.  Now sorted in ascending order
      return
      end
 
