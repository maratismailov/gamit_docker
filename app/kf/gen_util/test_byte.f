      program test_byte

      implicit none 
 
*     Outputs a grid of byte values
 
      byte xy(16384)
 
*         dim   - Dimension of array
*       i,j     - Loop counters
*       ierr    - IOSTAT error
*       size    - Size of grid
 
      integer*4 i,j, ierr, size, iel
 
*      val      - Value to put in grid
 
      real*8 val
 
*   out_name    - Name of output
 
      character*128 out_name
 
      write(*,100)
 100  format('Program to generate a byte file',/,
     .       'Enter gridsize (<=128) '$)
      read(*,*) size

      write(*,140)
 140  format('Enter the name of output file '$)
      read(*,'(a)') out_name
      open(200, file=out_name, iostat=ierr, form='unformated',
     .     access='direct', recl=size**2)
      if( ierr.ne.0 ) then
          write(*,120) ierr, out_name
 120      format(' IOSTAT error ',i4,' opening ',a)
          stop
      end if
 
*     Generate grid
      do i = 1, size
          do j = 1, size
              iel = (i-1)*size + j
              val = sin(i/5.)*cos(j/5.)*110. + 128
              xy(iel) = val
          end do
      end do
 
*     Now write
      write(200,rec=1) (xy(i),i=1,size**2)
      close(200)
      end
 
 
 
