      program memtest

      implicit none 

*     Program to memory allocation using the methods used in globk.
*     To use this program you need:
*     memtest.f and mallocg.c 
*     To install:
*     gcc -c mallocg.c 
*     g77 -o memtest memtest.f mallocg.o
*     memtest
* SAMPLE OUTPUT with memory size passed.
* dell4[277] memtest 128
* MEMTEST: Allocating   128 Mega-I*4 of memory
* Runstring: memtest <Mega I*4 words>
* NOTE: On a 32-bit machine, you can not go past 512
* Memory allocated
* Testing memory allocation: First 10 words
* Indx  Element    Value
*   1          1          1
*   2          2          2
*   3          3          3
*   4          4          4
*   5          5          5
*   6          6          6
*   7          7          7
*   8          8          8
*   9          9          9
*  10         10         10
* Testing memory allocation: Last 10 words
* Indx  Element    Value
*   1  134217719  134217719
*   2  134217720  134217720
*   3  134217721  134217721
*   4  134217722  134217722
*   5  134217723  134217723
*   6  134217724  134217724
*   7  134217725  134217725
*   8  134217726  134217726
*   9  134217727  134217727
*  10  134217728  134217728

* PROGRAM VARIABLES

      integer*4 vma_data(1)  ! Array that will be used to hold data
     .,         memsize_mI4  ! Number of Mega-I*4 words to allocate
     .,         i,j          ! Loop counter to write values into memory
     .,         jerr         ! IOSTAT error reading requested size to allocate
                             ! If there is an error 10 MI4 values are allocated
     .,         numwords, sizeword  ! Number of words of sizeword to be allocated
                             ! (Structure is used so that more than 2^31 bytes of
                             ! memory can be allocated while passing values as
                             ! integer*4)

      integer*8 istart_vma   ! Index in vma_data that corresponds to address 
                             ! with memory allocated
     .,         memassign    ! Integer*4 functiont to returns element number
                             ! in vma_data where memory can be allocated


      character*128 arg      ! First argument passed to program.

****  get the amount of memory to allocate
      call getarg(1,arg)

*     See if geterg returned the program name
      if( index(arg,'memtest').gt.0 ) then
         call getarg(2,arg)
      end if
      read(arg,*,iostat=jerr) memsize_mI4
      if( jerr.ne.0 ) then
          memsize_mI4 = 10
      endif

      write(*,100) memsize_mI4
 100  format('MEMTEST: Allocating ',i5,' Mega-I*4 words of memory',/,
     .       'Runstring: memtest <Mega I*4 words>',/,
     .       'NOTE: On a 32-bit machine, you can not go past 512')
****  Check to see if value too large
      if( memsize_mI4.ge.512 ) then
         write(*,120) memsize_mi4
 120     format('Selected memory ',i5,' is larger than allowed ',
     .          'for 32-bit machine; trying anyway')
*        Compute a new word size so that we pass without overflowing
*        integer*4
         sizeword = 1024   ! 1024 integer*4 words to be allocated
         numwords = (memsize_mI4*1024)
      else
*        Small enough value that we pass as integer*4
         sizeword = 1
         numwords = (memsize_mI4*1024)*1024
      end if

*     New version of memory allocation.  Pass the address of the
*     array to be partitioned 
*     Try to allocate the memory
      istart_vma  = memassign(numwords, sizeword, loc(vma_data)) 

      if( istart_vma .eq.0 ) then
         write(*,150) numwords*sizeword
 150     format('Unable to allocate ',i12,' I*4 words of memory')
         stop 'Unable to allocate memory'
      else
         write(*,160) 
 160     format('Memory allocated')
      endif

****  Print out some test data
      write(*,200) 
 200  format('Testing memory allocation: First 10 words ',/,
     .       'Indx  Element    Value ')
      do i = 1, 10
         vma_data(istart_vma+i-1) = i
         write(*,210) i,i,vma_data(istart_vma+i-1)
 210     format(i3,1x,I10,1x,i10)
      enddo

      write(*,250) 
 250  format('Testing memory allocation: Last 10 words ',/,
     .       'Indx  Element    Value ')

      do j = 1,10
         i = j + numwords*sizeword - 10
         vma_data(istart_vma+i-1) = i
         write(*,210) j, i,vma_data(istart_vma+i-1)
      enddo

      end


     
