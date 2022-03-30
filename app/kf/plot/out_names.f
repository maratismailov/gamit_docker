CTITLE OUT_NAMES
 
      subroutine out_names(iout, type, names, num)
c
c
c     Subroutine to out put a list of character strings
c
c Variables
c ---------
c iout -- the output lu number
c type -- the type of quantites being output
c names -- the list of names to be output
c num   -- the number of names of be output
c
      integer*4 iout, num
 
c
      character*(*) type, names(1)
 
c
c Local variables and functions
c -----------------------------
c i -- loop counter
c trimlen -- HP utility for length of string
c max_length -- maximum length of string
c
      integer*4 i, trimlen, max_length
 
c
c.... Output the description line
      write(iout,100)  num, type
  100 format(" There are ",i4,1x,a,1x," available:")
c
c.... Find the longest string
      max_length = 0
      do i = 1, num
          max_length = max( max_length,trimlen(names(i)) )
      end do
c
c.... Now output the list
      if( max_length.le. 10 ) then
          write(iout,200) (names(i)(1:max_length),i=1,num)
  200     format(5(1x,a,3x))
*                            ! one per line
      else
          write(iout,'(1x,a)') (trim(names(i)),i=1,num)
      end if
c
      return
      end
 
