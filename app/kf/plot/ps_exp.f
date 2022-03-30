      program ps_exp
 
*     This program will take a postscript file generated
*     with the ps.mono input to ctrans and stretch it so that it
*     filles the complete page.
 
*     This routine reads from std input and outputs to standard
*   output.
 
*            line   - Line read from PS file
 
      character*40 line
 
*           last_char   - Last character of line
 
      character*1 last_char
 
*    ierr           - IOSTAT error in read
*   jerr            - IOSTAT error decoding
*   trimlen         - Length of string emulation
*   lenl            - Length of line read
*   ix,iy           - Cooridinates of point
 
 
      integer*4 ierr, jerr, trimlen, lenl, ix,iy
 
****  Start reading until we reach end of the file
 
      ierr = 0
      do while ( ierr.eq.0 )
 
          read(*,'(a)', iostat=ierr ) line
 
*         if no error see if we should mod line.
*                                 ! OK Check if out
          if ( ierr.eq.0 ) then
 
*             Check for scale
*                                              ! Fix y scale to 0.160
              if( line.eq.'.125 .125 s' ) then
                  line = '.135 .175 s'
              end if
 
*             See if last character is m (move) or line (l)
              lenl = trimlen(line)
              if( line(lenl:lenl).eq.'m' .or.
*                                                     ! Try to raed position
     .            line(lenl:lenl).eq.'l'      ) then
                  read(line,*,iostat=jerr) ix, iy
*                                         ! Change y value and write line
                  if( jerr.eq.0 ) then
                      iy = iy - 900 
                      ix = ix - 200
                      last_char = line(lenl:lenl)
                      write(line,100) ix,iy, last_char
  100                 format(1x,i4,1x,i4,1x,a1)
                      lenl = 12
                  end if
              end if
 
*             Now write line back out
              write(*,'(a)') line(1:lenl)
*                     ! No error on read
          end if
*                     ! Loop until EOF
      end do
 
****  Thats all
      end
 
