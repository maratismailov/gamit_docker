      program ps_rot
 
*     This program will take a postscript file generated
*     with the ps.mono input to ctrans and rotate it and rescale
*     by a user specified amount (given in runstring as % increase)
 
*     This routine reads from std input and outputs to standard
*     output.
 
*            line   - Line read from PS file
 
 
      character*40 line
 
*           last_char   - Last character of line
 
 
      character*1 last_char
 
*    ierr           - IOSTAT error in read
*   jerr            - IOSTAT error decoding
*   trimlen         - Length of string emulation
*   lenl            - Length of line read
*   ix,iy           - Cooridinates of point
*   ixo,iyo         - Original coordinates of point
*   xoff, yoff      - Origin offsets to the applied in x and y.
*                     Should keep plot centered on the paper.
*   rcpar           - Reads runstring
 
      integer*4 ierr, jerr, trimlen, lenl, ix,iy, ixo, iyo, rcpar,
     .          xoff, yoff
 
*    dscale  - Scale change in percent
*    scale   - Scale to be passed to Postscript
 
 
      real*4 dscale, scale
 
****  See if scale passed in runstring
 
      lenl = rcpar(1,line)
      if( lenl.gt.0 ) then
          read(line,*, iostat=ierr ) dscale
          scale = (1.d0 + dscale/100)*0.125
      else
          scale = 0.125
      end if

***** Set the offsets
      xoff = 2600 * 0.125 / scale     
      yoff = 3265 * 0.125 / scale     
 
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
                  write(line,100) scale, scale
  100             format(f5.3,1x,f5.3,' s')
              end if
 
*             See if last character is m (move) or line (l)
              lenl = trimlen(line)
              if( line(lenl:lenl).eq.'m' .or.
*                                                     ! Try to read position
     .            line(lenl:lenl).eq.'l'      ) then
                  read(line,*,iostat=jerr) ixo, iyo
*                                         ! Change y value and write line
                  if( jerr.eq.0 ) then
                      iy =   (ixo-2600) + yoff 
                      ix =  -(iyo-3265) + xoff
                      last_char = line(lenl:lenl)
                      write(line,200) ix,iy, last_char
  200                 format(1x,i4,1x,i4,1x,a1)
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
 
