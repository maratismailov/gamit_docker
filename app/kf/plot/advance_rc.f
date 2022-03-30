
CTITLE ADVANCE_RC
 
      subroutine advance_rc( unit, irow, icol, page, file, cdate, 
     .                       mxr,  mxc )
 
*     This routine will advance the row and columns numbers and
*     check to see if we should output a new header
 
* PASSED VARIABLES
 
*   unit        - Unit number for output
*   irow, icol  - Current row and columns numbers
*   page        - Current page number
*   mxr, mxc    - number of rows and colmns
 
      integer*4 unit, irow, icol, page, mxr, mxc
 
*   file        - name of the bak file
*   cdate       - Current time encoded as character string
 
      character*(*) file, cdate
 
* LOCAL VARIABLE
 
*   lent        - Length of string
*   trimlen     - Gets length of string
*   ir          - Row counted from top of page instead of bottom
 
      integer*4 lent, trimlen, ir
 
*   view(4)     - View needed for this row and col position.
 
      real*4 view(4)
 
****  Increment row, and see if we need to go to next column
      irow = irow + 1
      if( irow.gt.mxr ) then
          irow = 1
          icol = icol + 1
          if( icol.gt.mxc ) then
              icol = 1
              page = page + 1
 
*             Now write header
              if( page.ne.1 ) then
*                                             ! Pause and Cause a new page
                  write( unit,* ) ' Id '
                  write( unit,* ) ' Erase '
              end if
*                                             ! Do preamble
              write( unit, 100)
 100          format(' View  0.05 0.95 0.95 0.99',/,
     .               ' Scale 0 1 0 1 ',/,
     .               ' Font 2',/,
     .               ' char 2 2 ')
              lent = max(1,trimlen(file))
              write(unit, 120 ) file(1:lent)
 120          format(' label 0.01 0.5 1 0 "',a,'"')
 
              write(unit, 140 ) page
 140          format(' label 0.5 0.5 1 0 "',i4,'"')
              lent = max(1,trimlen(cdate))
              write(unit, 160 ) cdate(1:lent)
 160          format(' label 0.75 0.5 1 0 "',a,'"')
 
****          Now reset the text quantities
              write(unit,200) max(0.9, 3.0/mxc)
 200          format(' font 1',/,
     .               ' char ', f6.2)
*                         ! new page
          end if
*                         ! Incrementing row and column
      end if
 
****  Now compute the view for this plot
      ir = mxr+1 -irow
      view(1) = (icol-1)*0.99/mxc + 0.100
      view(2) = icol*0.99/mxc 
 
      view(3) = (ir-1)*0.95/mxr + 0.05
      view(4) = ir*0.95/mxr
 
      write(unit, 300) view
 300  format(' View ',4f8.4)
 
****  Thats all
      return
      end
 
