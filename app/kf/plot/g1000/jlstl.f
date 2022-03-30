CTITLE JLSTL
 
      subroutine jlstl( line_type )
 
*     This routine uses the dashchar package if line style is greater
*     than 1.  Other wise the standard routine routines are used.  We
*     allow 8 line styles, and the thiknes of the line is determined
*     by the multiples of 8.
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   line_type   - Upto 8 styles, with values greater than eight
*               - genertating thicher lines.
 
 
      integer*4 line_type
 
* LOCAL VARIABLES
 
 
 
      integer*4 ipat(10)
 
      data ipat / o'177777', o'125252', o'146314', o'161616',
     .            o'177400', o'170350', o'170617', o'177770' ,
     .            o'101010',  o'000100'  /
 
***** First see if we are simple line
 
*                                 ! Value between 0 and 8
      jdash = mod(line_type,10)
      jwidth = int((line_type-1)/10) + 1
 
*     Set the dash pattern for jdash greater than 1
      if( jdash.gt.1 ) then
          call dashdb(ipat(jdash))
      end if
 
****  Now set the line width
      if( jwidth.gt.0 ) then 
          call setusv('LW', jwidth*500)
*         Mark size is now set in JCSIZ - character size rouitne
*         call setusv('MS', jwidth*1000)
      end if
 
***** Thats all
      return
      end
 
