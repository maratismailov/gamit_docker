CTITLE DUSE
 
      Logical Function Duse( data, n)
 
      implicit none 

 
*     Returns true if data type n (n=1 for delays, n=2 for rates) is
*     being used.  Routine assumes that there are 3 delay data types
*     (set with bits 1-3) and one rate type (bit 4 of data)
*                                 11:10 PM MON., 16 Feb., 1987
*
*   data    - bit mapped data type
*   n       - data type, 1 for delay and 2 for rate
 
      integer*4 data, n, cand
 
****  START, for delays mask the first 3 bits
 
      duse = .false.
      if( n.eq.1 ) then
          if( cand(data,7).ne.0 ) duse = .true.
      end if
 
      if( n.eq.2 ) then
          if( cand(data,8).ne.0 ) duse = .true.
      end if
 
***** Thats all
      return
      end
 
