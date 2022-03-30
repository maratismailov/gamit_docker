CTITLE LOGLU
 
 
      integer*4 function loglu(id)
 
      implicit none

*     routine to emulate the loglu function from the HP1000.
*     This function simply returns 6.
 
*   id      - Dummy argument which is not used.
 
      integer*4 id
                   
c     add this dummy statment to avoid a gfortran compiler warning
      id = 0       

      loglu = 6
 
*     Thats all
      return
      end
 
