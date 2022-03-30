CTITLE PROPER_GLB_RUN
 
      subroutine proper_glb_run ( abort )

      implicit none 
 
 
*     Reports the proper runstring for globk.  If abort is 1 then
*     the program will stop.
 
*   abort       - +1 to abort the program
*   id          - Dummy argument
*   iout        - Output device
*   loglu       - User LU number function
 
      integer*4 abort, id, iout, loglu
 
*     Get user LU
      call proper_runstring('globk.hlp','globk',0) 
      iout = loglu( id )
C     write(iout, 100)
C 100 format(//' Incorrect runstring for GLOBK:',/,
C    .         ' CI> GLOBK,<crt>,<prt>,<log>,list,markov,'
C    .         ' <commom>,<sort>,<direct>',/,
C    .         ' <> parameters are optional.  See GLOBK::HELP for',
C    .         ' details'/)
 
      if( abort.eq.1 ) then
          stop ' GLOBK terminating: Incorrect runstring'
      else
          write(iout,200)
  200     format(/' Minor error in runstring --',
     .            ' Processing continuing'/)
      end if
 
***** Thats all
      return
      end
 
