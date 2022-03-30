      Subroutine WRITE_SOLN( constraints,free_fix,phase_obs )

c     do the error analysis and write out the solution to the q-, o-, and/or h-file
c        King 930921

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      real*8 covkep(maxcov)

      character*4 free_fix, phase_obs
      character*5 constraints
      character*62 divider

      logical debug/.false./

      data divider
     ./'--------------------------------------------------------------'/
                    
c       print the header for the main solution

      if( free_fix.eq.'free' ) then
        if( logprt ) write(6,10) 
        write(10,10) 
   10   format(/,'----------------------------------------------------'
     .        ,/,'    **** Summary of biases-free solution ****'    ,/
     .          ,'----------------------------------------------------')
       elseif( free_fix.eq.'fixd' ) then
        if( logprt ) write(6,11) 
        write(10,11)  
   11   format(/,'----------------------------------------------------'
     .        ,/,'    **** Summary of biases-fixed solution ****'    ,/
     .          ,'----------------------------------------------------')
       else  
        if( logprt ) write(6,12) 
        write(10,12)
   12   format(//,'*******problem printing solution header',/)  
       endif

                                                               
c       do the error analysis
 
       if(debug) print *,'WRITE_SOLN sclerr ',sclerr           
                        
       call lsqerr( covkep,constraints,free_fix,phase_obs )
  
c       if the inversion worked write out the solution
       call lsqdo1( covkep,free_fix )  

c       indicate end-of-solution in q-file
       write(10,'(//,a,a5,a,a4,a,a4,/,a)')  'End of ',constraints
     .   ,' solution with ',phase_obs
     .   ,' observable and ambiguities ',free_fix,divider
                              
       return
       end

