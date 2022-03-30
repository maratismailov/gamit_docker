C       Test the subroutine get_rxnames

       implicit none

       include '../includes/makex.h'

       character*4 site
       character*80 rxfiles(10)
       
       integer*4 iwknstart,iwknstop

       real*8 sowstart,sowstop

       urinex = 1
       site = 'burg'            
                   
       frinex = ' '
       frinex(1:36) = '/data3c/simon/emed96/rinex/bulgaria/'
       frinex(37:72) = ' /data3c/simon/emed96/rinex/georgia/'
       print *,'In test frinex: ',frinex
       iwknstart = 871
       iwknstop  = 871
       sowstart = 86400.
       sowstop = 172800.  
       
       call get_rxfiles( site,iwknstart,sowstart,iwknstop,sowstop
     .                 , rxfiles )
       
       print *,'Back in test'
       print *,'site ',site
       print *,'rxfiles ',rxfiles

       stop
       end
 
       
