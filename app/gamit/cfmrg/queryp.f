Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      Subroutine queryp( nsite,ntpart
     .                 , igcchk,iatchk,icrchk,iorchk,iantchk,ieochk
     .                 , nepch,zenmod,numzen,idtzen
     .                 , gradmod,numgrad,idtgrad )
c
c     Query the user for coordinate, clock, atmosphere, and orbital partials
c     Values returned as '1' if yes, '2' if no
c
      implicit none
c
      include '../includes/dimpar.h'
c
      integer*4 ntpart,nsite,ioerr,nepch
     .      , igcchk,iatchk,icrchk,iorchk,ieochk,iantchk
     .      , nzen(maxsit),numzen,idtzen(maxatm),numwnd
     .      , ngrad(maxsit),numgrad,idtgrad(maxgrad),i

       character*1 ans,lowerc
       character*3 zenmod,gradmod
       character*256 message
c
c.....coordinate partials
       read(5,'(a)') ans
       if( lowerc(ans).eq.'y') then
         icrchk = 1
       elseif( lowerc(ans).eq.'n') then
         icrchk = 2
       else
         call report_stat('FATAL','CFMRG','queryp',ans,
     .   'Invalid input for coord partial availability: ',0)
       endif

c.....station clock partials (currently automatically used)
c
      igcchk=1
c
c.....atmospheric zenith angle parameters
c        
c      There are always atmospheric partials on the C-file (this should be hardwired)
       read(5,'(a)') ans
       if( lowerc(ans).eq.'y') then
         iatchk = 1
       elseif( lowerc(ans).eq.'n') then
         iatchk = 2
       else
         call report_stat('FATAL','CFMRG','queryp',ans,
     .   'Error, deciding if atm partial available: ',0)
       endif
c     if yes, get number of zenith delay parameters per station and set
c     the model and tabular epochs
      do i=1,nsite
        nzen(i) = 0
      enddo
      if(iatchk.eq.1) then
        read(5,*) (nzen(i),i=1,nsite)
      endif
      do i = 2,nsite
         if( nzen(i).ne.nzen(i-1) ) then
           write(message,'(a,i3,a,i3,a)') 'Numzen(',i,') ne numzen('
     .            ,i-1,')--must be same for all sites'
           call report_stat('FATAL','CFMRG','queryp',' ',message,0)
          endif
       enddo
       numzen = nzen(1)
       if( numzen.gt.maxatm )then
          write(message,'(a,i3,a,i3,a)') 'Number of zenith parameters ('
     .         ,numzen,') exceeds maxatm (',maxatm,')'
           call report_stat('FATAL','CFMRG','queryp',' ',message,0)
       endif         
       if( numzen*nsite.gt.2500 ) then
          write(message,'(a,i3,a,i3,a)') 'Number of zenith parameters ('
     .         ,numzen,') times number of sites (',nsite
     .     ,') cannot exceed 2500'
           call report_stat('FATAL','CFMRG','queryp',' ',message,0)
       endif
c      eventually allow different numbers for each site
       if( numzen.gt.1 ) then
          zenmod = 'PWL'
       else
          zenmod = 'CON'
       endif
       call uppers(zenmod)
       if( zenmod.eq.'CON' ) then
          numwnd = ( nepch+numzen-1)  / numzen
       elseif ( zenmod.eq.'PWL' ) then
         if( numzen.gt.1 ) then
            numwnd = ( nepch+numzen-2) / (numzen-1)
         else
            numwnd = 0.
         endif
       else
         call report_stat('FATAL','CFMRG','queryp',zenmod,
     .   'Error, unknown zenith delay model, zenmod: ',0)
         return
       endif
       idtzen(1) = 1
       if( numzen.gt.1 ) then
          do i = 2,numzen
            idtzen(i) = idtzen(i-1) + numwnd
          enddo
       endif


c.....orbital partials

       read(5,'(a)') ans   
       if( lowerc(ans).eq.'y') then
         iorchk = 1
         if( ntpart.gt.15 ) then
           ieochk = 1
         else
           ieochk = 2
         endif
       elseif( lowerc(ans).eq.'n') then
         iorchk = 2
         ieochk = 2
       else
         call report_stat('FATAL','CFMRG','queryp',ans,
     .   'Error, deciding if orbital partials available: ',0)
       endif

      if( iorchk.eq.1 .and. ntpart.lt.15 ) then 
         write(message,'(a,i3,a)')   
     .   'Orbit partials requested but ntpart (=',ntpart,') < 15'
         call report_stat('WARNING','CFMRG','queryp',' ',message,0)
         iorchk=2
      endif

c.....satellite antenna offset partials

       read(5,'(a)',iostat=ioerr) ans
       if(ioerr.eq.0 ) then
         if( lowerc(ans).eq.'y') then
           iantchk = 1 
         elseif ( lowerc(ans).eq.'n') then
           iantchk = 2                    
         else
           call report_stat('FATAL','CFMRG','queryp',ans,
     .    'Error, deciding if SV antenna partials available: ',0)
         endif
       else 
         call report_stat('WARNING','CFMRG','queryp',ans,
     .   'Unexpected end in batch file, assuming no SV antenna partials'
     .   ,0)
         iantchk = 2
       endif
c      possible add some more checks on npart and norb


c......atmospheric gradient parameters

       read(5,'(a)',iostat=ioerr) ans 
       if(ioerr.eq.0 ) then  
         if( lowerc(ans).eq.'y') then
c          if yes, get number of gradient parameters per station and set
c          the model and tabular epochs
           do i=1,nsite
             ngrad(i) = 0
           enddo
           read(5,*) (ngrad(i),i=1,nsite)
           do i = 2,nsite
             if( ngrad(i).ne.ngrad(i-1) ) then
                write(message,'(a,i3,a,i3,a)') 'Ngrad(',i,') ne ngrad('
     .                     ,i-1,')--must be same for all sites'
                call report_stat('FATAL','CFMRG','queryp',' ',message,0)   
c               eventually allow different numbers for each site
             endif
           enddo
           numgrad = ngrad(1) 
           if( numgrad.gt.1 ) then 
             if( numgrad.gt.maxgrad/2 )then
                write(message,'(a,i3,a,i3,a)') 
     .            'Number of gradient parameters * 2 ('
     .             ,2*numgrad,') exceeds maxgrad (',maxgrad,')'
                call report_stat('FATAL','CFMRG','queryp',' ',message,0)
             endif
             gradmod = 'PWL'  
             numwnd = ( nepch+numgrad-2) / (numgrad-1)  
             idtgrad(1) = 1
             do i = 2,numgrad
               idtgrad(i) = idtgrad(i-1) + numwnd
             enddo
           else
              gradmod = 'CON'
           endif   
         elseif( lowerc(ans).eq.'n') then
           numgrad = 1
           gradmod = 'CON'
         else
           call report_stat('FATAL','CFMRG','queryp',ans,
     .'Error deciding if # of multiple gradient parameters to be read'
     .       ,0)
         endif 
       else 
         call report_stat('WARNING','CFMRG','queryp',ans,
     .'Unexpected end in batch file, assuming 1 gradient parameter'
     .   ,0)
         numgrad = 1
         gradmod = 'CON'
       endif

      return
      end
