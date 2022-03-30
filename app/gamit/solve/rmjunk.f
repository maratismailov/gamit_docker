c
      Subroutine RMJUNK( jst,free,isite,mode )
c
c     remove site or satellite parameters if there is no
c     effective observation related to this site or satellite.
c     mode = 1: check site and sat observation for all sessions
c     mode = 2: remove atmospheric parameter if in one session
c               there is no observation for a site
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
c
c     nsite : number of stations this session  (in common /block2/)
c     jst    : number of observations each site this session
c     iseen  :
c     free   : array of parameter flags (=0 off, =1  on)
c     isite  : index for this site

      integer free(maxprm),jst(maxsit),jzen,ncoord,nclock
     .      , isite,mode,jgrad,i2,ilast,i1,i,j,k
                 
       logical debug/.false./
   
       if(debug) then
         print *,'RMJUNK mode nlive free ',mode,nlive, ntpart 
         write(*,'("FREE ",100(i5,1x,i1))')  (i,free(i),i=1,ntpart)
       endif 

      if( mode.eq.1 ) then

c     if there is any site no observation
c
      if ( sitpar ) then
         do 50 i = 1,nsite
            iusest(i) = 1
            if (jst(i).eq.0) iusest(i) = 0
            if (iusest(i).eq.0) then
               do 20 j = 1,3
                  i1 = (i-1)*3+j
                  if (free(i1).eq.0) goto 20
                  free(i1) = 0
                  nlive = nlive-1
                  msig = msig-1
 20            continue
            endif
 50      continue
      endif
                       
c     if there is any sat no observation
c
      if (satpar ) then
c         set number of parameters prior to this in the list
c         - includes coordinates (ncoord), zenith delays, gradients, and site clock parameters
          ncoord = nsite*3
          nzen = nsite*nzen 
c         if multiple zenith delays, add one for the average zenith delay parameter
          if(nzen.gt.1) nzen = nzen + nsite 
          ngrad = nsite*ngrad*2
          nclock = nsite 
          if( sitpar ) then
            ilast = ncoord + nzen + nclock + ngrad
          else
            ilast = nzen
          endif       
c         integrated orbital parameters
          do  i = 1,nsat
            iusesa(i) = 1
            if (iseen(i).eq.0) iusesa(i) = 0
            if(debug) print *,'RMJUNK isat iseen ',i,iseen(i) 
            if (iusesa(i).eq.0) then
               do  j = 1,norb
                  i2 = ilast+(i-1)*norb+j
                  if (free(i2).eq.1 ) then
                    free(i2) = 0    
                    if(debug) print *,'RMJUNK iorb ilast i2 free '
     .                        ,j,ilast,i2,free(i2)
                    nlive = nlive-1  
                    msig = msig-1
                  endif
               enddo
            endif
         enddo   
c        SV antenna offsets 
         ilast = ilast + nsat*norb  
         do  i = 1,nsat     
            if (iusesa(i).eq.0) then
               do  j = 1,3
                  i2 = ilast+(i-1)*3+j
                  if (free(i2).eq.1 ) then
                    free(i2) = 0
                    nlive = nlive-1  
                    msig = msig-1
                  endif
               enddo
            endif
         enddo   
      endif   

c     remove atmospheric zenith delay and gradient parameters - called by LSQUAR for each 
c     site only if no observations for this site this session

      else if (mode.eq.2) then

          ncoord = nsite*3
          if ( .not.sitpar ) ncoord = 0 
          jzen = 0 
          if ( zenpar ) then  
* MOD TAH 190615: Fixed bug and replaced nzen with jzen which is the 
*         quantity computed in the loop.   However since loop starts
*         at zero and increments by 1 each time. jzen is simply isite-1
*         so replaced the loop as well.
c           first the average zenith delay
C           do i = 1, isite-1
C            jzen = jzen + 1 
C           enddo 
C           ilast = ncoord + nzen  
*           Updated code:
            jzen = isite - 1
            ilast = ncoord + jzen

            if( free(ilast+1).eq.1 ) then
              free(ilast+1) = 0
              nlive = nlive - 1    
              msig  = msig  - 1  
            endif
c           now the piece-wise linear adjustments
            if( nzen.gt.1 ) then 
              jzen = 0
              do i = 1, isite-1
                jzen = jzen + nzen
              enddo  
              ilast = ncoord +  nsite + jzen  
              do k = 1,nzen 
                if( free(ilast+k).eq.1 ) then
                  free(ilast+k) = 0
                  nlive = nlive - 1    
                  msig  = msig  - 1  
                endif
              enddo
            endif
          endif

c         calculate the index for the gradient parameters for this site and reduce 
c         the number of live parameters for any that would have been adjusted   
          if ( gradpar ) then 
            jgrad = 0  
c           N/S gradient
            if ( zenpar ) jzen = (nzen+1)*nsite 
            do i = 1, isite-1
              jgrad = jgrad + ngrad
            enddo
            ilast = ncoord + jzen + jgrad 
            do k = 1,ngrad 
              if( free(ilast+k).eq.1 ) then
                free(ilast+k) = 0
                nlive = nlive - 1    
                msig  = msig  - 1  
              endif
            enddo
c           E/W gradient 
c             get # of N/S gradients + E/W gradients for other stations this session
            jgrad = ngrad*nsite
            do i = 1, isite-1
              jgrad = jgrad + ngrad
            enddo
            ilast = ncoord + jzen + jgrad 
            do k = 1,ngrad 
              if( free(ilast+k).eq.1 ) then
                free(ilast+k) = 0
                nlive = nlive - 1    
                msig  = msig  - 1  
              endif
            enddo
          endif 

      endif


         if(debug) then
           print *,'RMJUNK end nlive free',nlive 
           write(6,'((100i3,/))')  (free(i),i=1,1000)
          endif 


      return
      end


