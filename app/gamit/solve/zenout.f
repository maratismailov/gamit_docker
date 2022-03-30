Copyright 1993 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      Subroutine ZENOUT( free_fix,id  )

c     Write a summary to the O-file of the zenith-delay estimates,
c     for convenient grep-ing, scanning, and plotting.
c     R. King and T. Herring   Oct 1993

c    **Note: This subroutine works only if all sites have the same
c            number of zenith-delay parameters
      
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'                          
      include 'parameters.h' 

      logical zenslot

      real*8 sod,delt,seconds,tstart,adj,sig,post,cov,rho

      integer*4 izstart,izstop
     .       , avg_indx(maxsit),avg_ilive(maxsit),indxa,ilivea,jel,icov
     .       , idate(maxatm,5),iyr,idoy,id(maxsit),nztot,jdstart
     .       , jd,indx,ilive,is,i,j,k

      character* 1 bfixcd,upperc
      character* 3 lowerc
      character* 4 free_fix,sitcod(maxsit)
      character* 7 zenkey
                                     
c     Function
      integer*4 julday 

c     determine the parameter slots for the zenith delays

      zenslot = .false.
      izstart = 0
      izstop =  0
c     check if zenith delay
      do i=1,ntpart
         if( islot1(i).ge.21501.and.islot1(i).le.24000 ) then
           if(.not.zenslot) izstart = i
           zenslot = .true.
         else
           if(zenslot) then
             izstop = i-1
             goto 5
           endif
         endif
      enddo

c     if no multiple zenith parameters, skip the output

      if( izstart.eq.0 .and.izstop.eq.0 ) return

c     compute the total number and number/site of zenith delays

    5 nztot = izstop - (izstart-1)
      nzen  = nztot/nsite
      if( mod(nztot,nsite).ne.0 )
     .   call report_stat('WARNING','SOLVE','zenout',' '
     .       , 'Zenith parameters not same for all sites    ',0)

c     set the keyword
        zenkey = 'ATM_ZEN'

c     set the bias-fixing code letter
        if ( free_fix.eq.'fixd' ) then
            write(bfixcd,'(a1)') upperc('x')
        else
            write(bfixcd,'(a1)') upperc('r')
        endif

c   get the 4-character site codes from c-file name
      do i=1,nsite
        is = id(i)
        if (is.lt.1) is=1
        sitcod(i) = obfiln(is)(2:5)
        call uppers(sitcod(i))
      enddo

c    get the date (y m d h m)  --works now only for single sessions

       jdstart = julday( it0(1),it0(2),it0(3) )
       tstart  = 3600.d0*t00(1) + 60.d0*t00(2) + t00(3)
c      round to nearest minute
       tstart = 60.d0*dnint(tstart/60.d0)

       do i = 1,nzen
c         get time (sec) from start of tabular epoch for piecewise-linear model
          delt= (idtzen(i) - idtzen(1))*inter
c         for constant model, time is mean between two epoch
          if( lowerc(zenmod).eq.'CON' .and.nzen.gt.1 ) then
c            assume evenly spaced epochs
             delt = delt + (idtzen(2)-idtzen(1))*inter/2.
          endif
          jd = jdstart
          sod  = tstart
          call timinc(jd,sod,delt)
c         round to nearest minute
          sod = 60.d0*dnint(sod/60.d0)
          call dayjul(jd,iyr,idoy)
          idate(i,1) = iyr
          call monday(idoy,idate(i,2),idate(i,3),iyr)
          call ds2hms(iyr,idoy,sod,idate(i,4),idate(i,5),seconds)
       enddo
    
                                                  
c     find the slots for the average values to combine with the PWL estimates
               
      indxa = 3*nsite                                   
      ilivea = 0
      do i=1,indxa 
        if( free(i).gt.0 ) ilivea = ilivea + 1
      enddo 
      do i=1,maxsit
        avg_indx(i) = 0
        avg_ilive(i) = 0
      enddo
      do i=1,nsite       
        indxa = indxa + 1
        avg_indx(i) = indxa  
        if( free(indxa).gt.0 ) then
          ilivea = ilivea + 1
          avg_ilive(i) = ilivea 
        endif
      enddo
c      print *,'avg_index ',(avg_indx(i),i=1,10)
c      print *,'avg_ilive ',(avg_ilive(i),i=1,10)     
                                                                          

c     write the estimates to the o-file
                    
      write(15,'(/,a7,1x,a1,a,a16,a)') zenkey,bfixcd,'  q-file: ',qfiln
     .      ,'  Combination of avg and PWL values ' 
      indx = izstart - 1
      ilive = 0
      do i=1,indx
        if( free(i).gt.0 ) ilive = ilive + 1
      enddo     
c      print *,'izstart ilive ',izstart,ilive
      do i = 1,nsite  
       indxa = avg_indx(i)
       do j=1,nzen    
        indx = indx + 1
        if(indx.gt.maxprm) call report_stat('FATAL','SOLVE','zenout'
     .       ,' ','indx > maxprm ',0)  
         if( free(indx).gt.0 ) then  
          ilive = ilive + 1
          ilivea = avg_ilive(i)
          adj = adjust(indxa) + adjust(indx)  
          icov = jel(ilivea,ilive)                        
          cov = a(icov)
c          print *,'ilivea,ilive,icov,cov ',ilivea,ilive,icov,cov 
          rho = cov/sigma(ilivea)/sigma(ilive)
c          print *,'indxa indx rho ',indxa,indx,rho
          sig = sqrt(sigma(ilivea)**2 + sigma(ilive)**2 + cov)
          post = postvl(indxa) + adjust(indx)
c          print *,'postvla post '
c     .         ,postvl(indxa),post
          write(15,20) zenkey,bfixcd,sitcod(i),i,(idate(j,k),k=1,5)
c**     .               , adjust(indx),sigma(ilive),postvl(indx)  
     .                  ,adj,sig,post
   20     format(a7,1x,a,1x,a4,1x,i2,1x,i4,4i3,f8.4,' +- ',f8.4,2x,f8.4)
        endif
       enddo
      enddo

      return
      end
