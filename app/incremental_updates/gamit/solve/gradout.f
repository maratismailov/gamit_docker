Copyright 1998 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      Subroutine GRADOUT( free_fix,id  )

c     Write a summary to the O-file of the gradient-delay estimates,
c     for convenient grep-ing, scanning, and plotting.
c     R. King Sept 1998 from sb zenout (R. King and T. Herring Oct 1993)

c    **Note: This subroutine works only if all sites have the same
c            number of gradient-delay parameters
      
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      logical gradslot

      real*8 sod,delt,seconds,tstart,m2cyc,at10deg

      integer*4 igstart,igstop,idate(maxatm,5),iyr,idoy
     .        , id(maxsit),ngtot,jdstart,jd,indx,ilive,is,i,j,k

      character* 1 bfixcd,upperc
      character* 3 lowerc
      character* 4 free_fix,sitcod(maxsit)
      character* 7 labelns,labelew 
      character* 8 gradkey
                                                                      
c     Function
      integer*4 julday 

c     compute constants to convert values 
      m2cyc = 1/(299792458.d0/1.57542d9)
      at10deg = 1/(sin(10.d0/180.d0*pi)*tan(10.d0/180.d0*pi)+0.003)

c     determine the parameter slots for the gradient delays

      gradslot = .false.
      igstart = 0
      igstop =  0
c     check if gradient delay
      do i=1,ntpart
         if( islot1(i).ge.24001.and.islot1(i).le.29000 ) then
           if(.not.gradslot) igstart = i
           gradslot = .true.
         else
           if(gradslot) then
             igstop = i-1
             goto 5
           endif
         endif
      enddo

c     if no multiple gradient parameters, skip the output

      if( igstart.eq.0 .and.igstop.eq.0 ) return

c     compute the total number and number/site of each type (NS/EW) of gradient delays

    5 ngtot = igstop - (igstart-1)
      ngrad  = ngtot/nsite/2
      if( mod(ngtot,nsite).ne.0 )
     .   call report_stat('WARNING','SOLVE','gradout',' '
     .       , '# gradient parameters not same for all sites    ',0)

c     set the keywords 
        gradkey = 'ATM_GRAD'
        labelns = 'NS_GRAD'
        labelew = 'EW_GRAD'

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
        sitcod(i) = cfiln(is)(2:5)
        call uppers(sitcod(i))
      enddo

c    get the date (y m d h m)  --works now only for single sessions

       jdstart = julday( it0(1),it0(2),it0(3) )
       tstart  = 3600.d0*t00(1) + 60.d0*t00(2) + t00(3)
c      round to nearest minute
       tstart = 60.d0*dnint(tstart/60.d0)

       do i = 1,ngrad
c         get time (sec) from start of tabular epoch for piecewise-linear model
          delt= (idtgrad(i) - idtgrad(1))*inter
c         for constant model, time is mean between two epoch
          if( lowerc(gradmod).eq.'CON' .and.ngrad.gt.1 ) then
c            assume evenly spaced epochs
             delt = delt + (idtgrad(2)-idtgrad(1))*inter/2.
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

c     write the estimates to the o-file

      write(15,10) gradkey, bfixcd,qfiln
   10 format(/,a8,1x,a,'  q-file: ',a16)
      indx = igstart     
      ilive = 0
      do i=1,igstart
        if(free(i).gt.0) ilive = ilive + 1   
      enddo     
c     write the N/S gradients
      do i = 1,nsite
         do j=1,ngrad
           if(indx.gt.maxprm) call report_stat('FATAL','SOLVE','gradout'
     .          ,' ','indx > maxprm ',0)   
           if( free(indx).gt.0 ) then
             write(15,20) labelns,bfixcd,sitcod(i),i,(idate(j,k),k=1,5)
     .                  , adjust(indx) * at10deg/m2cyc
     .                  , sigma(ilive) * at10deg/m2cyc
     .                  , postvl(indx) * at10deg/m2cyc
             ilive = ilive + 1
           endif
           indx = indx + 1 
         enddo
       enddo  
c     write the E/W gradients
      do i = 1,nsite
         do j=1,ngrad
           if(indx.gt.maxprm) call report_stat('FATAL','SOLVE','gradout'
     .          ,' ','indx > maxprm ',0)    
           if( free(indx).gt.0 ) then
             write(15,20) labelew,bfixcd,sitcod(i),i,(idate(j,k),k=1,5)  
     .                  , adjust(indx) * at10deg/m2cyc
     .                  , sigma(ilive) * at10deg/m2cyc
     .                  , postvl(indx) * at10deg/m2cyc
             ilive = ilive + 1
           endif
           indx = indx + 1 
         enddo
      enddo  
   20 format(a7,1x,a,1x,a4,1x,i2,2x,i4,4i3,f8.4,' +- ',f8.4,2x,f8.4)  

      return
      end
