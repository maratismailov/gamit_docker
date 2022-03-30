      Subroutine eclout
c     Print out a time-ordered summary of eclipses at the end of each integration
c     Print out a PRN/TIME ordered yaw file at the end of each integration
c     R. King 14 Feb 95
c     S. McClusky 01 April 95; R King 28 Feb 97

      implicit none

      include '../includes/dimpar.h'   
      include '../includes/global.h'
      include '../includes/arc.h'

      character*1 eclipse,noeclipse
      integer*4 im1,id1,iyr1,idoy1,ihr1,imin1,im2,id2,iyr2,idoy2,ihr2
     .         ,imin2,yr,mo,day,doy,hr,min,indxs(100),indxf(100),js
     .         ,yaw_entries_old,insn,frqchn,svnstart(5),svnstop(5),i,j,k

      real*8 pjd1,pjd2,stime,ftime,etime,sec
     .      ,e_start(100),e_finish(100)
                       
      character*1 ast1,ast2,e_type(100)
      character*256 message    
      
      logical delecl(100)

      eclipse = 'E'
      noeclipse = ' '
      k = 0
      yaw_entries_old = yaw_entries
      stime = jdb + (tb/86400.d0)
      ftime = jdf + (tf/86400.d0)
      etime = jde + (te/86400.d0)

c     If satellite is eclipsing write out eclipse summary and yaw file entries
      do i=1,neclipse
         if( neclipse.gt.100) then 
           write(message,'(a,i4,a)' ) 'Too many eclipses (',neclipse,')'
           call report_stat('FATAL','ARC','eclout',' ',message,0)
         endif

c        If end time still zero, either the eclipse ended outside of integration 
c        span or the eclipse was within the backward-forward integration startup.
         if( eclipse_end(i).eq.0.d0 ) then               
           if( dabs(eclipse_start(i)-etime).lt.20.d0*diint ) then  
c            integration start case (startup is 20 steps), set end = start 
             eclipse_end(i) = eclipse_start(i)
c             print *,'in startup i eclipse_end ',i,eclipse_end(i)
           else
c            eclipse started or ended outside of span, set end to span limits
             if( eclipse_start(i).lt.etime ) then
               eclipse_end(i) = stime
             else
               eclipse_end(i) = ftime
             endif  
c             print *,'outside range i eclipse_end ',i, eclipse_end(i)
           endif
         endif

c        Reorder start and end time to account for backward and forward integraton.
         if( eclipse_start(i).gt.eclipse_end(i) ) then
           pjd1 = eclipse_end(i)
           pjd2 = eclipse_start(i)
         else
           pjd1 = eclipse_start(i)
           pjd2 = eclipse_end(i)
         endif   
c
c        Save all begin and end eclipse times to be sorted and output later.
         e_start(i) = pjd1
         e_finish(i) = pjd2    
         e_type(i) = eclipse_type(i)

      enddo   
c      print *,'After reorder neclipse ',neclipse
c      do i=1, neclipse
c       write(*,'(i3,2f20.12)') i,e_start(i),e_finish(i)
c      enddo


c     If satellite eclipsing then sort start and stop times, check eclipse times make sense.
      if ( neclipse .gt. 0 ) then
c
        write(iarh,'(1x)')
c
c       Sort eclipse start and finish times
        call indexx(neclipse,e_start,indxs)
        call indexx(neclipse,e_finish,indxf)   
c        print *,'neclipse, e_start indxs ',neclipse,e_start,indxs
c        print *,'neclipse  e_finish indxf ',neclipse,e_finish,indxf 
                                  
           
c*** rwk 050223: I don't think this works correctly, and is probably unncessary, comment for now
c       If integration started in shadow, some eclipse start and stop times will be
c       incorrect. Remove incorrect start and stop times.
c**        js = 0     
c**        do i = 1,neclipse
c**          if( i .eq. 1 ) then
c**            js = js+1
c**            e_start(js) = e_s(indxs(i)) 
c**          elseif( i .le. neclipse .and.
c**     .      e_s(indxs(i)) .gt. e_start(js)+0.15d0 ) then
c**            js = js+1
c**            e_start(js) = e_s(indxs(i)) 
c**          endif
c**        enddo    
c**        print *,'js e_start ',js,e_start
c**        jf = 0
c**        do i = neclipse,1,-1
c**          if( i .eq. neclipse ) then
c**            jf = jf+1
c**            e_finish(jf) = e_f(indxf(i))  
c**         elseif( i .ge. 1 .and.
c**     .      e_f(indxf(i)) .lt. e_finish(jf)-0.15d0) then
c**            jf = jf+1
c**            e_finish(jf) = e_f(indxf(i))
c**          endif
c**        enddo
c**        print *,' jf e_finish ',jf,e_finish 
c**
c**      print *,'After checking js jf ',js,jf
c**      do i=1, 7
c**       write(*,'(i3,2f20.10)') i,e_start(i),e_finish(i)
c**      enddo

 
c** rwk 050223: Alternative attempt at getting rid of start-up eclipses

c     First reorder by start time
        do i=1,neclipse
          do j=1,neclipse
            if(indxs(j).eq.i) then
               eclipse_start(i) = e_start(j)
               eclipse_end(i)  = e_finish(j)
               eclipse_type(i) = e_type(j)
            endif
          enddo
        enddo     

c        print *,'After reorder by start time neclipse',neclipse
c        do i=1,7
c         write(*,'(i3,2f20.10)') i,eclipse_start(i),eclipse_end(i)
c        enddo


c     Now remove any short eclipse segments that are within longer ones (slop = 2 min)
        do i=1,neclipse    
          do j=1,neclipse   
            delecl(j) = .false.
            if( (eclipse_start(j)+.0014d0).ge.eclipse_start(i) .and.
     .          (eclipse_end(j)-.0014d0).lt.eclipse_start(i) ) then
               delecl(j) = .true.
            endif
          enddo
        enddo

c        print *,'After check for removal '
c        do i=1, 7
c         write(*,'(i3,2f20.10,1x,l1)') 
c     .       i,eclipse_start(i),eclipse_end(i),delecl(i)
c        enddo

c     Reduce the number to only the valid ones
        js=0
        do i=1,neclipse
          if( .not.delecl(i) ) then     
            js = js + 1
            e_start(js) = eclipse_start(i)
            e_finish(js) = eclipse_end(i)
            e_type(js) = eclipse_type(i)
          endif
        enddo
         
c
c     Loop over eclipse start times outputing to archive and yaw files.  
                   
        do i = 1,js

c     Only write out the ones contained within the arc requested   
c          print *,'i,e_start,ftime stime ',i,e_start(i),ftime,stime
          if((e_start(i).le.ftime.and.e_start(i).ge.stime))then
             k=k+1   

c          If the SV does not start in eclipse, write the start time to the yaw file
            if ( stime .ne. e_start(i) .and. i .eq. 1 ) then  
              call dayjul(jdb,yr,doy)
              call monday (doy,mo,day,yr)
              call ds2hms(yr,doy,tb,hr,min,sec)
c*rwk no longer needed since yawrate and bias stored in common 
c     skipping here avoids a prn change over a day boundary
c              call svnav_read( -1,yr,doy,hr,min,gnss,iprn,insn
c     .                       , frqchn,antbody,sbmass,bias,yawrate
c     .                       , svnstart,svnstop )
              write(iyawtmp,100)iprn,yawrate,bias,yr,mo,day,hr,min
     .             ,noeclipse,eclipse_beta(i)
100           format(1x,i2,1x,f7.3,2x,a1,2x,i4,2(1x,i2),2x,i2,1x,i2  
c                  previously did not write yaw or beta angles
c    .              ,2x,a1)   
c                  now write beta but not yet yaw
     .              ,2x,a1,10x,f7.1)
c                   format including yaw angle (orbits/read_yaw.f)
c    .              ,2x,a1,3x,2f7.1)
              yaw_entries = yaw_entries + 1
            endif

c         Check to see if noon turn occurs between start of arc and first Earth eclipse.
c         If yes print noon turn line
            if ( k .eq. 1 .and. e_start(i)-0.25d0 .gt. stime .and.
     .           e_type(i).eq.'E' ) then 
              call pjdhms(e_start(i)-0.25d0,iyr1,idoy1,im1,id1,ihr1
     .                   ,imin1)
c              call svnav_read( -1,iyr1,idoy1,ihr1,imin1,gnss,iprn,insn
c     .                       , frqchn,antbody,sbmass,bias,yawrate
c     .                       , svnstart,svnstop )
              write(iyawtmp,100) iprn,yawrate,bias,iyr1,im1,id1,ihr1
     .                         , imin1,noeclipse,eclipse_beta(i)
              yaw_entries = yaw_entries + 1
            endif

c         Write eclipse line to yaw file
            call pjdhms( e_start(i),iyr1,idoy1,im1,id1,ihr1,imin1 )
c            call svnav_read( -1,iyr1,idoy1,ihr1,imin1,gnss,iprn,insn
c     .                     , frqchn,antbody,sbmass,bias,yawrate 
c     .                     , svnstart,svnstop )
            write(iyawtmp,100) iprn,yawrate,bias,iyr1,im1,id1,ihr1
     .                       , imin1,e_type(i),eclipse_beta(i)
            yaw_entries = yaw_entries + 1

c         Write noon turn line to yaw file 
            if( e_type(i).eq.'E' ) then
              call pjdhms(e_start(i)+0.25d0,iyr1,idoy1,im1,id1,ihr1
     .                   , imin1)
c              call svnav_read( -1,iyr1,idoy1,ihr1,imin1,gnss,iprn,insn
c     .                       , frqchn,antbody,sbmass,bias,yawrate 
c     .                       , svnstart,svnstop ) 
              write(iyawtmp,100)iprn,yawrate,bias,iyr1,im1,id1,ihr1
     .                         , imin1,noeclipse,eclipse_beta(i)
              yaw_entries = yaw_entries + 1  
            endif
c
c         Write a summary line to the archive file    
            call pjdhms( e_start(i),iyr1,idoy1,im1,id1,ihr1,imin1 )
            call pjdhms( e_finish(i),iyr2,idoy2,im2,id2,ihr2,imin2 ) 
c           put asterisk after time if same as start or end of span  
            ast1 = ' '
            ast2 = ' ' 
            if( e_start(i).le.stime ) ast1 = '*'
            if( e_finish(i).ge.ftime ) ast2 = '*'
             write(iarh,150) iprn,iyr1,idoy1,ihr1,imin1,ast1
     .                     , iyr2,idoy2,ihr2,imin2,ast2
     .                     , eclipse_type(i)
150           format(1x,'PRN ',i2,' Eclipse interval: '
     .                 ,2(i4,i5,2i3,a1,1x),1x,a1)    
c          endif on within span (k=k+1)
           endif 
c       end loop on eclipses
        enddo

c   Write out noninal yaw rate for non eclipsing satellites from start of arc
      elseif ( neclipse .le. 0 ) then
        call dayjul(jdb,yr,doy)
        call monday (doy,mo,day,yr)
        call ds2hms(yr,doy,tb,hr,min,sec)
c        call svnav_read( -1,yr,doy,hr,min,gnss,iprn,insn
c     .                  , frqchn,antbody,sbmass,bias,yawrate
c     .                  , svnstart,svnstop )
        write(iyawtmp,100)iprn,yawrate,bias,yr,mo,day,hr,min,noeclipse
            yaw_entries = yaw_entries + 1

c     endif on eclipses or not
      endif

c  kluge to make sure that at least one entry is made in the yaw file
      if( yaw_entries_old. eq. yaw_entries ) then
        call dayjul(jdb,yr,doy)
        call monday (doy,mo,day,yr)
        call ds2hms(yr,doy,tb,hr,min,sec)
c        call svnav_read( -1,yr,doy,hr,min,gnss,iprn,insn
c     .                 , frqchn,antbody,sbmass,bias,yawrate
c     .                 , svnstart,svnstop )
        write(iyawtmp,100)iprn,yawrate,bias,yr,mo,day,hr,min,noeclipse
            yaw_entries = yaw_entries + 1
      endif

      return
      end

