      Program AUTECL
c
c Purpose: To generate eclipse edit lines to be written to standard output
c
c S. McClusky simon@chandler.mit.edu April 1995
c R. King rwk@chandler.mit.edu  January 2014: Modify to read the binary y-file generated
c    by yawtab with event flags, rather than arcout.DDD. In this version, edit lines
c    are written for all SVs on the t-/y-files, even if they don't appear in session.info
c    (simpler and no harm done).  

c input: binary y-file
c        eclipse [+ time after eclipse]: delete eclipse data plus attitude recovery period data
c        noeclipse [+ time after eclipse]: just delete recovery period data
c
c output: file to be appended to autcln.cmd file
c
      implicit none
c
      include '../includes/dimpar.h'
c
      character*1 eclipse_out                    
      character*5 block
      character*80 arg(4),line
              
      logical debug,fcheck,endsv

      integer luy,lses,lout,ioerr,nversn
     .      , after_e,check_flg
     .      , iyr,idoy,is,inter,nepoch,nprn
     .      , ihr,imin,imo,iday,ievent(maxsat),ievent_last(maxsat)
     .      , nsat,itsat(maxsat),isvn(maxsat),iblk(maxsat)
     .      , sess_prn(maxsat),ep_start,ep_end,i
  
      real*8 pjd,attit(maxsat) 
     .     , t_start,t_end,sess_start,sess_end
     .     , start_night(maxsat),end_night(maxsat)
     .     , start_post(maxsat),end_post(maxsat)
     .     , start_noon(maxsat),end_noon(maxsat)

      character*5 ablk                       
      character*16 tfile,yfile  
      character*20 svantbody(maxsat)
      character*47 comment 
      character*256 message

      logical night, noon, post 

c     function
      integer*4 julday 

c     Skip if a previous step has failed
      if( fcheck('GAMIT.fatal') )  
     .  call report_stat('FATAL','AUTECL','utils/autecl',' '
     .                  ,'GAMIT.fatal exists: AUTECL not executed',0)
                
 
c  Logical unit numbers for; y-file (luy), session.info (lses), std out (lout)
c
      luy = 10
      lses = 11
      lout = 6

c  Initialize to keep G77 from complaining
      sess_start = 0.d0 
      sess_end = 0.d0   
      comment = ' '


c  Read input off command line if given
c
      call rcpar(1,arg(1))
      yfile = arg(1)
      call rcpar(2,arg(2))
      read(arg(2),'(i3)') idoy
      call rcpar(3,arg(3))
      eclipse_out = arg(3)(1:1)
      call rcpar(4,arg(4))
      read(arg(4),'(i3)') after_e
cd      print *,'AUTECL idoy eclipse_out aftr_e ',idoy,eclipse_out,after_e
 
      if (arg(1).eq." ") then
c       interactive
        write(*,10)
10    format(/,'########################################################
     .######',/,'                    Program AUTECL              '
     .,/,'Purpose: To generate autocln edit lines for removing eclipse'
     .,/,'         data from gamit solution.'
     .,/,'Usage:(1) Interactive (just answer the questions below)'
     .,/,'      (2) Command Line Input ( yfile  doy  Y or N  time)'
     .,/,' yfile = binary y-file from yawtab, e.g. yigsft.168 '
     .,/,' doy = day number of year'
     .,/,' Y/N = [Y] delete or [N] not delete eclipse period data'
     .,/,' time = Minutes of post eclipse attitude recovery data edited'
     ./,'##############################################################'
     .,/)
        print*, ('Enter y-file name: ')
        read(*, '(a)') yfile
        print*, ('Day number [doy]: ')
        read(*, '(i3)') idoy
        print*, ('Delete eclipse period data [Y/N]: ')
        read(*, '(a)') eclipse_out
        print*, ('Minutes of post eclipse attitude recovery data cut: ')
        read(*, '(i3)') after_e
      endif
                  

c  Set the editing controls 
 
      night = .false.
      post = .true.
      noon = .false.     
      if( eclipse_out.eq.'Y'.or.eclipse_out.eq.'y' ) night = .true.


c  Open the required y-file, session.info, and autcln.cmd input/output files.

      open(unit=luy,file=yfile,form='unformatted',access='sequential'
     .           ,status='old',iostat=ioerr)  
      if( ioerr.ne.0 )  call report_stat('FATAL','AUTECL','utils/autecl'
     .      ,yfile,'Error opening input y-file ',ioerr)
      open(unit=lses,file='session.info',status='old',iostat=ioerr)
      if( ioerr.ne.0 ) call report_stat('FATAL','AUTECL','utils/autecl'
     .    ,'session.info','Error opening session.info file',ioerr)
      write(lout,'(a)') 
     .     '## START OF YAW-MODEL DATA EDITS  ## '

                 
c  Read session.info to get observation interval, session start and duration

c     set check_flg to tell rsesfo to match only the day-of-year                                                      
      debug = .false.
      check_flg = 1
      call rsesfo( lses,debug,check_flg,iyr,idoy,is
     .         , ihr,imin,inter,nepoch,nprn,sess_prn )
cd     print *,'Aft RSESFO: y doy h m int epc ',iyr,idoy
cd    .          ,ihr,imin,inter,nepoch
      call check_y2k(iyr)    
      call monday(idoy,imo,iday,iyr)
cd      print *,'iyr imo iday ihr imin ',iyr,imo,iday,ihr,imin
      sess_start = julday(imo,iday,iyr)
     .             + ihr/24.d0 + imin/1440.d0
      sess_end = sess_start + ((inter*nepoch)/86400.d0)    
cd      print *,'sess_start sess_end',sess_start,sess_end

             
c  Read the y-file header to get the SV PRN, SVN, and BLK numbers
      
c     first read just the version number so that we have the correct header format
      read(luy,iostat=ioerr) nversn
      if( ioerr.ne.0 ) then  
        call report_stat('FATAL','AUTECL','utils/autecl',yfile
     .        ,'Error reading version number on the y-file',ioerr)
      elseif( nversn.eq.1051 ) then 
        rewind(luy)     
        read(luy,iostat=ioerr) nversn,tfile,t_start,t_end,inter,nepoch
     .         , nsat,(itsat(i),isvn(i),iblk(i),i=1,nsat)
      elseif( nversn.eq.1061) then 
        rewind(luy)
        read(luy) nversn,tfile,t_start,t_end,inter,nepoch
     .      ,nsat,(itsat(i),isvn(i),svantbody(i),i=1,nsat)
      else  
        write(message,'(a,i6,a)') 
     .     'Invalid y-file version number (',nversn,')'
        call report_stat('FATAL','AUTECL','utils/autecl',yfile
     .        ,message,ioerr)
      endif

c Read the y-file to get the satellite yaw events
c   keyed by ievent:  0=nominal yaw  1=night  -1=night recovery  2=noon turn                     
c   set the start time to zero at the completion of each event
                                                                      
      ioerr = 0                          
      do i=1,maxsat
        ievent(i) = 0 
        ievent_last(i) = 0 
        start_night(i) = 0.d0
        end_night(i) = 0.d0   
        start_post(i) = 0.d0
        end_post(i) = 0.d0
        start_noon(i) = 0.d0
        end_noon(i) = 0.d0
      enddo 
           
cd      print *,'night post noon ',night,post,noon
      do while (ioerr .eq. 0 .and. pjd.le.sess_end )
         
        read(luy,iostat=ioerr) pjd,(attit(i),ievent(i),i=1,nsat)
cd        do i=1,nsat
cd          if( ievent(i).ne.0 ) print *,'pjd i ievent '
cd     .           ,pjd,itsat(i),ievent(i)
cd        enddo                                                            

        do i=1,nsat 
          if( night ) then
            if( ievent_last(i).eq.0  .and. ievent(i).eq.1 ) then
               start_night(i) = pjd
cd               print *,'start_night ',pjd,itsat(i)
            elseif( ievent_last(i).eq.1 .and.
     .             (ievent(i).eq.-1.or.ievent(i).eq.0) ) then
               end_night(i) = pjd          
cd                 print *,'end_night ',pjd,itsat(i)
            endif
          endif
          if( post ) then
            if( (ievent_last(i).eq.0.or.ievent_last(i).eq.1)  .and. 
     .            ievent(i).eq.-1 ) then
              start_post(i) = pjd
cd              print *,'prn start_post ',itsat(i),start_post(i)
            elseif( ievent_last(i).eq.-1 .and.ievent(i).eq.0 ) then
              end_post(i) = pjd       
cd              print *,'prn end_post ',itsat(i),end_post(i)
            endif
          endif
          if( noon ) then
            if( ievent_last(i).eq.0  .and. ievent(i).eq.2 ) then
              start_noon(i) = pjd    
cd              print *,'start_noon ',pjd,itsat(i)
            elseif( ievent_last(i).eq.2 .and.ievent(i).eq.0 ) then
              end_noon(i) = pjd                 
cd              print *,'end_noon ',pjd,itsat(i)
            endif         
          endif   
          if( end_night(i).ne.0.d0 ) then
            call get_epoch(sess_start,start_night(i),inter,ep_start)         
cd            print *,'prn sess_start start_night ep_start '
cd     .             ,itsat(i),sess_start,start_night(i),ep_start 
            call get_epoch(sess_start,end_night(i),inter,ep_end)  
cd            print *,'prn sess_start end_night  ep_end '
cd     .             ,itsat(i),sess_start,end_night(i),ep_end 
            if( ep_end.gt.0  ) then
              if( ep_start.le.0 ) ep_start = 1 
              if( nversn.eq.1051 ) then 
                call block_name(iblk(i),ablk)
                comment = 
     .            '  ! Blk '//ablk//' eclipse or night turn            '
              elseif( nversn.eq.1061) then
                comment = '  ! '//svantbody(i)//' eclipse or night turn'
              endif
              write(lout,'(a,i2,2i6,a)') ' edit_site_sv  all '
     .                ,itsat(i),ep_start,ep_end,comment
            endif
            start_night(i) = 0.d0
            end_night(i) = 0.d0
          endif         
          if( end_post(i).ne.0.d0 ) then
            call get_epoch(sess_start,start_post(i),inter,ep_start) 
cd            print *,'prn sess_start start_post  ep_start '
cd     .             ,itsat(i),sess_start,start_post(i),ep_start
            call get_epoch(sess_start,end_post(i),inter,ep_end)    
cd            print *,'prn sess_start end_post ep_end '
cd     .             ,itsat(i),sess_start,end_post(i),ep_end 
           if( ep_end.gt.0 ) then    
              if( ep_start.le.0 ) ep_start = 1          
              if( nversn.eq.1051 ) then 
                call block_name(iblk(i),ablk)
                comment = 
     .            '  ! Blk '//ablk//' post-eclipse or turn             '
              elseif( nversn.eq.1061) then 
                comment = '  ! '//svantbody(i)//' post-eclipse or turn '
              endif
              write(lout,'(a,i2,2i6,a)') ' edit_site_sv  all '
     .                ,itsat(i),ep_start,ep_end,comment
            endif
            start_post(i) = 0.d0
            end_post(i) = 0.d0
          endif
          if( end_noon(i).ne.0.d0 ) then
            call get_epoch(sess_start,start_noon(i),inter,ep_start)
cd            print *,'prn sess_start start_noon inter ep_start '
cd     .             ,itsat(i),sess_start,start_noon(i),inter,ep_start
            call get_epoch(sess_start,end_noon(i),inter,ep_end)  
            if( ep_end.gt.0 ) then  
              if( ep_start.le.0 ) ep_start = 1
              if( nversn.eq.1051 ) then 
                call block_name(iblk(i),ablk)
                comment =
     .            '  ! Blk '//ablk//' noon turn                        '
              elseif( nversn.eq.1061 ) then
               comment = '  ! '//svantbody(i)//' noon turn             '
              endif
              write(lout,'(a,i2,2i6,a)') ' edit_site_sv  all '
     .               ,itsat(i),ep_start,ep_end,comment
            endif
            start_noon(i) = 0.d0
            end_noon(i) = 0.d0 
          endif   
        ievent_last(i) = ievent(i)
c       end loop on SVs at this epoch
cd        if( pjd.gt.2455730.3 ) stop 
        enddo           
                                      
c     end loop through file
      enddo                                   

      end

c************************************************************************************************8
       
      Subroutine get_epoch(t_start,t,inter,iepoch)
                                           
      implicit none
      integer*4 inter,iepoch
      real*8 t_start,t,epoch
      
      epoch = (t - t_start)*86400.d0/inter
      iepoch = idnint (epoch)
cd      print *,'epoch iepoch ',epoch,iepoch

      return
      end


