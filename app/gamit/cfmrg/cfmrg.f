Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      Program CFMRG
c
c.....Combine C-files header information into an M-file to be used in SOLVE

      implicit none

      include '../includes/dimpar.h'

      character*1  explct
      character*1  nbatch
      character*1  batch
      character*1  lowerc,upperc
      character*1  stbias(maxsit)
      character*3  zenmod,gradmod
      character*4  usite(maxsit)
      character*5  is2names(maxprm)
      character*12 sitnams(maxsit),aform
      character*16 cfname(maxsit)
      character*16 filnam
      character*16 fname(maxsit)
      character*16 obfiln(maxsit)
      character*16 tfiln(maxtfl)
      character*17 iscr
      character*20 alabel(maxprm),clabel(maxlab)
      integer is2nums(maxprm)
c      character*5  buf5  - used only for interactive
      character*4  buf4
      character*4  endchk
      character*3 satnam(maxsat)

      integer*4 idms(maxprm),inter(maxsit),ioerr    
     .      , iatchk,icrchk,igcchk,iorchk,iantchk,ieochk
     .      , norb,nparam,nrfile,nsat,itsat,nsite,ntfile,ntpart
     .      , iscrn,iterm,iprnt,iuc,ium,kpart,mpos,mtpart,nepch
     .      , ncfile,isat(maxsat),isprn(maxsat),icnt,isite
     .      , totsat(maxsat),numzen,numgrad
     .      , islot1(maxprm),islot2(maxprm),it0(3,maxsit)
     .      , nlabel,islotc(maxlab),idtzen(maxatm),idtgrad(maxgrad)
     .      , nskip,i,j,k

      logical ldual(maxsit),l1_only
      logical newday
      logical nowrt2
      logical fcheck

      real*8  preval(maxprm)
      real*8  satwgt(maxsat)
      real*8  stawgt(maxsit)
      real*8  t00(3,maxsit)
      real*8  elvcut(maxsit)
c
      real*8 aprcrd,aprclk,aprorb,aprant,aprval,adjust,apreop,aprzen
     .     , aprgrd
      common/mercom/aprval(maxprm),adjust(maxprm),aprcrd(3),
     $              aprclk(maxsit,3),aprorb(maxorb,maxsat),
     $              aprant(3,maxsat),apreop(6),aprzen(maxsit),
     $              aprgrd(maxsit,2)

      logical debug/.false./

c initialize
c
      data fname/maxsit*'                '/
      data usite/maxsit*'    '/
c     why are these weights = 5.?  change to 1.   (not currently used in SOLVE)
c     data satwgt/maxsat*5.d0/
c     data stawgt/maxsit*5.d0/
      data satwgt/maxsat*1.d0/
      data stawgt/maxsit*1.d0/
      data idms/maxprm*0/
      data nsite,nrfile/2*0/ 
      data ncfile/0/
      data iterm,iscrn,iprnt,iuc,ium/5,6,1,11,13/
      data explct,endchk,batch/'e',' end','b'/
      data nskip/1/
   
c.....skip CFMRG if a previous step has failed
      if( fcheck('GAMIT.fatal')) 
     . call report_stat('FATAL','CFMRG','cfmrg',' '
     .                ,'GAMIT.fatal exists: CFMRG not executed',0)
  
c     open the print file
      open (unit=iprnt,file='cfmrg.out',form='formatted'
     .     , status='unknown',iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL','CFMRG','cfmrg','cfmrg.out',
     .  'Error opening cfmrg print file: ',ioerr)
      endif

c     write the version number, date, and operator
      call cversn(iprnt)
             
c     write to screen the location of the complete output
      call report_stat('STATUS','CFMRG','cfmrg',' ',
     .  'Parameter summary written to file cfmrg.out',0)

c     case insensitive variables
      call uppers(explct)
      call uppers(endchk)
      call uppers(batch)

c.....decide operation status
c
c**   hide interactive option
c**   write(iscrn,'(a)') '(i)nteractive or (b)atch mode ? (a1)'
      read (iterm,'(a1)') nbatch
      nbatch = upperc(nbatch)
c
c.....initialize arrays
      do i=1,maxprm
        adjust(i) = 0.d0
        aprval(i) = 0.d0
      enddo     
      do i=1,3
        do j=1,maxsit
            it0(i,j) = 0
            t00(i,j) = 0.d0
            aprclk(j,i) = 0.d0
        enddo
      enddo
      do i=1,maxtfl
        tfiln(i) = ' '
      enddo                                                                    

c-----------------------------------------------------------------------
c    get site codes to analyze
c-----------------------------------------------------------------------

   10 read(5,'(a4)',iostat=ioerr) buf4
      if (ioerr .ne. 0) then
        call report_stat('FATAL','CFMRG','cfmrg',' ',
     .  'Error reading the CFMRG batch file (line 1) ',ioerr)
      endif
      if (buf4.ne.'    ' ) then
          nsite = nsite + 1
          call uppers(buf4)
          usite(nsite) = buf4
          goto 10
       endif
      if (nsite.gt.maxsit) then
        call report_stat('FATAL','CFMRG','cfmrg',' ',
     .  'Error, number of C-files (nsite) greater than (maxsit) ',0)
      endif
c
c-----------------------------------------------------------------------
c     get complete satellite list
c-----------------------------------------------------------------------
      do 110 icnt=1,maxsat
       satnam(icnt)='   '
  110 continue
      write(aform,113) maxsat
  113 format('(',i2,'a3)')
      read (iterm,aform) (satnam(icnt),icnt=1,maxsat)
      write(iprnt,'(1x)')
      write(iprnt,aform) (satnam(icnt),icnt=1,maxsat)
      write(iprnt,'(/)')
      itsat=0
      do 111 icnt=1,maxsat
      if(satnam(icnt).eq.'   ') go to 112
      itsat=itsat+1
      read(satnam(icnt),'(i3)') totsat(icnt)
  111 continue
  112 continue
c      write(6,*) (totsat(icnt),icnt=1,maxsat)

c-----------------------------------------------------------------------
c     read c-files for each day, set a priori data
c-----------------------------------------------------------------------

c.....   cfile loop
c
         do 100 i = 1, nsite+1
            iscr = '                 '
            read(5,'(a17)') iscr
            call uppers(iscr)
            if (iscr(1:4).eq.'    ') goto 200
c           with straight read:
            if (iscr(1:3).eq.endchk(2:4)) goto 210
            mpos=index(iscr,'.')

c.....      site code order

            write(iprnt,'(/)')
            do 50 isite=1,nsite
               if (usite(isite).eq.iscr(mpos-5:mpos-2)) then
c
c.....            open and read the c-file headers
c
                  cfname(isite) = iscr(1:16)
                  call lowers(cfname(isite))
c                  write(iscrn,'(1x)')
                  call getcfs ( iprnt,iuc 
     .                        , cfname(isite),ioerr 
     .                        , ldual(isite),nparam,nsat,ntpart 
     .                        , nepch,isprn,obfiln(isite),tfiln(1)
     .                        , sitnams(isite),elvcut(isite),preval 
     .                        , nlabel,clabel,islotc,norb 
     .                        , inter(isite) 
     .                        , it0(1,isite),t00(1,isite) )   
                  do icnt=1,nsat
                   isat(icnt)=isprn(icnt)
                  enddo 
                  if(debug) print *
     .               ,'CFMRG aft GETCFS isite obfiln tfiln(1)'
     .                                 ,isite,obfiln(isite),tfiln(1)
c
                  if (ioerr.ne.0) then
                     if (lowerc(nbatch).eq.lowerc(batch)) then
                        call report_stat('FATAL','CFMRG','cfmrg'
     .                    ,cfname(isite),'Error opening C-file: ',ioerr)
                     else
                        goto 100
                     endif
                  endif
                  fname(isite) = iscr(1:16)
                  call lowers(fname(isite))
                  ncfile=ncfile+1
                  nrfile=nrfile+1
c
c.....            set apriori site clock, coordinate, zenith delay, and satellite params
c
                  call aprior (isite,nsat,isprn,
     $                         norb,sitnams(isite),preval,itsat,totsat)
c
                  goto 100
               endif
   50       continue   
            call report_stat('WARNING','CFMRG','cfmrg',fname(isite),
     .      'Error, C-file not in station code list: ',0)
  100    continue
  200 continue
      call report_stat('FATAL','CFMRG','cfmrg',' ',
     .'Error, Sessions greater than maxtfl ',0)
c
  210 continue   

c
c-----------------------------------------------------------------------
c     select bias type (implicit or explicit)
c-----------------------------------------------------------------------
c
      call queryb (nsite,stbias)
c
c-----------------------------------------------------------------------
c     open output m-file
c-----------------------------------------------------------------------
c
      read(iterm, '(a16)') filnam
      call lowers (filnam)
      call mopens (filnam,'unknown',ium,ioerr)
      if (ioerr.ne.0) then
         call report_stat('FATAL','CFMRG','cfmrg',filnam,
     .   'Error opening M-file: ',ioerr)
      endif
c
c-----------------------------------------------------------------------
c     queries for coordinate, clock, atmosphere, and orbital/eop partials
c-----------------------------------------------------------------------
c
      call queryp ( nsite,ntpart
     $            , igcchk,iatchk,icrchk,iorchk,iantchk,ieochk
     $            , nepch,zenmod,numzen,idtzen
     $            , gradmod,numgrad,idtgrad )
c
c-----------------------------------------------------------------------
c     fill partial slots for all assigned partials (islot1), labels
c-----------------------------------------------------------------------
c
      call fills1 (igcchk,iatchk,
     $             icrchk,iorchk,iantchk,ieochk,
     $             nsite,itsat,norb,isat,
     $             ldual,usite,fname,stbias,clabel,
     $             mtpart,idms,islot1,alabel,totsat,nsat,
     $             numzen,numgrad,l1_only)
c                 
cd      print *,'CFMRG aft FILLS1 tfname ',tfname 
c**   divert summary to file.
      if( l1_only) write (iprnt,'(/,1x,3a,/)') 'Note: For simplicity in'
     . , ' SOLVE, the L2 bias slots are assigned in CFMRG even for L1-'
     . , 'only C-files'
      write (iprnt,'(1x,i5,1x,i5,1x,a20,1x,d22.15)')
     $                      (k,islot1(k),alabel(k),aprval(k),k=1,mtpart)
      write(iprnt,'(/,a,i2)')    ' Number of orbital parameters: ',norb
      write(iprnt,'(/,a,1x,a3,i4)')
     .          'Zenith model, number of parameters: '
     .          ,zenmod,numzen
      if( numzen.gt.1 )
     .    write(iprnt,'(a,30i4)')  'Epochs of tabular points:  '
     .                         , (idtzen(i),i=1,numzen) 
      write(iprnt,'(/,a,1x,a3,i4)')
     .          'Gradient model, number of parameters: '
     .          ,gradmod,numgrad
      if( numgrad.gt.1 )
     .    write(iprnt,'(a,30i4)')  'Epochs of tabular points:  '
     .                         , (idtgrad(i),i=1,numgrad)
      write(iprnt,'(a,a16)')   ' T-file:   ',tfiln(1) 
      write(iprnt,'(a,i4)')   ' No. sites: ',ncfile

c-----------------------------------------------------------------------
c     write the m-file
c-----------------------------------------------------------------------

c.... write the main m-file header

      ntfile = 1       
      call writm1 ( ium
     .            , nepch,mtpart
     .            , alabel,idms,islot1,aprval,adjust
     .            , itsat,totsat
     .            , nsite,sitnams
     .            , nrfile,obfiln
     .            , ntfile,tfiln,norb)

cd      print *,'CFMRG aft WRITM1 tfname ',tfname 

      write( iprnt,'(/,a,/,a,/,a)') 'C-file subheader:'
     .   ,' isite    k    islot2 is2names is2nums'
     .   ,' -----   ---   ------ -------- -------'
         nowrt2 = .true.
         do 500 isite=1,nsite

c.....         write the m-file session subheader once after finding cfile
c
               if (nowrt2) then
                  nowrt2 = .false. 
cd            print *,'CFRG writing m rec 2 nepch inter nskip ncfile '
cd     .                                   ,nepch,inter,nskip,ncfile 
                  call writm2 (ium,
     $                         it0(1,isite),t00(1,isite),
     $                         nepch,inter(isite),nskip,
     $                         ncfile,cfname,stawgt,
     $                         nsat,satwgt )
               endif
c
c....          fill islot2 (cfile-specific parameter slots)
c
               call fills2 ( icrchk,iatchk,igcchk,iorchk,iantchk,ieochk
     $                     , isite,itsat,nsite,norb
     $                     , ldual,stbias,kpart,islot2,nsat
     $                     , numzen,numgrad,is2names,is2nums)

c
c.....         write the third mfile header (c file subheader)
c                        
cd                print *,'isite gradmod,numgrad kpart '
cd     .          ,isite,gradmod,numgrad,kpart
               call writm3 (ium,
     $                      it0(1,isite),t00(1,isite),
     $                      kpart,islot2,
     $                      elvcut(isite),zenmod,numzen,idtzen,
     $                      gradmod,numgrad,idtgrad )
                      
      write (iprnt,'(1x,i5,1x,i5,2x,i5,4x,a5,3x,i3)')
     $       (isite,k,islot2(k),is2names(k),is2nums(k),k=1,kpart)
  500    continue
  510 continue
      close(unit=ium)
      call report_stat('STATUS','CFMRG','cfmrg',' ',
     .'Normal stop in CFMRG' ,0)
c
c
      end
