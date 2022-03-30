      subroutine writec
     . (namec0,namec1,namev,idelc,imode,isite,imsite,lumf,versn,jstat)

c     Write a new C-file with the cleaned data,
c     but copy the partials from the old one.
c     Update both the observed and the 0-C values

c     better storage structure May 89 K. Feigl
c     Correct M-file reading for omitted stations Kurt 91
c     allow deletion of old C-file: Shimon 91
c     rename output C-file name: Kurt June 91

c     variables:
c     namec0          old C-file name
c     namec1          new C-file name
c     namev           name of V-file for history
c     imode           0 for prefit, 1 for postfit
c     idelc           0 to delete old C-file (namec0)
c                     1 to keep old C-file and write new (namec1)
c                     2 to delete old C-file and rename output to namec0.
c                       (this gives the illusion of overwriting).
c     isite          the CVIEW or SCANDD index number of this site
c     imsite         the M-file index number of this site
c     lumf           unit for M-file, assumed open
c     versn          version of program writing C-file
                            

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer maxpar
      parameter (maxpar = 13+maxsat*maxorb)
c      3 crds + 3 clks + atm + 6 eop = 13

c     standard error flag decoders
      character*4 errcod

      logical bflag_diff(maxsat),pushbias(maxsat),lloel,lbias,lgood
      real*8 addadj,dt,domc
      real*8 preval(maxpar),offarp(3),offsl1(3),offsl2(3)
      real*8 antdaz  ! Antenna aligment from True N (deg).
C MOD TAH 200126: Changed to svantdx(3,2,maxsat) (added L1/L2 index)
      real*8 svantdx(3,2,maxsat)
      real*8 save(maxsav),sec,sod
      real*8 obsv(maxdat),omcs(maxdat),elev,tmpart(maxlab)
      real*8 azim
      integer ierfl,nabove(maxsat)
      integer imode,id,ii,jj,im,iy,msat,nsat,iprn,is
      integer nepoch
      integer lambda(maxsat,maxdat)
      integer ntext,ihr,min,iyr,nlabel,i,nsave,isite,imsite,nparts
      integer npart,ndat
      integer nparam,ndats,kdoy,idoy,j,mtime,iepoch,norb
      integer*4 ioerr
      integer*4 irunt(6),ihnsec
      integer  isnr(maxdat),mprn(maxsat),dattyp(maxdat)
      real*8 t00(3)
      integer it0(3)
      integer jslot(maxlab),islot2(maxprm),kpart,lumf
      character*16 obfiln,tfiln
      character*20 rlabel(maxlab)
      character*80 text(maxtxt)
      character*16 namev,namec0,namec1,uname,namect
      character*80 versn
      integer*4 jsprn(maxsat)
      integer idelc,jstat

      integer li,lo,lv,nblen
      real*8 slip1(maxsat),slip2(maxsat),slip1b(maxsat),slip2b(maxsat)
      real*8 yadd1,yadd2
      real*8 rnhalf,rnqrtr

c     variables for new c-file format oct. 90, nov 91, apr 92, may 95, aug 10
c      integer*4 okmet
      integer okmet
      integer      avlmet
      integer*4    jde,jdr
      integer*4    isessn
      integer*4    ietide,isptide,nclock,data_flag   
      integer*4    niextra,iextra(maxext),nrextra,ncextra
      real*4       ampl1,ampl2,atmlod(3)
      real*8  te,tr,ut1,ut1dot,xp,xpdot,yp,ypdot,psi,psidot,eps,epsdot
      real*8  rclock,pres,temp,relhum,zendel
      real*8  elvcut_cfiles(maxsit),clock(4)
      real*8  azmdot,elvdot,atmdel,tau
      real*8  svcepc,svcl1
      real*8  drate
      real*8  spare(maxspr)    
      real*8  atmlavg(3),hydrolavg(3),rextra(maxext) 
      real*8  latr_sph,lonr,radius
      character*1  skd
      character*3  rcvrsw
      character*4  antmod,svantmod(maxsat),dryzen,wetzen,drymap,wetmap
     .          ,  ionsrc,sitecd
      character*5  frame,precmod,nutmod,gravmod,sprmod,eradmod,antradmod
      character*6  magfield 
      character*10 antmod_snx,svantmod_snx(maxsat)
      character*8        etidemod,atmlmod,otidemod,atmtide,hydrolmod
     .,                  speopmod,cextra(maxext) 
      character*16 jfiln
      character*20 rctype,rcvnum,anttyp,antnum,svantbody(maxsat)
      character*32 sitnamt
      character*256 message
      real*4 swver
      integer nslip,nspare
      integer*2    islip(maxcsb),islpst(maxcsb)
      integer*4 coord_index,atm_index1,atm_index2,clock_index
     .        , orb_index(maxsat),svant_index(maxsat),eop_index
     .        , grad_index(2)
      real*8 grad_parts(2)
      real*8  l1z,l1n,l1e,l2z,l2n,l2e,antaz
      real*8 pi

c     initialise pi
      pi = 4.0d0 * datan(1.0d0)

c     initialise bflag_diff and pushbias arrays
      do i = 1,maxsat
        bflag_diff(i) = .FALSE.
        pushbias(i) = .FALSE.
        nabove(i) = 1
      enddo

c     initialize lambda arrays to avoid problems
      do i = 1,maxsat
         do j = 1,maxdat
            lambda(i,j) = 0
         enddo
      enddo

c     initialize record of slips
      do 15 i=1,maxsat
         slip1(i) = 0.d0
         slip2(i) = 0.d0
         slip1b(i) = 0.d0
         slip2b(i) = 0.d0
 15   continue

c     open a V-file to record cycle slips
      lv = 15
      open (unit   = lv,
     .      file   = namev,
     .      form   = 'formatted',
c    .      status = 'append',
c    .      status = 'append', **Won't work on HP
     .      status = 'unknown',
     .      position = 'append',
     .      iostat = ioerr,
     .      err    = 100)

      if (ioerr .ne. 0) then
         write (6,*) 'WRITEC: Could not open v-file. Using scratch.'
         open (unit = lv,status='scratch')
      endif

      go to 101

  100 continue
      if (ioerr .ne. 0) then
         call report_stat('WARNING','CVIEW','writec',namev,
     .   'Error, opening V-file, continuing',ioerr)
c         stop 'in WRITEC.'  Do not stop but continue
      endif

101   continue

      li = 11
      call copens (namec0,'OLD',li,ioerr)

c     Temporary name for writing in place
      if (idelc .eq. 2) then
         namect = 'tmp.writec'
      else
         namect = namec1
      endif

      lo = 12
      call copens (namect,'UNKNOWN',lo,ioerr)

      if (idelc .eq. 2) then
         write(6,543) namec0(1:nblen(namec0))
  543    format(' overwriting: ',a,$)
      else
         write(6,544) namec0(1:nblen(namec0)),namect(1:nblen(namect))
  544    format(' input c-file: ',a,' output c-file: ',a,$)
      endif

c     copy the header block 1 from input (unit li) to output (unit lo)
      call readc1 (li,ntext,text)

c     add a line to the textual history in the header
      ntext = min0(ntext+1,maxtxt)
      call getdat(irunt(1),irunt(2),irunt(3))
      call gettim(irunt(4),irunt(5),irunt(6),ihnsec )
      call getusr(uname)
      write (text(ntext),2000)
     .versn,uname(1:nblen(uname)),(irunt(i),i=1,6)
 2000 format (a40,1x,'Run by ',a,1x,'on',1x,
     .   i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2)

      call writc1 (lo,ntext,text)

      call readc2 (li
     .,            sitnamt,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,            npart,norb,gnss,nsat,jsprn,fL1,fL2
     .,            ndat,dattyp,lambda
     .,            skd,nepoch,inter,ircint,mtime,isessn
     .,            iy,im,id,ihr,min,sec
     .,            offarp,offsl1,offsl2, antdaz, svantdx
     .,            obfiln,tfiln,jfiln
     .,            frame,precmod,nutmod,gravmod,sprmod,eradmod,antradmod
     .,            ietide,isptide,speopmod
     .,            etidemod,otidemod,atmtide,atmlmod,hydrolmod  
     .,            atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap  
     .,            ionsrc,magfield
     .,            antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,            elvcut_cfiles(isite),nclock,clock
     .,            jde,te,jdr,tr
     .,            ut1,ut1dot,xp,xpdot,yp,ypdot
     .,            psi,psidot,eps,epsdot
     .,            avlmet
     .,            nslip,islip,islpst  
     .,            niextra,iextra,nrextra,rextra,ncextra,cextra)

c Kludge to fix cfile header bug introduced by autcln.
c convert cfile header cutoff from rad to deg if < 1....
      if ( elvcut_cfiles(isite) .lt. 1.d0 ) then
         elvcut_cfiles(isite) = elvcut_cfiles(isite)*180.d0/pi
      endif
      call writc2 (lo
     .,            sitnamt,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,            npart,norb,gnss,nsat,jsprn,fL1,fL2
     .,            ndat,dattyp,lambda
     .,            skd,nepoch,inter,ircint,mtime,isessn
     .,            iy,im,id,ihr,min,sec
     .,            offarp,offsl1,offsl2,svantdx
     .,            obfiln,tfiln,jfiln
     .,            frame,precmod,nutmod,gravmod,sprmod,eradmod,antradmod
     .,            ietide,isptide,speopmod
     .,            etidemod,otidemod,atmtide,atmlmod,hydrolmod  
     .,            atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap
     .,            ionsrc,magfield
     .,            antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,            elvcut_cfiles(isite),nclock,clock
     .,            jde,te,jdr,tr
     .,            ut1,ut1dot,xp,xpdot,yp,ypdot
     .,            psi,psidot,eps,epsdot
     .,            avlmet
     .,            nslip,islip,islpst     
     .,            niextra,iextra,nrextra,rextra,ncextra,cextra)
      call readc3 (li,nlabel,nparam,jslot,idms,rlabel,preval)
      call writc3 (lo,nlabel,nparam,jslot,idms,rlabel,preval)

      if(imode.eq.1) then
         call readm3 (lumf
     $             ,   it0,t00
     $             ,   kpart,islot2
     $             ,   elvcut(isite),zenmod,numzen,idtzen
     .             ,   gradmod,numgrad,idtgrad )

c     get mapping from parameter slots to C-file partial slots
         call pointers( imsite,coord_index,atm_index1,atm_index2
     .                , clock_index,orb_index,svant_index,eop_index
     .                , grad_index )
      endif

c     if no partials are being used,
c     get 'tt00' and 'it0' from c-file, not m-file.
      if (imode.eq.0) then
         iit0(1)= im
         iit0(2)= id
         iit0(3)= iy
         tt00(1)= ihr
         tt00(2)= min
         tt00(3)= sec
      endif
c
      do 50 ii = 1,nepoch
         call readc4 (li
     .,               msat,mprn
     .,               iepoch,iyr,kdoy,sod,rclock
     .,               okmet,zendel,pres,temp,relhum,atmlod 
     .,               sitecd,latr_sph,lonr,radius
     .,               l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,               nsave,save  )

         call writc4 (lo
     .,               msat,mprn
     .,              iepoch,iyr,kdoy,sod,rclock
     .,               okmet,zendel,pres,temp,relhum,atmlod
     .,               sitecd,latr_sph,lonr,radius
     .,               l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,               nsave,save )


c        seconds past initial epoch
         if (iepoch .ne. ii) then
            write(message,116)ii,iepoch
  116       format ('Epoch counter mismatch.',2(i5,1x))
            call report_stat('STATUS','CVIEW','writec',' ',message,0)
         else
            if (kdoy .gt. 0) then
               dt = 86400.0d0*dble(kdoy-idoy(iy,im,id))
     .            + sod - (3600.0d0*ihr+60.0d0*min+sec)
            else
               dt = dble(iepoch - 1) * dble(inter)
            endif
         endif

c        add receiver clock offset
c        RCLOCK = observed - true
         dt = dt - rclock

c        read the data, one sat at a time
         do 20 jj=1,msat

c           new C-file format
            call readc5 (li
     .,                  iprn
     .,                  elev,elvdot
     .,                  azim,azmdot
     .,                  atmdel
     .,                  svcepc,svcl1
     .,                  tau,drate
     .,                  ierfl,data_flag
     .,                  ndats
     .,                  obsv
     .,                  omcs
     .,                  isnr,ampl1,ampl2
     .,                  nspare
     .,                  spare
     .,                  nparts
     .,                  tmpart)



c           which is the correct PRN (sat id #)?
            is = 1
            do 16 j = 1,nsat
               if (iprn .eq. jsprn(is)) then
                  goto 17
               else
                  is = is + 1
               endif
 16         continue     
            write(message,'(a,i3)') 'Cannot find PRN ',jsprn(is)
            call report_stat('FATAL','CVIEW','writec',' ',message,0)
 17         continue
c
c           compare new and old error flags,
c           record changes in error flags to V-file,
c           must remove low elev flags put in in readc to account for
c           solve elevation cutoff.....

c   check to see if solve cutoff > cfiles cutoff
            if ( elvcut(isite) .gt. elvcut_cfiles(isite) ) then
              if ( elev .lt. elvcut(isite)*pi/180.d0 ) then
                nabove(is) = 0
                if ( ierfl .eq. ierr(iepoch,is,isite) ) then
                      ierr(iepoch,is,isite) = ierfl
                elseif ( lloel(ierr(iepoch,is,isite)) ) then
                   if ( pushbias(is) ) then
                      ierr(iepoch,is,isite) = 10
                      pushbias(is) = .FALSE.
                   else
                     ierr(iepoch,is,isite) = ierfl
                   endif
                elseif ( lbias(ierfl) ) then
                   pushbias(is) = .TRUE.
                   bflag_diff(is) = .TRUE.
                else
                   bflag_diff(is) = .TRUE.
                endif
              else
                nabove(is) = nabove(is) + 1
                if ( .not. bflag_diff(is) .and. nabove(is) .eq. 1 ) then
                   if ( lbias(ierr(iepoch,is,isite)) ) then
                      ierr(iepoch,is,isite) = ierfl
                      bflag_diff(is) = .FALSE.
                      pushbias(is) = .FALSE.
                   elseif ( .not. lgood(ierr(iepoch,is,isite)) ) then
                      if ( ierfl .eq. ierr(iepoch,is,isite) ) then
                        bflag_diff(is) = .FALSE.
                        pushbias(is) = .FALSE.
                        nabove(is) = 0
                      else
                        bflag_diff(is) = .FALSE.
                        pushbias(is) = .FALSE.
                      endif
                   else
                      bflag_diff(is) = .FALSE.
                      pushbias(is) = .FALSE.
                   endif
                else
                   bflag_diff(is) = .FALSE.
                   pushbias(is) = .FALSE.
                endif
              endif
            endif

            if (ierfl .ne. ierr(iepoch,is,isite)) then
                write (lv,12) namec0(2:6),iprn,iepoch,
     .          errcod(ierfl),errcod(ierr(iepoch,is,isite))
 12             format (1x,a5,1x,i3,1x,i4,2(1x,20x    ),2(1x,a4))
            endif

c           new error flags
            ierfl = ierr(iepoch,is,isite)

c    compute the atmospheric gradient partials
            call atm_grad_part(elev,azim,grad_parts)

c           find change due to editing
c           for this, we recreate the original observable
            if (imode .eq. 1) then
c               domc=addadj( is,isite,dt,iepoch,tpart )
               domc=addadj( is,dt,iepoch,norb,nparts,tmpart
     .                    , coord_index,atm_index1,atm_index2
     .                    , clock_index,orb_index(is),svant_index(is)
     .                    , eop_index,grad_index,grad_parts )
            else
               domc = 0.0d0
            endif

c           allow slips of magnitude 1/(2*lambda) cycle
            if (lambda(is,1) .ne. 0) then
             if (iabs(lambda(is,1)) .eq. 2) then
                  yadd1 = rnqrtr(yl1(ii,is,isite) - (omcs(1)-domc))
               else
                  yadd1 = rnhalf(yl1(ii,is,isite) - (omcs(1)-domc))
               endif
               omcs(1) = omcs(1) + yadd1
               obsv(1) = obsv(1) + yadd1
            else
               yadd1 = 0.0d0
            endif
            if (lambda(is,2) .ne. 0) then
               if (iabs(lambda(is,2)) .eq. 2) then
                  yadd2 = rnqrtr(yl2(ii,is,isite) 
     .                    - (omcs(2)-gear(jsat1)*domc))
               else
                  yadd2 = rnhalf(yl2(ii,is,isite) 
     .                    - (omcs(2)-gear(jsat1)*domc))
               endif
               omcs(2) = omcs(2) + yadd2
               obsv(2) = obsv(2) + yadd2
            else
               yadd2 = 0.0d0
            endif

c           the pseudoranges in obsv(3) and obsv(4) remain unchanged

c           new C-file format
            call writc5 (lo
     .,                  iprn
     .,                  elev,elvdot
     .,                  azim,azmdot
     .,                  atmdel
     .,                  svcepc,svcl1
     .,                  tau,drate
     .,                  ierfl,data_flag
     .,                  ndats
     .,                  obsv
     .,                  omcs
     .,                  isnr,ampl1,ampl2
     .,                  nspare
     .,                  spare
     .,                  nparts
     .,                  tmpart)

c           record the edited change to the V-file
            slip1(is) = yadd1
            slip2(is) = yadd2

c           If the change is not the same as the change at the
c           last epoch, assume that this is a cycle slip.
            if (dabs(slip1(is)-slip1b(is)) .gt. 0.10d0 .or.
     .          dabs(slip2(is)-slip2b(is)) .gt. 0.10d0) then
                write (lv,810) namec0(2:6),jsprn(is),iepoch,yadd1,yadd2
 810            format (1x,a5,1x,i3,1x,i4,2(1x,1pe20.4))
            endif

c           buffer the values of the slips
            slip1b(is) = slip1(is)
            slip2b(is) = slip2(is)
   20    continue
         if (mod(ii,100) .eq. 1) write (6,'(''.'',$)')
   50 continue

      close (lo)

c     deal with C-file names
      if (idelc .eq. 0) then
c        option to delete old c-file
         close (li, iostat=ioerr, status = 'delete')
         if (ioerr .eq. 0) then
            print *,' deleted: ',namec0(1:nblen(namec0))
         else
            print *,' ERROR!'
            print *,'WRITEC: error deleting ',namec0(1:nblen(namec0))
            call ferror (ioerr,6)
            close (li)
         endif
      else if (idelc .eq. 2) then
c        option to overwrite old C-file
         close (li, iostat=ioerr, status = 'delete')
         if (ioerr .ne. 0) then
            print *,' ERROR!'
            print *,'WRITEC: error deleting ',namec0(1:nblen(namec0))
            call ferror (ioerr,6)
            close (li)
         endif
c        Rename output file with name of input file,
c        creating the illusion of having worked "in place"
         call irename (namect,namec0,jstat)
         if (ioerr .eq. 0 .and. jstat. eq. 0) print *,' OK'
      else
         close (li)
         print *,' OK'
      endif


      close (lv)

      return
      end
