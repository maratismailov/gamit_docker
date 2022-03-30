
      subroutine readc
     .   (lunit,lumf,filnam,imode,isite,imsite,nepoch,nsat,ndats)
c
c     read a C-file called filnam on logical unit lunit
c
c     input:
c            lunit    logical file unit
c            lumf     m-file logical file unit
c            filnam  file name
c            imode   = 0 pre-fit, no partials
c                    = 1 post-fit with partials
c            isite   the CVIEW or SCANDD index number of this site
c            imsite  the M-file index number of this site
c     output
c            nepoch total number of epochs
c            nsat   total number of satellites
c            ndats   number of different observables
c                   e.g. 4 for L1,L2,P1 and P2
c
c
c     cleaned up by K. Feigl May 89
c     roots:
c        readc - R.W. King

      implicit none
c
      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

c     standard error flag decoders
      logical lmarg,lgood,lbias,pushbias(maxsat)

      integer maxpar
      parameter (maxpar = 13+maxsat*maxorb)
c      3 crds + 1 clks + 1 atm + 6 eop + 2 grad = 13
c
      real*8 addadj,dt,domc
     .     , preval(maxpar),offarp(3),offsl1(3),offsl2(3)
C MOD TAH 200126: Changed to svantdx(3,2,maxsat) (added L1/L2 index)
     .     , svantdx(3,2,maxsat)
     .     , rextra(maxext),save(maxsav),sec,sod
     .     , obsv(maxdat),omcs(maxdat),elev,tmpart(maxlab)
     .     , r1,r2,azim,grad_parts(2)
     .     , t00(3),latr_sph,lonr,radius
      integer*4 ierfl,imode,id,ii,jj,im,iy,msat,nsat,iprn,is
     .        , nepoch,lambda(maxsat,maxdat)
     .        , lunit,ntext,ihr,min,iyr,nlabel,i,nsave,isite,imsite
     .        , nparts,npart,ndat,nparam,ndats
     .        , jdoy,idoy,j,mtime,iepoch,norb,ioerr
     .        , isnr(maxdat),mprn(maxsat),dattyp(maxdat)
     .        , coord_index,atm_index1,atm_index2,clock_index
     .        , orb_index(maxsat),svant_index(maxsat),eop_index
     .        , grad_index(2),it0(3)
     .        , jslot(maxlab),islot2(maxprm),kpart,lumf
     .        , len,rcpar 
      integer*2 ii2              
      character*4 sitecd 
      character*32 sitnam
      character*16 filnam,obfiln,tfiln
      character*20 rlabel(maxlab)
      character*80 text(maxtxt),prog_name
      character*256 message

c     variables for new c-file format oct. 90, nov 91, apr 92, may 95, feb 03
      integer*4 okmet
      integer      avlmet
      integer*4    jde,jdr
      integer*4    isessn
      integer*4    ietide,isptide,nclock,data_flag 
      integer*4    niextra,iextra(maxext),nrextra,ncextra
      real*4       ampl1,ampl2,atmlod(3)
      real*8  te,tr,ut1,ut1dot,xp,xpdot,yp,ypdot,phi,phidot,eps,epsdot
      real*8  rclock,pres,temp,relhum,zendel
      real*8  elvcut_cfiles(maxsit),clock(4)
      real*8  azmdot,elvdot,nadang
      real*8  atmdel,tau,drate
      real*8  svcepc,svcl1
      real*8  spare(maxspr)  
      real*8  atmlavg(3),hydrolavg(3)
      real*8 antdaz  ! Antenna aligment from True N (deg).
      character*1  skd
      character*3  rcvrsw
      character*4  antmod,svantmod(maxsat),dryzen,wetzen,drymap,wetmap
     .            ,ionsrc
      character*5  frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
      character*6  magfield
      character*8  etidemod,atmlmod,otidemod,atmtide,hydrolmod
     .,            speopmod,cextra(maxext) 
      character*10 antmod_snx,svantmod_snx(maxsat)
      character*16 jfiln
      character*20 rctype,rcvnum,anttyp,antnum,svantbody(maxsat)
      character*32 sitnamt
      real*4 swver
      integer nslip,nspare
      integer*2    islip(maxcsb),islpst(maxcsb)
      integer jsprn(maxsat)
      real*8  l1z,l1n,l1e,l2z,l2n,l2e,antaz
      real*8 pi
                 
c     functions to return I*2 values (necessary for gcc compiler)
      integer*2 min02, max02

c     initialise pi
      pi = 4.0d0 * datan(1.0d0)

c     function to calculate multipath path delay
CKF   real*8 multip


c     get program name calling readc
      len = rcpar(0,prog_name)

c     initialize lambda arrays to avoid problems
      do i = 1,maxsat
         do j = 1,maxdat
            lambda(i,j) = 0
         enddo
      enddo

c     initialize start and end indexes
      do 3 i = 1,maxsat
         kk0(i,isite) = maxepc
         kk1(i,isite) = 1
  3   continue

c     initialise pushbias array
      do 4 i = 1,maxsat
         pushbias(i) = .FALSE.
  4   continue

      call report_stat('STATUS',prog_name,'readc',filnam
     .                ,'Opening C-file: ',0)

      call copens (filnam,'OLD',lunit,ioerr)

      call readc1 (lunit,ntext,text)
      write (6,'(a80)') (text(i),i=1,ntext)

c     new format
      call readc2 (lunit
     .,            sitnamt,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,            npart,norb,gnss,nsat,jsprn,fL1,fL2
     .,            ndat,dattyp,lambda
     .,            skd,nepoch,inter,ircint,mtime,isessn
     .,            iy,im,id,ihr,min,sec
     .,            offarp,offsl1,offsl2,antdaz,svantdx
     .,            obfiln,tfiln,jfiln
     .,            frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod   
     .,            ietide,isptide,speopmod
     .,            etidemod,otidemod,atmtide,atmlmod,hydrolmod 
     .,            atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap
     .,            ionsrc,magfield
     .,            antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,            elvcut_cfiles(isite),nclock,clock
     .,            jde,te,jdr,tr
     .,            ut1,ut1dot,xp,xpdot,yp,ypdot
     .,            phi,phidot,eps,epsdot
     .,            avlmet
     .,            nslip,islip,islpst  
     .,            niextra,iextra,nrextra,rextra,ncextra,cextra)

      SITNAM = SITNAMT(1:12)

c     copy lambda array
c     note conversion from integer*4 to integer*2
c     copy lambda array
      do 135 i=1,nsat
         do 133 j=1,ndat
            lambds(isite,i,j) = lambda(i,j)
 133     continue
 135  continue

c     determines the number of epoches between good data points
      if((inter .ge. ircint) .or. (inter .le. 0)) then
           inext = 1
      else
           inext = ircint/inter
      endif
      jnext(isite) = inext

c     copy isprn array
      do 137 i=1,nsat
         isprn(i)=jsprn(i)
  137 continue

      if (ndat .gt. maxdat) then
         write (message,117) ndats,maxdat
  117    format ('File contains more data types (',i5,') than ',
     .   'program dimensions MAXDAT: (',i5,')')
         close (lunit)
         call report_stat('FATAL',prog_name,'readc',filnam,message,0)
      else if (ndat .le. 2) then
         write (message,118)
  118    format ('C-file contains no pseudoranges')
         call report_stat('WARNING',prog_name,'readc',filnam,message,0)
      endif

      if (nepoch .gt. maxepc) then
         write (message,110) nepoch,maxepc
  110    format ('File contains more epochs (',i5,') than ',
     .   'program dimensions MAXEPC: (',i5,')')
         close (lunit)
         call report_stat('FATAL',prog_name,'readc',filnam,message,0)
      endif

      if (nepoch .gt. 32768) then
         write (message,112) nepoch
  112    format ('File contains more epochs (',i5,') than ',
     .   'Integer*2 allows. HELP!!!: 32768' )
         close (lunit)
         call report_stat('FATAL',prog_name,'readc',filnam,message,0)
      endif

      if (nsat .gt. maxsat) then
         write (message,115) nsat,maxsat
  115    format ('File contains more sats (',i5,') than ',
     .   'program dimensions MAXSAT: (',i5,')')
         close (lunit)
         call report_stat('FATAL',prog_name,'readc',filnam,message,0)
      endif

c     unchanged from old to new formats
      call readc3 (lunit,nlabel,nparam,jslot,idms,rlabel,preval)

c     get partials if we are in postfit mode
      if(imode.eq.1) then
         call readm3 (lumf
     $             ,  it0,t00
     $             ,  kpart,islot2
     .             ,  elvcut(isite),zenmod,numzen,idtzen
     .             ,  gradmod,numgrad,idtgrad )

         if( numzen.gt.maxatm ) then
           write(message,130)numzen,maxatm
 130       format( 'File contains more zenith-delay parameters (',i2,
     .     ') than program dimension MAXATM: (',i2,')')
           call report_stat('FATAL',prog_name,'readc',filnam,message,0)
         endif

c        get mapping from parameter slots to C-file partial slots
         call pointers( imsite,coord_index,atm_index1,atm_index2
     .                , clock_index,orb_index,svant_index,eop_index
     .                , grad_index )
c         if( imsite.eq.1 )
c     .       print *,'READC mtpart islot ',mtpart,(islot(i),i=1,mtpart)
c         print *,'READC imsite indices: crd atm clk svant eop orb '
c     .  ,imsite,coord_index,atm_index1,atm_index2,clock_index,svant_index
c     .  , eop_index,orb_index

       endif

       if ( imode .eq. 1 ) then
c        print solve and cfile elevation angle cutoffs to screen
         write(6,150)elvcut_cfiles(isite),elvcut(isite)
 150     format(1x,'C-file elevation angle cutoff: ',f7.3,
     .         ' Solve elevation angle cutoff: ',f7.3)
       else
c        print cfile elevation angle cutoff to screen
         write(6,155)elvcut_cfiles(isite)
 155     format(1x,'C-file elevation angle cutoff: ',f7.3)
       endif

      write(6,201) sitnam,nsat,nepoch
 201  format (1x,a12,4x,i2,' sats ',i5,' epochs')
      write(6,202) (isprn(j),j=1,nsat)
 202  format (1x,'PRNs   :',32(i2,1x,:))
      do 195 i=1,ndat
         write(6,203) (lambds(isite,j,i),j=1,nsat)
 195  continue
 203  format (1x,'lambdas:',32(i2,1x,:))

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

c     Current CFMRG and SOLVE assume the same number of zenith delays for all sites
c     read epoch record
      do 50 ii = 1,nepoch
c        new format
         call readc4 (lunit
     .,               msat,mprn
     .,               iepoch,iyr,jdoy,sod,rclock
     .,               okmet,zendel,pres,temp,relhum,atmlod
     .,               sitecd,latr_sph,lonr,radius
     .,               l1z,l1n,l1e,l2z,l2n,l2e,antaz 
     .,               nsave,save )

c        seconds past initial epoch
         if (iepoch .ne. ii) then
            write(message,116)ii,iepoch
  116       format ('Epoch counter mismatch.',2(i5,1x))
            call report_stat('FATAL',prog_name,'readc',filnam,message,0)
         else
            if (jdoy .gt. 0) then
               dt = 86400.0d0*dble(jdoy-idoy(iy,im,id))
     .            + sod - (3600.0d0*ihr+60.0d0*min+sec)
            else
               dt = dble(iepoch - 1) * dble(inter)
            endif
         endif

c        add receiver clock offset
c        RCLOCK = observed - true
         tag(ii,isite) = dt - rclock

c        initialize the error flags to avoid problems
         do 13 jj=1,nsat
            ierr(ii,jj,isite) = ignone
 13      continue

c        read the data, one sat at a time
         do 20 jj=1,msat
c           new C-file format
            call readc5 (lunit
     .,                  iprn
     .,                  elev,elvdot,azim,azmdot,nadang
     .,                  atmdel, svcepc,svcl1
     .,                  tau,drate
     .,                  ierfl,data_flag
     .,                  ndats,obsv,omcs,isnr,ampl1,ampl2
     .,                  nspare,spare
     .,                  nparts,tmpart)


c           which is the correct PRN (sat id)?
            is = 1
            do 16 j = 1,nsat
               if (iprn .eq. isprn(is)) then
                  goto 17
               else
                  is = is + 1
               endif
 16         continue
 17         continue
            if (is .gt. nsat) then
               write(message,'(a,i3)') 'PRN not found: ',iprn
               call report_stat('FATAL',prog_name,'readc',' ',message,0)
            endif

            if (lgood(ierfl)) tag(ii,isite) = dt

c           update start and end pointers
            if (lgood(ierfl) .or. lmarg(ierfl) .or. lbias(ierfl)) then
               ii2 = ii
               kk0(is,isite) = min02(kk0(is,isite),ii2)
               kk1(is,isite) = max02(kk1(is,isite),ii2)
            endif

c           error flag and bias pushing
            if ( elvcut(isite) .gt. elvcut_cfiles(isite) ) then
              if ( elev .lt. elvcut(isite)*pi/180.d0 ) then
                if ( lbias(ierfl) ) then
                  pushbias(is) = .TRUE.
                  ierr(ii,is,isite) = 4
                elseif ( lgood(ierfl) ) then
                  ierr(ii,is,isite) = 4
                else
                  ierr(ii,is,isite) = ierfl
                endif
              else
               if (pushbias(is).and.(lgood(ierfl).or.lbias(ierfl))) then
                  ierr(ii,is,isite) = 10
                  pushbias(is) = .FALSE.
                else
                  ierr(ii,is,isite) = ierfl
                endif
              endif
            else
              ierr(ii,is,isite) = ierfl
            endif

            ela (ii,is,isite) = sngl(elev)
            aza (ii,is,isite) = sngl(azim)
            clk (ii,is,isite) = rclock
 
c    compute the atmospheric gradient partials
            call atm_grad_part(elev,azim,grad_parts)

            if (imode. eq. 1) then
c              postfit mode: add adjustment
               domc=addadj( is,dt,iepoch,norb,nparts,tmpart
     .                    , coord_index,atm_index1,atm_index2
     .                    , clock_index,orb_index(is),svant_index(is)
     .                    , eop_index,grad_index,grad_parts ) 
c               if( iepoch.eq.78 ) then 
c                 print *,'clock_index is orb_index(is),svant_index(is) '
c     .                  , clock_index,is,orb_index(is),svant_index(is) 
c                 print *,'domc ',domc
c               endif
            else
               domc = 0.0d0
            endif

            if (lambda(is,1) .ne. 0) then
               yl1(ii,is,isite)=omcs(1)-domc
            else
               yl1(ii,is,isite)=0.0d0
            endif 
            if (lambda(is,2) .ne. 0) then
crwk 151230  yl2(ii,is,isite)=omcs(2)-domc*gear
             yl2(ii,is,isite)=omcs(2)-domc*fL2(is)/fL1(is) 
            else
               yl2(ii,is,isite)=0.0d0
            endif
            if (lambda(is,3) .ne. 0) then
ckf            pr1(ii,is,isite)=omcs(3)
               pr1(ii,is,isite)=omcs(3)-domc
               r1 = (obsv(3) - omcs(3))
            else
               pr1(ii,is,isite)=0.0d0
               r1 = 0.0d0
            endif
            if (lambda(is,4) .ne. 0) then
ckf            pr2(ii,is,isite)=omcs(4)
crwk 151230    pr2(ii,is,isite)=omcs(4)-domc*gear
               pr2(ii,is,isite)=omcs(4)-domc*fL2(is)/fL1(is)
               r2 = (obsv(4) - omcs(4))
            else
               pr2(ii,is,isite)=0.0d0
               r2 = 0.0d0
            endif

c              don't bother with the multipath calculation
ckf            if (lambda(is,1) .ne. 0 .and. lambda(is,3) .ne. 0
ckf     .          .and. lgood(ierfl)) then
ckf               rd1(ii,is,isite)
ckf     .            =sngl(multip(elev,r1,offsl1(1),clight/frq1))
ckf            else
ckf               rd1(ii,is,isite) = 0.0d0
ckf            endif
ckf            if (lambda(is,2) .ne. 0 .and. lambda(is,4) .ne. 0
ckf     .          .and. lgood(ierfl)) then
ckf               rd2(ii,is,isite)
ckf     .            =sngl(multip(elev,r2,offsl2(1),clight/frq2))
ckf            else
ckf               rd2(ii,is,isite) = 0.0d0
ckf            endif


CD             if (mod(ii,100).eq.0.0) then
CD             write (6,*) 'READC: data ',ii,is,isite,
CD     .            ierr(ii,is,isite),
CD     .            yl1 (ii,is,isite),
CD     .            yl2 (ii,is,isite),
CD     .            pr1 (ii,is,isite),
CD     .            pr2 (ii,is,isite),
CD     .            clk (ii,is,isite)
CD             endif
   20    continue
   50 continue


      write (6,'(//,5(1x,a7))') 'Channel','PRN','First','Last','Ndata'
      do 305 i = 1,nsat
        write (6,'(5(1x,i7))')  i,isprn(i),kk0(i,isite),kk1(i,isite),
     .                          max(int(kk1(i,isite)-kk0(i,isite)),0)
 305  continue

      close (lunit)
  100 continue
      if (ioerr .ne. 0) then
         call report_stat('FATAL',prog_name,'readc',filnam,
     .   'Error, reading C-file: ',ioerr)
      endif

CD     do 295 i=1,ndats
CD        write(6,203) (lambds(isite,j,i),j=1,nsat)
CD 295 continue

  300 return

      end

