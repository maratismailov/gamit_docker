       subroutine jreadc(lunit,filnam,nepoch,inter,nsat,isprn,
     .                   jd0,sod0,pl1,ttag,ierrf,
     .                   clkepc,clkrat,clkacc,ifil)
c
c     read a C-file called filnam on logical unit lunit
c
c     Input:
c            lunit   logical file unit
c            filnam file name
c     Output
c            nepoch total number of epochs
c            nsat   total number of sattelites
c            ndats  number of different observables
c                   e.g. 4 for L1,L2,P1 and P2
c            pl1    L1 phase (observed-calculated in cycles)
c            jd0    PEP JD (UTC) of initial epoch
c            sod0   second of UTC day of initial epoch
c            ttag   seconds after intial epoch (jd0,sod0)
c            ierrf  GAMIT error flag
c
c     cleaned up by K. Feigl May 89
c     roots:
c     readc - R.W. King
c
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/errflg.h'

c     the L1, O minus C phase data in cycles
      real*8 pl1 (maxepc,maxsat,maxsit)
c     the calculated delay in seconds
c**      real*8 tau (maxepc,maxsat,maxsit)
c     the accompaning error flag
      integer*4 ierrf (maxepc,maxsat,maxsit)
c     and the accompanying time tag
      real*8 ttag (maxepc)

c     Receiver clock polynomial coefficients: 0th,1st and 2nd order.
c     units should be s,s/s,1/s respectively
      real*8 clkepc,clkrat,clkacc

      integer maxpar
      parameter (maxpar=7+maxsat*maxorb)
c
      real*8 dt
      real*8 latr_sph,lonr,radius
      real*8 preval(maxpar),offarp(3),offsl1(3),offsl2(3)
      real*8 rextra(maxext),save(maxsav)
      real*8 obsv(maxdat),omcs(maxdat),tmpart(maxlab)
      real*8 azim,elev,zendel
      real*8 tau0,sod,sod0,sec,rclock
      real*8 svcepc,svcl1
      integer*4 inter,julday,jd0,jdoy,jdoy0
      integer*4 ierfl,nepoch,nsat
      integer*4 id,ii,jj,im,iy,msat,iprn,is
      integer ifil
      integer lambda(maxsat,maxdat)
      integer lunit,ntext,ihr,min,iyr,nlabel,i,nsave,nparts
      integer npart,ndat,norb
      integer nparam,ndats,j,mtime,iepoch
      integer niextra,iextra(maxext),nrextra,ncextra
      integer*4 ioerr
      integer isnr(maxdat),mprn(maxsat),dattyp(maxdat)
      integer jslot(maxlab),idms(maxlab)
      character*32 sitnam
      character*16 filnam,obfiln,tfiln
      character*20 rlabel(maxlab)
      character*80 text(maxtxt)
      character*256 message
      integer*4 isprn(maxsat)

c     variables for new c-file format oct. 90; May 95
      integer*4 okmet
      integer      avlmet,ircint,isessn
      integer*4    jde,jdr
      integer*4 ietide,isptide,nclock,data_flag
      real*4  ampl1,ampl2,atmlod(3)
      real*8  te,tr,ut1,ut1dot,xp,xpdot,yp,ypdot,psi,psidot,eps,epsdot
      real*8  pres,temp,relhum
      real*8  azmdot,elvdot
     .      , atmdel,drate
C MOD TAH 200126: Changed to svantdx(3,2,maxsat) (added L1/L2 index)
      real*8  atmlavg(3),hydrolavg(3),svantdx(3,2,maxsat)
      real*8 elvcut,clock(4)
      real*8  spare(maxspr)
      real*8  l1z,l1n,l1e,l2z,l2n,l2e,antaz
      real*8  fL1(maxsat),fL2(maxsat)
      character*3  rcvrsw
      character*4  antmod,svantmod,dryzen,wetzen,drymap,wetmap
     .          ,  sitecd
      character*5  frame,precmod,nutmod,gravmod,sprmod   
     .,            otidemod,atmtide,hydrolmod
      character*6  etidemod
      character*8  speopmod,atmlmod,cextra(maxext) 
      character*10 antmod_snx,svantmod_snx(maxsat)
      character*16 jfiln
      character*20 rctype,rcvnum,anttyp,antnum,svantbody(maxsat)
      character*32 sitnamt
      real*4 swver
      integer nslip,nspare
      integer*2    islip(maxcsb),islpst(maxcsb)
      character*1 skd,gnss

c      real*8 freql1

c     L1 frequency in Hz
c      DATA FREQL1/1.57542D9/
      dt = 0.d0

c     function to calculate multipath path delay
CKF   real*8 multip

      write (6,5) filnam
  5   format (1x,'Opening: ',a16)

      call copens (filnam,'OLD',lunit,ioerr)

      call readc1 (lunit,ntext,text)
      write (6,'(a80)') (text(i),i=1,ntext)

C MOD TAH 200126: Replaced out-dated read with latest version from 
C     readc.f
      call readc2 (lunit
     .,            SITNAMT,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,            npart,norb,gnss,nsat,isprn,fL1,fL2
     .,            ndat,dattyp,lambda
     .,            skd,nepoch,inter,ircint,mtime,isessn
     .,            iy,im,id,ihr,min,sec
     .,            offarp,offsl1,offsl2,svantdx
     .,            obfiln,tfiln,jfiln
     .,            frame,precmod,nutmod,gravmod,sprmod     
     .,            ietide,isptide,speopmod
     .,            speopmod,etidemod,otidemod,atmtide,atmlmod,hydrolmod     
     .,            atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap
     .,            antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,            elvcut,nclock,clock
     .,            jde,te,jdr,tr
     .,            ut1,ut1dot,xp,xpdot,yp,ypdot
     .,            psi,psidot,eps,epsdot
     .,            avlmet
     .,            nslip,islip,islpst     
     .,            niextra,iextra,nrextra,rextra,ncextra,cextra)

      SITNAM = SITNAMT(1:12)

      if (ndat .gt. 4) then
         write (message,117) ndats,4
  117    format (1x,'File contains more data types (',i5,') than ',
     .   'program dimensions allow (',i5,')')
         close (lunit)
         call report_stat('FATAL','MAKEJ','jreadc',filnam,message,0)
      else if (ndat .le. 2) then
         call report_stat('WARNING','MAKEJ','jreadc',' ',
     .   'File contains no pseudoranges ',0)
      endif

      if (nepoch .gt. maxepc) then
         write (message,110) nepoch,maxepc
  110    format ('File contains more epochs (',i5,') than ',
     .   'program dimensions allow (',i5,')')
         close (lunit)
         call report_stat('FATAL','MAKEJ','jreadc',filnam,message,0)
      endif

      if (nsat .gt. maxsat) then
         write (message,115) nsat,maxsat
  115    format ('File contains more sats (',i5,') than ',
     .   'program dimensions allow (',i5,')')
         close (lunit)
         call report_stat('FATAL','MAKEJ','jreadc',filnam,message,0)
      endif

      call readc3 (lunit,nlabel,nparam,jslot,idms,rlabel,preval)

      clkepc = preval(4)
      clkrat = preval(5)
      clkacc = preval(6)

      write(6,201) sitnam,nsat,nepoch
 201  format (1x,a12,4x,i2,' sats ',i5,' epochs')
      write(6,202) (isprn(j),j=1,nsat)
 202  format (1x,'PRNs   :',32(i2,1x,:))

c     set up intial epoch
      jd0 = julday (im,id,iy)
      sod0 = 3600.0d0*ihr + 60.0d0*min + sec
      jdoy0 = jd0 - julday(1,1,iy) + 1

c     read epoch record
      do 50 ii = 1,nepoch
         call readc4 (lunit
     .,               msat,mprn
     .,               iepoch,iyr,jdoy,sod,rclock
     .,               okmet,zendel,pres,temp,relhum,atmlod
     .,               sitecd,latr_sph,lonr,radius
     .,               l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,               nsave,save  )

c        seconds past initial epoch
         if (iepoch .ne. ii) then
            write(message,'(a,2(i5,1x))') 'Epoch counter mismatch '
     .             ,ii,iepoch
            call report_stat('FATAL','MAKEJ','jreadc',filnam,message,0)
         else
            if (jdoy .gt. 0) then
               dt = 86400.0d0*dble(jdoy-jdoy0) + sod - sod0
            else
               dt = dble(iepoch - 1) * dble(inter)
            endif
         endif

c        and add receiver clock correction
         ttag(ii) = dt - rclock

c        initialize the error flags to avoid problems
         do 13 jj=1,nsat
            ierrf(ii,jj,ifil) = ignone
 13      continue

c        read the data, one sat at a time
         do 20 jj=1,msat
            call readc5 (lunit
     .,                  iprn
     .,                  elev,elvdot
     .,                  azim,azmdot
     .,                  atmdel
     .,                  svcepc,svcl1
     .,                  tau0,drate
     .,                  ierfl,data_flag
     .,                  ndats
     .,                  obsv
     .,                  omcs
     .,                  isnr,ampl1,ampl2
     .,                  nspare
     .,                  spare
     .,                  nparts
     .,                  tmpart)

c           which is the correct PRN (sat id)?
            is = 0
            do 16 j = 1,maxsat
               if (iprn .eq. isprn(j)) then
                  is = j
               endif
 16         continue

            if (is .eq. 0) then  
              write(message,'(a,i3)') 'Error, PRN not found: ',iprn
              call report_stat('FATAL','MAKEJ','jreadc',' ',message,0)
            endif

            ierrf(ii,is,ifil) = ierfl

c           we want phase in Doppler convention
c           so that its time derivative is opposite
c           of that of the range. KF 910202
            if (lambda(is,1) .ne. 0) then
ckf 910202     pl1(ii,is,ifil) = omcs(1)
               pl1(ii,is,ifil) = -omcs(1)
c**               tau(ii,is,ifil) = tau0
            else
               pl1(ii,is,ifil)=0.0d0
c**               tau(ii,is,ifil) = 0.d0
            endif
   20    continue
   50 continue

      close (lunit)
      if (ioerr .ne. 0) then
         call report_stat('FATAL','MAKEJ','jreadc',filnam,
     .   'Error reading C-file: ',ioerr)
      endif

      return

      end

