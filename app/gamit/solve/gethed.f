      Subroutine GETHED( lunit,icf )
c
c     Get header for the site from the C-file

c     Note:  The session number on C-file header records observation sessions
c            on the day; this is not used by SOLVE and is different from
c            the SOLVE session number for solutions

* MOD TAH 200219: Updated ietide generation from itides value for pole itide
*     secular/mean pole.
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'  
      include 'parameters.h'
      include 'models.h'
                    
      logical kbit
      character*1 skd,gnss
      character*3 rcvrsw 
      character*8 cextra(maxext)
      character*32 sitnm
      character*16 jfiln
      character*20 rctype,rcvnum,anttyp,antnum
      character*80 text(maxtxt)
      character*256 message
      integer*4          dattyp(maxdat),ntext
     .,                  ircint,iobssn
     .,                  lambda(maxsat,maxdat),nslip,avlmet
      integer*4          lunit,icf,lnpart,im,ir,lnsat,lndat
     .,                  llpoch,linter,mtime,iy,id,ihr,min,nclock
     .,                  ii,itides,mzen,mgrad
     .,                  niextra,iextra(maxext),nrextra,ncextra
      integer*2          islip(maxcsb),islpst(maxcsb)
      real*8             offarp(3),offsl1(3),offsl2(3),rextra(maxext)
     .,                  clock(4),elvcut_cfile(maxsit)
     .,                  elvcut_mfile(maxsit),sec
     .,                  fL1(maxsat),fL2(maxsat)
      real*8 antdaz  ! Antenna aligment from True N (deg).
      real*4  swver
          
      logical debug/.false./

c     initialise some character variables
      rctype = ' '
      rcvnum = ' '
      anttyp = ' '
      antnum = ' '
      rcvrsw = ' '
      swver = 0.d0
 
c     Read the C-file header
                  
        if(debug) print *,'GETHED reading C-file icf lunit ',icf,lunit
        call readc1(lunit,ntext,text)
c       save the MODEL version and run time for H-file
        im = index(text(ntext-1),'Model')
        ir = index(text(ntext),'Run')
        cfdate = ' '
        if( im.gt.0 )  cfdate(1:35) = text(ntext-1)(im:im+35)
        if( ir.gt.0 )  cfdate(36:75) = text(ntext)(ir:ir+40)  
        call readc2 (lunit
     .,         sitnm,rctype,rcvnum,rcvrsw,swver,anttyp,antnum
     .,         lnpart,norb
     .,         gnss,lnsat,isprn,fL1,fL2
     .,         lndat,dattyp,lambda
     .,         skd,llpoch,linter,ircint,mtime,iobssn
     .,         iy,im,id,ihr,min,sec
     .,         offarp,offsl1,offsl2, antdaz, svantdx
     .,         obfiln,tfiln,jfiln
     .,         frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
     .,         itides,isptide,speopmod
     .,         etidemod,otidemod,atmtide,atmlmod,hydrolmod  
     .,         atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap      
     .,         ionsrc,magfield
     .,         antmod_snx(icf),antmod(icf)
     .,         svantbody,svantmod_snx,svantmod
     .,         elvcut_cfile(icf),nclock,clock
c     .,         elvcutx,nclock,clock
     .,         jde,te,jdr,tr
     .,         ut1,ut1dot,xp,xpdot,yp,ypdot
     .,         psi,psidot,eps,epsdot
     .,         avlmet
     .,         nslip,islip,islpst  
     .,         niextra,iextra,nrextra,rextra,ncextra,cextra )  

cd        print *,'GETHED icf elvcut_cfile(icf) '
cd    .          ,        icf,elvcut_cfile(icf)
cd        print *,' iblk ',iblk
cd        print *,' svantmod_snx ',svantmod_snx
cd        print *,' svantmod     ',svantmod
cd         stop

c       Make sure character variables do not contain nulls
        if (ichar(rctype(1:1)) .eq. 0 ) rctype = ' '
        if (ichar(rcvnum(1:1)) .eq. 0 ) rcvnum = ' '
        if (ichar(rcvrsw(1:1)) .eq. 0 ) rcvrsw = ' '
        if (ichar(anttyp(1:1)) .eq. 0 ) anttyp = ' '
        if (ichar(antnum(1:1)) .eq. 0 ) antnum = ' '

c       MODEL version and run time
        im = index(text(ntext-1),'Model')
        ir = index(text(ntext),'Run')
        cfdate = ' '
        if( im.gt.0 )  cfdate(1:35) = text(ntext-1)(im:im+32)
        if( ir.gt.0 )  cfdate(36:75) = text(ntext)(ir:ir+40)

c        copy lambda array
         do ii=1,nsat
            lwave(icf,ii,1) = iabs(lambda(ii,1))
            lwave(icf,ii,2) = iabs(lambda(ii,2))
         enddo
                                 
c        write the frequency information to the q-file and set the L1/L2 ratio
c        (write only once, not for each c-file)
         if( icf.eq.1 ) then                 
           write(10,'(/,a,a1,a,2(1x,f8.3),a)') 
     .      'Observation frequencies for GNSS ',gnss,': ' 
           if( gnss.eq.'R') then  
             call report_stat('WARNING','SOLVE','gethed',' '
     .        ,'gear ratio not set correctly for Glonass',0)
             do ii=1,nsat
               write(10,'(i3,a,i4,2(1x,f8.3),a)') 
     .             ii,' PRN',isprn(ii),fL1(1)/1.d6,fL2(1)/1.d6,' MHz'
             enddo
* MOD TAH 180311: Added calculation of gear for GLONASS
             gear = fL2(1)/fL1(1) 
             write(10,'("Setting gear to ",F6.3," using frequencies ",
     .           2(f10.3)," Mhz")') gear, fL1(1)/1.d6,fL2(1)/1.d6

           else
             write(10,'(2x,2(f10.3),a)') fL1(1)/1.d6,fL2(1)/1.d6,' MHz'
             gear = fL2(1)/fL1(1) 
           endif  
         endif                          
        

c        copy receiver and antenna information and write to the q-file
            rcvr_type (icf) = rctype
            rcvr_sn   (icf) = rcvnum
            rcvr_sw   (icf) = rcvrsw
            rcvr_swver(icf) = swver
            ant_type  (icf) = anttyp
            ant_sn    (icf) = antnum
         do ii = 1,3
            anoff(icf,ii,1) = offarp(ii)
            anoff(icf,ii,2) = offsl1(ii)
            anoff(icf,ii,3) = offsl2(ii)
         enddo   
* MOD TAH 200205: Save the antenna azimith offset from c-file
         s_antdaz(icf) = antdaz

         if( icf.eq.1 ) then
           if( logprt ) write(6,'(/,a,7x,2a)') ' Station information'
     .        ,'Receiver            SwVer   Antenna'
     .        ,'             Ht to ARP'
           write(10,'(/,a,7x,2a)') ' Station information'
     .        ,'Receiver            SwVer   Antenna'
     .        ,'             Ht to ARP'
         endif
         if( logprt ) 
     .   write(6,'(i3,2x,a4,2x,a12,2x,a20,1x,f5.2,2x,a20,2x,f7.4)') 
     .        icf,cfiln(icf)(2:5),sitnam(icf)
     .       ,rctype,swver,anttyp,offarp(1)
         write(10,'(i3,2x,a4,2x,a12,2x,a20,1x,f5.2,2x,a20,2x,f7.4)')
     .        icf,cfiln(icf)(2:5),sitnam(icf)
     .       ,rctype,swver,anttyp,offarp(1)
 
c        store sampling intervals
         intsam(icf) = ircint

c        convert the 'tides applied' variable to descriptions of each type of tide
         ietide = 0           
         iotide = 0           
         iatide = 0
         if( kbit(itides,1) ) ietide = ietide + 1
         if( kbit(itides,2) ) ietide = ietide + 2
* MOD TAH 200219: Changed the pole-tide foe IERS20.  Work backwards through the
*        pole-tide birs to be constistent wit implemenation.
C        if( kbit(itides,3) ) ietide = ietide + 4
c*  This changed rwk/tah 110525 : TAH 200219: Removed: Done above.
c*         if( kbit(itides,5) ) ietide = ietide + 16       
C        if( kbit(itides,5) ) ietide = ietide + 8
*        Stick to the same scheme where only SE-type tides are recorder here 
*        The bit mapping is as below
*           Bit                                             Value
*             1 = solid earth tides                             1
*             2 = frequency dependant K1 model                  2
*             3 = Pole tide applied to zero mean pole           4
*             4 = Pole tide applied to IERS2010 mean pole       8
*             5 = Pole tide applied to IERS2020 mean pole      16
*   The other models are recorded in different variables.  In GLOBK
*   these get combined back together and in 32-bit gamit_mod variables 

         if( kbit(itides,7) ) then
            call sbit(ietide,5,1)
         elseif( kbit(itides,5) ) then
            call sbit(ietide,4,1)
         elseif( kbit(itides,3) ) then
            call sbit(ietide,3,1)
         endif

c        temporarily, at least: set separate bits for SP and LP o-tides
         if( kbit(itides,4) ) iotide = iotide + 3   
         if( kbit(itides,6) ) iatide = iatide + 1

c        replaces readc3
         read(lunit)

c---- Read merge-observation (c-file) header from the m-file  

c        (1st 2 records of m-file read in read_bfl)
         call readm3 (11,
     .                itor(1,icf),tor(1,icf),
     .                npartm(icf),islot2(1,icf),  
     .                elvcut_mfile(icf),zenmod,mzen,
     .                idtzen,gradmod,mgrad,idtgrad )  
                 
       if(debug) print *,'GETHED icf npartm elvcut_mfile zenmod mzen '
     .                   ,icf,npartm(icf),elvcut_mfile(icf),zenmod,mzen
                                               
c      Confirm that the number of zenith and gradient parameters on record 3 
c      of the m-file matches the number inferred from islot1 on record 1
         if(mzen.ne.nzen) then
           write(message,'(a,i2,a,i2,a)') 'mzen on m-file rec 3 ('
     .       ,mzen,' not = nzen from islot1 (',nzen,')'
            call report_stat('FATAL','SOLVE','gethed',' ',message,0)
         endif
         if(mgrad.ne.ngrad) then
           write(message,'(a,i2,a,i2,a)') 'mgrad on m-file rec 3 ('
     .       ,mgrad,' not = ngrad from islot1 (',ngrad,')'
            call report_stat('FATAL','SOLVE','gethed',' ',message,0)
         endif
c
c      There are now three different arrays of elevation cut angles:
c         elvcut_cfile(maxsit) -- values assigned in MODEL/AUTCLN for each site 
c         elvcut_mfile(maxsit) -- values assigned in CFMRG from cfile headers
c         elvcut_solve(maxsit) -- values assigned in SOLVE read from solve batch file
c
c       elevcut(maxsit) -- values assigned in SOLVE for each site (same over sessions).
c       elevcut(maxsit) is assigned the value that describes the elevation extent of data used
c       from each station in the solution. IE elevcut(i) cannot be lower than the low elevation
c       cutoff in the cfile (elvcut_cfile(i)), but it can be higher as defined by the solve
c       elevation angle cutoff (elvcut_solve(i)). The values in elevcut are written to
c       the mfile (in update.f) for use by CVIEW.  We've added the value to the C-file record (#3)
c       of the M-file but we do not yet make use of it in CVIEW.  Thus, now the C-file value and
c       the M-file values of elevation angle cutoff can be different, with the solution elevation
c       angle elevcut being written to the M-File
c
        if (elvcut_solve(icf) .gt. elvcut_cfile(icf)) then
           elevcut(icf) = elvcut_solve(icf)
        else
           elevcut(icf) = elvcut_cfile(icf)
        endif

c        write(6,'(/,1x,2i5)') nsat,npartm(icf)
c        write(6,'(12i5)') (islot2(ji,icf),ji=1,npartm(icf))

      return
      end



