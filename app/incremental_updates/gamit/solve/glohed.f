Copyright 1993 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

       SUBROUTINE GLOHED( numkey,oldrms,iparm )

c     Form the head block of the H-file
c     See lsqprt.f for documention of versions

* MOD TAH 200219: Changed ietide etc output format for IERS20 pole tide.

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'            
      include 'parameters.h'
      include 'models.h'

      integer *4 iold,isol,jsat,isite,ielem,ieop
     .         , julday,iparm(maxprm),numkey,lutfil,lunit,idoy4
     .         , iyr1,imo1,idy1,ihr1,min1,iyr2,imo2,idy2,ihr2,min2
     .         , iyr3,imo3,idy3,ihr3,min3,iyr4,imo4,idy4,ihr4,min4
     .         , itmp,idiag,ia,i1,ij4,ii,j0,j1,ib,kk,i,j,k  
              
      real*8 fac,ut1dum,utcoff,taiutc,tmpp,var,sec1,sec2,sec3,sec4
     .     , oldrms,temp(maxprm)

c     dummy
      integer*4 nsessn 

      character*1  symbl
      character*4  scode,stname(maxsit)
      character*5  tframe,tprecmod
      character*20 glab
      character*80 comt(3)
      character*82 header
      character*256 message    

*     tfile header info
      integer*4 ioerr
      integer ite(3),itb(3),itstp(3),nepchs
      real*8  tee(3),tbb(3),tstp(3),sdelt,xpole,ypole
      character*80 head

      fac = 2.0d0*95.146835d0

c        Read T-file to get ephemeris ICs  

      call lowers(tfiln)
      lutfil = lunit()
      open (lutfil, file=tfiln, form='unformatted', status='old'
     .     , err=5, iostat=ioerr)
5     if (ioerr.ne.0) then
         call report_stat('WARNING','SOLVE','glohed',tfiln
     .                   , 'Problem opening t-file',ioerr) 
         do i = 1, 3
            ite(i) = 0
            tee(i) = 0.d0
         enddo               
         tframe = '     '
      else
        read(lutfil) head,ite,tee,itb,tbb,itstp,tstp,sdelt,nepchs,
     .            ut1dum,xpole,ypole
        call check_y2k(ite(3))
        read(lutfil) comt
        tframe = '     '
        if( comt(1)(21:21).ne.'E' ) then
           if( index(comt(1)(30:40),'5').ne.0)  then
             tframe='B1950'
           elseif( index(comt(1)(30:40),'2').ne.0)  then
             tframe='J2000'
           elseif(comt(1)(35:40).eq.'     ' ) then
             tframe='B1950'
           else
             write(message,'(a,1x,a5)')
     .             'Unknown inertial frame in tfile header: '
     .             ,comt(1)(36:40)
             call report_stat('FATAL','SOLVE','glohed',' ',message,0)
           endif
           tprecmod = comt(1)(42:46)
          if(tprecmod.eq.'     '.and.tframe.eq.'B1950') tprecmod='IAU68'
          if(tprecmod.eq.'     '.and.tframe.eq.'J2000') tprecmod='IAU76'
        endif
      endif
      if( ioerr.eq.0 ) close(lutfil)

c     Check T-file inertial frame vs C-file from /modcom/ (read in GETHED)

      if( frame.eq.'     ' ) then
        frame = tframe
        precmod = tprecmod
      else
        if( frame.ne.tframe .or. precmod.ne.tprecmod ) then
           write(message,'(a,a5,1x,a5,a,a5,1x,a5,2x,a)')
     .     'Inertial frames from C-file (',frame,precmod
     .   , ') and T-file (',tframe,tprecmod,') do not match'
           call report_stat('WARNING','SOLVE','glohed',' ',message,0)
        endif
      endif
     

c       Continue writing H-file (first 6 lines written in LSQPRT and READ_BFL)
   
c     for ocean loading, 1=short-period, 2=long period 
c       iotide set from c-file 'itides' in gethed.f     
c     for atmospheric loading, 1=on
      if( atmlmod(2:2).eq.' ' ) then
        iaload = 0
      else
       iaload = 1
      endif
c     for atmopsheric tides, 1=S1/S2
c       iatide set from c-file 'itides' in gethed.f
      if( hydrolmod(1:1).eq.' ') then 
        ihload = 0
      else
        ihload = 1
      endif
* MOD TAH 2002019: changed format from 1x,i1 to i3 to accomate
*     large values but still keep columns aligend
C     write(21,'(1x,a,6(a,a8,1x,i2))') 
      write(21,'(1x,a,6(a,a8,   i3))') 
c      write(21,'(1x,a,a,a8,i2)')
     .       'Models:','  SP E-rot ',speopmod,isptide
     . ,'  SE-tide ',etidemod,ietide,'  O-load ',otidemod,iotide
     . ,'  Atm-load ',atmlmod,iaload,' Atm-tide ',atmtide,iatide
     . ,'  Hydrol-load ',hydrolmod,ihload
      write(21,'(1x,5(a,a4),a,a6)') 
     .    'Atm models:  DryZen ',dryzen,'  WetZen ',wetzen
     .    ,'  DryMap ',drymap,'  WetMap ',wetmap
     .    ,'  IonSrc ',ionsrc,'  MagFld ',magfield
      header(1:49)='       user  soln  diff  phase constraints biases'  
      header(50:82)='         parameters      h-file'
      write(21,'(a,/,a,15(1x,a5))')
     .    header,' keys:',(keyword(I),i=1,numkey)
c
      if ( err_mod(1) .eq. 'baseline' ) then
         write (21,10) aphi*fac,bphi*fac
 10      format (1x,'Assumed standard deviation of measurement error: ',
     .           f7.1,' mm  + ',f5.2,' ppm')
      else
         do i = 1,nsite
            i1 = (i-1)*3 +1
            stname(i) = rlabel(i1)(1:4)
         enddo
         write (21,12) (stname(i),i=1,nsite)
 12      format(1x,'Assumed stn obs errors  ',50(a4,1x))
         write (21,13) (err_mod(i),i=1,nsite )
 13      format(1x,'Assumed elev model:     ',50(a4,1x))
         write (21,14) (sit_err(i)*fac,i=1,nsite )
 14      format(1x,'Assumed con std dev mm: ',50(f4.1,1x))
         write (21,15) (sit_elv(i)*fac,i=1,nsite )
 15      format(1x,'Assumed elv std dev mm: ',50(f4.1,1x))
         write (21,16) (isprn(i),i=1,nsat )
 16      format(1x,'Assumed sat obs errors  ',50(i5,1x))
         write (21,17) (sat_err(i)*fac,i=1,nsat )
 17      format(1x,'Assumed std dev mm:     ',50(f5.1,1x))
      endif

      if (l2flag.ge.2)
     .write(21,20) akappa*fac,bkappa*fac
 20   format(1x,'Assumed ionosphere error: '
     1        ,f7.1,' mm  + ',f7.2,' ppm')
                    
c**  rwk 190510: Now always one session but keep the line for backward compatibility
      nsessn = 1 
      write (21,30) nsessn,nsite,nsat
 30   format(1x,'Number of sessions: ',i3,/,
     .       1x,'Number of stations: ',i3,/,
     .       1x,'Number of svs: ',i3)

      write (21,'(1x,a)') 'Ephemeris file: '
      write(21,'(1x,i3,1x,a16,1x,a82)') i,tfiln,head
       
      write (21,'(1x,3a)')
     .       'Name of track stations      Obs '  
     .,      '   Receiver type        Serial number       Rcvr swver'
     .,      '  Ant type       Dome  Serial number'
      do i = 1,nsite
         i1 = (i-1)*3 +1
         scode = rlabel(i1)(1:4)
         write (21,50) i,scode,sitnam(i),jusit(i)
     .               , rcvr_type(i),rcvr_sn(i),rcvr_sw(i),rcvr_swver(i)
     .               , ant_type(i),ant_sn(i)
      enddo
 50   format(1x,i3,2x,a4,2x,a12,2x,i7,3x,a20,1x,a20,1x,a3,1x,f5.2
     .      ,2x,a20,1x,a20) 
      write (21,'(1x,''Data files: '',2x,a75)') cfdate
* MOD TAH 200205: Add the azimuth offset for the antenna. AntDAZ
      write (21,'(6x,a,7x,a,9x,a,7x,a,7x,a,4x,a,a)')
     .     'code','C-file','Offset(U,N,E for ARP'
     .   , 'Offset(U,N,E for L1)','Offset(U,N,E for L2)'
     .   , 'AntDAZ   Ant mod        Elev cut  Num Zen '
     .   , 'Atm load (N E U mm)  Hydro load (N E U mm)'
      do i = 1,nsite
        scode = cfiln(i)(2:5)  
        call uppers(scode)
        write (21,60) i,scode,cfiln(i),((anoff(i,j,k),j = 1,3),k=1,3)
     .               ,s_antdaz(i), antmod_snx(i),antmod(i)
     .               ,elevcut(i),nzen
     .               ,(atmlavg(i,j),j=1,3),(hydrolavg(i,j),j=1,3)  
      enddo
* MOD TAH 200205: Added 1x,F5.1,1x, for AntDAZ column.
 60   format(1x,i3,2x,a4,2x,a16,3(3(2x,f7.4)),1x,F5.0,1x
     .  ,3x,a10,1x,a4,3x,f4.1,4x,i3,2x,3f7.2,1x,3f7.2)
                                     
      write (21,'(1x,a)') 'Satellites:'
* MOD TAH 200126: Mod to output L1/L2 values of satellite PCO (svantdx)
      write (21,70) 
  70  format('Chan  PRN  Blk                     Obs      AntMod   ',
     .       '      L1  DX    DY      DZ     L2  DX    DY     DZ') 

      do i=1,nsat  
        write (21,'(1x,i2,2x,i4,2x,a20,2x,i6,2x,a10,1x,a4,2(1x,3f8.4))') 
     .    i,isprn(i),svantbody(i),iseen(i),svantmod_snx(i),svantmod(i)
     .      ,(svantdx(:,j,i),j=1,2)
      enddo 
      iyr1 = itor(3,1) 
      call check_y2k(iyr1)
      imo1 = itor(1,1)
      idy1 = itor(2,1)
      ihr1 = int(tor(1,1))
      min1 = int(tor(2,1))
      sec1 = tor(3,1)
      utcoff = taiutc( julday(imo1,idy1,iyr1) ) - 19.d0
      iyr2 = itor(3,1) 
      call check_y2k(iyr2)
      imo2 = itor(1,1)
      idy2 = itor(2,1)
      ihr2 = int(tor(1,1))
      min2 = int(tor(2,1))
      sec2 = tor(3,1) + dble(inter)*dble(nepoch-1)
c
      tmpp = sec2/60.0d0
      itmp = int(tmpp)
      sec2 = sec2 - dble(itmp)*60.0d0
      min2 = min2 + itmp
      itmp = min2/60
      min2 = min2 - itmp*60
      ihr2 = ihr2 + itmp
      itmp = ihr2/24
      ihr2 = ihr2 - itmp*24
      idy2 = idy2 + itmp
c
      iyr3 = ite(3) 
      imo3 = ite(1)
      idy3 = ite(2)
      ihr3 = int(tee(1))
      min3 = int(tee(2))
      sec3 = tee(3)
c
      call dayjul( jdr,iyr4,idoy4 )
      call monday( idoy4,imo4,idy4,iyr4 )
      call ds2hms( iyr4,idoy4,tr,ihr4,min4,sec4 )
c
      write(21,'(a,i7,3x,a,i4,a,i5,4x,a,i5,a,i3)')
     .   ' Double-difference observations: ',nobs
     .,   'Epoch numbers ',istart,' to',iend
     .,   'interval: ',inter, ' s   decimation: ',idecim
c
      write(21,'(a,4i4,2x,i3,2x,f7.3,a,f4.1,a)')
     .   ' Start time: ', iyr1, imo1, idy1, ihr1, min1, sec1
     .,  '  GPST    (GPST-UTC=', utcoff, ')'
c
      write(21,'(a,4i4,2x,i3,2x,f7.3)')
     .   ' End time  : ', iyr2, imo2, idy2, ihr2, min2, sec2
c
      write(21,'(a,4i4,2x,i3,2x,f7.3,3x,7(a,a5,a,a5,a,a5))')
     .   ' ICs time  : ', iyr3, imo3, idy3, ihr3, min3, sec3
     .  ,' Frame: ',frame,'   Precession: ',precmod
     .  ,'  SRP model: ',srpmod,'  Nutation: ',nutmod
     .  ,'  Gravity model: ',gravmod
     .  ,'  Earth-radiation model: ',eradmod
     .  ,'  Antenna-radiation model: ',antradmod

c ----New block added Nov 1991

c ----Key word and time change 11 June 1992

c*    write(21,'(a,8x,i8,2x,f10.3)')
c*   .   ' IC epoch(day & second) : ', jde, te
      write(21,'(a)')
     .  ' Earth orientation data'

      write(21,'(a,4i4,2x,i3,2x,f7.3)')
     .   ' EOP time  : ', iyr4, imo4, idy4, ihr4, min4, sec4

c*    write(21,'(a,i8,2x,f10.3)')
c*   .   ' E-rotation epoch(day & second) : ', jdr, tr

      write(21,'(a,2f12.6)')
     .   ' UT1(sec) & rate(sec/day) :         ', ut1, ut1dot

      write(21,'(a,2f12.6)')
     .   ' X pole(asec) & rate(asec/day) :    ', xp, xpdot

      write(21,'(a,2f12.6)')
     .   ' Y pole(asec) & rate(asec/day) :    ', yp, ypdot

      write(21,'(a,2f12.6)')
     .   ' Delta-psi(asec) & rate(asec/day) : ', psi,psidot

      write(21,'(a,2f12.6)')
     .   ' Delta-eps(asec) & rate(asec/day) : ', eps,epsdot

c----end new block

      write(21,'(a,i5,3x,a,i5)')
     .   ' Total parameters: ', ntpart
     .,   'live parameters: ', nlive

      write(21,'(a,e12.5,4x,a,e12.5)')
     .   ' Prefit nrms: ',oldrms
     .,   'Postfit nrms:',sclerr

c
c output parameters of sites and orbits
c
      write(21,240)
 240  format (/,8x,'label (units)',17x,'a priori',16x,'adjustment')
      do 260 i = 1,lpart
         ia = islot1(i)
c             skip clocks, multiple zenith delays, and atmosphere gradient parameters
         if ( ia.ge.401.and.ia.le.500   .or.
     .        ia.ge.21501.and.ia.le.29000 ) goto 260
         glab = rlabel(i)
         if (ia .gt. 0 .and. ia.le. 200) then
            glab(16:18) = 'rad'
         endif
         k = free(i)
         symbl = ' '
         if (k.ne.0) symbl = '*'
         write(21,250)
     .     i,symbl,glab,preval(i),adjust(i)
 260  continue
 250  format (1x,i4,a1,a20,1x,f23.16,4x,d23.16)

c
c     print out sub-covariance matrix:
c           no biases, no atmospheric parameters, no clock
c
      write(21,'(1x,/,1x,a)') 'Covariance matrix: '
      iold = 0
      do 100 i = 1,msig-1
         idiag = i*(i+1)/2
         i1 = iparm(i)
         ia = islot1(i1)
c        skip clocks, multiple zenith delays, and atmosphere gradient parameters
         if ( ia.ge.401.and.ia.le.500    .or.
     .        ia.ge.21501.and.ia.le.29000 ) goto 95
         j0 = 0
         do 90 j = 1,i
            j1 = iparm(j)
            ib = islot1(j1)
c           skip clocks, multiple zenith delays, and atmosphere gradient parameters
            if ( ib.ge.401.and.ib.le.500    .or.
     .           ib.ge.21501.and.ib.le.29000 ) goto 90
            j0 = j0 + 1
            ij4 = iold + j
            temp(j0) = a(ij4)
 90      continue
         write(21,110) i1,(temp(j),j=1,j0)
 110     format (i4,'. ',500(5(1x,d23.16),:,/,6x))
 95      iold = idiag
 100  continue

c      print out a priori covariances (non-zero elements)

      write (21,'(1x,a)') 'A priori covariance matrix:'
c     get solution type (1 = tight; 2 = loose)
      if( keyword(12)(3:3).eq.'C' ) then
        isol = 1
      elseif( keyword(12)(3:3).eq.'L' ) then
        isol = 2
      else
        call report_stat('FATAL','SOLVE','glohed',' '
     .       ,'Keyword C/L for solution missing',0)
        return 
      endif
      do 200 i = 1,msig-1
        i1 = iparm(i)
        ia = islot1(i1)
c       station coordinates
        if( ia.le.300 ) then
           isite = mod(ia,100)
           if( ia.le.100 ) then
c            latitude
             if( isol.eq.1 ) var = stat_apr (isite,1)
             if( isol.eq.2 ) var = stat_apr2(isite,1)
c            convert from km to radians using radius
             var = ( var/preval(i1+2) )**2
           elseif( ia.gt.100.and.ia.le.200 ) then
c            longitude
             if( isol.eq.1 ) var = stat_apr (isite,2)
             if( isol.eq.2 ) var = stat_apr2(isite,2)
c            convert from km to radians using radius and lat
             var = ( var/preval(i1+1)/dcos(preval(i1-1)) )**2
           elseif( ia.gt.200.and.ia.le.300 ) then
c            radius
             if( isol.eq.1 ) var = stat_apr (isite,3)**2
             if( isol.eq.2 ) var = stat_apr2(isite,3)**2
           endif
           write(21,'(1x,i4,1x,i4,2x,d23.16)') i1,i1,var    


c       average zenith delay  
        elseif( ia.gt.300 .and.ia.le.400 ) then 
          isite = mod(ia,100)
          if( isol.eq.1 ) var = zen_apr (isite)**2
          if( isol.eq.2 ) var = zen_apr2(isite)**2
           write(21,'(1x,i4,1x,i4,2x,d23.16)') i1,i1,var    

c       satellite parameters
        elseif( ia.gt.500 .and.ia.le.2700 ) then
           jsat = mod(ia,100)
           ielem = ia/100 - 4
           if( ielem.le.6 .and. satwt_type.eq.'KEP' ) then
c             input was Keplerian, so Cartesian covariance is non-diagoonal
              kk = (ielem*ielem - ielem)/2
              do ii = 1,ielem
                 kk = kk + 1
                 if( isol.eq.1 ) var = covorbx (kk,jsat)
                 if( isol.eq.2 ) var = covorbx2(kk,jsat)
                 j1= i1 + ii - ielem
                 write(21,'(1x,i4,1x,i4,2x,d23.16)') i1,j1,var
               enddo
           else
              if( isol.eq.1 ) var = sat_apr (jsat,ielem)**2
              if( isol.eq.2 ) var = sat_apr2(jsat,ielem)**2
              write(21,'(1x,i4,1x,i4,2x,d23.16)') i1,i1,var
           endif
c          Earth-rotation parameters
        elseif( ia.gt.80000 .and.ia.le.80006 ) then
           ieop = ia - 80000
           if( isol.eq.1 ) var = eop_apr (ieop)**2
           if( isol.eq.2 ) var = eop_apr2(ieop)**2
           write(21,'(1x,i4,1x,i4,2x,d23.16)') i1,i1,var
        endif
 200  continue

      continue
      RETURN
      END
