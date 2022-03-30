Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994. All rights reserved.
 
      Subroutine wrthed 

       implicit none

      include '../includes/dimpar.h'     
      include '../includes/units.h'  
      include '../includes/global.h'
      include '../../libraries/includes/const_param.h'
      include '../includes/model.h'

c     ***internal to this routine ***
      integer*4    idms(maxlab),islot(maxlab),idoy,nclock
     .          ,  iy,im,id,ihr,min,i,j 
      character*1 skd 
     
       real*8 clock(4),sec   


      character*8 otidemod1
      character*32 sitnamt    
      character*256 message   


      data preval/maxprm*0.0d0/   

c**rwk 190826: Temporary to avoid changing the c-file format for 10.71, write 
c              only a single-frequency value  for the SV antenna PCO 
C MOD TAH 200126: Code below not needed and commented out.
C     real*8 svantdxx(3,maxsat) 
C     do j=1,maxsat
C       do i=1,3
C         svantdxx(i,j) = svantdx(i,1,j)
C       enddo
C    enddo \
C MOD TAH 200126: Added range L1 to L2 conversion factors
*     pcf1 =  fL1**2/(fL1**2-fL2**2)
*     pcf2 = -fL2**2/(fL1**2-fL2**2)   
      real*8 pcf1, pcf2


c     In CFMRG and SOLVE, the parameters are given unique id #s (islot values)
c     according to the parameter type and the number of stations, satellites,
c     zenith delays, biases, etc., as documented in subroutine cfmrg/fills1.f.
c     For purposes of identification of parameter types on the C-file, we 
c     write into islot on record 3 a truncated version of these id #s, derived
c     simply by using the first Id # of the group (i.e., for station 1 or satellite 1).
c     Thus, for example, station latitude (islot = 1-100) becomes 1, orbital and
c     satellite direct radiation pressure (islot = 1101-1200) becomes 1101.
c     This scheme is slightly different from that original used in GAMIT, but
c     we are free to change now because the C-file islot numbers have not been 
c     recently used by either CFMRG or SOLVE; rather the order of partials on the C-file
c     has been assumed and the correct partials selected by counting stations, satellites,
c     and the number of orbital parameters.  From henceforth, however, CFMRG will 
c     check the C-file islot values in determining how to set up the M-file islots
c     for SOLVE.   rwk 980904
                                        

      do i=1,maxlab
        islot(i) = 0
        rlabel(i) = ' '
      enddo
           
c----- Assign labels, slots, and apriori values for station parameters

      if( npart.lt.5 ) then
        call report_stat('FATAL','MODEL','wrthed',' ',
     1 'npart less than 5, something is wrong with npart' ,0)
      else
        nlabel = 5
c         station latitude
        islot(1) = 1
        rlabel( 1) =   'SITE GEOC LAT  dms  '     
        preval(1) = latr_sph
c         station longitude
        islot(2) =101
        rlabel( 2) =   'SITE GEOC LONG dms  ' 
        preval(2) = lonr
c         station radius
        islot(3) =201
        rlabel( 3) =   'SITE RADIUS    km   '
        preval(3) = radius/1.d3
        idms(1) = 1
        idms(2) = 1
        do  i = 3,maxlab
         idms(i) = 0
        enddo     
c         atmosphere
        islot(4) =301
        rlabel( 4) =   'SITE ATMZEN  m      ' 
c       change units from light-seconds to meters
        preval(4)=zendel0*vel_light
c         station clock epoch  
        islot(5) =401
        rlabel( 5) =   'SITE CLOCK-# EP sec '  
c       station clock - only epoch usually estimated; model values separately on C-file
c                       these previously in slots 4-6
        preval(5) = clkepc   
c       a priori values for rate and acceleration, not estimated no longer stored 
c       in preval slots 6 & 7
        nclock = 4
        clock(1) = clkepc
        clock(2) = clkrat
        clock(3) = clkacc
        clock(4) = clkcub 
      endif      


c ---- Assign labels, slots, and apriori values for integrated orbital parameters
               
c     Current restriction:  SV antenna offsets and EOPs included among parameters only if
c                           orbits are also.  

c     See comments in setup.f on the definition of norbpart
c     Minimum number of parameters is 3 coordinates, zenith delay, and clock 
      nparam = 5
c     Add orbit parameters 
      if( norbpart.gt.0 ) then 
        do j=1,nchan
          do i=1,norbpart
            nparam = nparam + 1
            preval(nparam) = saticst(i,j)
          enddo
        enddo     
        if( icsnamt(1).ne.'    ' .and. icsnamt(1).ne.'X   ' ) then
          write(message,'(a,15(a4,1x))') 'Unknown orbital IC type'
     .                  ,(icsnamt(i),i=1,norbpart)
          call report_stat('WARNING','MODEL','wrthed',' ',message,0)
        endif   
        nlabel = 14
        rlabel(6) = 'ORBIT X    km       '
        rlabel(7) = 'ORBIT Y    km       '
        rlabel(8) = 'ORBIT Z    km       '
        rlabel(9) = 'ORBIT Xdot km/s     '
        rlabel(10)= 'ORBIT Ydot km/s     '
        rlabel(11)= 'ORBIT Zdot km/s     '
        islot(6)  =501
        islot(7)  =601
        islot(8)  =701
        islot(9)  =801
        islot(10) =901
        islot(11) =1001
c          ECOM1, ECOM2, ECOMC, and UCLR models all use at least 9 parameters
        if( srpmod.eq.'ECOM1'.or.srpmod.eq.'BERNE'.or.
     .      srpmod.eq.'ECOM2'.or.srpmod.eq.'ECOMC'.or.
     .      srpmod.eq.'UCLR1'.or.srpmod.eq.'UCLR2' ) then     
          nlabel = 18 
          islot(12) =1101
          islot(13) =1201
          islot(14) =1301  
          islot(15) = 1401
          islot(16) = 1501
          islot(17) = 1601
          islot(18) = 1701
          islot(19) = 1801
          islot(20) = 1901   
          write(rlabel(12),'("RAD PRES DIRECT     ")')
          write(rlabel(13),'("Y AXIS BIAS         ")')
          write(rlabel(14),'("B AXIS BIAS         ")') 
          write(rlabel(15),'("COS U DIRECT        ")')
          write(rlabel(16),'("SIN U DIRECT        ")')
          write(rlabel(17),'("COS U Y AXIS        ")')
          write(rlabel(18),'("SIN U Y AXIS        ")')
          write(rlabel(19),'("COS U B AXIS        ")')
          write(rlabel(20),'("SIN U B AXIS        ")') 
          islot(15) = 1401
          islot(16) = 1501
          islot(17) = 1601
          islot(18) = 1701
          islot(19) = 1801
          islot(20) = 1901   
          if( srpmod.eq.'ECOM2'.or.srpmod.eq.'ECOMC') then 
            write(rlabel(21),'("COS 2U DIRECT       ")')
            write(rlabel(22),'("SIN 2U DIRECT       ")')
            write(rlabel(23),'("COS 4U DIRECT       ")')
            write(rlabel(24),'("SIN 4U DIRECT       ")')
            islot(21) = 2001
            islot(22) = 2101
            islot(23) = 2201 
            islot(24) = 2301 
          endif
        endif  

c----- Assign labels, slots, and  and apriori values for SV antenna offsets

        nlabel = 5+norbpart+3
        if( nlabel.gt.maxlab ) then
          write(message,'(a,i2,a,i2,a)') 'Number of labels (',nlabel
     .        ,') exceeds maxlab (',maxlab,')'
          call report_stat('FATAL','MODEL','wrthed',' ',message,0)
        endif  
        rlabel(5+norbpart+1) = 'SVANT X AXIS        '
        rlabel(5+norbpart+2) = 'SVANT Y AXIS        '
        rlabel(5+norbpart+3) = 'SVANT Z AXIS        '
        islot(5+norbpart+1) = 2501
        islot(5+norbpart+2) = 2601
        islot(5+norbpart+3) = 2701  
        do j=1,nchan
          do i=1,3
            nparam = nparam + 1
* MOD TAH 200126: Change preval to LC combination since only LC update
*           can be estimated.  Compute J dependent for Glonass incase
*           pcf1 =  fL1**2/(fL1**2-fL2**2)
*           pcf2 = -fL2**2/(fL1**2-fL2**2)   
*           L1/L2 values are releases
C           preval(nparam) = svantdxx(i,j)
            pcf1 =  fL1(j)**2/(fL1(j)**2-fL2(j)**2)
            pcf2 = -fL2(j)**2/(fL1(j)**2-fL2(j)**2)
            preval(nparam) = pcf1*svantdx(i,1,j) + pcf2*svantdx(i,2,j) 
          enddo
        enddo

c----- Assign labels,  slots, and apriori values for EOP parameters

        nlabel = 5+norbpart+3+6
        if( nlabel.gt.maxlab ) then
          write(message,'(a,i2,a,i2,a)') 'Number of labels (',nlabel
     .        ,') exceeds maxlab (',maxlab,')'
          call report_stat('FATAL','MODEL','wrthed',' ',message,0)
        endif  
        rlabel(5+norbpart+3+1) =   'X POLE       arcs   '
        rlabel(5+norbpart+3+2) =   'X POLE RATE  arc/d  '
        rlabel(5+norbpart+3+3) =   'Y POLE       arcs   '
        rlabel(5+norbpart+3+4) =   'Y POLE RATE  arcs/d '
        rlabel(5+norbpart+3+5) =   'UT1-TAI      sec    '
        rlabel(5+norbpart+3+6) =   'UT1-TAI RATE sec/d  '
        islot(5+norbpart+3+1) = 8001
        islot(5+norbpart+3+2) = 8002
        islot(5+norbpart+3+3) = 8003
        islot(5+norbpart+3+4) = 8004
        islot(5+norbpart+3+5) = 8005
        islot(5+norbpart+3+6) = 8006  
        preval(nparam+1) = xp
        preval(nparam+2) = xpdot
        preval(nparam+3) = yp
        preval(nparam+4) = ypdot
        preval(nparam+5) = ut1
        preval(nparam+6) = ut1dot
        nparam = nparam + 6
c     end check for orbit parameters
      endif

c        Set remaining dimensions and dummy variables

      sitnamt(1:16) = sitnam
      do i = 17, 32
         sitnamt(i:i) = ' '
      enddo
                            
c     Ocean loading is always applied in a CoM system (if input
c     model is CE, the CM correction (CMC) is applied, so remove
c     the 'E' (or 'M') designator in the 8th character

      otidemod1 = otidemod
      otidemod1(8:8) = ' ' 

c         Convert the times

c     write the output times in GPST
      mtime = 2
      call dayjul(jd0,iy,idoy)
      call monday(idoy,im,id,iy)
      call ds2hms(iy,idoy,t0,ihr,min,sec)

c     ***write the c-file headers***

      call writc1 (iuc
     .,            ntext,text)
           
      call writc2 (iuc
     .,            sitnamt,rctype,rcvrsn,rcvrswx,swver,anttyp,antsn
     .,            npart,norbpart,gnss,nchan,ischan,fL1,fL2
     .,            ndat,dattyp,lambda
     .,            skd,nepoch,inter,ircint,mtime,isessn
     .,            iy,im,id,ihr,min,sec
     .,            offarp,offl1,offl2, antdaz, svantdx
     .,            obfiln,tfiln,jfiln
     .,            frame,precmod,nutmod,gravmod,srpmod 
     .,            eradmod,antradmod
     .,            ietide,isptide,speopmod
     .,            etidemod,otidemod,atidemod,atmlmod,hydrolmod  
     .,            atmlavg,hydrolavg,dryzen,wetzen,drymap,wetmap  
     .,             ionsrc,magfield
     .,            antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,            elvcut,nclock,clock
     .,            jdet,tet,jdmidpt,tmidpt
     .,            ut1,ut1dot,xp,xpdot,yp,ypdot
     .,            psi,psidot,eps,epsdot
     .,            avlmet
     .,            nslip,islip,islpst 
     .,            niextra,iextra,nrextra,rextra,ncextra,cextra)

      call writc3 (iuc,nlabel,nparam,
     .             islot,idms,rlabel,preval)



C        Print the C-File header information

          write(iprnt,47)
47    format(/,'---------------------------------------------------',/,
     1       1x,'** OUTPUT C-FILE HEADER INFORMATION **',/)
      WRITE(IPRNT,'(A80)') (TEXT(I),I=1,NTEXT)
      WRITE(IPRNT,'(//,1X,A,//)') 'Record 2 :'
      WRITE(IPRNT,'(1X,2A,/)') 'Site: ',SITNAM
      write(iprnt,'(1x,a,2(1x,a20))') 'Rcvr, SN: ',rctype,rcvrsn
      write(iprnt,'(1x,a,2(1x,a20))') 'Ant,  SN: ',anttyp,antsn
      WRITE(IPRNT,'(1X,A,1x,a3,f6.2/)') 'Rcvr, SW: ',rcvrswx,swver  
      write(iprnt,'(1x,a,1x,2i4,2x,a1,i4)') 'npart norbpart gnss nchan:'
     .           ,npart,norbpart,gnss,nchan 
      WRITE(IPRNT,'(1X,a,3f8.4,A,3F8.4,A,3F8.4,A,/)')
     1     'Antenna offsets (Up North East) ARP:',offarp,'  L1:',offl1
     2,    '  L2:',offl2,'  meters'
      WRITE(IPRNT,'(1X,A,a,3x,a,I2,3X,A,I4,3X,A,I4,3x,a,i4,/)')
     .    'gnss=',gnss,'mtime=',mtime,'nepoch=',nepoch,
     .    'inter=',inter,'ircint=',ircint
      WRITE(IPRNT,'(1X,A,I5,2I3,2X,2I3,F7.3,a,i2,/)')
     1     'Start time:',IY,ID,IM,IHR,MIN,SEC,' (GPST)  Session:',isessn
      write(iprnt,'(1x,a,10x,a,6i3)') ' Chan  PRN','Data type: '
     .                   ,(dattyp(i),i=1,ndat)
      write(iprnt,'(13x,a,2x,a,3x,a)') 'L1 Freq',' L2 Freq'
     .                  ,'Lambda'
      do i=1,nchan
        write(iprnt,'(4x,i2,2x,a1,i2,1x,2f10.3,6i3)')
     .                    i,gnss,ischan(i),fl1(i)*1.d-6,fL2(i)*1.d-6
     .                 , (lambda(i,j),j=1,ndat)
      enddo
      WRITE(IPRNT,'(/,1X,a,a16,a,a16,a,a16)')
     1   'Raw data file: ',obfiln,' T-file: ',tfiln,'  J-file: ',jfiln
      WRITE(IPRNT,110) nslip,(islip(i),islpst(i),i=1,nslip)
110   FORMAT(/,1X,'Bias flags (',i4,') : ',14i4,/,24x,14i4)
      write(iprnt,'(/,1x,a)') 
     .   'Models:  Frame Preces  Nutat  Grav  RadPres ERad AntRad'  
      write(iprnt,'(10x,7(a5,2x))') frame,precmod,nutmod,gravmod,srpmod
     .   ,eradmod,antradmod
      write(iprnt,'(/,9x,a)')
     .  ' SP-EOP   E-tides   O-load      Atm-load  Atm-tide Hydrol-load'
      write(iprnt,'(10x,6(a8,2x))')
     .            speopmod,etidemod,otidemod,atidemod,atmlmod,hydrolmod
      write(iprnt,'(/,2(1x,a,a4,1x,a4))') 
     .  'Zenith-delay models: ',dryzen,wetzen,
     .       '   Mapping functions: ',drymap,wetmap 
      write(iprnt,'(/,1x,a4,1x,a6)') 'Ionosphere models: '
     .         ,ionsrc,magfield
      write(iprnt,'(/,1x,a,a10,1x,a4)') 'Rcvr AntMod ',antmod_snx,antmod
* MOD TAH 200126: Removed output sice already written above.
C     write(iprnt,'(/,a,/,2a)') ' SV Antmod'
C    .             ,   '  PN AntBody                 AntMod           '
C    .             ,' Offsets (X Y Z)'
C     do i=1,nchan               
C       write(iprnt,'(1x,i3,1x,a20,a10,1x,a4,7f10.4)') ischan(i)
C    .   ,svantbody(i),svantmod_snx(i),svantmod(i),(svantdx(:,j,i),j=1,2)
C     enddo
      write(iprnt,114) jdet,tet
114   format(/,1x,'Ephemeris epoch (PEP JD, sec GPST) : ',i8,f12.3)    
      write(iprnt,115) jdmidpt,tmidpt,ut1,ut1dot,iuttyp
     .     ,xp,xpdot,yp,ypdot,psi,psidot,eps,epsdot
115   format(/,1x,'Earth rotation' 
     .    ,/,3x,'  Epoch (PEP JD, sec GPST)   :',i8,f12.3
     .    ,/,3x,'  UT1-TAI, UT1dot (s, s/d)   :',2f10.5,'  UT1 type:',i2
     .    ,/,3x,'  xp,  xpdot (arcs, arcs/d)  :',2f10.5
     .    ,/,3x,'  yp,  ypdot (arcs, arcs/d)  :',2f10.5
     .    ,/,3x,'  psi, psidot (arcs, arcs/d) :',2f10.5
     .    ,/,3x,'  eps, epsdot (arcs, arcs/d) :',2f10.5 )
      write(iprnt,120) avlmet
120   format(/,1x,'Met data (series binary flag) :',i3)
      write(iprnt,121) preval(4)
121   format(1x,'Zenith delay (m) :',f7.4)
      WRITE(IPRNT,'(3X,A,I2,3X,A,I2,3X,A,I2)')
     1  'npart=',npart,'norb=',norbpart,'nlabel=',nlabel     
      write(iprnt,'(1x,a,1x,20i4)') 
     .    'niextra iextra',(iextra(i),i=1,niextra)
      write(iprnt,'(1x,a,1x,20d12.3)') 
     .    'nrextra rextra',(rextra(i),i=1,nrextra)
      write(iprnt,'(1x,a,1x,20(1x,a8))') 
     .    'ncextra cextra',(cextra(i),i=1,ncextra)
      WRITE(IPRNT,'(//,1X,A,/)') 'Record 3 :'
      WRITE(IPRNT,150) (rlabel(I),islot(I),idms(I),I=1,NLABEL)
150    FORMAT(1X,'rlabel, islot, idms =',A20,I5,I3,33(/,26X,A20,I5,I3))
      WRITE(IPRNT,160) (PREVAL(I),I=1,NPARAM)
160    FORMAT(1X,'preval=',3D22.14,/,(12X,3D22.14))

      RETURN
      END
