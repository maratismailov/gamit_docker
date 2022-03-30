Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 2017.  All rights reserved.

      Subroutine read_otides

c     Read from the IERS fes2004C2.dat file coefficients through degree 
c     and order 12 (ARC hard-wire) for the prograde and retrograde terms 
c     for 18 waves.  Decode the entries by Doodson Number but save the 
c     wave name for readablity.  Units for fes2004C2.dat are 1.e-12.

c     R. King 13 March 2017 
 
      implicit none  

      include '../includes/dimpar.h'
      include '../includes/arc.h'

      integer*4 ioerr,ideg,iord,iarg,Darg(18,6),i,j,k
      logical eoh,eof
      real*8 cplus,splus,cminus,sminus
      character*120 header
* MOD TAH 200503: Added additional book keeping on model used
      character*16 otmodel  ! Name of model based on 'Ocean tide model line'
        
c prograde (+,p) and retrograde (-,m) Stokes harmonic coefficients 
c from the tables for each of 18 Doodson waves (o[zc][pm]; coefficients
c summed across all waves, evaluated in sbfn.f -  arc.h
c     (same degree and order as fixed harmonics)
c      common/otides/ozcp(18,11),ozsp(18,11),ozcm(18,11),ozsm(18,11)
c     .              otcp(18,77),otsp(18,77),otcm(18,77),otsm(18,77)
c     .            , czcotid(11),czsotid(11),ctcotid(77),ctsotid(77)
c     .            , otid_doodson(18) 
c      real*8 ozcp,ozsp,ozcm,ozsm,otcp,otsp,otcm,otsm
c     .       czcotid,czsotid,ctcotid,ctsotid
c      character*7 otid_doodson

c Using the mapping in egm08.f but adding sin zonals, we have
c  Zonals
c   deg 2-12 --> oz[cs][pm](1-11)
c  Tesserals
c  2, 1-2  ot[cs][pm](1-2)
c  3, 1-3  ot[cs][pm](3-5)
c  4, 1-3  ot[cs][pm](6-9)
c  5, 1-5  ot[cs][pm](10-14)
c  6, 1-6  ot[cs][pm](15-20)
c  7, 1-7  ot[cs][pm](21-27)
c  8, 1-8  ot[cs][pm](28-35)
c  9  1-9  ot[cs][pm](36-44)
c  10 1-10 ot[cs][pm](45-54)
c  11 1-11 ot[cs][pm](55-65)
c  12 1-12 ot[cs][pm](66-77)
c    k = i(i-1)/2 + j - 1

      character*3 waves(18),arg
      data waves /'Om1','Om2','Sa ','Ssa','Mm ','Mf ','Mtm','Msq'
     .           ,'Q1 ','O1 ','P1 ','K1 ','2N2','N2 ','M2 ','S2 '
     .           ,'K2 ','M4 '/
      character*7  Darga
      logical debug/.false./
        
c  Set the Doodson numbers (in common/otides) for the FES2004 model 
c    See doodson_angle.f for documentation of Doodson numbers and their
c    relationship to Brown's fundamental arguements

      otid_doodson(1)  = ' 55.565'
      otid_doodson(2)  = ' 55.575'
      otid_doodson(3)  = ' 56.554'
      otid_doodson(4)  = ' 57.555'
      otid_doodson(5)  = ' 65.455'
      otid_doodson(6)  = ' 75.555'
      otid_doodson(7)  = ' 85.455'
      otid_doodson(8)  = ' 93.555'
      otid_doodson(9)  = '135.655'
      otid_doodson(10) = '145.555'
      otid_doodson(11) = '163.555'
      otid_doodson(12) = '165.555'
      otid_doodson(13) = '235.755'
      otid_doodson(14) = '245.655'
      otid_doodson(15) = '255.555'
      otid_doodson(16) = '273.555'
      otid_doodson(17) = '275.555'
      otid_doodson(18) = '455.555'

      do i=1,18
        do k=1,nczon1
          ozcp(i,k) = 0.d0
          ozsp(i,k) = 0.d0
        enddo 
        do k=1,nctes2
          otcp(i,k) = 0.d0
          otsp(i,k) = 0.d0
          otcm(i,k) = 0.d0
          otsm(i,k) = 0.d0
        enddo
      enddo
                
c     skip the header comments
      eoh =.false.
      otmodel = 'UNKNOWN'
      do while (.not.eoh) 
        read(iotide,'(a)',iostat=ioerr) header
        if(ioerr.ne.0) call report_stat('FATAL','ARC','read_otides'
     .      ,otidfname,'Error reading header of FES2004 otides',ioerr)
        if(header(1:7).eq.'Doodson' ) eoh = .true.
* MOD TAH 200503: Save model name if found.  Use of * for read will
*       get first word only.
        if(header(1:17).eq.'Ocean tide model:' ) then
           read(header(19:),*) otmodel
        endif
      enddo       
      eof = .false.                             
      do while(.not.eof )                                                 
        read(iotide,'(a7,1x,a3,2i4,2f10.0,2x,2f10.0)',iostat=ioerr) 
     .       Darga,arg,ideg,iord,cplus,splus,cminus,sminus
C       write(*,'(a,1x,a,1x,2I4)') darga,arg,ideg,iord
        if(ioerr.eq.-1 ) then
           eof = .true.
        elseif(ioerr.ne.0 ) then
           call report_stat('FATAL','ARC','read_otides',otidfname, 
     .        'Error reading tidal waves for ' // trim(otmodel),ioerr)
        else     
c         get the argument index   
          iarg = 0    ! Initialize incase this line is not found.     
          do i=1,18
            if(arg.eq.waves(i)) iarg=i
          enddo                                     
c         put the coefficients into linear arrays and convert the units
* MOD TAH 200503: Only save values if the argument iarg is found.
*         (FES2014 has more tidal lines than FES2004).
          if( iarg.gt.0 ) then 
             if( iord.eq.0.and.ideg.gt.1.and.ideg.le.nczon1 ) then       
                ozcp(iarg,ideg-1) = cplus
                ozsp(iarg,ideg-1) = splus  
C               write(*,'(3i3,2f8.2)') ideg,iord,iarg, 
C    .                       ozcp(iarg,i),ozsp(iarg,i)
             elseif( iord.gt.0.and.iord.le.nctess.and.
     .               ideg.gt.1.and.ideg.le.nczon1 ) then
                k= ideg*(ideg-1)/2 + iord -1  
                otcp(iarg,k) = cplus
                otsp(iarg,k) = splus
                otcm(iarg,k) = cminus
                otsm(iarg,k) = sminus
C               write(*,'(4i3,4f8.2)') ideg,iord,iarg,k, 
C    .             otcp(iarg,k),otsp(iarg,k),otcm(iarg,k),otsm(iarg,k)             
             endif
           endif
        endif
      enddo       
                              
c DEBUG      
      if(debug) then
c      write out the coefficients by wave and harmonic degree/order
       write(*,'(a,1x,a)') 'Values from table',otmodel
       do iarg =1,18
         write(*,'(2x,a3,1x,a7)') waves(iarg),otid_doodson(iarg)
         do i = 1,nczon1 
           ideg = i+1                        
           iord = 0   
           write(*,'(2i3,2f8.2)') ideg,iord,ozcp(iarg,i),ozsp(iarg,i)
           do iord = 1,i+1   
             k= ideg*(ideg-1)/2 + iord -1 
             write(*,'(2i3,4f8.2)') ideg,iord,otcp(iarg,k),otsp(iarg,k)
     .                                       ,otcm(iarg,k),otsm(iarg,k)             
           enddo
         enddo
       enddo
      endif                 

c     Apply the units conversion (do this here rather than earlier so that we can
c     print out the coefficients in the original units for debugging)
      do iarg = 1,18        
        do i=1,nczon1 
          ozcp(iarg,i) = ozcp(iarg,i) * 1.d-11
          ozsp(iarg,i) = ozsp(iarg,i) * 1.d-11
        enddo
        do i = 1,nctes2   
          otcp(iarg,i) = otcp(iarg,i) * 1.d-11
          otsp(iarg,i) = otsp(iarg,i) * 1.d-11
          otcm(iarg,i) = otcm(iarg,i) * 1.d-11
          otsm(iarg,i) = otsm(iarg,i) * 1.d-11
        enddo
      enddo

* MOD TAH 200503: Report model found.
      call report_stat('STATUS','ARC','read_otides',otidfname, 
     .        'MODEL ' // trim(otmodel),0)

cd      print *,'leaving read_otides '
      return
      end 
  

                     
