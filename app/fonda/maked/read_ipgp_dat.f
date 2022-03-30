      subroutine read_ipgp_dat(ifil1,ifil2)
c
c     read IPGP data file and create FONDA format data file
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      integer ifil1,ifil2

      character*8 name1,name2
      character*30 note
      character*80 line
      integer icode1,icode2
      integer iyear,iday,month,id1,id2,idonold,idon,idon2
      integer i,j
      integer ide,imi,itp
      integer*4 julday
      integer in,nexp
      integer ierr
      logical old1,old2,ok
c
      small = 1.0d-5

c  Input format is like this, acording to IDON
c  IDON
c -1  Date: year month day (comment)
C- 1- Mark to Mark distance  ex:  1 11 12 1000.000   0.005  0.   0.
C- 2- Distance reduite ( D0 )
c- 3- Distance axe optique - Rep. ou cible
c            ex:  3 11 12 1000.000   0.005  1.52 1.43
c- 4- Denivellees repere-repere:  ex:  4 11 12   10.002   0.001  0.   0.
c- 5 ou 9 - Angles horizontaux
C            ex:  9 11 12    0.      0.0001 0.   0.
C                 5 11 13  155.0321  0.001  0.   0.
c: IDON=9 indique la direction devant erre prise comme refence,
c          mais doit etre la premire de la serie
C- 6- Angles zenitaux (sur rep., theod. ou cible): angle corrige 
c          de refraction a partir visees dir. et inv. ou calcul
C            ex:  6 11 13   98.1234  0.0002 1.45 0.05
c- 7 - Difference de coordonnees locales Est: dE
c- 8 - Difference de coordonnees locales Nord: dN
C- 12 - DX en coordonnees geocentriques (issu de GPS)
C- 13 - DY en coordonnees geocentriques (issu de GPS)
C- 14 - DZ en coordonnees geocentriques (issu de GPS)
c---------------------------------------------------------------------------       

      idonold = -1

c     First, sort out all working sites and experiment number.
c     site id
      note = '{nsit (number of sites)'
      write (ifil2, 30) nsit,note
      write (ifil2,'(a8)') (sname(i),i=1,nsit)
c
c     count number of experiments
      nexp = 0
c     skip header
      read (ifil1,'(a)',iostat=ierr) line
      do 50 while (ierr .eq. 0) 
         read (ifil1,*,iostat=ierr) idon
         if (idon .eq. 0) then
c           end of experiment
c           next line should contain a date
            read (ifil1,*,err=130,end=130,iostat=ierr) iyear,month,iday
            idonold = idon
         else if   ( idon.eq.1.or.
     .               idon.eq.2.or.
     .               idon.eq.3.or.
     .               idon.eq.5.or.
     .               idon.eq.9) then
            if (idon.ne.idonold) then
               in = 0
               ok = .true.
               do while (ok)
                  read (ifil1,*,end=200) idon2
                  in = in +1
                  if (idon2 .eq. 0) ok = .false.
                  if(idon2.eq.1.or.
     .               idon2.eq.2.or.
     .               idon2.eq.3.or.
     .               idon2.eq.5.or.
     .               idon2.eq.9) then
                     if ((idon2.ne.idon.and.idon2.ne.5)
     .                    .or.idon2.eq.9) then
                        ok = .false.
                     endif
                  endif
               enddo
 200           continue 
               do i = 1,in
                  backspace (ifil1)
               enddo
               if (idon.eq.1 .or. idon.eq.2 .or. idon.eq.3) then
                  nexp = nexp + 1
               else if (idon.eq.9) then
                  nexp = nexp + 1
               endif
            endif
            idonold = idon
         endif
 50   continue
 130  continue 
      rewind (ifil1)
      note = '{number of experiments            '
      print *,'Found nexp = ',nexp
      write (ifil2,30) nexp,note
      
c     try to get date from header
      read (ifil1,*,end=900,iostat=ierr) iyear,month,iday
      time1 = 1900.0d0 + julday(month,iday,iyear,1)/365.2422d0 

      idonold = -1
      do 90 while (ierr .eq. 0)
 85      read (ifil1,*,iostat=ierr)
     .        IDON,icode1,icode2,DATI,ERRI,HSTA,HVIS
         write (name1,'(i8.8)') icode1
         write (name2,'(i8.8)') icode2

         if (idon .eq. 0 .and. ierr .eq. 0) then
c           end of experiment
c           next line should contain a date
            read (ifil1,*,err=900,end=900,iostat=ierr) 
     .      iyear,month,iday
            time1 = 1900.0d0 + julday(month,iday,iyear,1)/365.2422d0 
            idonold = idon
         else if   ( idon.eq.1.or.
     .               idon.eq.2.or.
     .               idon.eq.3.or.
     .               idon.eq.5.or.
     .               idon.eq.9) then
            if (idon.ne.idonold) then
               ngp = 1
               in=0
               ok = .true.
               do while (ok)
                  in = in+1
                  read (ifil1,*,end=110) idon2
                  if (idon2 .eq. 0) ok = .false.
                  if(idon2.eq.1.or.
     .               idon2.eq.2.or.
     .               idon2.eq.3.or.
     .               idon2.eq.5.or.
     .               idon2.eq.9) then
                     if ((idon2.ne.idon.and.idon2.ne.5)
     .                    .or.idon2.eq.9) then
                        ok = .false.
                     else
                        ngp = ngp+1
                     endif
                  endif
               enddo
 110           continue 
               do i = 1,in
                  backspace (ifil1)
               enddo
               if (idon.eq.1 .or. idon.eq.2 .or. idon.eq.3) then
                  iexp = iexp + 1
                  note = '{exp. index, obs. type, number'
                  itp = 4
                  write (ifil2,140) iexp,itp,ngp,note
               else if (idon.eq.9) then
                  iexp = iexp + 1
                  note = '{exp. index, obs. type, number'
                  itp = 2
                  write (ifil2,140) iexp,itp,ngp,note
                  ngp = 1
               endif
            endif
            idonold = idon
         endif
   
         if   ( idon.eq.1.or.
     .          idon.eq.2.or.
     .          idon.eq.3.or.
     .          idon.eq.5.or.
     .          idon.eq.9) then

c           identify two sites
            old1 = .false.
            old2 = .false.
            do 70 j = 1,nsit
               if (name1.eq.sname(j)) then
                  old1 = .true.
                  id1 = j
               endif
               if (name2.eq.sname(j)) then
                  old2 = .true.
                  id2 = j
               endif
               if (old1.and.old2) goto 100
 70         continue
            if (.not.old1) then
               print *,'READ_IPGP_DAT: no coords for ', name1
            endif
            if (.not.old2) then
               print *,'READ_IPGP_DAT: no coords for ', name2
            endif

c           come here after finding sites
 100        continue
         endif

c        here is where we choose the data
         if (idon .eq. 1) then
c           mark to mark distance
            robs = dati
c           fonda uncertainty is in mm!
            err = erri * 1000.0d0
            write (ifil2,80) time1,robs,err,name1,name2
         else if (idon .eq. 2) then
            print *,'READ_IPG_DAT: skipping idon = ',idon
         else if (idon .eq. 3) then
c           instrument to instrument distance
c           normale au point A : P0
            CALL CNORM(X(id1),Y(id1),Z(id1),DN1,DN2,DN3)    
            XAP= X(id1)+HSTA*DN1
            YAP= Y(id1)+HSTA*DN2
            ZAP= Z(id1)+HSTA*DN3
            CALL CNORM(X(id2),Y(id2),Z(id2),DN1,DN2,DN3)    
            XPV= X(id2)+HVIS*DN1
            YPV= Y(id2)+HVIS*DN2
            ZPV= Z(id2)+HVIS*DN3
            DX=X(id2)-X(id1)
            DY=Y(id2)-Y(id1)
            DZ=Z(id2)-Z(id1)
            DC=dsqrt(DX**2+DY**2+DZ**2)
            DXX=XPV-XAP
            DYY=YPV-YAP
            DZZ=ZPV-ZAP
            DDC=dsqrt(DXX**2+DYY**2+DZZ**2)
            robs=dati*DDC/DC
c           fonda uncertainty is in mm!
            err = erri * 1000.d0
            write (ifil2,80) time1,robs,err,name1,name2
c           should match format number 80 in outdat.f
 80         format (1x,f10.4,f14.5,f14.8,2x,a,2x,a)
         else if (idon .eq. 5 .or. idon .eq. 9) then
c           first direction in a series
            if (idon .eq. 9) then
               dirzero = dati
c              print *,'Setting dirzero to ',dirzero
            endif
            dati = dati - dirzero
            if (dati .ge. 400.d0) dati = dati - 400.0d0
            if (dati .lt.   0.d0) dati = dati + 400.0d0
c           convert grads to rads
            robs = dati * pi/200.0d0
            if (robs .lt. 0) robs = robs + 2.0d0*pi
c           convert grads to seconds
            err =  erri * 0.9d0 * 3600.0d0
            call rtodms(robs,ide,imi,s1,1)
            write (ifil2,120) time1,ide,imi,s1,err,name1,name2
c           should match format number 120 in outdat.f
 120        format (1x,f10.4,2i5,f14.8,f9.4,2x,a,2x,a)
         else if (idon .eq. 0) then
            print *,'Found start of: ',iyear,month,iday
         else
            print *,'READ_IPGP_DAT: skipping idon = ',idon
         endif
  90  continue

  30  format (i5,25x,a30)
c     should match format number 40 in outdat.f
 140  format (1x,3i5,14x,a30)
 900  continue
      return
      end

c=============================================
      SUBROUTINE CNORM(XX,YY,ZZ,DN1,DN2,DN3)
C---------------------------------------------
c----Calcul des composantes du vecteur normal l'ellipsoide
c        au pt Po
c
      IMPLICIT REAL*8 (A-H,L,O-Z)

      include 'maked.fti'

c     JCR uses different constants than dong
ckf      A=CC/DSQRT(1.D0+EP2)
ckf      B=A*A/CC

c     give JCR his constants from Dong's
      a = radius
      b = radius * (1.0d0 - 1.0d0/finv)
      ep2 = (a*a - b*b)/(b*b)

c     print *,'CNORM: a,b,finv ',a,b,finv

c     Begin JCR's code.
      P= DSQRT(XX**2+YY**2)
      PB=B*P
      ZA=ZZ*A
      TETA= DATAN2(ZA,PB)
      LBD= DATAN2(YY,XX)
      ST=DSIN(TETA)
      CT=DCOS(TETA)
C
      FX=ZZ+EP2*B*ST*ST*ST
      FY=P-EP2*A*CT*CT*CT/(1.d0+EP2)
      FI= DATAN2(FX,FY)
c
      DN1=dcos(FI)*dcos(LBD)
      DN2=dcos(FI)*dsin(LBD)
      DN3=dsin(FI)
c
      RETURN
      END
