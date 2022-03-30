      Subroutine FILPAR ( istat, jsat, ipcnt  )

c     Fill the partials submatrix using partials from the C-files
c     Current version written by P. de Jonge and modified for additional parameters by R. King
c     Sept 1998 - April 2019 

c     ------ -------- --------------------------------------------------
                   
C-File order   M-file order             Slot value
c----------    ----------------------   ----------
c  1            'SC GEOC. LAT. (DMS) ',     1-100
c  2            'SC LONG. (DMS)      ',  101-200
c  3            'SC RADIUS (KM)      ',  201-300
c  4            'SC ATMOS-n          ',  301-400  (average zenith delay)
c  5            'SC CLOCK-n EP SECS  ',  401-500  
c------------------------------------------------------------------------------------
c  6            'ORBIT ELEMENT 1 PNnn',  501-600
c  7            'ORBIT ELEMENT 2 PNnn',  601-700
c  8            'ORBIT ELEMENT 3 PNnn',  701-800
c  9            'ORBIT ELEMENT 4 PNnn',  801-900
c 10            'ORBIT ELEMENT 5 PNnn',  901-1000
c 11            'ORBIT ELEMENT 6 PNnn',  1001-1100
c 12            'RAD PRES DIRECT PNnn',  1101-1200
c 13            'Y AXIS BIAS     PNnn',  1201-1300
c 14            'B AXIS BIAS     PNnn',  1301-1400    (or X AXIS BIAS or Z AXIS BIAS)  
c-------------------------------------------------------------------------------------
c 15            'COS U DIRECT    PNnn',  1401-1500  |
c 16            'SIN U DIRECT    PNnn',  1501-1600  |
c 17            'COS U Y         PNnn',  1601-1700  |  optional (ECOM1/BERNE or ECOM2-extended models)
c 18            'SIN U Y         PNnn',  1701-1800  |
c 19            'COS U B         PNnn',  1801-1900  |
c 20            'SIN U B         PNnn',  1901-2000  |
c-------------------------------------------------------------------------------------   
c 21            'COS 2U DIRECT   PNnn',  2001-2100  |
c 22            'SIN 2U DIRECT   PNnn',  2101-2200  |  optional (ECOM2 model)
c 23            'COS 4U DIRECT   PNnn',  2201-2300  |                                       
c 24            'SIN 4U DIRECT   PNnn',  2301-2400  |

c-------------------------------------------------------------------------------------
c The following C-file order assumes ECOM1/BERNE or ECOM2 
c Subtract 4 for ECOM1; subtract 10 3-parameter models; subtract 23 if no orbits
c-------------------------------------------------------------------------------------
c 25            'SVANT X AXIS    PNnn',  2501-2600
c 26            'SVANT Y AXIS    PNnn',  2601-2700
c 27            'SVANT Z AXIS    PNnn',  2701-2800  
c------------------------------------------------------------------------------------
c 28            'PNnnmmmmk BIAS L1   ',  2901-7400 
c 29            'PNnnmmmmk BIAS L2-L1',  7401-11900
c------------------------------------------------------------------------------------
c                nn is the nn'th zenith delay or gradient for that site
c                s  is the session number               
c 30            'SC ATM ZEN  nn s    ', 21501-24000 Multiple zenith delays 
c 31            'SC N/S GRAD nn s    ', 24001-26500 Multiple gradient parameters
c 32            'SC E/W GRAD nn s    ', 26501-29000 Multiple gradient parameters
c-------------------------------------------------------------------------------------
c 33            'X POLE (ARCS)       ', 80001-80001
c 34            'X POLE RATE (ARCS/D)', 80002-80002
c 35            'Y POLE (ARCS)       ', 80003-80003
c 36            'Y POLE RATE (ARCS/D)', 80004-80004
c 37            'UT1-TAI (SEC)       ', 80005-80005
c 38            'UT1-TAI RATE (SEC/D)', 80006-80006
c---------------------------------------------------------------------------     

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      double precision
     +  tmpart
      integer 
     +  jzen, jgrad, istat, jsat, ipcnt
      integer 
     +   iparm, ki, ki1, kl, jjsat, keop, iold
c      integer - for debug -- add to calling arguments
c     +  ml
      character 
     +   kzen*3, kgrad*2, message*80
  
      logical debug/.false./
cd      integer*4 i

      iold=ipcnt+1 
      if(debug.and.iepoch.eq.2031) 
     .   print *,'FILPAR istat npartm ',istat,npartm(istat)
      do iparm = 1,npartm(istat)   
c        npartm is the number of partials, including biases, designated by the m-file

         ki = islot2(iparm,istat)
         ki1 = islot1(ki)
         kl = ki1/100+1

         if (kl.ge.6 .and. kl.le.27) then
            ! For integrated orbit partials and SV antenna offsets
            jjsat = Mod(ki1,100)
            if (jsat.eq.jjsat) then
               tmpart = tpart (kl, istat, jsat)
            else
               goto 10
            endif 

         else if (kl.ge.29 .and. kl.le.119) then
            ! bias parameter - skip (a lot)
            goto 10                        

         else if (kl.ge.215 .and. kl.le.240) then
            if (zenest) then
               ! For multiple zenith delay slots (25 times)
               ! Find which zenith delay it is
               kzen=rlabel(ki)(16:18)
cd               print *,'FILPAR kl ki rlabel kzen ',kl,ki,rlabel(ki),kzen
               read (kzen,'(i3)') jzen
               call PZENTH (istat, jsat, jzen,  tmpart)
            else
               tmpart = 0.d0
            endif 

         else if ( kl.ge.241 .and. kl.le.265 ) then
            if (gradest) then
               ! For multiple N/S gradient parameters
               ! Find which gradient delay it is
               kgrad=rlabel(ki)(17:18)
               read(kgrad,'(i2)') jgrad
               call PGRAD ( istat,jsat,jgrad,1,tmpart )
            else
               tmpart = 0.d0
            endif 

         else if ( kl.ge.266 .and. kl.le.290 ) then
            if (gradest) then
               ! For multiple E/W gradient parameters
               ! Find which gradient delay it is
               kgrad=rlabel(ki)(17:18)
               read(kgrad,'(i2)') jgrad
               call PGRAD ( istat,jsat,jgrad,2,tmpart )
            else
               tmpart = 0.d0
            endif 
              
         else if (kl.eq.801) then
            if (eoppar .and. npartc.ge.20) then
               ! For new earth rotation partials, since all are put 
               ! in 80001-80006 slots (5 times)

               ! npartc is the number of partial derivative types on the C-file,
               ! set in model/setup, written in model/cfout, and read in solve/getcda.
               ! npartc =  5 => lat lon rad atm stn_clk (no orbital or global parameters)
               !        = 20 => 5 +  9 orbital + 6 EOP
               !        = 26 => 5 + 15 orbital + 6 EOP

               keop = npartc - 6
               if (ki1.eq.80001) then
                  tmpart = tpart (keop+1, istat, jsat)
               else if (ki1.eq.80002) then
                  tmpart = tpart (keop+2, istat, jsat)
               else if (ki1.eq.80003) then
                  tmpart = tpart (keop+3, istat, jsat)
               else if (ki1.eq.80004) then
                  tmpart = tpart (keop+4, istat, jsat)
               else if (ki1.eq.80005) then
                  tmpart = tpart (keop+5, istat, jsat)
               else if (ki1.eq.80006) then
                  tmpart = tpart (keop+6, istat, jsat)
               endif
            else
               tmpart = 0.d0
            endif

         else if (kl.lt.6) then
            ! coordinates and the like (4 times)
            tmpart = tpart (kl, istat, jsat)

         else  
           write(message,'(a,4i6)') 'Unknown parameter istat ki ki1 kl '
     .          ,istat,ki,ki1,kl
           call report_stat('FATAL','SOLVE','filpar',' ',message,0)
         endif 

         ipcnt = ipcnt+1
         ipntc(ipcnt) = ki
         c(ipcnt) = tmpart       
         if(debug.and.iepoch.eq.2031) 
     .    print *,'FILPAR istat jsat iparm ki ki1 kl ipcnt ipntc tmpart'
     .            , iparm,istat,jsat,ipcnt,ki,ki1,kl,ipntc(ipcnt),tmpart

 10      continue
      end do
    
cd      if(debug) then 
cd        write (*,*) ' FILPAR iold ipcnt ipntc ',iold,ipcnt
cd        write (*,'(5i16)') (ipntc(i),i=iold,ipcnt)
cd      write (*,'(5g16.6)') (c(i),i=iold,ipcnt)
cd        write (*,*) 
cd      endif 

      return
      end

