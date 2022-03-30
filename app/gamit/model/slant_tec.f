      Subroutine slant_tec (zenith,radius,shellht1,shellht2,MSLMhght,
     .  alpha,VTEC,STEC,yearad,doy,IRC,debug)
C Subroutine to calculate slant TEC from the vertical TEC interpolated from the IONEX file
C In :	zenith		angle z between vertical and ray path/degrees
C	radius		average radius of earth/km, use 6371km
C	hght		height of single layer ionosphere/km (CODE use 450km)
C	alpha		constant for modified Single Layer Mapping Function
C	VTEC		vertical Total Electron Content interpolated from IONEX file/TECU
C Out:	STEC		slant TEC calculated using MSLM function/TECU
C	IRC		error integer

C MLSM Mapping function F(z) = 1/cos(z'), where sin(z') =R(sin(alpha*z)/(R+MLSM_H)
C Force declaration of variables
      implicit none

C Passed variables
      Real*8 zenith, radius,  alpha, VTEC, STEC,shellht1,shellht2,
     .       MSLMhght
      Integer*4 yearad,doy,IRC,debug

C Other variables
      Real*8 sinzd, Fzen, pi,radzenith,hght
      Integer*4 mapfn
C Set error variable & initial mapping function value
      IRC=1  
      mapfn = 0   
C Debug: Test values only!!
C      doy = 252
C      yearad=2001
C      VTEC=50
C      zenith=85
C      Print*, 'yearad',yearad,'doy',doy


C Convert zenith angle to radians
	pi=4D0*atan(1D0)
C     Debug
C      Print*,'MODEL\slant_tec pi is ',pi
            radzenith=zenith*pi/180D0
C Debug      
C      Print*,'MODEL\slant_tec zenith: ',zenith
C      Print*,'MODEL\slant_tec zenith/rads: ',radzenith
C------------------------------------------------------------------
C Decide if 1) 1/cos(z) or 2) MLSM mapping function should be used
C CODE used MLSM mapping fn to create the IONEX files starting on day 252 2001.
C MLSM Mapping function F(z) = 1/cos(z'), where sin(z') =R(sin(alpha*z)/(R+MLSM_H)
      If (yearad.gt.2001) then
         mapfn = 2
      Elseif (yearad.eq.2001 .AND. doy.gt.251) then
         mapfn = 2
      Else
         mapfn = 1
      Endif
C debug
      If (debug.ge.3) then
        Print*,'MODEL/slant_tec: mapfn',mapfn
      Endif
C-----------------------------------------------------------------
C Use correct mapping function to calculate Fzen
      If (mapfn.eq.1) then
C                Fzen=1.d0/cos(radzenith)  wrong mapping fn (from IONEX file)
C     Use the Single Layer Mapping function
C         Set height of iono shell appropriately
          If (yearad.lt.1998) then
               hght = shellht1
C            Print*,'pre,shellht1',shellht1
          Else if (yearad.eq.1998.and.doy.lt.087) then
               hght= shellht1
C            Print*,'justpre,shellht1',shellht1
          Else
               hght=shellht2
C             Print*,'post,shellht2',shellht2
          Endif
C         Calculate sin(z')
          sinzd=(radius/(radius+hght))*sin(radzenith)
C         Calculate F(z)
          Fzen=1.0D0/sqrt(1.0D0-(sinzd*sinzd))
C     Debug
            if(debug.ge.2) then
      Print*,'MODEL\slant_tec sin(z)`', sinzd,'radius',radius,
     .'hght',hght, '   without alpha'
            endif

      Elseif (mapfn.eq.2) then
C     Use the Modified Single Layer Mapping function
C         Set height of iono shell appropriately
          hght=MSLMhght
C        Calculate sin(z')
        sinzd=(radius/(radius+hght))*sin(alpha*radzenith)
C        Calculate F(z)
        Fzen=1.0D0/sqrt(1.0D0-(sinzd*sinzd))
C Debug
        if(debug.ge.2) then
      Print*,'MODEL\slant_tec sin(z)`', sinzd,'radius',radius,
     .'hght',hght,'alpha',alpha
        endif

      Else
        call report_stat('FATAL','MODEL','slant_tec',' ',
     .  'ionospheric mapping function problem',0)
      Endif
C------------------------------------------------------------------
C Calculate slant TEC
      STEC=VTEC*Fzen
C Debug
        if(debug.ge.2) then
       Print*,'MODEL\slant_tec: STEC',STEC,'VTEC', VTEC
        endif
C Attempt to check value is realistic
      If (STEC.le.(VTEC-1.0)) then
         IRC=0
         Print*,'MODEL\slant_tec  STEC is less than VTEC'
	Elseif (VTEC.le.0.00001) then
	IRC=0
	Print*, 'MODEL\slant_tec: Very small VTEC', VTEC
      Endif

      End
