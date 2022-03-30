       Subroutine UCLSRPeditIIA (angb,along_x,along_y,along_z)
C EJP Jan 2011 Code adapted from a UCL model described below (edited from fortran 90 to fortran 77+, 
C duplicated calculations removed): 
! This program provides a simple example of how to use the currently publicly available suite of UCL  
! non-conservative force modelling functions for GPSIIA satellites in the FORTRAN 90 language.
! The suite currently includes the following forces: Solar Radiation Pressure (SRP), Thermal Re-Radiation 
! (TRR) and Antenna Thrust (AT)
!
! There are currently two types of GPSII satellites, this model should ONLY be used with block IIA satellites,
! i.e. those block II satellites with SVN>=23. The model's creators have attempted to ensure that the functions 
! will operate correctly with the IIA satellites currently on orbit, but the suite's authors assume no 
! responsibility for errors resulting from the use of this code or any code derived from it.
!
! The model and function suite have been produced at University College London by:
!
! Dr Ant Sibthorpe, Dr Marek Ziebart and Alex Zoutsos.
! Dept. of Civil, Environmental and Geomatic Engineering
! University College London [UCL] 
! Gower Street
! London
! UK
! WC1E 6BT
! +44 (0)207 679 3638 
!
! May 2008
!
! See www.cege.ucl.ac.uk for further contact details.
!
! See the end of this file for expected output to compare your implementation against!

!       declare variables
      implicit none
!       Constant which sets UCL Fourier coefficient vector size, total number of a and b coefficients
      INTEGER*4 UCL_ARRAY_SIZE
!       Declare vectors to hold UCL IIA NAP Fourier coefficients
      REAL*8 IIA_X(73), IIA_Y(73), IIA_Z(73)
!       Mass of the SVN (in Kg) whose SRP/TRR/AT accelearations you are evaluating, default IIA mass is 972.9Kg
      REAL*8 Mass
C       angb:   sun-satellite-earth angle (radians)
      Real*8 angb
C     logical to only call the first time
      Logical alpha
C     X,Y,Z spacecraft body fixed accelerations in km s-2
      Real*8 along_x,along_y,along_z,EvalFourier

      save alpha,IIA_X,IIA_Y,IIA_Z
      data alpha /.true./

C     Set up constants
      Parameter (UCL_ARRAY_SIZE = 73)
      Mass = 972.9d0

!     Populate the UCL Fourier coefficient vectors - first call only
      if (alpha) then
        Print*, "ARC/UCLSRPeditIIA UCLSRP block 3 first call"
        call Fill_GPSIIA_Arrays(UCL_ARRAY_SIZE, IIA_X, IIA_Y, IIA_Z)
        alpha=.false.
C     Debug
Cd      Print*, 'alpha',alpha
Cd      Print*, 'IIA_X',IIA_X
Cd      Print*, 'IIA_Y',IIA_Y
Cd      Print*, 'IIA_Z',IIA_Z
      end if
Cd      Print*, 'IIA_X',IIA_X
Cd      Print*, 'IIA_Y',IIA_Y
Cd      Print*, 'IIA_Z',IIA_Z

!     Evaluate Fourier series for IIA satellite (returns SRP accels in km s-2)
      along_x = EvalFourier(UCL_ARRAY_SIZE, IIA_X, angb)
      along_y = EvalFourier(UCL_ARRAY_SIZE, IIA_Y, angb)
      along_z = EvalFourier(UCL_ARRAY_SIZE, IIA_Z, angb)

	! Mass scaling
        along_x = along_x*(880.0d0/Mass)
	along_y = along_y*(880.0d0/Mass)
	along_z = along_z*(880.0d0/Mass)

Cd      Print*, "along_x",along_x	
Cd      Print*, "along_y",along_y
Cd      Print*, "along_z",along_z

      Return
      END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to fill the UCL GPSIIA Fourier coefficient arrays which are passed as arguments.
! Inputs:
!	UCL_ARRAY_SIZE - constant INTEGER detailing the total number of a and b coefficients
! The following vectors all have length (UCL_ARRAY_SIZE)
!	IIA_X - REAL vector to hold the X GPSIIA coefficients
!	IIA_Y - REAL vector to hold the Y GPSIIA coefficients
!	IIA_Z - REAL vector to hold the Z GPSIIA coefficients
! Outputs:
!	None explicitly, the passed vectors are filled with the required coefficient values
!
! This function need be called only once during your program, unless for some reason you
! decide to alter the values in the vectors (NOT recommended).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Fill_GPSIIA_Arrays(UCL_ARRAY_SIZE,IIA_X,IIA_Y,IIA_Z)

       IMPLICIT NONE
	! Declare variables
C	INTEGER :: UCL_ARRAY_SIZE ! Constant which sets UCL Fourier coefficient vector size
C	REAL*8 (UCL_ARRAY_SIZE)  IIA_X,IIA_Y,IIA_Z
       Real*8  IIA_X(73),IIA_Y(73),IIA_Z(73)
       Real*8  D_X(73),D_Y(73),D_Z(73)
       Integer*4 UCL_ARRAY_SIZE,i
!     Fill vectors - First value is a0, followed by a1...an [cosine terms], then b1...bn [sine terms]
C DATA (A(I),I=1,4)/4*0.0/,A(20)/-1.0/ 

      data D_X / 1.13742775913152d-010,-4.36156370787797d-011,
     . -6.56633236215792d-011,-1.06733547064433d-011,
     . -4.92861461132937d-012, 1.06765186373103d-010,
     . 1.41619425497468d-011,-1.31547743272176d-011,
     . 3.75290900117057d-012,-3.78101194360702d-011, 
     . 8.50309007526584d-013,1.45300929026223d-011,
     . 9.71724825673089d-013,3.41551553277815d-012,
     . -3.39355439522452d-012,-5.29375376038608d-012,
     . 3.19975169915582d-012,-8.53875459808516d-012,
     . -2.29342228839055d-012,9.58566722504006d-012,
     . -3.79369792474613d-013,-8.08555279040173d-012,
     . 3.31145138283496d-012,5.55912489267685d-012,
     . -2.54325029179415d-012,-8.53897886886031d-012,
     . 3.14129841308580d-012,2.89278539307637d-012,
     . -3.94370683083198d-012,5.23319177743290d-013, 
     . 2.57968578057132d-012,-4.25658848597599d-012,
     . -2.91404854297686d-012,-3.00177026674080d-013,
     . 3.86937586891595d-012,-1.34034272966455d-012,
     . -4.30253775305435d-012,-9.95872106516941d-008,
     . -5.09653667061830d-012,9.51099847674255d-011, 
     .  6.66561013331401d-011,7.74832114107248d-010,
     . -1.17138031574701d-010,-2.00786163499211d-010,
     . -1.33472078689119d-010,2.72415767490814d-010,
     . -5.26340243202333d-011,-4.26319548342257d-011,
     . -3.60221931280744d-011,1.15850541999874d-010, 
     . -1.98095356700602d-012,-4.51673426726555d-011,
     . -2.60545435455372d-011,9.11578763178333d-011,
     . -2.19828933429935d-011,-2.73844138254457d-011,
     . -1.78396156452063d-011,7.05758503007920d-011,
     . -1.93421038329100d-011,-1.58748163636388d-011,
     . -1.10575807982788d-011,4.14439483298088d-011,
     . -1.28850972212297d-011,-1.32147550046877d-011,
     . -4.56308308626010d-012,3.30133218019322d-011,
     . -6.27362224404010d-012,-2.64719426221200d-011,
     . -2.37370559742605d-012,4.19197584355166d-011, 
     . -4.14298615271804d-012,-4.09403937099301d-011,
     . -5.46623548286160d-023 /

      data D_Y / -2.88477737476970d-013,-1.65540302309771d-011,
     . -4.90944228800647d-012,-3.29454560184814d-012,
     . 2.00953328718567d-012,5.89828041872537d-012,
     . 1.28742563385268d-012,1.05688486534863d-012,
     . -3.24482599431020d-012,1.41079104001917d-012, 
     . 2.39219941192254d-012,2.41833183273243d-013,
     . -1.56987081121230d-012,1.35785580591530d-012,
     . 2.06266274250845d-012,-2.99630000842186d-013,
     . -1.90619757288763d-013,1.98458772839972d-012,
     . 1.69643229446440d-013,-2.69769451552149d-014, 
     . 6.41736909606966d-013,1.23919332935740d-012,
     . -2.77137163147756d-013,6.51964001809384d-013,
     . 2.46904954814202d-013,8.16826507505550d-013,
     . 3.55643821396458d-013,6.44458541773549d-013,
     . 5.95063488977099d-014,1.37860515717556d-012,
     . 5.26274564524086d-013,8.21858997089367d-013,
     . 3.43158625556963d-013,1.32660959519075d-012,
     . 2.84125255136575d-014,1.08950119127110d-012,
     . 3.97941673842014d-014,1.85343209657004d-012,
     . -1.12790392842564d-011,-1.12556906101797d-011,
     . -9.02815910415205d-012,3.03772534734829d-012,
     . 2.16192181626127d-013,-2.36215722740622d-012,
     . -3.67813933524960d-012,3.96160894125812d-013,
     . -1.25875906969005d-012,7.41204029952882d-013,
     . -2.84077869988806d-012,-1.03481739736461d-012,
     . 5.74831734588647d-013,3.26444349587658d-013,
     . -1.94589704735720d-012,3.03806507666250d-013,
     . -7.86949139543321d-014,-1.30545362743699d-012,
     . -1.33296032370070d-012,7.55658175412787d-013,
     . -6.37908797083073d-013,-3.63182610258572d-013, 
     . 7.91207484585023d-014,-1.22029596623807d-013,
     . -1.18282841053089d-012,-9.08522474517731d-015,
     . -5.97346568894574d-013,-1.24113787868241d-013,
     . -2.38596775330194d-013,-1.09329528613636d-013,
     . -5.51914521145069d-013,-4.72892975089711d-013,
     . 5.58739654797244d-015,1.41530726231829d-013,
     . 5.36974066363256d-026 /
      data D_Z / 4.64850076755965d-010,-1.02164805827643d-007,
     . 1.14946534395330d-010,-7.48797596412440d-011,
     . -3.82353743789959d-010,-3.74570466082087d-010,
     . 1.34876188999835d-010,-2.64652269045367d-011,
     . -5.90732510280964d-011,9.71003829044845d-011, 
     . 3.25123015853716d-012,3.29416924471246d-011,
     . -6.94307512002836d-012,2.05535518169029d-011,
     . -9.21305815860011d-012,2.19014861352533d-011,
     . -7.17435374704896d-012,3.90820290143745d-011,
     . -1.57588512058770d-011,3.54215825252667d-011,
     . -1.55582702525711d-011,2.60175937742926d-011,
     . -1.41018795500849d-011,1.94265137082134d-011,
     . -9.93816803269634d-012,3.30874035477319d-011,
     . -1.67135163223982d-011,2.06838898364183d-011,
     . -5.92456272946695d-012,2.11360778090008d-011,
     . -1.99641112522905d-011,3.32495533388403d-011,
     . -1.98312950434205d-011,2.15407118497605d-011,
     . -9.09387540716478d-012,2.79012267033183d-011,
     . -1.59925356574976d-011,-2.46172078663837d-011,
     . -4.32872641323069d-011,-1.49501479246936d-010,
     . -5.20505680619302d-011,-4.34512396709160d-011,
     . 2.21492089196324d-012,8.62770329806625d-011,
     . 1.54030327904695d-011,2.07106199140511d-011,
     . -1.16771560935035d-012,-5.50376641610658d-012,
     . -6.18034699159949d-012,1.47872102222678d-011,
     . 7.14625583504367d-012,1.70544546841010d-011,
     . -3.92198031890481d-012,-7.96659216516498d-012,
     . -2.19185772311288d-012,5.73296642644366d-012,
     . 9.30429247903458d-013,3.32484524346559d-012,
     . -1.78021348418644d-012,5.42669987784157d-012,
     . 3.25938880256126d-012,5.01213430832117d-012,
     . -2.24237612199282d-013,-3.15795332319993d-012,
     . -2.52274827262122d-013,4.73666053369799d-012,
     . 2.45238311232140d-012,3.79470777342276d-012,
     . -1.80688149899745d-012,-2.11333456884934d-012, 
     . 3.63709158484317d-013,1.83929114484771d-012,
     . -5.96191106043032d-022 /                 
      do i=1,73
         IIA_X(i)=D_X(i)
         IIA_Y(i)=D_Y(i)
         IIA_Z(i)=D_Z(i)
      enddo       
      RETURN
      END 
! End of Fill_UCL_GPSIIA_Arrays subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function to evaluate the UCL Fourier series for SRP
! You should not need to call this function directly, as it is called by Get_UCL_GPSIIA_Accelerations()
! Inputs:
!	UCL_ARRAY_SIZE - Constant INTEGER detailing the total number of a and b coefficients
!	UCL_IIA_COEFFS - One of the REAL vectors containing UCL Fourier coefficients
!	EPS - REAL value of the Earth-Probe-Sun angle in radians
! Outputs:
! 	REAL value equal to the Fourier coefficients computed for the independent argument of EPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION EvalFourier(UCL_ARRAY_SIZE, UCL_IIA_COEFFS, EPS)

      IMPLICIT NONE
	! Declare variables
	REAL*8 Evalfourier ! The variable to be evaluated
	INTEGER*4 UCL_ARRAY_SIZE ! Constant which sets UCL Fourier coefficient vector size
	REAL*8 UCL_IIA_COEFFS(UCL_ARRAY_SIZE) ! Fourier coefficient vector
	REAL*8 EPS ! Variable for the EPS angle
	
	INTEGER*4 counter,AStop,BStart ! Simple counters
	REAL*8 harmonic
	REAL*8 sine_coeff, cosine_coeff
	
        AStop = 1+((UCL_ARRAY_SIZE-1)/2) ! a(1)...a(n) run from (2) to AStop in the Fourier Series vector
	BStart = AStop-1 ! b(1)...b(n) run from BStart+2 to BStart+AStop in the Fourier Series vector
	
	EvalFourier=UCL_IIA_COEFFS(1)*0.5 ! a(0)is the first coefficient in the Fourier Series vector
	
      DO counter=2,AStop ! Default increment of +1 accepted
	 	harmonic=counter-1 ! Cast INTEGER counter to REAL, -1 accounts for fact that counter starts from 2
		cosine_coeff= UCL_IIA_COEFFS(counter) ! a(n)
		sine_coeff = UCL_IIA_COEFFS(BStart+counter) ! b(n)
	EvalFourier = EvalFourier+(cosine_coeff*dCOS(EPS*harmonic))
	EvalFourier = EvalFourier+(sine_coeff*dSIN(EPS*harmonic))
      END DO

	! Scale acceleration into Km/s^2
	EvalFourier = EvalFourier*0.001

      RETURN
      END FUNCTION
! End of Evaluate_UCL_SRP_Fourier subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Expected output for a range of different satellite (& Sun) positions.
! Use these to test your implementation against.  They sample a full
! range of EPS angles; the final set of values give an example of the 
! accelerations for an arbitrary position.

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	26378.137000,	0.000000 /)
!  Km/s^2:  0.0000000E+00  0.9223632E-10  0.0000000E+00
!   m/s^2:  0.0000000E+00  0.9223631E-07  0.0000000E+00

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	22844.136747,	13189.068500 /)
!  Km/s^2:  0.4105931E-13  0.9096618E-10  0.1308964E-11
!   m/s^2:  0.4105931E-10  0.9096618E-07  0.1308964E-08 

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	13189.068500,	22844.136747 /)
!  Km/s^2:  0.9113176E-14  0.9060482E-10  0.5200465E-12
!   m/s^2:  0.9113175E-11  0.9060482E-07  0.5200465E-09

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	0.000000,	26378.137000 /)
!  Km/s^2: -0.1384050E-13  0.8793327E-10  0.6423167E-12
!   m/s^2: -0.1384050E-10  0.8793327E-07  0.6423167E-09

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	-13189.068500,	22844.136747 /)
!  Km/s^2: -0.1023300E-15  0.9093967E-10 -0.1109729E-11
!   m/s^2: -0.1023300E-12  0.9093966E-07 -0.1109729E-08

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	-22844.136747,	13189.068500 /)
!  Km/s^2: -0.2050739E-13  0.9096619E-10 -0.1036570E-11
!   m/s^2: -0.2050739E-10  0.9096620E-07 -0.1036570E-08

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	-26378.137000,	0.000000 /)
!  Km/s^2:  0.0000000E+00  0.9163920E-10  0.0000000E+00
!   m/s^2:  0.0000000E+00  0.9163919E-07  0.0000000E+00

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	-22844.136747,	-13189.068500 /)
!  Km/s^2:  0.2050739E-13  0.9096619E-10  0.1036570E-11
!   m/s^2:  0.2050739E-10  0.9096620E-07  0.1036570E-08

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	-13189.068500,	-22844.136747 /)
!  Km/s^2:  0.1023300E-15  0.9093967E-10  0.1109729E-11
!   m/s^2:  0.1023300E-12  0.9093966E-07  0.1109729E-08

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	0.000000,	-26378.137000 /)
!  Km/s^2:  0.1384050E-13  0.8793327E-10 -0.6423167E-12
!   m/s^2:  0.1384050E-10  0.8793327E-07 -0.6423167E-09

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	13189.068500,	-22844.136747 /)
!  Km/s^2: -0.9113176E-14  0.9060482E-10 -0.5200465E-12
!   m/s^2: -0.9113175E-11  0.9060482E-07 -0.5200465E-09

! ECI_SUN=(/ 0.0, -150000000.0, 0.0 /) 
! ECI_GPS=(/ 0.000000,	22844.136747,	-13189.068500 /)
!  Km/s^2: -0.4105931E-13  0.9096618E-10 -0.1308964E-11
!   m/s^2: -0.4105931E-10  0.9096618E-07 -0.1308964E-08

! ECI_SUN=(/ 143025970.23676714, -36214128.793360971, -15699944.830049932 /) 
! ECI_GPS=(/ 7941.359692, -24367.938402, -7661.956963 /) 
!  Km/s^2: -0.8927984E-10  0.2377457E-10  0.1014883E-10
!   m/s^2: -0.8927984E-07  0.2377457E-07  0.1014883E-07

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
