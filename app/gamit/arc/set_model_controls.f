      Subroutine set_model_controls

c     Check for availability and consistency of the force models;
c     assign the parameters for the radiation-pressure and gravity
c     models.     
c     R. King  20 May 2012

      implicit none
                                   
      include '../includes/dimpar.h' 
      include '../includes/global.h'
      include '../includes/arc.h'
                    
      integer*4 i,i1
      character*5 upperc 
      character*256 message
      logical new/.true./,debug/.false./  

c     local for zero-ing out tesserals: values are the index of the 1st degree higher than requested
      integer*4 tessindx(11)/1,3,6,10,15,21,28,35,45,55,66/

c  Load the conversion factors
 
      call dtwopi

           
c  Set the solar radiation pressure model

      call uppers(srpmod)
      if ( srpmod .eq. 'BERNE' .or. srpmod .eq. 'ECOM1' ) then
        modrad = 1
        nics = 15
      elseif ( srpmod .eq. 'ECOM2'.or.srpmod.eq.'ECOMC' ) then
        modrad = 2   
        nics = 19     
      elseif( srpmod. eq. 'UCLR1' ) then
        modrad = 7
        nics = 15      
      elseif( srpmod. eq. 'UCLR2' ) then
        modrad = 8
        nics = 15 
      else
         write(message,'(a,a5)' ) 'Invalid radiation-pressure model: '
     .        ,srpmod
         call report_stat('FATAL','ARC','arc',' ',message,0)
      endif    

c  Set up the gravity field 

c     set gm for harmonics same as gm by default
      gmhfct = 1.d0
cd      print *,'SET MODEL CONTROLS gravmod ',gravmod
      if (upperc(gravmod).eq.'IGS92') then
        nctess = 8 
        nczone = 8 
        call igs92
cd        write (iarh,'(1x,a)') 'IERS92/IGS standards for model constants'
      elseif (upperc(gravmod).eq.'EGM96') then
cd        print *,'calling EGM96 '
        nctess = 10
        nczone = 10         
        call egm96   
cd        print *,'after EGM96 cchar(1-2) ',cchar(1),cchar(2)
cd        print *,'nczone nctess nczon1 nctes1  ',
cd     .           nczone,nctess,nczon1,nctes1
cd        write (iarh,'(1x,a)') 
cd     .                 'EGM96/IERS2000 standards for model constants'   
      elseif (upperc(gravmod).eq.'EGM08'.or.
     .        upperc(gravmod).eq.'EGR08' ) then
        call egm08
        nctess = 12
        nczone = 12
cd        write (iarh,'(1x,a)') 
cd     .                 'EGM08/IERS2010 standards for model constants'
cd        print *,'nczone nctess nczon1 nctes1 '
cd     .         , nczone,nctess,nczon1,nctes1
      else
        call report_stat('FATAL','ARC','set_model_constrols','gravmod',
     .  'Requested model constants not known:',0)
      endif     
      nctes1 = nctess-1
      nczon1 = nczone-1
      nctes2 = nczone*(nczone-1)/2 + nctess - 1  
      write(iarh,'(2a,i3,a,i3,a,i3,/)') 
     .   'Degree and order of harmonics used for '
     .  ,'static gravity ',gravdeg,'  solid-Earth tides ',etidedeg
     .  ,'  ocean tides ',otidedeg
c     unscale the harmonic coefficients
      new = .false.
      if( new ) then  
        if(debug) print *,'before SCLCOF1 cchar(1-2) ',cchar(1),cchar(2)
        call sclcof1(-1,0,nczone,nctess,czhar)
        call sclcof1(-1,1,nczone,nctess,cchar)
        call sclcof1(-1,1,nczone,nctess,cshar) 
        if( debug ) then      
          write(*,'(a)') ' Low order harmonics after unscaling '
          write(*,'(a,4d18.10)') ' Zonals (2-4) ',(czhar(i),i=1,4)
          write(*,'(a,2d18.10)') ' Cosine Tesserals 2,1 2,2 '
     .       ,(cchar(i),i=1,2)
          write(*,'(a,2d18.10)') ' Sine   Tesserals 2,1 2,2 '
     .       ,(cshar(i),i=1,2)
          write(*,'(a,2i4)') ' nczone nctess ',nczone,nctess
        endif
      else
        if(debug) print *,'before SCLCOF cchar(1-2) ',cchar(1),cchar(2)
        call sclcof    
        if( debug ) then
          write(*,'(a)') ' Low order harmonics after unscaling '
          write(*,'(a,4d18.10)') ' Zonals (2-4) ',(czhar(i),i=1,4)
          write(*,'(a,2d18.10)') ' Cosine Tesserals 2,1 2,2 '
     .         ,(cchar(i),i=1,2)
          write(*,'(a,2d18.10)') ' Sine   Tesserals 2,1 2,2 '
     .       ,(cshar(i),i=1,2)
          write(*,'(a,2i4)') ' nczone nctess ',nczone,nctess      
        endif
      endif 
c     if requested a lower degree than the input field, zero out the 
c     higher or coefficients
      if( gravdeg.lt.nczone ) then
        if( gravdeg.lt.2 ) then
          i1 = 1
        else
          i1 = gravdeg 
        endif  
        do i=i1,nczon1                
          czhar(i) = 0.d0
          czhar(i) = 0.d0
        enddo           
        if( gravdeg.lt.2 ) then
          i1 = 1
        else
          i1 = tessindx(otidedeg)
        endif
        do i=i1,nctes2
          cchar(i) = 0.d0
          cshar(i) = 0.d0
        enddo 
      endif

  
c Set the Love numbers and frequency-dependent terms for the solid-Earth tidal accerations
  
c     From IERS 2010 standards as revised 10 August 2012, Table 6.3, 
c     and coded in NGS ORB global.h:

c // Nominal values of solid Earth tide
const double k2mr[3] = {  0.30190,  0.29830,  0.30102 }; // Re k2m, m=0,1,2
const double k2mi[3] = { -0.00000, -0.00144, -0.00130 }; // Im k2m, m=0,1,2
const double k2mp[3] = { -0.00089, -0.00080, -0.00057 }; // k2m(+), m=0,1,2
const double k3m[4]  = { 0.093, 0.093, 0.093, 0.094 };   // k3m

      k2mr(1) = 0.30190d0
      k2mr(2) = 0.29830d0
      k2mr(3) = 0.30102d0
      k2mi(1) = -0.00000d0
      k2mi(2) = -0.00144d0
      k2mi(3) = -0.00130d0
      k2mp(1) = -0.00089d0
      k2mp(2) = -0.00080d0
      k2mp(3) = -0.00057d0
      k3m(1) =  0.093d0
      k3m(2) =  0.093d0
      k3m(3) =  0.093d0
      k3m(4) =  0.094d0        
      call iers_etides

      return
      end


