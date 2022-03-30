      Program TESTEPHRED

c     Temporary program to test new coding for the planetary ephermerides
c      rwk 180111

      implicit none
     
      include '../includes/dimpar.h'  
      include '../includes/units.h'
      include '../includes/global.h' 
      include '../includes/arc.h'

      integer*4 ivel,ispeed,i
      

      real*8 pjd,rvec(6)

      ibody = 1             
c     isun and ilun set to 33 and 34 in ephdrd 
      ivel = 0   
      ispeed = 1 
c     turn on velocities
      ivel = 1  
      ispeed = 0 
      frame = 'J2000' 

c     turn on venus and jupter
      lbody = 1 
   
c     start with the epoch used for the ocean tide tests in ARC  (2016 201 12 hrs)   
c  called SOLRED t ccor    2457589.5005924073        69042083.677939311       -124259339.94129406       -53868225.366969630        2250051.3260036074        1063265.6072976775        460992.45086404867     
       pjd = 2457589.5005924073d0 
c  use the date in the testsp3 directory 2018 1 17 12 hrs 
c       pjd = 2458136.5005924073d0  

c     Print the epoch
      write(*,'(a,f18.10)') 'Testing PJD ',pjd

c     Open and interpolate soltab. and luntab.
      call ephdrd(pjd)     
      call evrtcf  
      call solred(ispeed,pjd,rvec)
      write(*,'(a,6f16.2)') 'SOLRED Earth   ',(rvec(i),i=1,6)
      call lunred(ispeed,pjd,rvec)
      write(*,'(a,6f16.2)') 'LUNRED Moon    ',(rvec(i),i=1,6)
                              
c     Open and read the ephemeris 
c       nbody740.2020 (binary PEP) 
c       nbody740.2020.asc (ascii PEP) 
c       JPL.DE200  (binary JPL
c     and link the choice to nbody  

      call ephred(ibody,pjd,lbody,ivel)

c     Intepolate to the requested epoch
c     Earth
       call ephtrp(pjd,3,ivel,rvec)
       write(*,'(a,6f16.2)') 'EPHTRP Earth   ',(rvec(i),i=1,6)

c     Moon
       call ephtrp(pjd,10,ivel,rvec)
       write(*,'(a,6f16.2)') 'EPHTRP Moon    ',(rvec(i),i=1,6)

c     Venus     
       call ephtrp(pjd,2,ivel,rvec)
       write(*,'(a,6f16.2)')'EPHTRP Venus   ',(rvec(i),i=1,6)

c     Jupiter          
       call ephtrp(pjd,5,ivel,rvec)
       write(*,'(a,6f16.2)')'EPHTRP Jupiter ',(rvec(i),i=1,6)

      stop
      end 

     



