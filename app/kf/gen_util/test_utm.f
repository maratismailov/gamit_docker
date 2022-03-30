      program test_utm

      implicit none 

      include '../includes/utmut.h' 

      include '../includes/const_param.h'

      real*8 incrd(3), outcrd(3)
      integer*4 zone
      character*8 intype, indatum
      character*8 outtype, outdatum
      character*1 hemi
      logical OK

****  Test old calls
C     write(*,*) 'Geod Lat Long Ht '
C     read(*,*) incrd
C     incrd(1) = (90-incrd(1))*pi/180
C     incrd(2) = incrd(2)*pi/180
C     print *,'GEOD ',incrd
C     call GEOD_to_UTM(incrd, outcrd, zone, hemi, 'WGS84')
C     print *,outcrd, zone, hemi

***** Get the input and output systems
      write(*,*) 'Input type, datum and outputs? '
      read(*,*) intype, indatum, outtype, outdatum
      OK = .true.
      do while (OK )
         if( intype.ne.'UTM') then
             read(*,*) incrd
         else
             read(*,*) incrd, zone, hemi
         endif
         if( incrd(1).eq.0 ) OK = .false.
         if( OK ) then 
            if( intype.eq.'GEOD' ) then
               incrd(1) = (90.d0 - incrd(1))*pi/180.d0
               incrd(2) = incrd(2)*pi/180
            end if
            call geod_to_geod(incrd, outcrd, intype, outtype,
     .                           indatum, outdatum, zone, hemi)
            if( outtype.eq.'UTM' ) then
                write(*,140) outcrd, zone, hemi
 140            format(3F13.4,1x,i2,1x,a1)
            else if( outtype.eq.'GEOD' ) then
                write(*,130) (pi/2-outcrd(1))*180/pi,outcrd(2)*180/pi,
     .                        outcrd(3) 
 130            format(2F16.9,1x,F12.4)
            else
                write(*,120) outcrd
 120            format(3F15.4)
            endif
         endif 
      enddo



      end


