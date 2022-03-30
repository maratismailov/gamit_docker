      program test_wild

* Install with
* gfortran test_wild.f gen_util_lib.a ../../libraries/comlib/com_lib.a -o test_wild

      implicit none

      character*256 file, key, jdstr, newfile
      integer*4 lenrun

      integer*4 rcpar, i

      real*8 mjd 

      lenrun = rcpar(1,file)
      lenrun = rcpar(2,key)
      lenrun = rcpar(3,jdstr) 
      if( lenrun.gt.0 ) read(jdstr,*) mjd 
      if( lenrun.eq.0 ) then
         print *,'Usage: test_wild <file> <key> <MJD>'
         print *,
     .     'EXAMPLE: test_wild "@_<WWWW><D>.bak" mit19122.gdl 58653.6'
         print *,'   <file> is the name in globk command'
         print *,'   <key>  would be typically the globk .gdl file name'
         stop 'Incomplete runstring test_wild'
      endif
      write(*,100) trim(file), trim(key), mjd
 100  format('For Input file  : ',a,/,
     .       'And Key         : ',a,' MJD ',F10.2)
      call wild_card(file, key)
      newfile = file
      call sub_char(newfile,'.bak','.com')
      write(*,120) trim(file)
      call wild_date(file, mjd, 'I' )
      write(*,120) trim(file)
 120  format('Output filename : ',a)

      do i = 1,3
          mjd = mjd + i
          call wild_date(file, mjd, 'I' )
          write(*,140) MJD, trim(file)
 140      format('Iterated ',F10.2,' Name ',a)
          call wild_date(newfile, mjd, 'I' )
          write(*,140) MJD, trim(newfile)
      end do

       

      end

