      subroutine change_excl(ex_sat,ex_sit,nsit,nsat)

c  Subroutine to change the list of excluded sites/sats for 
c  reading vscan.out
c
c  P. Tregoning
c  27th November 1996
c
c  input/output:
c                nsit :  number of sites to exclude                  I*4
c                nsat :  number of satellites to exclude             I*4
c              ex_sit :  list of site numbers to be excluded         C*2(maxsit)
c              ex_sat :  list of satellite numbers to be excluded    C*2(maxsat)

      implicit none

      include '../includes/dimpar.h'

      integer nsat,nsit,i
      character*1 char1
      character*2 ex_sit(maxsit),ex_sat(maxsat)

      write(*,'(a)')' Add another site to be excluded ? (y) : '
      read(*,'(a)')char1
      if(char1.eq.'y'.or.char1.eq.'Y'.or.char1.eq.' ')then  
        char1 = 'y'
      endif

      do while (char1.eq.'y')
        nsit = nsit+1
        write(*,'(a,$)')' Enter site number : '
        read(*,'(i2)')i
        write(ex_sit(nsit),'(i2)')i
        write(*,'(a)')' Add another site to be excluded ? (y) : '
        read(*,'(a)')char1
        if(char1.eq.'y'.or.char1.eq.'Y'.or.char1.eq.' ')then  
          char1 = 'y'
        endif
      enddo

      write(*,'(a)')' Add another satellite to be excluded ? (y) : '
      read(*,'(a)')char1
      if(char1.eq.'y'.or.char1.eq.'Y'.or.char1.eq.' ')then  
        char1 = 'y'
      endif

      do while (char1.eq.'y')
        nsat = nsat+1
        write(*,'(a,$)')' Enter satellite number : '
        read(*,'(i2)')i
        write(ex_sat(nsat),'(i2)')i
        write(*,'(a)')' Add another satellite to be excluded ? (y) : '
        read(*,'(a)')char1
        if(char1.eq.'y'.or.char1.eq.'Y'.or.char1.eq.' ')then  
          char1 = 'y'
        endif
      enddo



      return
      end
