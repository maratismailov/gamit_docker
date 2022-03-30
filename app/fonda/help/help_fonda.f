      program help_fonda
c
c     guide book for FONDA learner
c
      character*50 comand
      character*6 name
      integer system,ierr,nfile
c
      comand(1:50) =
     .   'more                            ../help/fonda.help'
      ierr = system(comand)
      
 10   print*,'         *******************************     '
      print*,' package name: maked,solvem,diagno,fmodel,disply,design'
      Print*,' <Please enter package name, wrong name leads to exit>'
      read (5,'(a)') name
      
      if (name(1:5).eq.'MAKED'.or.name(1:5).eq.'maked') then
         comand(1:50) =
     .   'more                            ../help/maked.help'
         ierr = system(comand)
 110     print*,'         *******************************     '
         print*,' Please pick up a number to see the examples:'
         print*,' 1=driving file,  2=input file,    3=mapping file'
         print*,' 4=output file,   5=vmodel file,   6=priori file'
         print*,' 7=other file,    8=other package, 9=exit'
         read (5,*) nfile
         comand(1:31) =
     .   'more          ../example/maked/'
         if (nfile.eq.1) comand(32:50) = 'maked.ngs.drv      '
         if (nfile.eq.2) comand(32:50) = 'maked.ngs.in       '
         if (nfile.eq.3) comand(32:50) = 'maked.ngs.map      '
         if (nfile.eq.4) comand(32:50) = 'maked.ngs.out      '
         if (nfile.eq.5) comand(32:50) = 'vmodel.ngs         '
         if (nfile.eq.6) comand(32:50) = 'ngsdat.pri         '
         if (nfile.eq.7) goto 110
         if (nfile.eq.8) goto 10
         if (nfile.le.0.or.nfile.gt.8) goto 1000
         print*,'           ------------------------------  '
         ierr = system(comand)
         goto 110
      endif
      
      if (name.eq.'SOLVEM'.or.name.eq.'solvem') then
         comand(1:50) =
     .   'more                                ../solvem.help'
         ierr = system(comand)
 210     print*,'         *******************************     '
         print*,' Please pick up a number to see the examples:'
         print*,' 1=driving file,  2=input file,    3=mapping file'
         print*,' 4=solution file, 5=residual file, 6=strain file'
         print*,' 7=other file,    8=other package, 9=exit'
         read (5,*) nfile
         comand(1:30) =
     .   'more        ../example/solvem/'
         if (nfile.eq.1) comand(31:50) = 'solvem.ngs.drv      '
         if (nfile.eq.2) comand(31:50) = 'solvem.ngs.dat      '
         if (nfile.eq.3) comand(31:50) = 'solvem.ngs.map      '
         if (nfile.eq.4) comand(31:50) = 'solvem.ngs.out      '
         if (nfile.eq.5) comand(31:50) = 'solvem.ngs.res      '
         if (nfile.eq.6) comand(31:50) = 'solvem.ngs.str      '
         if (nfile.eq.7) goto 210
         if (nfile.eq.8) goto 10
         if (nfile.le.0.or.nfile.gt.8) goto 1000
         print*,'           ------------------------------  '
         ierr = system(comand)
         goto 210
      endif
      
      if (name.eq.'DIAGNO'.or.name.eq.'diagno') then
         comand(1:50) =
     .   'more                           ../help/diagno.help'
         ierr = system(comand)
         goto 10
      endif
      
      if (name.eq.'FMODEL'.or.name.eq.'fmodel') then
         comand(1:50) =
     .   'more                           ../help/fmodel.help'
         ierr = system(comand)
         goto 10
      endif
      
      if (name.eq.'DISPLY'.or.name.eq.'disply') then
         comand(1:50) =
     .   'more                           ../help/disply.help'
         ierr = system(comand)
         goto 10
      endif
      
      if (name.eq.'DESIGN'.or.name.eq.'design') then
         comand(1:50) =
     .   'more                           ../help/design.help'
         ierr = system(comand)
         goto 10
      endif

 1000 continue
      stop 
      end 

