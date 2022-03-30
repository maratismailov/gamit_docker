      program exclude

c  Program to read the vscan.out file from scandd and exclude the residuals
c  of a list of sites and satellites as requested by the user. The program
c  will read a qfile and find out which sites and satellites have no observations
c  and will suggest these to the user. The user is then allowed to add to this
c  list as required.
c
c  P. Tregoning
c  27th November, 1996

      implicit none
      
      include '../includes/dimpar.h'

      integer i,j,ierr,nsit,nsat,iqfile,iscrn
      character qfile*80, char1*1, line*85,ex_sit(maxsit)*2
     .          ,ex_sat(maxsat)*2,message*256
      CHARACTER*40 VERS
      logical use_q,flag   
     
      parameter (iscrn = 6)

      call uversn(iscrn,vers)

c  write a header of sorts
      write(iscrn,90)
90    format(/,'#######################################################'
     .,'#############',/,'                    Program VEXCLUDE         '
     .,//,'This utility will read a qfile and remove from the vscan.out'
     .,' file any',/,'sites and satellites which were excluded from the'
     .,' SOLVE solution. The',/,'input  vscan.out file will be '
     .,'overwritten  without the excluded sites ',/,'and satellites.'
     .,' The user can also add to this list interactively.'
     .,//,'############################################################'
     .,'########',//)


      write(*,'(a,$)')' Enter qfile to scan for excluded sites/sats : '
      read(*,'(a)')qfile

      open(unit=10,file=qfile,status='old',iostat=ierr)
      if(ierr.ne.0)then           
       print*,' have an error opening qfile'
        write(message,'(a,a,a)')'Error opening ',qfile
     .                      ,'Proceeding without it!'
        call report_stat('WARNING','VEXCLUDE','vexclude',' ',
     .  message,0)
        use_q = .false.
      else   
        iqfile = 10
        call read_q(iqfile,ex_sit,ex_sat,nsit,nsat) 
        use_q = .true.
      endif

c  open the vscan.out file, and some temporary files
      open(unit=20,file='vscan.out',status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(message,'(a)')'Error opening vscan.out - does it exist?'
          call report_stat('FATAL','VEXCLUDE','vexclude',' ',message,0)
          stop
        endif

      open(unit=30,file='tmp.pt',status='unknown')
      open(unit=40,status='scratch')
                      
c  now we have the sites/sats to exclude as per qfile. 
c  Now list them and ask user to modify.

      write(*,'(/,a)')' Sites to be excluded : '
      if(nsit.gt.0)then
        write(*,'(1x,a2)')(ex_sit(j),j=1,nsit)
      endif

      write(*,'(/,a)')' Satellites to be excluded : '
      if(nsat.gt.0)then   
        write(*,'(1x,a2)')(ex_sat(j),j=1,nsat)
      endif

      write(*,'(/,a)')' Do you wish to add to this (n): '
      read(*,'(a1)')char1

      if(char1.ne.'n'.and.char1.ne.'N'.and.char1.ne.' ') then
        call change_excl(ex_sat,ex_sit,nsit,nsat)
      endif


c  now have list of things to exclude. Now read vscan.out and rewrite it but
c  exclude the sites and satellite entries which are not required.

      do while (ierr.eq.0)

        read(20,'(a)',iostat=ierr)line
        if(ierr.eq.0)then

c  very crudely at first, if there is a character match in the
c  particular column for either site(s) or satellite(s) then don't
c  write out the line!
          flag = .true.
          if(nsat.gt.0)then  
            do i=1,nsat    
              if(line(3:4).eq.ex_sat(i).or.line(9:10).eq.ex_sat(i))then    
                flag = .false.
              endif
            enddo
          endif  
 
          if(nsit.gt.0)then
            do i=1,nsit
              if(line(15:16).eq.ex_sit(i).or.line(21:22)
     .                 .eq.ex_sit(i))then
                flag = .false.
              endif
            enddo
          endif 

          if(flag)then
             write(30,'(a)')line
          endif

        endif
      enddo

c  finally transfer the temp file to vscan.out
      ierr = 0   
      rewind(20)
      rewind(30) 
      do while (ierr.eq.0)
        read(30,'(a)',iostat=ierr)line
        if(ierr.eq.0)then
          write(20,'(a)')line
        endif
      enddo


      stop ' Normal stop in VEXCLUDE. Vscan.out has been overwritten.'
      end
