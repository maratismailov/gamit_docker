c   Utility to make a station.info file interactively
c
c    R. King 27 April 2001; revised for simpler, rigid format 20 January 2020
c
      implicit none

c       station.info variables
             
      character*1 ans                   
      character*4 sitcod
      character*5 htcod,radome  
      character*6 rcvcod,antcod
      character*15 anttyp
      character*16 sname
      character*20 rcvsn,antsn 
      real*4 anth,antn,ante,swver 
      integer*4 start(5),stop(5)
                             
c       function
      logical fcheck,end

c       other variables
             
      character*16 uname 
      character*32 fname

      integer*4 ustnfo,irunt(6),ihnsec,ioerr,i
      data ustnfo/10/
         
      write(*,*) ' '            
      write(*,*) 'Create a station.info.[user] file interactively' 
      write(*,*) ' ' 
      write(*,*)  'Use GAMIT 6-character codes for receiver and antenna'
      write(*,*)  '  (lowercase ok) ' 
      write(*,*) ' ' 
      call getdat(irunt(1),irunt(2),irunt(3))
      call gettim(irunt(4),irunt(5),irunt(6),ihnsec )
      call getusr(uname)
             
c Open the file         
      fname = 'station.info.'//uname 
      if( fcheck(fname) ) then
        write(*,*) 'File ',fname,' exists, overwrite?'
        read(*,*) ans
        if( ans.eq.'y' ) then 
          open(unit=ustnfo,file=fname,status='old') 
        else
          stop
        endif  
      else
        open(unit=ustnfo,file=fname,status='new')  
      endif 

c Write the header
                     
      write(ustnfo,'(a,a16,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)')
     .    '*  Station.info created manually using make_stnfo by ',uname
     .    ,irunt(1),"-",irunt(2),"-",irunt(3),":",irunt(4),":",irunt(5) 
      write(ustnfo,'(a)') '*'
      write(ustnfo,'(3a)')
     .  '*SITE  Station Name      Session Start      Session Stop     '
     . ,'  Ant Ht   HtCod  Ant N    Ant E    RcvCod  SwVer  '
     . ,'Receiver SN           AntCod  Dome   Antenna SN'
    
c Read in the entries
                          
      end = .false.
      do while (.not.end)   

c       reset the defaults
        sitcod = ' '
        sname = ' '
        start(1) = 1900
        start(2) = 1     
        stop(1) = 9999
        stop(2) = 999
        do i=3,5
          start(i) = 0
         stop(i) = 0
        enddo
        anth = 0.
        antn = 0.
        ante = 0.
        htcod = ' '
        rcvcod = ' '
        swver = 0.0
        rcvsn = ' '
        antcod = ' '   
        radome = ' '
        antsn = ' '
        write(*,*) 'Enter 4-char sitcod (space to quit)'
        read(*,'(a4)',iostat=ioerr) sitcod    
        if( sitcod(1:1).eq.' ' .or. ioerr.ne.0 ) then
          end = .true.  
          goto 999
        endif      
        call uppers(sitcod)
        write(*,*) 'Enter station name (up to 16 chacters)'
        read(*,'(a16)') sname      
        write(*,*) 'Enter start year and day-of-year'
        read(*,*) start(1),start(2)
        write(*,*) 'Enter start hr min sec '
        read(*,*) start(3),start(4),start(5)
        write(*,*) 'Enter stop year and day-of-year'
        read(*,*) stop(1),stop(2)
        write(*,*) 'Enter stop hr min sec '
        read(*,*) stop(3),stop(4),stop(5) 
        write(*,*) 'Enter antenna height (m)'
        read(*,*) anth         
        write(*,*) 
     .     'Enter the 5-char type of height measurement (CR for DHARP)'
        read(*,'(a5)',iostat=ioerr) htcod     
        if( htcod.eq.'     ' ) htcod = 'DHARP'
        if( ioerr.eq.0 ) call uppers(htcod)
        write(*,*) 'Enter antenna offsets N E (m) '
        read(*,*,iostat=ioerr) antn,ante
        write(*,*,iostat=ioerr) 'Enter 6-char receiver code'
   1    read(*,*,iostat=ioerr) rcvcod       
        if( ioerr.ne.0 ) then
          write(*,*) 'Receiver code must be entered' 
          goto 1
        endif
        call uppers(rcvcod)
        write(*,*) 'Enter receiver serial number (CR to omit)'
        read(*,'(a20)',iostat=ioerr) rcvsn           
        write(*,'(a)') 'Enter GAMIT firmware code (CR for 0.0)'
        read(*,'(f10.0)',iostat=ioerr) swver  
        write(*,*) 'Enter 6-char antenna code'
    2   read(*,*,iostat=ioerr) antcod       
        if( ioerr.ne.0 ) then
          write(*,*) 'Antenna code must be entered '
          goto 2 
        endif
        call uppers(antcod)
        write(*,*) 'Enter 5-char radome code (CR for NONE)'
        read(*,'(a5)',iostat=ioerr) radome 
        if( radome.eq.'     ' ) radome = 'NONE '  
        call uppers(radome)
        write(*,*) 'Enter antenna serial number (CR to omit)'
        read(*,'(a20)',iostat=ioerr) antsn     
        
c       Write the entry   

        write(ustnfo,'(1x,a4,2x,a16,2(2x,i4,1x,i3,3i3),f9.4,2x,a5,2f9.4
     .      ,2x,a6,2x,f5.2,2x,a20,2x,a6,2x,a5,2x,a20)',iostat=ioerr)   
     .   sitcod,sname,start,stop,anth,htcod,antn,ante,rcvcod,swver,rcvsn
     .  ,antcod,radome,antsn
        if( ioerr.ne.0 ) then 
           write(*,'(a)') 'Error writing values to station.info '
        endif
    
      enddo

 999  stop
      end


     


            
      

