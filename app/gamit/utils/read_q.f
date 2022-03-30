      subroutine read_q(iqfile,ex_sit,ex_sat,nsit,nsat)

c  Subroutine to make the list of excluded sites/sats to be excluded. It
c  reads the qfile and determines which sites/sats have no observations.
c
c  P. Tregoning
c  27th November 1996
c
c  input/output: 
c              iqfile :  unit number of qfile (= 10)                 I*4
c                nsit :  number of sites to exclude                  I*4
c                nsat :  number of satellites to exclude             I*4
c              ex_sit :  list of site numbers to be excluded         C*2(maxsit)
c              ex_sat :  list of satellite numbers to be excluded    C*2(maxsat)

      implicit none

      include '../includes/dimpar.h'

      integer nsat,nsit,iqfile,satnum,stanum,ierr,i,pos
      character*2 ex_sit(maxsit),ex_sat(maxsat)
      character*256 line
              
      stanum = 0
      satnum = 0

c read the qfile until we get to the section of the table of ddiffs per site/sat
      do while (ierr.eq.0)                         
        read(iqfile,'(a)',iostat=ierr)line
        if(ierr.eq.0)then

c get the number of stations
          if(line(2:18).eq.'Tracking stations')then  
            do while (line(2:3).ne.'  ')  
              stanum = stanum + 1
              read(iqfile,'(a)')line 
            enddo
          endif

          if(line(2:20).eq.'Satellites observed')then  
            do while (line(2:3).ne.'  ')  
              satnum = satnum + 1
              read(iqfile,'(a)')line   
            enddo
          endif
              
c find the list of sites/sats used
          if(line(2:14).eq.'Stations used')then  
            read(iqfile,'(a)')line
            read(iqfile,'(a)')line    
            do i=1,stanum - 1    
              pos = 1 + 4*i   
              if(line(pos:pos).eq.'N') then
                nsit = nsit+1
                write(ex_sit(nsit),'(i2)')i  
              endif
            enddo
          endif

          if(line(2:16).eq.'Satellites used')then
            read(iqfile,'(a)')line
            read(iqfile,'(a)')line
            do i=1,satnum - 1   
              pos = 1 + 4*i
              if(line(pos:pos).eq.'N') then
                nsat = nsat+1
                write(ex_sat(nsat),'(i2)')i
              endif
            enddo
          endif


c and finally put a break so that the whole qfile is not read

          if(line(2:26).eq.'Assumed ionosphere error')then  
            return
          endif

        endif
      enddo

      end
