      subroutine satatt(iyaw,time,nsat,yatt,ievent,iepoch)

c  subroutine to read one line of the y table and return all the attitudes
c  of the satellites in the array yatt
c
c  IN:
c      iyaw    :  unit number of yaw table
c      time    :  julian time of observation
c      nsat    :  number of satellites in tfile/y table
c
c  OUT:
c      yatt    :  attitudes of satellites at time of observation
c
c P Tregoning
c 1 December 1997
c Mod 27 Feb 1998: fix bug in interpolating between tabular values

      implicit none
      include '../includes/dimpar.h'

      integer*4 iyaw,nsat,ievent(nsat),i
c      for debug
     .   ,iepoch
      real*8 yatt(maxsat),time,ytime,ytime_old,yatt_old(maxsat),factor
     .       ,temp
      logical match
      character*256 message

      save ytime_old,yatt_old
      match = .false.
          
      do while(.not.match)
c  read the next entry in the yaw table
        read(iyaw) ytime,(yatt(i),ievent(i),i=1,nsat)           
cd         print*, " debug readyaw",yatt(i),nsat,iyaw,ytime,i


c  see whether this time is within 1 second of observation time  
        if(dabs(ytime-time).lt.1.d0/86400.d0)then
          match = .true.
        elseif(ytime-time.gt.1.d0/86400.d0)then 

c  in this case just do a linear interpolation between the current values of attitudes
c  and the previous values
          factor = (time-ytime_old)/(ytime-ytime_old) 
          do i = 1,nsat    
c convert attitudes to 0 - 360 range if at the [180,-180] boundary  
            if(dabs(yatt(i)-yatt_old(i)).gt.200)then  
              if(yatt(i).lt.0.d0)yatt(i) = yatt(i) + 360.d0
              if(yatt_old(i).lt.0.d0)yatt_old(i)= yatt_old(i) + 360.d0  
            endif

            temp = yatt(i)     
            yatt(i) = yatt(i) - factor*(yatt(i) - yatt_old(i))
            yatt_old(i) = temp
            if(yatt(i).gt.180.d0)yatt(i) = yatt(i) - 360.d0
            if(yatt_old(i).gt.180.d0)yatt_old(i) = yatt_old(i) - 360.d0
          enddo 
          match = .true. 
        elseif(time-ytime.gt.1.d0/86400.d0)then    
          do i = 1,nsat
            yatt_old(i) = yatt(i)
          enddo
        else

          write(message,'(a,a)')'Time is greater than ytime - missed'
     .          ,' match in y table'   
          call report_stat('WARNING','MODEL','satatt',' ',message,0)
          write(message,'(a,a)')'Contact pault@rses.anu.edu.au'
          call report_stat('FATAL','MODEL','satatt',' ',message,0)
        endif
        ytime_old = ytime
      enddo
    
      return
      end
