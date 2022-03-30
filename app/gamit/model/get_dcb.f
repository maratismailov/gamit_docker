      Subroutine get_dcb( iudcb,jd,numdcb,prndcb,dcb ) 

c     Read the file of differential code biases as a function of time 
c     and SV and return the values 

      implicit none                    

      include '../includes/dimpar.h'



c     Input
c     -----   
c      iudcb  unit number for dcb.dat
c      jd     PEP Julian Day of observation

c     Output
c     ------                  
c       numdcb          number of entries in dcb array 
c       prndcb(maxsat)  PRN #s for dcb array
c       dcb(maxsat)     corrections to C1 or P2' in ns
                            

      integer*4 iudcb,numdcb,jd,prndcb(32),ioerr
     .        , jdoy,jyear,jmon,jday,jdi,julday
      
      real*8 dcb(32)   

      character*80 line

      logical end_of_header,first_found,done
                              

c  Find the line before the first first epoch
           
      end_of_header = .false.
      do while (.not.end_of_header ) 
        read( iudcb,'(a)',iostat=ioerr) line
c        print *,'header ',line
        if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','get_dcb'
     .      ,' ','Error reading header of dcb.dat file',ioerr)
        if( line(1:4).eq.'*---' ) then
           end_of_header = .true.
        endif
      enddo          

c   Loop through entries in the table to find the latest epoch before the observations

c       Four cases:
c         1) Reading first epoch line: decode
c         2) Reading DCBs:  store and continue 
c         3) Encountered second epoch: decode, test, and continue
c              if before obs epoch. store and continue
c              if after obs epoch, exit 
c         4) Encountered EOF after one epoch: exit
             
      numdcb = 0                
      first_found = .true.
      done = .false.  
      do while( .not.done )

        read( iudcb,'(a)',iostat=ioerr) line  
c          print *,'read line ',line
        
        if( ioerr.eq.-1 ) then
          if( first_found .and.numdcb.gt.0 ) then 
            done = .true.          
c            print *,'eof after first_found, exit'
          else
            call report_stat('FATAL','MODEL','get_dcb'
     .        ,'dcb.dat','Unexpected EOF on DCB file',ioerr)
          endif
        elseif( ioerr.ne.0 ) then
           call report_stat('FATAL','MODEL','get_dcb'
     .        ,'dcb.dat','Error reading DCB file',ioerr)
        elseif(  line(2:6).eq.'Epoch') then  
            read(line(8:15),'(2i4)',iostat=ioerr) jyear,jdoy
            if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','get_dcb'
     .        ,' ','Error reading epoch from dcb.dat file',ioerr) 
            call monday(jdoy,jmon,jday,jyear)   
c            print *,'jdoy jmon jday jyear ',jdoy,jmon,jday,jyear
            jdi = julday(jmon,jday,jyear)  
            first_found = .true. 
c            print *,'read jdi ',jdi
c           reset the dcb count only if this epoch is to be used
            if( jdi.le.jd ) numdcb = 0
        else
c         assume this is a DCB line
          if( jd.ge.jdi ) then
c           tabular epoch preceeds obs epoch: store the value
            numdcb = numdcb + 1 
c            print *,'numdcb ',numdcb   
            if( numdcb.gt.32 ) call report_stat('FATAL','MODEL'
     .                     ,'get_dcb',' ','numdcb > 32' ,0)
            read(line,'(1x,i2,f9.3)',iostat=ioerr) 
     .               prndcb(numdcb),dcb(numdcb)   
c            print *,'read dcb ',numdcb
            if( ioerr.ne.0 )  call report_stat('FATAL','MODEL'
     .          ,'get_dcb',' ','Error decoding DCB line' ,ioerr)
          else
c           tabular epoch follow obs epoch: exit
            done = .true.  
c            print *,'New epoch > obs epoch, exit '
          endif
        endif     
      
      enddo
                
      return
      end



