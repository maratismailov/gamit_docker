      Subroutine SCANSP3( nsvs,issat,nepchs,jd,t,x,trun
     .                  , nsvt,itsat,nepcht )

c     Scan an SP3 file to determine which SVs can be used in a t-file,
c     by flagging SVs with no entries.
c      If trun = F, remove the flagged SVs entirely (final orbits)
c      if trun = T, truncate the t-file when any SV is not longer valid (ultrarapids or predictions)

      implicit none               

      include '../includes/dimpar.h'

      integer*4 maxepc
      parameter (maxepc=300)   
  
      logical trun

      integer*4 nsvs,issat(maxorb),nepchs
c        input SP3 number of SVs, array of PRN numbers, epochs

      integer*4 jd
      real*8 t(maxepc),x(3,maxsat,maxepc)
c        input times and coordinates from SP3 file
                       
      integer*4 nsvt,itsat(maxorb),nepcht
c        output number of SVs, array of PRN numbers, epochs for t-file
                                                        

      integer*4 badsvs(maxsat),badepc(maxsat)
c       index in issat array for SVs with missing coordinates;
c       epoch at which coordinates first missing
                
c       Local
      integer*4 isbad,badsv,is,i,j,k                                     
      real*8 test
      character*256 message
                             
c       Function
      logical checksvs

c       Initialization
                    
      nepcht = nepchs
      nsvt = nsvs
      do j=1,nsvt
        itsat(j) = issat(j)
      enddo 
      isbad = 0 
      do j=1,maxorb
        badsvs(j) = 0
        badepc(j) = 0
      enddo    
                                                 
c       Scan the array looking for bad coordinates
           
      do k=1,nepchs     
       do j=1,nsvs
         test = dsqrt(x(1,j,k)**2+x(2,j,k)**2+x(3,j,k)**2)  
         if( test.le.15.d3 .or. test.gt.30.d3 ) then
            if( .not.checksvs(j,maxorb,badsvs) ) then
              isbad = isbad + 1 
              badsvs(isbad) = j
              badepc(isbad) = k    
            endif
          endif
        enddo
      enddo   

c        If any bad SVs found, remove the SV or truncate the file
                      
      if( isbad.ne.0 ) then

c       first remove any SVs bad from the beginning 
        do is=1,isbad
          if( badepc(is).eq.1 ) then    
            do j=badsvs(1),nsvt
              itsat(j)=itsat(j+1)
              do k=1,nepcht
                do i=1,3
                  x(i,j,k) = x(i,j+1,k)
                enddo
              enddo
            enddo
            write(message,'(a,i3,a)') 'PRN ',itsat(badsvs(is))
     .       ,' bad at first epoch not written to t-file '
            call report_stat('WARNING','SP3TOT','orbits/scansp3',' '
     .                      , message,0)
            nsvt = nsvt -1 
          endif
        enddo      

c       then if trun = F, remove any SVs bad at any point      
        if( .not.trun ) then 
          do is=1,isbad  
            if( badepc(is).ne.1 ) then 
              write(message,'(a,i3,a,i4,a)') 'PRN ',itsat(badsvs(is))
     .         ,' bad at epoch ',badepc(is),' not written to t-file '
              call report_stat('WARNING','SP3TOT','orbits/scansp3',' '
     .                         , message,0)
              do j=badsvs(is),nsvt   
                itsat(j)=itsat(j+1)
                do k=1,nepcht
                  do i=1,3
                    x(i,j,k) = x(i,j+1,k)
                  enddo
                enddo
              enddo
              nsvt = nsvt -1 
            endif
          enddo   
c       or if trun = T truncate the file if any SV bad at a later epoch
        else  
c         find the earliest epoch of bad SVs
          nepcht = nepchs
          do i=1,isbad
            if( badepc(i).lt.nepcht) then
              nepcht = badepc(i) - 1
              badsv = badsvs(i) 
            endif
          enddo   
          write(message,'(a,i4,a,i3)') 'T-file stopped at epoch ',nepcht
     .        ,' for bad PRN ',issat(badsv)
          call report_stat('WARNING','SP3TOT','orbits/scansp3',' '
     .       ,message,0)
        endif     
                                                       
      
      endif
      return
      end


c--------------------------------------------------------
      Function checksvs(indx,maxorb,badsvs)

      implicit none

      logical checksvs             
       
      integer*4 maxorb,indx,badsvs(maxorb),i
        
      checksvs = .false.
      do i=1,maxorb
        if(badsvs(i).eq.indx) checksvs = .true.
      enddo                           
      return
      end
      
      
