Copyright (c) Massachusetts Institute of Technology,1986, 1996. All rights reserved.

      Subroutine GET_SIT_APR    

c     Get a priori constraints for station coordinates

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      real*4 temp(3) 

      character*120 wcmd
      character*16 code
      character*4 snam,upperc,buf4   

      integer i,j,i1,i2,i0,ib,ic,type,ifile,in
      integer lcmd,count_arg,lift_arg


c  Set the defaults to 100 m, but if no input for any station, keep sitwgt = .false.

      do i=1,nsite
         do j=1,3
           stat_apr(i,j) = 100.d0
           stat_apr2(i,j) = 100.d0
         enddo
      enddo

c   Get a priori for tightly constrained solution    

c     Read first the batch file, then the n-file
                
      i2 = 0
      do 160 ifile = 1,2
          
        if( ifile.eq.1 ) then
          in = 5
        else
          in = 13  
          if( nfiln(1:1).eq.' ') goto 160  
        endif               
        call getcmd(in,'aprior',wcmd,lcmd,2)             
        do 150 i0 = 1,1000
           call getcmd(in,'tight_apr_co',wcmd,lcmd,3)   
           if (lcmd.le.0) goto 160  
           if( ifile.eq.2.and.i0.eq.1 ) 
     ,         call report_stat('WARNING','SOLVE','get_sit_apr',' ',
     .  'Coordinate constraints loosened: see APTOL in prefit q-file',0)
           i2 = i2+1
c          decompose the command line
           ic = count_arg(wcmd)
c          not enough arguments
           if (ic.le.1) goto 150
c          pointer to
           ib = lift_arg(wcmd,code,1)
           if (ib.le.0) goto 150
           type = 1
           if (upperc(code(1:4)).eq.'ALL_') type = 2
           if (type.eq.1) snam = upperc(code(1:4))
           read (wcmd(ib+1:lcmd),*) (temp(i),i=1,3) 
           if ( temp(1).eq.0.d0.or.temp(2).eq.0.d0.or.temp(3).eq.0.d0 )
     .        call report_stat('FATAL','SOLVE','get_sit_apr',' ',
     .         'Coordinate tight constraints are zero in batch file',0)
           do 140 i=1,nsite
              i1=3*(i-1)+1
              if(free(i1).eq.0) go to 140
              if (type.eq.1) then
                buf4=upperc(rlabel(i1)(1:4))
                if (buf4.ne.snam) goto 140
              endif
              do j=1,3
                stat_apr(i,j)=temp(j)
              enddo        
  140      continue
  150   continue
  160 continue
      if (i2.gt.0) then
        sitwgt = .true.
        if( logprt ) write(6,161)
        write(10,161)
  161   format(/,
     .    4x,'A priori coordinate errors in meters',/,
     .    'Station                       Latitude  Longitude  Radius',/)
        do 165 i=1,nsite
          if( logprt ) write( 6,162) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                 ,(stat_apr(i,j),j=1,3)
          write(10,162) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                 ,(stat_apr(i,j),j=1,3)
  162     format(i3,2x,a4,1x,a12,3x,3(2x,f9.4))
c         all units are in meters (convert to kms)
          do  j=1,3
             stat_apr(i,j)=stat_apr(i,j)/1000.d0
          enddo
  165   continue
      endif

c      Get a priori for loosely constrained solution

      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 180 i0 = 1,1000
         call getcmd(5,'loose_apr_co',wcmd,lcmd,3)
         if (lcmd.gt.0) then
           i2 = i2+1
c          decompose the command line
           ic = count_arg(wcmd)
c          not enough arguments
           if (ic.le.1) goto 180
c          pointer to
           ib = lift_arg(wcmd,code,1)
           if (ib.le.0) goto 180
           type = 1
           if (upperc(code(1:4)).eq.'ALL_') type = 2
           if (type.eq.1) snam = upperc(code(1:4))
           read (wcmd(ib+1:lcmd),*) (temp(i),i=1,3)    
           if ( temp(1).eq.0.d0.or.temp(2).eq.0.d0.or.temp(3).eq.0.d0 )
     .        call report_stat('FATAL','SOLVE','get_sit_apr',' ',
     .         'Coordinate loose constraints are zero in batch file',0)
           do 170 i=1,nsite
              i1=3*(i-1)+1
              if(free(i1).eq.0) go to 170
              if (type.eq.1) then
                 buf4=upperc(rlabel(i1)(1:4))
                 if (buf4.ne.snam) goto 170
              endif  
c             all units are in meters (convert to kms)
              do  j=1,3
                stat_apr2(i,j)=temp(j)
              enddo 
 170       continue 
         endif   

c########################################################################################
c The following loop is not necessary, and screws up the stat_apr2 array. Causing incorrect
c loose constraints to be applied. Next 3 lines of code must be commented out. TAH/SCM 99/04/02 
c        all units are in meters (convert to kms)
C         do  j=1,3
C            stat_apr2(i,j)=temp(j)/1000.d0
C        enddo
c########################################################################################

 180  continue  

c     all units are in meters (convert to kms)
      do i=1,nsite
        do  j=1,3
           stat_apr2(i,j)=stat_apr2(i,j)/1000.d0
        enddo
      enddo

      return
      end
