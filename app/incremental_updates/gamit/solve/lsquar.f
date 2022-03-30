Copyright 1995 Massachusetts Institute of Technology and The Regents of the
c              University of California, San Diego.  All Rights Reserved.

      Subroutine LSQUAR

c Perform the least squares solution
c
c flow chart :
c
c     LSQUAR
c        |                             |
c        +---> (initialization)        |
c        |                             |
c        +--->  WSTAT                  |        } station constraints
c        |                             |
c        +--->  WSAT                   |        } satellite constraints
c        |                             |
c        +---> RMJUNK                  |        } remove junk sit/sat para.
c        |                             |
c        +---> (recording whole normal matrix)
c        |                             |
c        +---> QHEAD    <--------------+        } head output
c        |
c        +---> GETSLT                           } least square solution
c

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
             
      integer*4 icf,ifile,mbias0,last_nonbias,last_bias
     .        , ioerr,ierinv,iirow,ia1,ij4,ii1,ii2,ji,ist,ii,jj,j
                                                              

     
      real*8 r2lc
 
      character*16 lowerc
      character*80 message
      
      logical debug/.false./
      integer*4 i 

* TEST CODE DECLARATION
C     integer*4 nb1_start  ! Start of L1 bias parameters
C    .,         nb1_end    ! End of L1 Bias parameters


c     option 2 doesn't need epoch loop and session loop
c
c     if (iforce.eq.2) then
c        if (l2flag.eq.1) call wlmode(free,1)
c        goto 60
c     endif
c                 c     initialization

c     normal matrix
      ij4 = ntpart*(ntpart+1)/2
      call zero1d(1,ij4,a)
c lpart is =0  at this point, initialize instead by maximum dimension
c      ij4 = lpart*(lpart+1)/2
      ij4 = maxcom*(maxcom+1)/2
      call zero1d(1,ij4,alc)
c     right hand term
      call zero1d(1,ntpart,b)
      call zero1d(1,ntpart,blc)
c     satellite observation counter
      call zero1i(1,maxsat,iseen)
c     site observation counter
      call zero1i(1,nsite,jusit)

      nones = 0
      nobs = 0
      nd1obs = 0
      nd2obs = 0
      r2sum = 0.d0

c     apply station coordinate/satellite constraints

      if( sitwgt ) call wstat
      if(debug) print *,'LSQUAR aft wstat A ',(a(i),i=1,3)

      if( satwgt ) call wsat
      if(debug) print *,'LSQUAR aft wsat A ',(a(i),i=1,3)

* TEST CODE TAH 050428:
*     Get the parameter number range for L1 biases.  These will have
*     apiori weight defined. 7400 is end of L1; 11900 end of L2-L1
C     nb1_start = 0
C     do i = 1, ntpart     
c RWK 190429: change slot range for new code
C        if( islot1(i).ge.2901 .and. islot1(i).le.11900 ) then
C           if( nb1_start.eq.0 ) nb1_start = i
C           nb1_end = i
C           ij4 = i*(i+1)/2
C           a(ij4) = 1.d-4    ! Apriori sigama is a 100 cycles
C        end if
C     end do      
     
cd     if(debug)  print *,'LSQUAR BL1 Bias parameters ',nb1_start,' to '
cd    .   ,nb1_end, 'Tot and com ',ntpart,maxcom
* END TEST CODE
c     c-file counter 


c----------------------------------------------------------------------      _
c

c     Initialize the effective usage index array
c
      do  ii = 1,nsite 
        do jj = 1,nsat
          iuse(ii,jj) = 0
        enddo
      enddo
         
c     option to delete (stations and/or satellites)

c           new batch-file style  
         call read_bfl_sess

c     options to weight zenith delay parameters

         if( zenest ) call wzen 
         if(debug) print *,'LSQUAR aft wzen A ',(a(i),i=1,3)

c     apply n/s (1) and e/w (2) atmospheric gradient constraints

         if( gradest ) then
            call wgrad
         endif            
         if(debug) print *,'LSQUAR aft wgrad A ',(a(i),i=1,3)

               
c     apply polar motion and ut1 constraints

         if( eopest ) call weop    
         if(debug) print *,'LSQUAR aft weop A ',(a(i),i=1,3)


c     get header series information from c- and m-files

         call report_stat('STATUS','SOLVE','lsquar',' '
     .                   , 'Reading C-file headers',0)
         ifile = 30
         do icf = 1,nsite
            ifile = ifile+1
            call copens (lowerc(cfiln(icf)),'old',ifile,ioerr)
            if (ioerr.ne.0) then
               call report_stat('FATAL','SOLVE','lsquar',cfiln(icf)
     .                          ,'Unable to open C-file',ioerr)
            else
              write(message,'(a,a16,a,i3)') 'Opened ',cfiln(icf)
     .          ,'on unit ',ifile
              call report_stat('STATUS','SOLVE','lsquar',' '
     .                    ,message,0) 
            endif
            if( logprt ) write(6,'(/,a,a16)') 
     .        'Opening c-file >: ',cfiln(icf)
            call gethed( ifile,icf,islot2,isprn )
         enddo
                     
cd         if(debug) then 
cd           print *,'LSQUAR aft gethed cfiln ',(cfiln(i),i=1,nsite)
cd            print *, 'ntpart rlabel ',ntpart 
cd           do i=1,ntpart
cd             print *,i,rlabel(i)
cd           enddo 
cd         endif 
         
c     current code does not allow atmospheric (spatial) constraints if
c     there is more than one zenith delay parameter per site.

         if( iatcon.eq.1 .and. nzen.gt.1 ) then
            write(message,'(3a)')
     .           'Atmospheric spatial constraints '
     .         , 'cannot be used with multiple zenith delays--'
     .         , 'constraints turned off'
c            call report_stat('WARNING','SOLVE','lsquar',' ',message,0)
            iatcon = 2
         endif

c     check bias parameter
c
         call bcheck( mbias0,last_nonbias )  
         if(debug) print *,'LSQUAR aft bcheck A ',(a(i),i=1,3)

                 
c     read data records from c-file and increment normal equations
c                          
          call normd( islot1,islot2,r2sum,coords,isprn,icf,rlabel ) 
          r2lc = bn22(1)

c     close the c-files 
c
      do ii = 1,nsite
        close(unit = ii+30)
      enddo

c     record r2sum on scratch file for loose solutions
            write(28) r2lc

         if (iband.gt.0) then
  
c     set up mapping D-operator
c                    
            if(debug) print *,'LDQUAR before DOPT nlive last_nonbias '
     .         ,nlive,last_nonbias 
            call report_stat('STATUS','SOLVE','lsquar',' '
     .          , 'Setting up mapping operator for bias parameters',0)
            call dopt( mbias0,last_nonbias)
            if(debug)  
     .        print *,'LSQUAR after DOPT nlive last_nonbias irowd'
     .             ,nlive,last_nonbias,irowd(1)
            ia1 = last_nonbias
            call report_stat('STATUS','SOLVE','lsquar',' '
     .              , 'Calculating new normal equation submatrices',0) 
            call dbias( mbias0,last_nonbias,last_bias )  
            if(debug) print *,
     .        'LSQUAR after DBIAS ia1 last_nonbias last_bias iband'
     .                     ,ia1,last_nonbias,last_bias,iband 

c     Store the WL-free solution  
c     (do this now because the AN22 matrix is used by FNDDBI)
                 
             
            iirow = l1bias
            write(28) iirow
c     write LC mode N21
            ij4 = iirow*lpart
            call tapdrv(28,iirow,lpart,ij4,clc,3)    
            if(debug) print *,'Writing CLC to 28 ',iirow,lpart,ij4 
c     write LC mode N22
            ij4 = iirow*(iirow+1)/2    
            call tapdrv(28,1,iirow,ij4,an22,2)          
            if(debug) print *,'Writing AN22 to 28 ','1',iirow,ij4 


c     Check if there are dependent bias parameters
                 
            call report_stat('STATUS','SOLVE','lsquar',' '
     .                , 'Finding and removing dependent biases',0)    
      if(debug) print *,'calling FNDDBI ia1 last_bias nlive '
     .                                 ,ia1,last_bias,nlive
cd            print *,'free ',(free(i),i=1,1400)
            call fnddbi(free,ia1,last_bias) 
      if(debug)  print *,'After FNDDBI nlive ',nlive
cd            print *,'free ',(free(i),i=1,1400)
             
c     Read and the WL-free solution to restore AN22, then weight the 
c     bias parameters and re-write
            
            rewind(28)  
            read(28) r2lc
            iirow = l1bias
            read(28) iirow
c     read LC mode N21
            ij4 = iirow*lpart   
            call tapdrv(28,iirow,lpart,ij4,clc,6)   
c            print *,'Reading CLC from 28 ',iirow,lpart,ij4 
c     read LC mode N22
            ij4 = iirow*(iirow+1)/2    
            call tapdrv(28,1,iirow,ij4,an22,5)          
c            print *,'Reading AN22 from 28 ','1',iirow,ij4 
c           print *,'LSQUAR before WBIAS'     
c           call printa(last_bias) 
c           call printan22(l1bias)  
            if( bias_apr.gt.0.d0 ) then    
              write(message,'(a,f6.1,a)') 
     .          'Applying a priori ',bias_apr,' cyc sigma on biases'
              call report_stat('STATUS','SOLVE','lsquar',' ',message,0)
              call wbias( ia1,last_bias,iband )   
              if(debug) print *,'LSQUAR aft WBIAS last_bias iband nlive'
     .            ,last_bias,iband,nlive
c             call printa(last_bias) 
c             call printan22(l1bias)                    
            endif
            rewind(28)
            write(28) r2lc 
            iirow = l1bias
            write(28) iirow
c     write LC mode N21
            ij4 = iirow*lpart
            call tapdrv(28,iirow,lpart,ij4,clc,3)   
c            print *,'Writing CLC to 28 ',iirow,lpart,ij4 
c     write LC mode N22
            ij4 = iirow*(iirow+1)/2    
            call tapdrv(28,1,iirow,ij4,an22,2)          
c            print *,'Writing AN22 to 28 ','1',iirow,ij4 

c      end if for iband
       endif
                  
c      Remove atmospheric parameters for site with zero observation
       do ii = 1,nsite
c rwk 070413: change this since multisession no longer supported
c            ii1 = istid(ii)     
            ii1 = ii
            ii2 = 0
            do ji = 1,nsat
               ii2 = ii2+iuse(ii,ji)
            enddo
            if (ii2.eq.0 .and. zenpar )
     .        call rmjunk(jusit,free,ii,2)
            jusit(ii1) = jusit(ii1)+ii2
            if(debug) print *,'LSQUAR aft rmjunk mode2 ii1 jusit nlive '
     .          ,ii1,jusit(ii1),nlive 
       enddo

* TEST CODE TAH 050428:
*     Get the parameter number range for L1 biases.  These will have
*     apiori weight defined 7400 is end of L1; 11900 end of L2-L1
C     nb1_start = 0
C     do i = 1, ntpart
C        if( islot1(i).ge.2901 .and. islot1(i).le.11900 ) then
C           if( nb1_start.eq.0 ) nb1_start = i
C           nb1_end = i
C           ij4 = i*(i+1)/2
C           print *,i,islot1(i),a(ij4)
C           a(ij4) = a(ij4)+ 1.d-4    ! Apriori sigama is a 100 cycles
C        end if
C     end do
cd      if(debug) print *,'BL1 Bias parameters ',nb1_start,' to ',nb1_end, 
cd     .       ' Tot and com ',ntpart,maxcom

c     close the M-file  
      close(unit=11)

c     change NTPART to save space
      ntpart = last_bias

c     if there is any site or sat no observation, remove corresponding parameters

      call rmjunk( jusit,free,ist,1 )
      if(debug)  print *,'LSQUAR after rmjunk mode1 nlive msig '
     .   ,nlive,msig  

c      call printa(last_bias)

c   Store the total bias-free solution
      rewind 27
      write(27) ntpart
      ij4  =  ntpart*(ntpart+1)/2
      call tapdrv(27,1,ntpart,ij4,a,2)
c      print *,'Writing A to 27 ','1',ntpart,ij4
      call tapdrv(27,1,ntpart,ntpart,b,1)      
      write (27) (free(j),j=1,ntpart) 
      write (27) r2sum
      write (27) nlive
      rewind 27  
      if(debug) print *,'LSQUAR storing bias-free nlive =',nlive


c     print effective observations

      call qhead1(1)

      nobs = nobs+nd2obs
      if(l2flag.lt.2) go to 999
      nobs = 2*nobs
      do 50 ii = 1,nsat
         iseen(ii) = 2*iseen(ii)
 50   continue
c
 999  continue
c     skip solution if single-differences (delay this step until
c     double differences)
c     if(indx.eq.1) go to 200
c
c  this a relic from old code:
c 60   continue
c
c     print parameter number information
c
      call qhead1(2)
c
c     solve normal equation to get least squares solution
c     and postfit goodness

      call report_stat('STATUS','SOLVE','lsquar',' '
     .          , 'Solving initial normal equations',0)
      call getslt( ierinv )
      call report_stat('STATUS','SOLVE','lsquar',' '
     .          , 'Finished solving initial normal equations',0)
                                                 
      if( ierinv.ne. 0)  call report_stat('WARNING','SOLVE','lsquar',' '
     .    ,'Bad inversion of normal equations',0)
      if(debug) print *,'LSQUAR after GETSLT last_bias ',last_bias 
c       call printa( last_bias )
 
c this a relic from single difference code
c 200  continue
c
      return
      end

c*** DEBUG routine

      subroutine printa ( n )

c     Print in rows of 2 x 20 the full current normal equations --diagonal element
c     and the adjacent off-diagonal element


      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
     
      integer*4 i,j,jj,id(10),id1(10),n
      
      j= 0
      do i=1,n 
        j = j +1
        if( j.gt.10 ) stop 10
        id(j) = i*(i+1)/2 
        id1(j) = id(j) -1   
        if( id(j).gt.maxnrm ) stop 20
        if( mod(i,10).eq.0)  then     
          write(*,'(i4,a,20e10.2)') i-9,'A   '
     .          ,(a(id(jj)),a(id1(jj)),jj=1,10)   
          j=0 
        endif
      enddo       
      return
      end


      subroutine printan22 ( n )

c     Print in rows of 2 x 20 the AN22 matrix --diagonal element
c     and the adjacent off-diagonal element


      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
     
      integer*4 i,j,jj,id(10),id1(10),n
      
      j= 0
      do i=1,n 
        j = j +1
        if( j.gt.10 ) stop 10
        id(j) = i*(i+1)/2 
        id1(j) = id(j) -1   
        if( id(j).gt.maxnrm ) stop 20
        if( mod(i,10).eq.0)  then     
          write(*,'(i4,a,20e10.2)') i-9,'AN22'
     .          ,(an22(id(jj)),an22(id1(jj)),jj=1,10)   
          j=0
        endif
      enddo       
      return
      end



       
      
