      Subroutine READ_BIASES( last_nonbias,lbias1,lbias1dd,nlen,ndex )
                               
c     Read double-difference bias combinations from an AUTCLN output
c     and set the pointers for SOLVE.    
c     R. King  26 January 2004     
                      
c     Input:
c        last_non_bias : parameter number of last non-bias parmaeter  
c     Output:
c        lbias1        : number of live double-difference combinations (=nsite*nsat - dead)
c        lbias1dd      : number of viable L1 double-difference biases (n-file w/o missing sites/sats)
c        nlen          : number of independent baselines used for biases
c        ndex          : baseline indices (points to full set of all combinations)

c     The logic for pointers is taken from the D. Dong routine DOPT.
                                   
      implicit none  
                                                          
      include '../includes/dimpar.h'
      include 'solve.h'  
      include 'parameters.h' 

c       Arguments
      integer*4 last_nonbias,lbias1dd

c       Local  

      integer*4 maxbas
      parameter ( maxbas=maxsit*(maxsit-1)/2 )
      integer*4 ndex(maxbas*2)               
      integer*4 ioerr
     .        , ivalue,lcmd,indx,isum,lbias1,nlen,ix1,ix2,i,j,k,jp 
      integer*4 isigmatmp(maxprm),ipnt2dtmp(maxobs)
   
     
      real*8 rvalue

      character*4 sitcod(maxsit),satcod(maxsat),btype,cvalue
     .    ,sit1(maxprm),sit2(maxprm),sat1(maxprm),sat2(maxprm)
      character*120 wcmd  
      character*128 line  
      character*256 message

      logical eof,match,debug/.false./

c       Functions

      integer*4 iarray,ic4array

               
      if(debug) print *,'READ_BIASES last_nonbias ',last_nonbias

c  Read the file  

c     N-file opened and unit assigned in get_err_apr, called by read_bfl_sess.
      rewind(unit=13,iostat=ioerr)    
      if( ioerr.ne.0 )  call report_stat('FATAL','SOLVE','read_biases'
     .        ,' ','Problem rewinding unit 13',ioerr)
c**   debug
c      read(13,'(a)') line
c      print *,'LINE 1',line 
c**   end debug  
c     Rewind and stop on key word        
      call getcmd(13,'biases',wcmd,lcmd,2)  
      if ( lcmd.ne.0 ) 
     .     call report_stat('FATAL','SOLVE','read_biases'
     .        ,nfiln,'Biases not found in  N-file',ioerr)  
      eof = .false.                                           
      k = 0
      do while (.not.eof ) 
        read(13,'(a)',iostat=ioerr) line      
        if( ioerr.eq.-1 ) then 
          eof = .true.
        elseif ( ioerr.ne.0 )  then
          call report_stat('FATAL','SOLVE','read_biases'
     .        ,nfiln,'Error reading bias line from N-file',ioerr) 
        else 
          indx = 2  
c         read the parameter label (should always be 'B1L2')
          call read_line(line,indx,'CH',ioerr,rvalue,btype)  
          if (ioerr.eq.-1 ) then
            eof = .true.
          elseif( ioerr.gt.0 .or. btype.ne.'B1L2' ) then
             call report_stat('FATAL','SOLVE','read_biases'
     .        ,nfiln,'Error reading btype from N-file',ioerr)  
          else  
            if( btype.ne.'B1L2' )  
     .         call report_stat('FATAL','SOLVE','read_biases'
     .         ,nfiln,'Unexpected parameter type on N-file',ioerr)   
c           found an entry, increment the counter and read the sites and sats
            k = k + 1
            call read_line(line,indx,'CH',ioerr,rvalue,sit1(k) )
            call read_line(line,indx,'CH',ioerr,rvalue,sit2(k) )  
            if ( ioerr.ne.0  )
     .         call report_stat('FATAL','SOLVE','read_biases'
     .           ,nfiln,'Error reading sites from N-file',ioerr)   
            call read_line(line,indx,'I4',ioerr,ivalue,cvalue )  
            sat1(k)(1:2) = 'PN'
            write(sat1(k)(3:4),'(i2)') ivalue       
c*            if( sat1(k)(3:3).eq.' ' ) sat1(k)(3:3) = '0'
            call read_line(line,indx,'I4',ioerr,ivalue,cvalue )   
            if ( ioerr.ne.0  )
     .        call report_stat('FATAL','SOLVE','read_biases'
     .         ,nfiln,'Error reading satellites from N-file',ioerr)   
            sat2(k)(1:2) = 'PN'
            write(sat2(k)(3:4),'(i2)') ivalue 
c*           if( sat2(k)(3:3).eq.' ' ) sat2(k)(3:3) = '0'   
c           read the integer ambiguity (normally zero), stored as real  
C            call read_line(line,indx,'I4',ioerr,ivalue,cvalue )
C            wlval(k) = dfloat(ivalue)
* MOD TAH 050713: Read WideLane as real*8 to support codeless receivers
            call read_line(line,indx,'R8',ioerr,wlval(k),cvalue )  
c           read the confidence level 
            call read_line(line,indx,'R8',ioerr,wlconf(k),cvalue ) 
c           read whether free or fixed from autcln
            call read_line(line,indx,'CH',ioerr,rvalue,rorx(k) )   
            if ( ioerr.ne.0  )
     .        call report_stat('FATAL','SOLVE','read_biases',nfiln,
     .       'Error reading bias, confidence or flag from N-file',ioerr) 
          endif
        endif                      
      enddo 
c     set number of biases read from n-file; save in common /acbiases for get_widelane
      numwl = k 
c     (number of viable DD biases set later to lbias1dd and saved by DOPT in common /bbii/ 
      if(debug) then 
        print *,'READ_BIASES numwl',numwl
        do i=1,numwl
         write(*,'(i3,4(1x,a4),1x,a1)')  
     .        i,sit1(i),sit2(i),sat1(i),sat2(i),rorx(i)
        enddo
      endif

c  Get the site/sat indices

c     temporary(?): recover the site, sat list from the global one-way bias labels
c** begin debug
c      print *,'last_nonbias nsat rlabel ',last_nonbias, nsat 
c      j= last_nonbias +1
c      k= last_nonbias + nsat
c      print *,'sat start stop ',j,k 
c      write(*,'(4(1x,a20))') (rlabel(i),i=j,k) 
c** end debug
      k = 0
      do i= last_nonbias + 1, last_nonbias + nsat
        k = k + 1
        satcod(k) = rlabel(i)(1:4) 
      enddo   
c** begin debug       
c      j= last_nonbias +1    
c      k = last_nonbias+nsat*nsite-nsat+1
c      print *,'sit start stop ',j,k 
c      write(*,'(4(1x,a20))') (rlabel(i),i=j,k)    
c** end debug
      k = 0         
      do i = last_nonbias+1, last_nonbias+nsat*nsite-nsat+1, nsat   
        k = k + 1
        sitcod(k) = rlabel(i)(5:8)
      enddo  
c      print *,'satcod ',(satcod(i),i=1,7)
c      print *,'sitcod ',(sitcod(i),i=1,7)                  
                   

c  See which site/sat combinations are not available  
     
       
c     Use the iuse array of included observations, then check this
c     against the biases assigned by AUTCLN.

c      print *,' READ_BIASES nsite nsat iuse at start: ',nsite,nsat
c      write(*,'(7i7)') ((iuse(i,j),j=1,nsat),i=1,nsite)  
      isum = 0         
      ix1 = 0  
      call zero1i(1,maxprm,isigma)
      do i =1,nsite
        do j=1,nsat  
          isum = isum + 1   
          isigma(isum) = 0 
          if( iuse(i,j).gt.0 ) then
            ix1 = ix1 + 1 
            isigma(isum) = ix1
            ipnt2d(ix1) = isum  
          endif
        enddo
      enddo  
      lbias1 = ix1  

      isum = 0         
      ix1 = 0  
      call zero1i(1,maxprm,isigmatmp)
      do i =1,nsite
        do j=1,nsat  
          isum = isum + 1
          match = .false.
          do k = 1, numwl  
c*           print *,'i sitcod j satcod k sit/sat ',
c*c     .        i,sitcod(i),j,satcod(j),k,sit1(k),sit2(k),sat1(k),sat2(k)
            if( (sitcod(i).eq.sit1(k).or.sitcod(i).eq.sit2(k)) .and.
     .          (satcod(j).eq.sat1(k).or.satcod(j).eq.sat2(k))
     .           ) then
              if( .not.match ) then
c                print *,'match i j k ix1 ',i,j,k,ix1
                ix1 = ix1 + 1 
                isigmatmp(isum) = ix1
                ipnt2dtmp(ix1) = isum      
                match = .true.
c                print *,'ix1 isigmatmp ipnt2dtmp '
c     .          ,ix1,isigmatmp(ix1),ipnt2dtmp(ix1)
              endif
            endif
          enddo
        enddo
      enddo       
c      print *,'ix1 lbias1 ',ix1,lbias1
      if (ix1.ne.lbias1 )  then
        if( ix1.eq.0 ) then
          call report_stat('FATAL','SOLVE','read_biases',' ',
     .   'Zero WL biases read from N-file, check processing sequence',0)
       else
          write(message,'(a,i4,a,i4,a)') 
     .    'Number of 1-way biases from AUTCLN (',ix1
     .  ,') does not match number from site/sat combinations with obs ('
     .  ,lbias1,')'   
         call report_stat('WARNING','SOLVE','read_biases',' ',message,0)
          print *,' isigma ', lbias1   
          write(*,'(10i7)') (isigma(i),i=1,nsite*nsat)  
         print *,' isigma AUTCLN ',ix1
         write(*,'(10i7)') (isigmatmp(i),i=1,nsite*nsat)  
        endif
      endif                
 

c**  print out what we have so far

c      print *,'READ_BIAS nsite sitcod ',nsite  
c      write(*,'(10a4)') (sitcod(i),i=1,nsite)
c      print *,' nsat satcod ',nsat
c      write(*,'(10a4)') (satcod(i),i=1,nsat)  
c      print *,' isigma '
c      write(*,'(10i7)') (isigma(i),i=1,50)
c      print *,' ipnt2d '
c      write(*,'(10i7)'),(ipnt2d(i),i=1,50)
                       

c  Fill the site/sat index arrays
       
c  ipntd (# L1 or L2 one-way biases), giving DD combinations in groups of four
c     e.g.  9  30  8  29   2  30 1 29  where
c      9 30 8 29 are the one-way bias indices for 
c      sit1/sat1  sit2/sat1   sit1/sat2  sit2/sat2
c     do this by matching the site/sat names in the lists   

c   In the following code, 'i' is the index over the DD biases (L1 or L2)
c   'j' is the index over one-way biases, pointing to the original parameter labels
c   'ipnt2d' gives the position in the original parameter array (relative to the 
c      last non-bias parameter of the live biases; i.e. ipnt2d(1-3) = 1 3 4 if 
c      1-way bias #2 is omitted because there were no observations to that site/sat 
c      combination.  
c   'k' is the index to the elements of the ipntd array, four times the size of ipnt2d.
                   

c      print *,' Fill ipntd--mbias ',mbias  
      k = 0     
      lbias1dd = numwl 
      do i = 1, numwl    

        miss_bias(i) = .false.
        
c       set sit1, sat1 index (1st of 4)
        do j = 1, mbias/2  
          jp = last_nonbias + j
c        print *,'i sit1 sat1 jp rlabel ',i,sit1(i),sat1(i),jp,rlabel(jp)
          if( sit1(i).eq.rlabel(jp)(5:8) .and.
     .        sat1(i)(3:4).eq.rlabel(jp)(3:4) ) then 
c        print *,'i sit1 sat1 jp rlabel ',i,sit1(i),sat1(i),jp,rlabel(jp)
            k=k+1 
            ipntd(k) = iarray(j,ipnt2d,mbias/2)  
c           if no observation for this site or sat, set flag to skip this DD bias
            if( ipntd(k).eq.0 ) miss_bias(i) = .true.   
c            print *,'k ipntd ',k,ipntd(k) 
          endif                          
        enddo   
c       set site2 sat1 index (2nd of 4)
        do j = 1, mbias/2   
          jp = last_nonbias + j
c        print *,'i sit2 sat1 jp rlabel ',i,sit2(i),sat1(i),jp,rlabel(jp)
          if( sit2(i).eq.rlabel(jp)(5:8) .and.
     .        sat1(i)(3:4).eq.rlabel(jp)(3:4) ) then   
c        print *,'i sit2 sat1 jp rlabel ',i,sit2(i),sat1(i),jp,rlabel(jp)
             k=k+1       
             ipntd(k) = iarray(j,ipnt2d,mbias/2)   
            if( ipntd(k).eq.0 ) miss_bias(i) =.true.
c             print *,'k ipntd ',k,ipntd(k) 
          endif
        enddo 
c       set site1 sat2 index (3rd of 4)
        do j = 1, mbias/2       
          jp = last_nonbias + j
c        print *,'i sit1 sat2 jp rlabel ',i,sit1(i),sat2(i),jp,rlabel(jp)
          if( sit1(i).eq.rlabel(jp)(5:8) .and.
     .        sat2(i)(3:4).eq.rlabel(jp)(3:4) ) then   
c        print *,'i sit1 sat2 jp rlabel ',i,sit1(i),sat2(i),jp,rlabel(jp)
             k=k+1 
             ipntd(k) = iarray(j,ipnt2d,mbias/2) 
             if( ipntd(k).eq.0 ) miss_bias(i) =.true.
c             print *,'k ipntd ',k,ipntd(k) 
          endif
        enddo 
c       set site2 sat2 index (4th of 4)
        do j = 1, mbias/2   
          jp = last_nonbias + j
c        print *,'i sit2 sat2 jp rlabel ',i,sit2(i),sat2(i),jp,rlabel(jp)
          if( sit2(i).eq.rlabel(jp)(5:8) .and.
     .        sat2(i)(3:4).eq.rlabel(jp)(3:4) ) then 
c        print *,'i sit2 sat2 jp rlabel ',i,sit2(i),sat2(i),jp,rlabel(jp)
             k=k+1
             ipntd(k) = iarray(j,ipnt2d,mbias/2)    
            if( ipntd(k).eq.0 ) miss_bias(i) =.true.
c             print *,'k ipntd ',k,ipntd(k) 
          endif
        enddo 

c       if bias missing (no obs on a site or sat), remove this quadruplet 
        if( miss_bias(i) ) then   
          lbias1dd = lbias1dd - 1
          do j=1,4
            ipntd(k+1-j) = 0  
          enddo
          k=k-4 
        endif
      enddo   
               
c      print *,' lbias1dd k ipntd ',lbias1dd,k
c      write(*,'(10i7)') (ipntd(i),i=1,lbias1dd*4)


c  Fill the baseline index arrays  


c     ndex   = list of all independent baselines, sorted by length 
c              e.g., 1 5  2 6   means BL 1 is site1/site5, BL 2 is site2/site6 
c              This array is filled in DOPT because the full list is needed
c              in forming the normal equations for the DD biases from the one-way
c              biases; that is, there seems to be an implicit order of independent
c              baselines that includes all stations.
c     ipntdt = pointers from the bias mapping array ipntd to the BLs in ndex
c              e.g.  1 3 means that the first four elements of ipntd refer to BL 1,
c              and the second four refer to BL 3.  Since some BLs are not formed
c              for biases, all values of ndex will not appear in ipntdt.  
                                               
      call zero1i(1,maxdd,ipntdt)   
      k = 0
      do i=1,numwl  
        do j=1,nlen*2,2   
          ix1 = ic4array(sit1(i),sitcod,nsite)
          ix2 = ic4array(sit2(i),sitcod,nsite)  
c          print *,'j ix1 ndex ix2 ndex1 '
c     .           ,j,ix1,ndex(j),ix2,ndex(j+1)
          if( ndex(j).eq.ix1.and.ndex(j+1).eq.ix2 .and.
     .       .not.miss_bias(i) ) then  
            k = k + 1
            ipntdt(k) = j/2 + 1    
c           shift the N-file values to skip missed biases
            rorx(k) = rorx(i)
            wlval(k) = wlval(i)
            wlconf(k) = wlconf(i)
          endif
        enddo
      enddo                                      
c      print *,' ipntdt '
c      write(*,'(10i7)'),(ipntdt(i),i=1,140)
               

c        irowd    operator matrix element per row counter
c          irowd(1)  = total number of non-zero element
c          irowd(2)  = 0 (always)
c          irowd(3+) = running count of the number of elements
c                      per row of D-matrix
c          irowd(2+) set in DOUBLE ??
        
      irowd(1) = lbias1dd
      irowd(2) = 0             
      irowd(3) = 0 
      do i=3,lbias1dd+2
         irowd(i)= irowd(i-1) + 4 
      enddo
c      print *,' irowd '
c      write(*,'(10i7)') (irowd(i),i=1,50)

c        irowdt counts the number of times each site/SV combination appears in the
c        double difference mapping pointers ipntd.  For the non-N-file case, this
c        is done in DOPT with a call to FILLD.  In both cases, DOPT replaces the
c        elements by their acculumated sums.    
c           irowdt(1) = initially 0, then after accumulated, equal to the number of 
c                       number of live elements = nsite * nsat  less any combinations 
c                       missing observations        
c           irowdt(2) = always zero
c           irowdt(3 to lbias1+1) = number of times a site/sat combinations appears
                  
      call zero1i(1,maxnd,irowdt)   
      do i=1,lbias1
        do j=1,lbias1dd*4   
c        print *,'i ipntd(j) irowdt(i+2)',i,ipntd(j),irowdt(i+2) 
          if( ipntd(j).eq.i ) then
            irowdt(i+2) = irowdt(i+2) + 1  
c            print *,'i jp irowdt(i+2)',i,jp,irowdt(i+2) 
          endif 
        enddo     
      enddo   
c      print *,' irowdt '
c      write(*,'(10i7)') (irowdt(i),i=1,50)  

c         dt is the double difference operator, given always by four elements,
c         (-1 1 1 -1 ) for each double-difference set, = 4 x #dd biases
c         The array dt, which points by column, is set in DOPT

      do i = 1 , lbias1dd*4-1, 4
        d(i)   = -1.d0
        d(i+1) =  1.d0
        d(i+2) =  1.d0
        d(i+3) = -1.d0 
      enddo             
c      print *,'d '  
c      write(*,'(10f7.3)') (d(i),i=1,150)    
     
      return
      end



