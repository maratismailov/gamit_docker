C
      SUBROUTINE LOOS12( iop )
C
C     loose constraints for other modes besides LC mode
C     iop = 3: bias free, loose constraint solution
C     iop = 4: bias fixed, loose constraint solution
c                      
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
c
      integer fre2(maxprm),iop
      logical zeroad
      character*256 message

      integer*4 nded,nlive0,ierinv,ifast,k0,k1,k2,k3,imb,imb1
     .        , ifix,i1,ia1,ia2,ij,ib,il,id,i,j
             
c**   For debug
c      integer*4 iold,idiag,j0,ij4
c      real*8 temp(maxprm)
c**

      real*8 rnew,goodft,rcond
 
      equivalence (isigma,fre2)


c      Biases-free solution

      if( iop.eq.3 ) then
     
c---- write a header to the screen and q-file
      write(message,'(a)') 
     .    'Performing L1 and/or L2 biases-free loose solution'
      call report_stat('STATUS','SOLVE','lcloos',' ',message,0)
      write(10,'(/,a,/)') message

c---- store neccesary values
      do i = 1,ntpart
         fre2(i) = free(i)
         adorg(i) = adjust(i)
      enddo
c---- copy normal matrix
      rewind 27
      read (27) ntpart
      k1 = ntpart*(ntpart+1)/2
      call tapdrv(27,1,ntpart,k1,a,5)
      call tapdrv(27,1,ntpart,ntpart,b,4)
c     do we want the original value of r2sum, or the value currently
c     stored in common/slvaux/ ?  Similarly, do we want the original
c     values of nlive and free, or the values stored in in solve.h
c     common /block2/ and explicit common /params/ ?
      read (27) (free(j),j=1,ntpart)
      read (27) r2sum
      read (27) nlive                    
c**   DEBUG
c      write(*,'(a,2i5,e10.3,i4)') 'LOOS12 ntpart k1 r2sum nlive '
c     .     ,ntpart,k1,r2sum,nlive
c      write(*,'(a,20i3)') 'free ',(free(i),i=1,ntpart)     
c      write(21,'(a)') ' Normal equations: '
c      iold = 0
c      do i = 1,nlive
c         idiag = i*(i+1)/2      
c         j0 = 0
c         do j=1,i
c           j0 = j0 + 1
c           ij4 = iold + j
c           temp(j0) = a(ij4)
c         enddo
c         write(*,'(i4,500(5(1x,d11.2),/,6x))') i,(temp(j),j=1,j0)
c         iold = idiag
c      enddo
c** end debug
          
c This now seems to work ok, at for L1-only 
c      write( message,'(a)')
c     .     'Untested loose-solution code for L1/L2 independent'
c      write(10,'(//)')
c      write(10,'(a)') message
c      call report_stat('WARNING','SOLVE','loos12',' ',message,0)
c     Note:  This code has been failing on a unit 27  EOF because free, r2sum,
c            and nlive were not written in LSQUAR.  I fixed LSQUAR but have
c            not tested the loose-solution result.   rwk  930929
c---- restore initial values
      i1 = lpart
      do i = 1,mbias
         i1 = i1+1
         if (free(i1).eq.1.and.idxb(i).lt.0) idxb(i) = iabs(idxb(i))
      enddo

c---- new station constraints
      if ( sitest  ) call lwstat

c---- new satellite constraints
      if ( satest ) call lwsat

c-----zenith-delay weighting
      if ( zenest )  call lwzen

c----- atmospheric gradient constraints. iflag = 1 N/S gradient, = 2 E/W gradient   
      if ( gradest ) then
          call lwgrad 
      endif
  
c-----new earth orientation parameter constraints
      if ( eopest )  call lweop
                            

c**   DEBUG
c      write(*,'(a,2i5,e10.3,i4)') 'Aft rewgt LOOS12 
c      write(*,'(a,20i3)') 'free ',(free(i),i=1,ntpart)     
c      write(21,'(a)') ' Normal equations: '
c      iold = 0
c      do i = 1,nlive
c         idiag = i*(i+1)/2      
c         j0 = 0
c         do j=1,i
c           j0 = j0 + 1
c           ij4 = iold + j
c           temp(j0) = a(ij4)
c         enddo
c         write(*,'(i4,500(5(1x,d11.2),/,6x))') i,(temp(j),j=1,j0)
c         iold = idiag
c      enddo
c** end debug 

c---- reorder normal matrix with live parameters
      zeroad = .true.
      nded = ntpart-nlive
c**** temporary change (good only when adjust = 0)
      ij = 0
      ib = 0
      do 80 i = 1,ntpart
         il = i*(i-1)/2
         if (free(i).eq.0) goto 80
         ib = ib+1
         b(ib) = b(i)
         do 90 j = 1,i
            id = il+j
            if (free(j).eq.0) goto 90
            ij = ij+1
            a(ij) = a(id)
 90      continue
 80   continue
      call nrmscl(adjust,nlive,0,1)
      call report_stat('STATUS','SOLVE','loos12',' '
     .  , 'Solving normal equation for loose constraints',0)
c      call vinv2(3,nlive,ierinv)
      call inver2(a,b,3,nlive,rcond,ierinv)
      call nrmscl(adjust,nlive,1,1)
      IF(IERINV.NE.0)
     .  call report_stat('FATAL','SOLVE','loos12',' '
     .        ,'Bad inverse: two many parameters adjusted',0)
      call zero1d(1,ntpart,adjust)

c---- derive solution
      call solve2( adjust,nded,rnew,zeroad)
c
c---- in the case of x2 = 0, all results have been obtained
       if ( zeroad ) goto 40
c
c---- ifast = 0 means do both chi2 and params
      ifast = 0
      call solve3( adjust,nded,rnew,ifast )
c
c---- compute postfit goodness of fit
 40   goodft = chi2/dble(nobs-nlive)
      goodft = dabs(goodft)
      sclerr = dsqrt(goodft)

c      Biases-fixed solution

      elseif (iop.eq.4 ) then
       
c---- write a header to the screen and q-file
      call report_stat('STATUS','SOLVE','lcloos',' '
     .          , 'Performing L1/L2 biases-fixed loose solution',0) 
      write(10,'(/,a)') 'Performing L1/L2 biases-fixed loose solution'


c---- replace biases by original fixed value
      k0 = idxb(1)
      if (k0.lt.0) k0 = -k0
      k1 = k0
      ifix = 1
      k3 = 0
         k2 = l1bias
         ia1 = k1
         ia2 = k1+k2-1
         k1 = k1+k2*2
         if (l2flag.eq.2) ia2 = ia2+k2
         imb1 = 0
         do 152 i = ia1,ia2
            k3 = k3+1
            if (free(i).eq.1.and.fre2(i).eq.0) then
               imb = msig+imb1
               nfix(ifix) = imb
               bdev(ifix) = adorg(i)-adjust(i)
               adjust(i) = adorg(i)
               free(i) = 0
               idxb(k3) = -iabs(idxb(k3))
               nlive0 = nlive
               nlive = nlive-1
               call bnew( nlive0,ifix )
            endif
            if (free(i).eq.1) imb1 = imb1 + 1
 152     continue
         k3 = k3+k2

      chi2 = sclerr**2*dble(nobs-nlive)

      endif

c-------------------------------------------------------------

      return
      end

