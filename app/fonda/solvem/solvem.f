      program solvem
c
c     this program solves observation equations and
c     derives solutions.
c                        Danan Dong at MIT 930710
c
ccc   line added by Gilbert FERHAT 25/03/94
ccc   see line beginning with ccc 
ccc   in order to have a sketch file (file #28)
ccc   in order to have a GMT file (file #29)
c
c     control keys: (key)
c       1. free all parameters (no rank deficiency)
c       2. fix some station or velocity or baseline orientation
c       3. fix network center and no rotation (both location and velocity)
c
c     files (idatf):
c       No     mode       function
c        8 --  input   driving file
c       10 --  input   observation data 
c       13 --  input   priori value 
c       14 --  input   sequential combining file list
c       15 --  input   earthquake correction file
c       16 --  input   refraction correction file
c       18 --  input   subnetwork list for strain rate calculation
c       20 -- output   solutions (coordinate+velocity)
c       21 -- output   velocity mapping
c       22 -- output   residuals
c       23 -- output   strain rate
c       24 -- output   event record file
c       25 -- output   loose constraint h-file
c       26 -- output   modified network coordinate file
ccc     28 -- output   GMT sketch file
ccc     29 -- output   GMT file
ccc     30 --  input   Lambert file 
c
c       31 --  input   sequential input subfile 
c       33 -- output   temporary file to store normal matrix and right hand term
c       34 --  input   inner coordinate site list
c       35 --  input   outer coordinate site list
c       36 --  input   model coordinate site list
c    
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
c     character*64 priori_coord,gmt_sketch_file
      character*16 version
      character*64 drvfil
      character*52 runstr
      integer system,iclarg,imnetf,imnetf2,get_path
      integer i1,i2,ii,ierr,idim,key,isum,nobs,k,npar,i
c      real*8 temp2(400000)
      integer ios,idatf
      integer kobs(maxsit)
      data fixcnt/0/
      data fixsit/maxfix*0/
      common/center/xs,ys,zs,slxy
     
c     set the version number
      call museum(version)
c
c     If the driving file name is missing, the program will give
c     the description only.
      i1 = iclarg(1,drvfil)
      if (i1.le.0) then
         i2 = get_path('FONDA_H',drvfil)
         runstr(1:25+i2) = 
     .   'head -16 ' // drvfil(1:i2) // 'help/solvem.help'
         ierr = system(runstr)
         stop
      endif

      call pline (6,80,'-',2)
      print*, '               SOLVEM: version ',version 
      call pline (6,80,'-',3)
      idim = 3
c
ccc   opening driving file = 8
      open (8,file=drvfil,status='old',iostat=ios)
      if (ios .eq. 0) then
         print *,' drive file name: ',drvfil(1:i1)
      else
         print *,' SOLVEM: error opening drive file: ',drvfil
         call ferror (ios,6)
         stop 'SOLVEM: 1'
      endif
c
c     get all file names and open files
      call readrv(key,1,chi)
ccc   job = 1 :  get input/ouput informations
c     output head part 1
      call wthead(20,key,2,1,version,drvfil)
ccc   20 = output solutions (coordinate+velocity) file
ccc   2  = mode : verbose 
ccc   1  = job : get i/o informations
c
c     get reference frame parameters
      call geotab(frame,1,radius,finv,tx,ty,tz)
      f = 1.0d0/finv
      e2 = 2.0d0*f-f*f
c     print*, 'IN SOLVEM' 
c     print*,'reference frame parameters'
c     print*, 'inv= ',inv
c     print*, 'f= ',f
c     print*, 'e2= ',e2
c 
c     get priori values
c     idatf is the file type - see list above
      idatf = 10
      if (iomode(5).gt.0) idatf = 14
      call getsit(idatf,13,idim)
c
ccc   make GMT sketch file
ccc   iomode(1) is for file #13 = priori value   
      if (iomode(1).gt.0) then
ccc       iomode(5) is for file #14 = sequential combining file list
          if (iomode(5).gt.0) then
              if (iomode(15).gt.0) then 
                   print *,' starting with GMT file'
                   call gmt_sketch(28)
              endif
          endif
ccc
      endif
      if (iomode(16).gt.0) then
         print *,'make a GMT file'
c        call gmt_file(priori_coord,gmt_sketch_file)
         call gmt_file
      endif
      close (29)
      close (13)
ccc
c     get earthquake influenced site list
      if (iomode(3).gt.0) then
         call get_quake_list(15)
      endif
c
c     get refraction corrections
      if (iomode(13).gt.0) then
         call get_refractn(16)
      endif
c
c     prechecking, identify "ill" sites and "ill" observations
      if (idatf.eq.10) call prechk(idatf,timemin,timemax)
      if (idatf.eq.14) call prechk_full(idatf)

      do i = 1,nsit
         kobs(i) = map(i) + map(maxsit+i)
      enddo
c
c     subtract the effective site from whole sites
      call subtra(isum,1)
c
c     initialization
      call initn
c
c     put loose constraints on normal matrix
      call lswght(1)
c     
c     construct the normal matrix
      if (idatf.eq.10) call frmnor(idatf,nobs)
      if (idatf.eq.14) call frmnor_seqe(idatf,nobs)
c
      if (iaux.gt.0.or.jaux.gt.0.or.jeaux.gt.0) call lswght(2)
c
c     get geometry center
c     call geocnt(xs,ys,zs,slxy,slxz,slyz,so2y,1)
c
c     subtract the effective submatrix from whole normal matrix.
 60   call subtra(isum,2)
      print*,' effective site : ',gdsit
c
c     solve the normal equation (by-passed if not output h-file)
      if (iomode(11).gt.0) then 
         call lsq_soln(npar)
c        output global h-file
         print*, ' output h-file'
         call hwrite(25,nobs,version,drvfil,kobs,timemin,timemax)
      endif
c
c     get apriori covariance matrix
      call readrv(key,2,chi)
c
      if (iomode(11).gt.0) then 
         rewind (33)
         read (33) npar
         do i = 1,npar
            i1 = i*(i-1)/2
            read (33) (anorm(i1+k),k=1,i)
         enddo
         read (33) (bnorm(k),k=1,npar)
         close (33)
      endif
      do i = 1,nlive
         i1 = i*(i+1)/2
         if (i.le.3) k = i*(i-1)/2
         if (i.gt.3) k = 3+(i-3)*3
         do 70 ii = 1,3
            if (ii.gt.i) goto 70
            if (i.le.3) i2 = k+ii
            if (i.gt.3) i2 = i1+ii-3
            anorm(i2) = anorm(i2)+aprm(k+ii)
 70      continue
      enddo
      print*,' update solution'
      call lsq_soln(npar)
      chi = chi2
c      goto 65
c
c     get constraint informations
 65   call readrv(key,3,chi)
cmk   is this the correct number of dof - no nlive is not correct 
cmk    in the case of constraints. Need to compute the redundancy
cmk    some other way
c     rmsn = dsqrt(chi/(nobs-nlive))
      write(*,'(2a,i5,i5,f8.2)') ' SOLVEM: after constraints: ',
     .   'nobs,nlive,chi:',nobs,nlive,chi
      close (8)
c
c     copy velocity solution
      print*,' copying velocity solution...'
      call cpsoln(1)
c
c     output results by standard format
      imnetf = 26
      imnetf2 = 37
      call output(20,21,imnetf,imnetf2,key,1,version)
c
c     get residual output
      if (iomode(8).gt.0.and.idatf.eq.10) then
            call residl(22,10,nobs)
      endif
      if (iomode(8).gt.0.and.idatf.eq.14) call res_list(22,14,nobs)
c
c     get strain rate solution
      if (iomode(9).gt.0) call calstr(23,3,iomode(9))
c
      if (iomode(2).gt.0) close (10)
      if (iomode(5).gt.0) close (14)
      close (20)
      close (21)
      if (iomode(8).gt.0) close (22)
      if (iomode(9).gt.0) close (23)
      if (iomode(3).gt.0) close (15)
      if (iomode(10).gt.0) close (24)
      if (iomode(11).gt.0) close (25)
      if (iomode(12).gt.0) close (26)
      goto 100
c
c     troubled stop
 1000 print*,'suicide due to ..'
      stop
c
c     normal stop
 100  print*,' Normal stop. Congratulations!'
      call pline (6,80,'-',4)
      stop
      end
