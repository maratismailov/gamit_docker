      subroutine readrv(key,job,chi)
c
c     get all control information from the driving file
c
ccc   modified by Gilbert FERHAT in order to add the file #28
ccc   =  GMT sketch file
ccc   version 21 03 94 
ccc   line added are beginnig with cccc
c
c     job = 1:   get i/o file name and general informations
c     job = 2:   get apriori coor. and velo. informations
c     job = 3:   get constraint informations
c
c     input/output mode: (iomode)
c       1. priori coordinate 
c       2. observation data
c       3. seismic deformation correction 
c       4. model constraint 
c       5. sequential combining
c       6. output solution
c       7. output velocity map
c       8. residuals
c       9. strain rate
c      10. event record
c      11. loose constraint h-file
c      12. modified network coordinate 
c      13. refraction correction
c      14. gamma rate
ccc    15. GMT sketch ?
ccc    16. GMT file
ccc    17. Lambert file
c      18. Velocity units 
c
c
c     solution combination mode: (smode)
c       1. Gauss-Markov model
c       2. Gauss-Helmert model
c       3. model coordinate approach
c       4. sequential updating model
c       5. Kalman filtering model
c     1. or 3. (?) are the only currently implemented!
c
c     parameter mode: (pmode)
c       1. geocentric (Cartesian)
c       2. geodetic (local Cartesian)
c       3. geodetic (local lon, lat and vertical similar to DYNAP)
c
c     control keys: (key)
c       1. free all parameters (no rank deficiency)
c       2. fix some station or velocity or baseline orientation
c       3. fix network center and no rotation (both location and velocity)
c
c     files:
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
c       26 -- output   modified network coordinate file uncorrelated 
c       27 -- output   gamma rate file
ccc     28 -- output   GMT sketch file
ccc     29 -- output   GMT file
ccc     30 -- input    a priori value LAMBERT coord.
c
c       31 --  input   sequential input subfile 
c       33 -- output   temporary file to store normal matrix and right hand term
c       34 --  input   inner coordinate site list
c       35 --  input   outer coordinate site list
c       36 --  input   model coordinate site list
c       37 --  output  modified network coordinate file correlated
 
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

c
c     character*64 priori_coord, gmt_sketch_file
      character*16 code, upperc
      character*120 wcmd
      integer job,i,i1,lcmd,iloop,id_sit,j,j1
      integer isx,ib,ic,ie,lift_arg
      integer match_name,key,mchkey,nblen
      dimension apr_val(6)
      common/center/xs,ys,zs,slxy
c     initialise array

      go to (50,100,150) job
c
c     begin set up of files
 50   continue

c     priori coordinate and velocity file
      call getcmd(8,'priori',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(1) = 1
         open (13,file=wcmd(1:lcmd),status='old',err=1000)
         print*, ' coordinate file name: ',wcmd(1:lcmd)    
         priori_coord = wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(1) = 0
      endif
c
c     input data file 
      call getcmd(8,'input',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(2) = 1
         open (10,file=wcmd(1:lcmd),status='old',err=1000)
         print*, ' input data file name: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(2) = 0
      endif
c
c     correction file name
      call getcmd(8,'correction',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(3) = 1
         open (15,file=wcmd(1:lcmd),status='old',err=1000)
         print*, ' correction file name: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(3) = 0
      endif
c
c     sequential combining file name
      call getcmd(8,'seque',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(5) = 1
         open (14,file=wcmd(1:lcmd),status='old',err=1000)
         print*, ' sequential file name: ',wcmd(1:lcmd)
      else
         iomode(5) = 0
      endif
c
c     output file 
      call getcmd(8,'output',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(6) = 1
         open (20,file=wcmd(1:lcmd),status='unknown',err=1000)
         print*, ' output file name: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(6) = 0
      endif

c     mapping file name
      call getcmd(8,'mapping',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(7) = 1
         open (21,file=wcmd(1:lcmd),status='unknown',err=1000)
         print*, ' mapping file name: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(7) = 0
      endif
c
c     residual file
      call getcmd(8,'residu',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(8) = 1
         open (22,file=wcmd(1:lcmd),status='unknown',err=1000)
         print*, ' residual file name: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(8) = 0
      endif
c
c     strain rate file
      call getcmd(8,'strain_file',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         ib = lift_arg(wcmd,code,1)
         open (23,file=code(1:ib),status='unknown',err=1000)
         print*, ' strain rate file name: ',code(1:ib)
         call blank(wcmd)
         call getcmd(8,'strain_list',wcmd,lcmd,1)
         if (lcmd.gt.0) then
            iomode(9) = 1
            ib = lift_arg(wcmd,code,1)
            open (18,file=code(1:ib),status='old',err=1000)
            print*, ' strain site list file: ',code(1:ib)
         else
            iomode(9) = 0
         endif
      else 
         iomode(9) = 0
      endif
c
ccc   GMT sketch file
c
ccc   10 lines added below   
      call getcmd(8,'GMT sketch file',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(15) = 1        
         open (28,file=wcmd(1:lcmd),status='unknown',err=1000)
         print*, ' GMT sketch file: ',wcmd(1:lcmd)
         gmt_sketch_file = wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(15) = 0
      endif
ccc
c
ccc   GMT  file
c
ccc   10 lines added below   
      call getcmd(8,'GMT file',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(16) = 1         
         open (29,file=wcmd(1:lcmd),status='unknown',err=1000)
         print*, ' GMT file: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(16) = 0
      endif
c
ccc   10 lines added below   
      call getcmd(8,'LAMBERT file',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(17) = 1         
         open (30,file=wcmd(1:lcmd),status='old',err=1000)
         print*, ' LAMBERT file: ',wcmd(1:lcmd)
         lambert = wcmd(1:lcmd) 
         call blank(wcmd)
         close(30)
      else
         iomode(17) = 0
      endif
ccc
c
c     implement option of input/output files with 
c        velocity units of m/yr rather than mm/yr
c     FONDA internal workings still the same though
      call getcmd(8,'Velocity',wcmd,lcmd,1)
      if ((lcmd.gt.0).and.(wcmd(1:4).eq.'m/yr')) then
         iomode(18) = 1
      else
         iomode(18) = 0
      endif
ccc


c     gamma file
      call getcmd(8,'gamma_file',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(14) = 1
         ib = lift_arg(wcmd,code,1)
         open (27,file=code(1:ib),status='unknown',err=1000)
         print*, ' gamma rate file name: ',code(1:ib)
         call blank(wcmd)
      else
         iomode(14) = 0
      endif
      

c
c     event record file
      call getcmd(8,'event',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(10) = 1
         open (24,file=wcmd(1:lcmd),status='unknown',err=1000)
         print*, ' event record file name: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(10) = 0
      endif
c
c     loose constraint h-file
      call getcmd(8,'global',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(11) = 1
         open (25,file=wcmd(1:lcmd),status='unknown',err=1000)
         print*, ' global h-file file name: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(11) = 0
      endif
c
c     modified network coordinate file
      call getcmd(8,'modif',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(12) = 1
         open (26,file=wcmd(1:lcmd),status='unknown',err=1000)
cmk      add file for uncorrelated coords
         open (37,file=(wcmd(1:lcmd)//'.cor'),status='unknown',err=1000)
         print*, ' modified network file (uncorrelated): ',wcmd(1:lcmd)
         print*, ' modified network file (correlated): ',
     .    wcmd(1:lcmd)//'.cor'
         call blank(wcmd)
      else
         iomode(12) = 0
      endif
c
c     refraction correction file
      call getcmd(8,'refra',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         iomode(13) = 1
         open (16,file=wcmd(1:lcmd),status='unknown',err=1000)
         print*, ' refraction correction file: ',wcmd(1:lcmd)
         call blank(wcmd)
      else
         iomode(13) = 0
      endif
c
c     read reference frame
      call getcmd(8,'frame',wcmd,lcmd,1)
      if (lcmd.gt.0) frame = wcmd(1:6)
      call blank(wcmd)
      print*,' reference frame: ',frame
c
c     input format code
c     fix to pmode = 1, fcode = 2      
      pmode = 1
      fcode = 2
c     call getcmd(8,'paramet',wcmd,lcmd,1)
c     if (lcmd.gt.0) then 
c        ib = lift_arg(wcmd,code,1)
c        if (mchkey(code,'geoc',ib,4).gt.0) pmode = 1
c        if (mchkey(code,'loca',ib,4).gt.0) pmode = 2
c        if (mchkey(code,'geod',ib,4).gt.0) pmode = 3
c        ib = lift_arg(wcmd,code,2)
c        if (mchkey(code,'char_4',ib,6).gt.0) fcode = 1
c        if (mchkey(code,'char_8',ib,6).gt.0) fcode = 2
c        if (mchkey(code,'char_4_8',ib,8).gt.0) fcode = 3
c     endif
c     print*,' parameter mode: ',pmode,fcode
c
c     get reference time
      call getcmd(8,'rtime',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         read(wcmd,*) code,rtime
         ib = nblen(code)
         if (mchkey(code,'sphe',ib,4).gt.0) comode = 1
         if (mchkey(code,'geod',ib,4).gt.0) comode = 2
         if (mchkey(code,'geoc',ib,4).gt.0) comode = 3
         if (mchkey(code,'topo',ib,4).gt.0) comode = 4
      endif
      print*,' rtime,comode:',rtime,comode
c
c     determine solution combination mode
cmk   currently 'Mark' is default and 'model' is only
cmk    other implemented option - check this!
      call getcmd(8,'soluti',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         ib = lift_arg(wcmd,code,1)
         if (mchkey(code,'Mark',ib,4).gt.0) smode = 1
         if (mchkey(code,'Helm',ib,4).gt.0) smode = 2
         if (mchkey(code,'model',ib,5).gt.0) smode = 3
         if (mchkey(code,'sequ',ib,4).gt.0) smode = 4
         if (mchkey(code,'Kalm',ib,4).gt.0) smode = 5
      endif 
      print*, ' solution combination mode:', smode
c
c     It is unnecessary to define minic and miniv
      minic = 1
      miniv = 1
c     assign minimum observation condition
c     call getcmd(8,'minimum',wcmd,lcmd,1)
c     if (lcmd.gt.0) read(wcmd,*) minic,miniv
c     print*, ' minimum observation condition:',minic,miniv
c
c     set up outlier criteria for angle, length, coordinate and velocity
      call getcmd(8,'criter',wcmd,lcmd,1)
      if (lcmd.gt.0) read(wcmd,*) cria,cril,cric,criv
      write(6,'(a,4f10.3)') ' outlier criteria: ',cria,cril,cric,criv
c
c     loose constraints for coordinate and velocity
c     unit:             meter        meter/year
      call getcmd(8,'loose',wcmd,lcmd,1)
      if (lcmd.gt.0) read(wcmd,*) wcoe,wcon,wcou,wvee,wven,wveu
      write(6,'(a,6f10.3)') ' loose constraint: ',
     .   wcoe,wcon,wcou,wvee,wven,wveu
c
c     get control key
      call getcmd(8,'key',wcmd,lcmd,1)
      if (lcmd.gt.0) read(wcmd,*) key
      print*, ' control key:', key

      goto 200
c
c     get apriori informations
c
c     determine frame of the apriori informations
 100  call getcmd(8,'end',wcmd,lcmd,2)
      call getcmd(8,'apr_frame',wcmd,lcmd,3)
      if (lcmd.gt.0) read(wcmd,*) code
      if (code(1:3).eq.'xyz') id_frame = 1
      if (code(1:3).eq.'enu') id_frame = 2
      call zero1d(1,nlive*6,aprm)
      do 130 iloop = 1,1000
         call getcmd(8,'apr_value',wcmd,lcmd,3)
         if (lcmd.le.0) goto 115
         ie = lift_arg(wcmd,code,1)
         code = upperc(code)
         if (ie.le.0) goto 130
         read(wcmd(ie+1:lcmd),*) (apr_val(i),i=1,6)
         if (code(1:3).eq.'ALL') then
            id_sit = 0
            call put_apr_wght(apr_val,id_sit)
         else
            i1 = nblen(code)
            id_sit = match_name(nsit,i1,sname,code(1:i1))
            if (id_sit.gt.0.and.id_sit.le.nsit) then
               call put_apr_wght(apr_val,id_sit)
            else
                print *,' REDRV: could not deal with ',code
            endif
         endif
 130  continue
c
c     put apriori information on the auxiliary parameters
 115  if (iaux.gt.0) then
         acov = 3.6d3
         do i = nlive-iaux+1,nlive
            if (i.le.3) i1 = i*(i+1)/2
            if (i.gt.3) i1 = (i-1)*3
c            if (iomode(11).gt.0) aprm(i1) = acov
            if (iomode(11).gt.0) aprm(i1) = 1.0d0/acov/acov
            if (iomode(11).le.0) aprm(i1) = 1.0d0/acov/acov
         enddo
      endif
c
c     apriori uncertainty for episodic parameters
      if (jeaux.gt.0) then
         ic = 0
         ie = nlive-iaux-jeaux
         do 125 ib = 1,iq_sit
            if (quake_use(ib).le.0) goto 125
            do j = 1,3
               ic = ic+1
               i = ie+ic
               if (i.le.3) i1 = i*(i+1)/2
               if (i.gt.3) i1 = (i-2)*3
               do j1 = 1,j
                  acov = gvm((ib-1)*6+j*(j-1)/2+j1)
                  aprm(i1+3+j1-j) = acov
               enddo
            enddo
 125     continue
      endif
c
c     get Markov parameters (not finished yet)
 120  continue
      goto 200
c
c     get constraint informations
 150  call getcmd(8,'command',wcmd,lcmd,2)
c     
c     get fixed parameters
      do 70 iloop = 1,1000
         call getcmd(8,'fix_para',wcmd,lcmd,3)
         if (lcmd.le.0) goto 110
         call put_fix(wcmd(1:lcmd),lcmd,chi)
 70   continue
c
c     link site
 110  call getcmd(8,'command',wcmd,lcmd,2)

c     this next section doesn't do anything yet
c     it needs to go into the output section 
c     somewhere
      if (iomode(6).gt.0) then
         call pline (20,80,'-',2)
         write(20,'(6x,a4,9x,3(a5,5x),8x,a4,8x,a5)')
     .   'type','site1','site2','site3','chi2','dchi2'
         call pline (20,80,'-',1)
      endif
      do 30 iloop = 1,1000
         call getcmd(8,'link_site',wcmd,lcmd,3)
         if (lcmd.le.0) goto 140
         call put_link(wcmd,lcmd,chi,1)
         
      if (code.eq.'excl') then
         nexc = isx
         if (nexc.gt.0) then
            read (8,*) (iesit(ie),ie = 1,nexc)
            do 160 ie = 1,nexc
               iexc(iesit(ie)) = 1
 160         continue
         endif
      endif
      if (code.eq.'exit'.or.code.eq.'    ') then
         goto 80
      endif
 30   continue
c
c     link direction
 140  call getcmd(8,'command',wcmd,lcmd,2)
      do 90 iloop = 1,1000
         call getcmd(8,'link_dire',wcmd,lcmd,3)
         if (lcmd.le.0) goto 170
         call put_link(wcmd(1:lcmd),lcmd,chi,2)
 90   continue  
c
c     link baseline length
 170  call getcmd(8,'command',wcmd,lcmd,2)
      do 180 iloop = 1,1000
         call getcmd(8,'link_leng',wcmd,lcmd,3)
         if (lcmd.le.0) goto 182
         call put_link(wcmd(1:lcmd),lcmd,chi,3)
 180  continue  
c
c     link velocity gradient
 182  call getcmd(8,'command',wcmd,lcmd,2)
      do 185 iloop = 1,1000
         call getcmd(8,'link_gradient',wcmd,lcmd,3)
         if (lcmd.le.0) goto 186
         call put_link(wcmd(1:lcmd),lcmd,chi,4)
 185  continue  
c
c     link episodic parameter
 186  call getcmd(8,'command',wcmd,lcmd,2)
      do 187 iloop = 1,1000
         call getcmd(8,'link_episo',wcmd,lcmd,3)
         if (lcmd.le.0) goto 190
         call put_link(wcmd(1:lcmd),lcmd,chi,5)
 187  continue  
c
c     assign network center movement
 190  call getcmd(8,'command',wcmd,lcmd,2)
      do 195 iloop = 1,1000
         call getcmd(8,'assign_cent',wcmd,lcmd,3)
         if (lcmd.le.0) goto 210
         call put_assign(wcmd(1:lcmd),lcmd,chi,1)
 195  continue  
c
c     assign network rotation
 210  call getcmd(8,'command',wcmd,lcmd,2)
      do 215 iloop = 1,1000
         call getcmd(8,'assign_rota',wcmd,lcmd,3)
         if (lcmd.le.0) goto 220
         call put_assign(wcmd(1:lcmd),lcmd,chi,2)
 215  continue  
c
c     assign site movement constraint
 220  call getcmd(8,'command',wcmd,lcmd,2)
      do 225 iloop = 1,1000
         call getcmd(8,'assign_sdir',wcmd,lcmd,3)
         if (lcmd.le.0) goto 260
         call put_assign(wcmd(1:lcmd),lcmd,chi,3)
 225  continue  
c
c     assign episodic parameter
 260  call getcmd(8,'command',wcmd,lcmd,2)
      do 265 iloop = 1,1000
         call getcmd(8,'assign_epi',wcmd,lcmd,3)
         if (lcmd.le.0) goto 230
         call put_assign(wcmd(1:lcmd),lcmd,chi,7)
 265  continue  
c
c     assign inner coordinate constraint
 230  call getcmd(8,'command',wcmd,lcmd,2)
      do 235 iloop = 1,1000
         call getcmd(8,'assign_inner',wcmd,lcmd,3)
         if (lcmd.le.0) goto 240
         call put_assign(wcmd(1:lcmd),lcmd,chi,4)
 235  continue  
c
c     assign outer coordinate constraint
 240  call getcmd(8,'command',wcmd,lcmd,2)
      do 245 iloop = 1,1000
         call getcmd(8,'assign_outer',wcmd,lcmd,3)
         if (lcmd.le.0) then
            if (iomode(6).gt.0) call pline (20,80,'-',3)
            goto 250
         endif
         call put_assign(wcmd(1:lcmd),lcmd,chi,5)
 245  continue  
c
c     assign model coordinate constraint
 250  call getcmd(8,'command',wcmd,lcmd,2)
      do 255 iloop = 1,1000
         call getcmd(8,'assign_model',wcmd,lcmd,3)
         if (lcmd.le.0) goto 200
         call put_assign(wcmd(1:lcmd),lcmd,chi,6)
 255  continue  

 80   icnum = ic
      ibnum = ib
c
c     fix velocity in dynap parameterization mode
      if (pmode.eq.3) then
         do i = 1,nsit
            i1 = (i-1)*6
            fix(i1+3) = 1
            fix(i1+4) = 1
            fix(i1+5) = 1
            fix(i1+6) = 1
         enddo 
      endif
c
      goto 200

c     troubled stop
 1000 print*,'suicide due to .. ',wcmd
      if (lcmd.lt.0) stop
c
 200  continue
      return
      end
