      subroutine getfil (asites,infile,jnfile,vfile,mfile,vers,
     .nobs,nsat,ncfls,ncfls1,luin,
     .imode,lumf,it0,t00,kpart,islot2,idyoyr,
     .more,edtbl,lmfile)

c     get the input file names and read them into memory

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

c     Hide a few variables.
      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

c     integer*2 icombs
      character*1 lowerc,a1,a2,cfv
      character*16 infile(maxsit),buff16,vfile,jnfile(maxsit),mfile
      character*4 asites(maxsit)
      character*80 vers,prog_name
      character*256 message
      logical more,lask,fcheck,lmfile,edtbl,toomny,ltemp,lclmf
      logical lcldoy,lclcfv
      character*3 icldoy
      integer idyoyr,iii
      integer*4 imode,ncfls,ncfls1,nwave,ndat,nobs,nsat,nband
      integer luin,lumf
      integer i,isite
      integer nblen,len,rcpar

      integer ii,ioerr
c     FUNCTION to return command line arguments
      integer iclarg

c     need this junk for M-file read:
      real*8 t00(3),elvcut_dum
      integer it0(3)
      integer islot2(maxprm),kpart
c
c     get the program name
      len = rcpar(0,prog_name)

c     assume we can write output (OK, unless we have RINEX or X-files)
      edtbl = .true.

c     assume we can look at all the c-files
      toomny = .false.
      ltemp = .false.

c     get M-file from command line if it is there
      lclmf = .false.
      lclcfv = .false.
      lcldoy = .false.
      ii = iclarg(1,mfile)
      if (ii .gt. 0) then
         lmfile = .true.
         lclmf  = .true.
c        check to see if day number is there, as well
         iii = iclarg(2,icldoy)
         if (iii .gt. 0 ) lcldoy = .true. 
c        and finally, check for a substitute c-file version letter (6th character)
         iii = iclarg(3,cfv)
         if (iii.gt.0 ) lclcfv = .true.
      else
         print *,'Enter RINEX, X-, C-, or M-file name(s), ! to finish'
 10      continue
            read (5,'(a)',end=20) buff16
            i = nblen(buff16)
            if (buff16(1:1) .eq. '!') then
               more = .false.
            else if (lowerc(buff16(1:1)) .eq. 'm' .and.
     .               lowerc(buff16(i:i)) .ne. 'o') then
               lmfile = .true.
               lclmf  = .false.
               mfile = buff16
               do while (.not. fcheck(mfile))
                  print *,mfile,' does not exist. Enter new name.'
                  read (*,*) mfile
               enddo
               more = .false.
            else
               imode = 0
               ncfls = ncfls+1
               call ljust(16,buff16)
c              check if C-file exists
               do while (.not. fcheck(buff16)) 
                  print *,buff16,' does not exist. Enter new name.'
                  read(5,'(a)') buff16   
                  if (buff16(1:1) .eq. '!') then
                     more = .false.
                  endif
               enddo
               infile(ncfls) = buff16(1:16)
            endif
         if (more) goto 10
 20      continue
      endif

c     OK, use M-file to get list of C-file names
      if (lmfile) then

c        Derive day number from file name
         ii = index (mfile,'.')
         read (mfile(ii+1:ii+3),'(i3)',iostat=ioerr) idyoyr

         call readm(lumf,mfile,idyoyr,ncfls,infile)


c        Error reading day number, ask for it
c        Also ask if more than one session in M-file.
c        Read day number from command line, if possible
         if (ioerr.ne.0 .or. idyoyr.gt.366 ) then
            if (lcldoy) then
               read (icldoy,'(i3)',iostat=ioerr) idyoyr
               if (ioerr.ne.0) then
                 write(message,'(a,i3)') 'Problems reading day number: '
     .                                     ,icldoy
                 call report_stat('WARNING',prog_name,'getfil',' '
     .                            ,message,0)
                endif
            else
               write (6,9000) 'Enter day of year to read.'
               read (*,*) idyoyr
            endif

*           re-read mfile with correct day
            close(lumf)
            call readm(lumf,mfile,idyoyr,ncfls,infile)
         endif

c        Don't ask stupid questions if M-file is on command line
         if (.not. lclmf) then
            write (6,9000) 'Use post-fit residuals?'
            if (lask()) then
               imode = 1
            else
               imode = 0
            endif
         else
c           assume post-fits are desired if M-file is on Command Line
            imode = 1
         endif
      endif


c     If C-files names come from M-file, edit the list.
      if (lmfile) then   
         if( .not.lclcfv ) then
            print *,mfile,' points to these C-files:'
         else
            print *,'C-files requested : '
         endif                    
         ncfls1 = 0
         do 900 i=1,ncfls
            call lowers(infile(i)) 
            if( lclcfv ) infile(i)(6:6) = lowerc(cfv)
            if (fcheck(infile(i))) then
               write (6,'(i2,1x,a)') i,infile(i)(1:nblen(infile(i))) 
               ncfls1 = ncfls1 + 1
            else
               write (6,'(i2,1x,a,a)') i,infile(i)(1:nblen(infile(i))),
     .         ' does not exist'
            endif
  900    continue     

         if( ncfls1.eq.0 ) then
           if( mfile(12:16).eq.'autcl' ) then   
              print *,' ' 
              print *,'C-files on ',mfile
     .          ,' not correct for autcln postfit'
              print *,'Try: ',prog_name(1:8), mfile,' ',idyoyr,'  c' 
              print *,' '
           endif
           call report_stat('FATAL',prog_name,'getfil',' '
     .            ,'No C-files available',0 )
         endif

         if (ncfls .gt. ncvsit) then
            write (message,901) ncfls,ncvsit
901         format ('Number of C-files (',i2,') exceeds dimension (',
     .      i2,') NCVSIT. I would strongly suggest you pick a subset.')
            call report_stat('WARNING',prog_name,'getfil',' ',message,0)
            toomny = .true.
         endif

c        allow user to delete some C-files from list in M-file
c        deleted files will have a space as the first char in name
c        Don't ask this question if M-file name is on command line
         if (.not. lclmf) then
            if (.not. toomny) then
               write (6,9000) 'Use all these stations?'
               ltemp = lask()
            endif
         else
            ltemp = .true.
         endif

         if (.not.ltemp .or. toomny) then
            do 915 i=1,ncfls
               write (6,903) i,infile(i)(1:nblen(infile(i)))
 903           format (1x,'Use C-file # ',i2,1x,a,'?',$)
               if (.not. lask()) then
                  infile(i)(1:1) = ' '
               endif
 915        continue
         endif

c        change series letter if desired
c        Don't ask stupid questions if M-file name is on command line
         if (.not. lclmf) then
            write (6,9000) 'Use this series of C-files?'
            if (.not.lask()) then
               write (6,9000) 'Enter letter, e.g. b for c????b.??? '
               read (*,'(a1)') a1
               do 925 i=1,ncfls
                  infile(i)(6:6) = a1
  925          continue
            endif
         endif
      endif

c     make sure all files exist, and copy site names
      write (6,'(1x,a)') 'You have selected these files: '
      ncfls1 = 0
      do 950 i=1,ncfls
         call lowers(infile(i))
         if (infile(i)(1:1) .ne. ' ') then
            if (fcheck(infile(i))) then
               print *,infile(i)
               ncfls1 = ncfls1 + 1
c              jnfile(ncfls1) = infile(i)
               a1 = infile(i)(1:1)
               a2 = infile(i)(nblen(infile(i)):nblen(infile(i)))
               if (lowerc(a2) .eq. 'o') then
                  asites(ncfls1) = infile(i)(1:4)
               else if (a1.eq.'c' .or. a1.eq.'x' .or. a1.eq.'f') then
                  asites(ncfls1) = infile(i)(2:5)
               else
                  asites(ncfls1) = infile(i)(1:4)
               endif
            else
               call report_stat('WARNING',prog_name,'getfil',infile(i),
     .         'File does not exist: ',0)
               infile(i)(1:1) = ' '
            endif
         endif
 950  continue

      if (ncfls1 .gt. ncvsit) then 
         call report_stat('FATAL',prog_name,'getfil',' ',
     .   'Silly, you picked too many C-files!',0)
      endif

c     read files and store data in arrays
c     note that isite is not the same as i, because i loops over
c     ALL files, while isite just handles those we are reading!
c      'i' must be used to index the adjustments on the m-file  (READC only)
c      'isite' must be used to store quantities for CVIEW (equivalent to 'i' for X-
c         and RINEX files since no m-file used.
      isite = 0
      do 1000 i=1,ncfls
         a1 = infile(i)(1:1)
         a2 = infile(i)(nblen(infile(i)):nblen(infile(i)))
         if (lowerc(a2) .eq. 'o') then
            write(6,980)
 980        format (//,80('_'),/)
            write(6,*) 'Assuming this is a RINEX file.'
            write(*,*) 'Enter GNSS code to request (G R C E J I)'
            read(5,'(a)') gnss
            call uppers(gnss)
            isite = isite + 1
            call rinex (luin,infile(i),isite,nobs,nsat)
            edtbl = .false.
         else if (lowerc(a1) .eq. 'c' .or. lowerc(a1) .eq. 'f') then
            write(6,980)
            isite = isite + 1
            call readc(luin,lumf,infile(i),imode,isite,i,nobs,nsat,ndat)
         else if (lowerc(a1) .eq. 'x') then
            write(6,980)
            isite = isite + 1
            call readx(luin,infile(i),imode,isite,nobs,nsat,nband,nwave)
            edtbl = .false.
         else if (a1 .eq. ' ') then
c           Tricky! We must read C-file header even if we have no
c           intention of reading the C-file!
            if (imode.eq.1) call readm3 ( lumf,it0,t00,kpart,islot2
     .                                ,elvcut_dum,zenmod,numzen,idtzen
     .                                ,gradmod,numgrad,idtgrad )
         endif

*        MOD TAH 160109: Map all the GLONASS omc values to the central
*        L1 (1602 Mhz) and L2 (1246 Mhz) freqencies so that doouble 
*        differences can be formed (integr ambiguities are lost).
*        Make sure we have values
         if( gnss(1:1).eq.'R' ) then
             call remapf( isite, nobs, nsat)
         endif

 1000 continue

 9000 format (1x,a,1x,$)
         
      return
      end

C TEMPORARY FOR TESTING: MOVE LATER

      subroutine remapf( isite, nobs, nsat )

c
c     Routine to remap GLONASS frequencies to one frequecy so that
c     doible differences can be computed.  Only needed for gnss==R
c
c     input:
c            isite   the CVIEW or SCANDD index number of this site
c            nobs    Number of epochs (from getfil)
c            nsat    Number of channles
c     output
c            None    (Values in common are changed).
c
c     gnss    'R' for glonass: Used from COMMON
c
c     Coder  T. Herring.

      implicit none
c
      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

* Passed
      integer*4 isite   ! Site nuuber
      integer*4 nobs    ! Number of epochs
      integer*4 nsat    ! Number of channels

* Local
      integer*4 ep   ! Epoch counter
      integer*4 ch   ! Channels counter

      real*8 fL1_R0, fL2_R0, fL3_R0    ! GLONASS 0 frequencies at L1, L2, L3

***** Verify this is GLONASS data
      if( gnss(1:1).ne.'R' ) RETURN

      print *,'Mapping GLONASS phase to common frequency'

      fL1_R0 = 1602.0d6  ! GLONASS L1 center Hz
      fL2_R0 = 1246.0d6  ! GLONASS L2 center Hz
      fL3_R0 = 1201.0d6  ! GLONASS L3 center Hz
      

*     OK Loop over epochs 
      do ep = 1, nobs
*        Loop over all possible channels that this epoch
*        (Fix later to used actual number of channels)
         do ch = 1, nsat
*           Scale fL1 to fL1_R0
            yl1(ep,ch,isite) = yl1(ep,ch,isite)*fL1_R0/fL1(ch)
            pr1(ep,ch,isite) = pr1(ep,ch,isite)*fL1_R0/fL1(ch)
            yl2(ep,ch,isite) = yl2(ep,ch,isite)*fL2_R0/fL2(ch)
            pr2(ep,ch,isite) = pr2(ep,ch,isite)*fL2_R0/fL2(ch)
         end do
      end do

****  Thanks all
      return
      end


  
