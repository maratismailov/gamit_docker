      Program CVIEW
c
c     Interactive editor for C-files
c
c     Kurt Feigl and Shimon Wdowinski
c
c     with contributions by:
c        Yehuda Bock
c        Mark Murray
c        Peter Morgan
C        Greg Beroza
c        Justin Revenaugh
c        Tom Herring

c
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

c     Hide a few variables.
      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

c     integer*2 icombs
      character*16 infile(maxsit),vfile,jnfile(maxsit),mfile
      character*4 asites(maxsit)
      character*80 vers
      character*256 message
      logical more,lask,lmfile,edtbl,interactive
      integer idyoyr
      integer*4 imode,ncfls,ncfls1,nobs,nsat,iscrn
      integer luin,lumf
c     length of LIST
      integer llist
c     0 to delete old C-files, 1 to keep
      integer idelc
      integer jstat

c     need this junk for M-file read:
      real*8 t00(3)
      integer it0(3)
      integer islot2(maxprm),kpart

c     file unit numbers
c     screen
      iscrn = 6
c     input file unit
      luin = 11
c     M-file unit
      lumf = 13

c     get program version and write status
      call cversn(vers)
      write(iscrn,'(a)')' '
      write(message,5)vers
    5 format('Started CVIEW ',a80)
      call report_stat('STATUS','CVIEW','cview',' ',message,0)

      vfile = 'vcview.out'
      ncfls=0
      more = .true.
      lmfile = .false.
      interactive = .true.

c     obtain the input file names and read the files
      call getfil (asites,infile,jnfile,vfile,mfile,vers,
     .nobs,nsat,ncfls,ncfls1,luin,
     .imode,lumf,it0,t00,kpart,islot2,idyoyr,
     .more,edtbl,lmfile)
           
c     assign the frequency combinations needed for editing
c     (fundamental fL1,fL2 read from c-file or set in subroutines rinex and readx called by getfil)
      call set_freqs(nsat)
               
c     warn about forming between-satellite differences for Glonass
      if( gnss.eq.'R' ) 
     .  write(6,'(2a)') '**WARNING: Between-satellite and double '
     .                 ,'differences will not be valid for Glonass FDMA'

c     Loop over editing process.
      more = .true.
      do 1200 while (more)

c        open LIST file and reads the list
         call read_list (llist)

c        the EDITOR
         call editor (asites,nobs,nsat,ncfls1,imode,llist)

c        write out C-files if possible.
         if (edtbl) then
            call finish (lumf,luin,ncfls,ncfls1,
     .      infile,jnfile,mfile,vfile,asites,vers,imode,
     .       idyoyr,idelc,interactive,jstat)
         endif

c        May NOT update C-files, once the original is disturbed,
c        because WRITEC needs it for comparison.
c        lost the original basis for
         if (  idelc.eq.0 .or. idelc.eq.2 .or. idelc.eq.3
     .       .and. jstat .eq.0 ) then
            more = .false.
         else
            write (6,1210)
 1210       format (//)
            write (6,9000) 'Do you wish to continue editing?'
            more = lask()
         endif
 1200 continue

      if (imode .eq. 1) then
         close(lumf)
      endif

 9000 format (1x,a,1x,$)
c
      stop
      end
