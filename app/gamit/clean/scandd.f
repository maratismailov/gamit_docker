      program scandd
c
c     Scan all double difference combinations
c     Report RMS and JUMPS
c
c
c     Shimon Wdowinski
c     May, 1991
c     Kurt Feigl
c     June, 1991
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
      logical more,lmfile,edtbl,fcheck
      integer idyoyr,imode,ncfls,ncfls1,nobs,nsat,luin,lumf
      integer nblen,i
      integer iv,ir,idd,iscrn,ioerr

c     need this junk for M-file read:
      real*8 t00(3)
      integer it0(3)
      integer islot2(maxprm),kpart
           
c     Exit if a previous step has failed
  
      if( fcheck('GAMIT.fatal') )
     .  call report_stat('FATAL','SCANDD','scandd',' '
     .                  ,'GAMIT.fatal exists: SCANDD not executed',0)
  
c     file unit numbers
c     screen
      iscrn = 6
c     input file unit
      luin = 11
c     M-file unit
      lumf = 13

      call cversn(vers)
      write(iscrn,'(a)')' '
      write(message,5)vers
    5 format('Started SCANDD ',a80)
      call report_stat('STATUS','SCANDD','scandd',' ',message,0)

c     output V-file
      vfile = 'vscan.out'
      iv = 47
      open (unit=iv,file=vfile,status='unknown',iostat=ioerr)
      if (ioerr .ne. 0 ) then
        call report_stat('FATAL','SCANDD','scandd',vfile,
     .  'Error opening output file: ',ioerr)
      endif

c     output RMS file
      ir = 48
      open(unit=ir,file='scan.rms',status='unknown',iostat=ioerr)
      if (ioerr .ne. 0 ) then
        call report_stat('FATAL','SCANDD','scandd','scan.rms',
     .  'Error opening output file: ',ioerr)
      endif

c     output RMS file
      idd = 49
      open(unit=idd,file='scan.dd',status='unknown',iostat=ioerr)
      if (ioerr .ne. 0 ) then
        call report_stat('FATAL','SCANDD','scandd','scan.dd',
     .  'Error opening output file: ',ioerr)
      endif

      ncfls=0
      more = .true.
      lmfile = .false.

c     obtain the input file names and read the files
      call getfil (asites,infile,jnfile,vfile,mfile,vers,
     .nobs,nsat,ncfls,ncfls1,luin,
     .imode,lumf,it0,t00,kpart,islot2,idyoyr,
     .more,edtbl,lmfile)

      do i = 1,ncfls1
         write (iv,'(i2,1x,a)') i,infile(i)(1:nblen(infile(i)))
      enddo

      call scannerdd (nobs,nsat,ncfls1,iv,ir,idd)

      close (iv)
      close (ir)

      write(6,8000)

      write(6,2100)
 2100 format(//,'Created scan.rms, scan.dd, and vscan.out',//,
     .      'Sort scan.rms by running the program "sorter"',//,
     .      'Sort vscan.out using the script "sortv"',//)

      call report_stat('STATUS','SCANDD','scandd',' ',
     .'Normal end of SCANDD',0)

 8000 format (/,80('_'),//)
 9000 format (1x,a,1x,$)
c
      stop
      end
