      program scanrms
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
      IMPLICIT REAL*8 (A-H,O-Z)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

c     Hide a few variables.
      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

c     integer*2 icombs
      character*16 infile(maxsit),vfile,jnfile(maxsit),sfile,mfile
      character*4 asites(maxsit)
      character*80 vers
      logical more,lmfile,edtbl
      integer idyoyr
      integer imode,ncfls,ncfls1,nobs,nsat
      integer luin,lumf
      integer i
      integer nblen
      integer iv,ir

c     need this junk for M-file read:
      real*8 t00(3)
      integer it0(3)
      integer islot2(maxprm),kpart

c#      call sversn(vers)
      call cversn(vers)

c     file unit numbers
c     input file unit
      luin = 11
c     M-file unit
      lumf = 13

c     output V-file
      vfile = 'vscan.out'
      iv = 47
      open (unit=iv,file=vfile,status='unknown')
      write (iv,*) vers

c     output RMS file
      ir = 49
      sfile = 'scan.rms'
      open(unit=ir,file='scan.rms',status='unknown')

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

      call scanner (nobs,nsat,ncfls1,iv,ir)

      close (iv)
      close (ir)

         write(6,8000)

         write(6,2100)
 2100    format(//,'Created scan.rms and vscan.out',//,
     .      'Sort scan.rms by running the program "sorter"',//,
     .      'Sort vscan.out using the script "sortv"',//)

      print *,'Normal end of SCAN'

 8000 format (/,80('_'),//)
 9000 format (1x,a,1x,$)
c
      stop 0
      end
