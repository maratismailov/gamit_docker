      program scanflg
c
c     Double difference cycle cleaner
c
c     Shimon Wdowinski
c     March, 1991
c
c     Roots:
c        CVIEW by Kurt Feigl

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
      character*16 infile(maxsit),vfile,jnfile(maxsit),mfile
      character*4 asites(maxsit)
      character*80 vers
      logical more,lmfile,interactive,edtbl
      integer idyoyr,idelc,jstat
      integer imode,ncfls,ncfls1,nobs,nsat
      integer luin,lumf

c     need this junk for M-file read:
      real*8 t00(3)
      integer it0(3)
      integer islot2(maxprm),kpart

c
      call cversn(vers)

c     input file unit
      luin = 11
c     M-file unit
      lumf = 13

      ncfls=0
      more = .true.
      lmfile = .false.

c     obtain the input file names and read the files
      call getfil (asites,infile,jnfile,vfile,mfile,vers,
     .nobs,nsat,ncfls,ncfls1,luin,
     .imode,lumf,it0,t00,kpart,islot2,idyoyr,
     .more,edtbl,lmfile)


      if (ncfls1 .gt. 1) then
c        ******** cleaning the data  **************
         open(unit=47,file='slip.epo',status='unknown')
         open(unit=48,file='slip.com',status='unknown')
         open(unit=49,file='slip.inf',status='unknown')

         call driver2 (nobs,nsat,ncfls1)

         close(47)
         close(48)
         close(49)

c        write out C-files if possible.
         if (edtbl) then
            call finish (lumf,luin,ncfls,ncfls1,
     .      infile,jnfile,mfile,vfile,asites,vers,imode,
     .       idyoyr,idelc,interactive,jstat)
         endif

         write(6,8000)
         write(6,2100)
 2100    format(//,' Unresolved silps are listed in the following',
     .                                                   ' file:',/,
     .           '        sort slip.epo > epo.srt  ',///)

      else
         write(6,2220)
 2220    format(//,' I need at least 2 c-files for SCANFLG',/)

      endif

      if (imode .eq. 1) then
         close(lumf)
      endif

 8000 format (/,80('_'),//)
 9000 format (1x,a,1x,$)
c
      stop
      end
