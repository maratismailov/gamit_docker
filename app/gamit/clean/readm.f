      subroutine readm(luin,mfile,idyoyr,ncfls,cfiles)
c
c     read an m-file named mfile on unit luin
c
c     return array with cfiles names
c     clean up of readm routine K. Feigl May 89
c     standard open and readm routines M. Murray July 89
c
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'

      character*16 mfile
      character*3 buf3
      character*80 prog_name
      character*256 message

      integer*4 ioerr,len,rcpar
      integer luin,nsite,nrfile,ntfile,ntsat
      integer ncfls,kpart,i,nepoch,idyoyr,idytst,nsat,mpos,nskip
c
c     main header arrays *
      character*20 alabel(maxprm)
      character*12 site(maxsit)
      integer*4 isat(maxsat),norb
      real*8 aprval(maxprm)
      character*16 rfname(maxsit),tfname(maxtfl)
c
      character*16 cfiles(maxsit),cftmp
      real*8 stawgt(maxsit),satwgt(maxsat)
c
c     c-file header arrays *
      integer islot2(maxprm)
c
c     mfile record 3
c-      character*3    zenmod,gradmod - in cview.h
c-      integer*4      numzen,numgrad - in cview.h
c-      integer*4      idtzen(maxatm),idtgrad(maxgrad/2) - in cview.h
      real*8         elvcut_dum
c
c     get the program name
      len = rcpar(0,prog_name)
c
c     open merge file
      call mopens (mfile,'old',luin,ioerr)
      if (ioerr.ne.0) then
        call report_stat('FATAL',prog_name,'readm',mfile,
     .  'Error opening M-file: ',ioerr)
      endif
c
      call report_stat('STATUS',prog_name,'readm',mfile,
     .'Reading M-file: ',0)

c     read main header
      call readm1 (luin,
     $             nepoch,mtpart,
     $             alabel,idms,islot,aprval,adjust,
     $             ntsat,isat,
     $             nsite,site,
     $             nrfile,rfname,
     $             ntfile,tfname,norb )
c
c     read the c-file list
   50 continue
      call readm2 (luin,
     $             iit0,tt00,
     $             nepoch,inter,nskip,
     $             ncfls,cfiles,stawgt,
     $             nsat,satwgt)

      cftmp=cfiles(1)
      mpos = index(cfiles(1),'.')
      write(unit=buf3,fmt='(a3)') cftmp(mpos+1:mpos+3)
      read(unit=buf3,fmt='(i3)') idytst
                  
c     normal finish
      return
      end
