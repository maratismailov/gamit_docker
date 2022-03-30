      subroutine finish (lumf,luin,ncfls,ncfls1,
     .infile,jnfile,mfile,vfile,asites,vers,imode,
     .idyoyr,idelc,interactive,jstat)

c     write the C-files, using M-file if necessary.
     
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

c     Hide a few variables.
      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

c     integer*2 icombs
      character*16 infile(maxsit),buff16,vfile,jnfile(maxsit),mfile
      character*4 asites(maxsit)
      character*80 vers
      integer idyoyr
      integer jstat
      integer*4 imode,ncfls,ncfls1
      integer luin,lumf
      integer i,isite,ipick
c     idelc           0 to delete old C-file (namec0)
c                     1 to keep old C-file and write new (namec1)
c                     2 to delete old C-file and rename output to namec0.
c                       (this gives the illusion of overwriting).
      integer idelc
c     .true. to ask questions
      logical interactive
c     .true to update C-files
      logical updatec

c     need this junk for M-file read:
      real*8 t00(3)
      real*8 elvcut_dum
      integer it0(3)
      integer islot2(maxprm),kpart

      if (interactive) then
         write(6,1700)
 1700    format(//)
         write(6,1710)
 1710    format (1x,'You may:'/,
     .   '  1 Write new C-files, incrementing series letter.',/,
     .   '  2 Write new C-files, incrementing series letter',
     .   ' and deleting old C-files.',/,
     .   '  3 Overwrite input C-files.',/,
     .   '  4 Not write any C-files.',/,
     .   '  5 Not write any C-files, but save vcview.out',//,
     .   ' Only 1 and 4 will let you continue editing.',///)
         call imenu (ipick,5)
         if (ipick .eq. 1) then
            updatec = .true.
            idelc = 1
         else if (ipick .eq. 2) then
            updatec = .true.
            idelc = 0
         else if (ipick .eq. 3) then
            updatec = .true.
            idelc = 2
         else if (ipick .eq. 4) then
            updatec = .false.
            idelc = 1
         else
            updatec = .false.
            idelc = 3
         endif
      else
         updatec = .true.
         idelc = 2
      endif

      if (updatec) then
         if (imode .eq. 1) then
            close(lumf)
c           read into jnfile, rather than infile because infile
c           has our edited list of stations.
            call readm(lumf,mfile,idyoyr,ncfls,jnfile)
         endif
         write(6,1720)
 1720    format(//,' Writing C-files:',/)
         isite = 0
         do 2000 i=1,ncfls
            if (infile(i)(1:1) .ne. ' ') then
c              construct new output C-file name
               call upnam1 (infile(i),buff16)
               call lowers(buff16)
               write(6,980)
c              and write to it
               isite = isite + 1
               call writec (infile(i),buff16,vfile,
     .              idelc,imode,isite,i,lumf,vers,jstat)
            else
c             Tricky! We must read C-file header even if we have no
c             intention of reading the C-file!
              if (imode.eq.1) then
                 call readm3 (lumf,it0,t00,kpart,islot2
     .             ,  elvcut_dum,zenmod,numzen,idtzen
     .             ,  gradmod,numgrad,idtgrad )
              endif
            endif
 2000    continue
      endif

      if( ipick .eq. 5 ) then
        if (imode .eq. 1) then
          close(lumf)
          call readm(lumf,mfile,idyoyr,ncfls,jnfile)
        endif
        write(6,2500)
 2500   format(//,' Writing VCVIEW file:',/)
        isite = 0
        do 3000 i = 1,ncfls
          if (infile(i)(1:1) .ne. ' ') then
            write(6,980)
            isite = isite + 1
            call write_vcview
     .          (infile(i),vfile,imode,isite,i,lumf,vers,jstat)
          else
            if (imode.eq.1) then
              call readm3 (lumf,it0,t00,kpart,islot2
     .             ,  elvcut_dum,zenmod,numzen,idtzen
     .             ,  gradmod,numgrad,idtgrad )
            endif
          endif
 3000   continue
      endif

 980  format (/,80('-'))
 9000 format (1x,a,1x,$)
c
      return
      end

