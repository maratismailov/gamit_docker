Copyright 1995 Massachusetts Institute of Technology and the University of
California at San Diego.  All Rights Reserved.

      Subroutine UPDATE(filnam)
c
c     Update M-file for C-file editor

      implicit none
c
      include '../includes/dimpar.h' 
      include 'solve.h'
      include 'parameters.h'
      include 'models.h'

c     main header
      integer*4    ntpartm
      real*8       temp(maxprm)
c
c     session header arrays
      integer*4     integ,nrfile
      real*8       swgt(maxsit),dwgt(maxsit)
      character*16 cfile(maxsit)
c
c     C-file header arrays
      integer*4 mzen,mgrad 
      real*8       elvcut_mfile(maxsit)

c     Local
      integer*4 ioerr,nskip,ntepoc,icf,ntfile 
      character*16 filnam,tmpnam,tfilnm(maxtfl)
c       need local m-file version of labels to avoid overwriting 
c       DD bias labels in parameters.h
      character*9 rlabelm(maxprm) 
                         


      if( logprt ) write(6,'(/,2a,/,2a)' ) 'Updating M-file : ',minf
     .                      , 'New M-file      : ',moutf
      if( iqflag.eq.1 )
     .write(10,'(/,2a,/,2a)' ) 'Updating m-file : ',minf
     .                       , 'New m-file      : ',moutf

c     read the original m-file

c     now closed in lsquar to avoid mysterious clobbering on Solaris 2 at GSFC
c      close(unit=11)
      tmpnam=minf
      call lowers(tmpnam)

      call mopens (tmpnam,'old',11,ioerr)
      if (ioerr.ne.0) then
         call report_stat('WARNING','SOLVE','update',tmpnam
     .          , 'Cannot open old M-file--no upodate',ioerr)
         return
      endif
c
c     main header
      call readm1 (11,
     $             ntepoc,ntpartm,
     $             rlabelm,idms,islot1,preval,temp,
     $             nsat,isprn,
     $             nsite,sitnam,
     $             nrfile,obfiln,
     $             ntfile,tfilnm,norb)
      
c
c     session header
         call readm2 (11,
     $                it0,t00,
     $                nepoch,integ,nskip,
     $                nsite,cfile,swgt,
     $                nsat,dwgt(1) )
c
c      c-file headers on the m-file
         do icf=1,nsite
            call readm3 (11,
     $                 itor(1,icf),tor(1,icf),
     $                 npartm(icf),islot2(1,icf),
     $                 elvcut_mfile(icf),zenmod,mzen,idtzen(1),
     $                 gradmod,mgrad,idtgrad(1))
         enddo
      close(11)

c  update the m-file with new adjustments
c
      filnam = moutf
      call lowers (filnam)
c
      call mopens (filnam,'unknown',12,ioerr)
      if (ioerr.ne.0) then
         call report_stat('WARNING','SOLVE','update',filnam
     .          , 'Cannot open new M-file--no upodate',ioerr)
         return
      endif

c
      call writm1 (12,
     $             ntepoc,ntpartm,
     $             rlabelm,idms,islot1,preval,adjust,
     $             nsat,isprn,
     $             nsite,sitnam,
     $             nrfile,obfiln,
     $             ntfile,tfilnm,norb)
c
         call writm2 (12,
     $                it0,t00,
     $                nepoch,integ,nskip,
     $                nsite,cfile,swgt,
     $                nsat,dwgt )
c
         do icf=1,nsite
            call writm3 (12,
     $                   itor(1,icf),tor(1,icf),
     $                   npartm(icf),islot2(1,icf),
     $                   elevcut(icf),zenmod,nzen,idtzen(1),
     $                   gradmod,ngrad,idtgrad(1) )
        enddo
      close(12)
c
cd     print*,'update ntpartm ntpart: ',ntpartm,ntpart
cd     do i = 1, ntpart
cd        print*,'update temp,adjust: ',i,temp(i),adjust(i)
cd     enddo
c
      return
      end
