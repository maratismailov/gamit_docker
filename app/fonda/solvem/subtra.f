      subroutine subtra(isum,mode)
c
c     mode=1: subtract effective site from whole sites
c     mode=2: subtract effective submatrix from whole normal matrix
c
c     itoj:  from effective to total index
c     jtoi:  from total to effective index
c     map :  from total to live parameter index
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      integer isum,mode,i,isit,i1,i2,j,npar,id,ibad
      integer iuse(maxsit),vuse(maxsit),nbase,vbase
      equivalence (iuse,map)
      equivalence (vuse,map(maxsit+1))

      small = 1.0d-15
      if (mode.eq.2) goto 60
c
c     discard sites with too few observationd
      isit = 0
      nvsit = 0
      vbase = miniv
      if (miniv.eq.0) vbase = 1
      do 20 i = 1,nsit
         jtoi(i) = 0
         i1 = (i-1)*6+1
         i2 = i*6
         if ((minic.gt.0.and.iuse(i).ge.minic).or.vuse(i).ge.vbase) then
            isit = isit+1
            itoj(isit) = i
            jtoi(i) = isit
            if (iuse(i).ge.minic.and.vuse(i).lt.vbase) nvsit = nvsit+1
            goto 20
         endif
         do j = i1,i2
            fix (j) = 1
         enddo
 20   continue
      gdsit = isit
            
      if (mode.eq.1) goto 100
c
c     copy whole normal matrix to temporary file 33
 60   npar = nsit*6+jaux+iaux	
      if (iomode(3).ge.1.and.iq_sit.gt.0) then
         nbase = nsit*6
         jeaux = 0
         ibad = 0
         do 40 i = 1,iq_sit
            if (quake_use(i).gt.0) goto 40
            fix(nbase+i*3-2) = 1
            fix(nbase+i*3-1) = 1
            fix(nbase+i*3) = 1
            ibad = ibad+1
 40      continue
         jeaux = (iq_sit-ibad)*3
         npar = npar+iq_sit*3
      endif

      nlive = 0
      isum = 0
      do 50 i = 1,npar
         map(i) = 0 
         if (fix(i).eq.1) goto 50
         nlive = nlive+1
         map(i) = nlive 
         do 10 j = 1,i
            if (fix(j).eq.1) goto 10
            id = i*(i-1)/2+j
            isum = isum+1
            anorm(isum) = anorm(id)
 10      continue
         bnorm(nlive) = bnorm(i)
 50   continue
c
c     temporary file in the case of output h-file
      if (iomode(11).gt.0) then
         open(unit=33,form='unformatted',status='scratch')
         write (33) nlive	
         do 30 i = 1,nlive
            i1 = i*(i-1)/2
            write (33) (anorm(i1+j),j=1,i)
 30      continue
         write (33) (bnorm(j),j=1,nlive)
      endif
c      print*,'mapping:',(map(i),i=1,npar)
 100  continue
      return
      end
c     

