      subroutine res_list(ifile,idatf,nobs)
c
c     calculate residuals with full covariance data
c
c     data type: (it)
c         1. astrometric azimuth
c         2. horizontal angle
c         3. horizontal direction
c         4. baseline length
c         5. zenith height
c         6. leveling
c
c        11. astrometric azimuth rate
c        12. horizontal angle rate
c        13. horizontal direction rate
c        14. baseline length rate
c        15. zenith height rate
c        16. leveling rate
c
c        21. 3-D geocentric coordinate
c        22. 3-D geocentric velocity
c        23. 3-D geodetic coordinate
c        24. 3-D geodetic velocity
c        25. 3-D geocentric baseline vector
c        26. 3-D geocentric baseline rate vector
c        27. 3-D geodetic baseline vector
c        28. 3-D geodetic baseline rate vector
c
c        31. 3-D spherical coordinate
c        32. 3-D spherical frame velocity
c        33. 3-D Cartesian coordinate
c        34. 3-D Cartesian frame velocity
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : mm/year
c        err      : mm
c        time     : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      character*6 unit(3)
      character*25 desit(41)
      character*80 subfile
      integer ifile,idatf,nobs,indsav,iazsav,i,iobs,ie
      integer iobs1,it,isubf,iefob,ib,ieff,ie1,k,imsit
      integer istsav,iauxsav,iex,lsit,ifrm
      dimension indsav(12),coesav(12)
      common/saved/azisav,coesav,obssav,esave,tsave,iazsav,indsav
      common/save2/istsav,iauxsav
      common /iframe/ifrm
c
      data (desit(i),i=1,6)/'astrometric azimuth      ',
     . 'horizontal angle         ','horizontal direction     ',
     . 'baseline length          ','zenith height            ',
     . 'leveling                 '/
      data (desit(i),i=11,16)/'astrometric azimuth rate ',
     . 'horizontal angle rate    ','horizontal direction rate',
     . 'baseline length rate     ','zenith height rate       ',
     . 'leveling rate            '/
      data (desit(i),i=21,24)/
     . '3-D geocentric coordinate','3-D geocentric velocity  ',
     . '3-D geodetic coordinate  ','3-D geodetic velocity    '/
      data (desit(i),i=25,28)/
     . '3-D geoc baseline vector ','3-D geoc bsln rate vector',
     . '3-D geod baseline vector ','3-D geod bsln rate vector'/
      data (desit(i),i=31,34)/
     . '3-D spherical coordinate ','3-D spherical velocity   ',
     . '3-D Cartesian coordinate ','3-D Cartesian velocity    '/
c
c     print title
      print*,' ------  Residual calculation begins ......'
      print*,'   Observation Type    # Obs.  Pre-rms  Post-rms'
      write (ifile,*) '*  ----------  Residual list  ----------'
c
      rewind (idatf)
      read (idatf,*) nsit
      do i = 1,nsit
         read (idatf,'(2x)') 
      enddo
    
      ieff = 0
      ifrm = 0
      rms1 = 0.0d0
      rms2 = 0.0d0
c     get subfile number
      read (idatf,*) isubf
      do 50 ie = 1,isubf
c        get subfile name
         read (idatf,'(a)') subfile
         open (41,file=subfile,status='old',err=50)
         read (41,*) lsit
         do i = 1,lsit
            read (41,'(2x)')
         enddo
         read (41,*) iexp
         iobs = 0
         do 40 iex = 1,iexp
c           read experiment index, obs. type and number
            read (41,*) ie1,it,iobs1
            k = iobs1
            if (it.gt.20) k = iobs1*3
            if (jaux.gt.0.and.(it.eq.31.or.it.eq.32)) ifrm = ifrm+1
            write (ifile,'(''* exp. : '',i3,5x,''type:'',
     .          a25,5x,''No:'',i6)')
     .          iex,desit(it),k
            prerms = 0.0d0
            pstrms = 0.0d0
            wrms  = 0.0d0
c
            if (it.gt.0.and.it.le.30) 
     .      call get_res_sgl(ifile,41,iobs1,it,prerms,pstrms,
     .         unit,wrms,iefob,iobs,imsit)
            if (it.gt.30.and.it.le.40) 
     .      call get_res_full(ifile,41,iobs1,it,prerms,pstrms,
     .         unit,wrms,iefob)
         if (iefob.le.0) goto 40
         rms1 = rms1+prerms
         rms2 = rms2+pstrms
         ib = iefob-imsit
         ieff = ieff+ib
         prerms = dsqrt(prerms/ib)
         pstrms = dsqrt(pstrms/ib)
         wrms = dsqrt(wrms/ib)

         print '(a25,i4,2f10.2)',desit(it),ib,prerms,pstrms
         write (ifile,*) '*.............. statistics ................'
         write (ifile,*) '*  observations:',ib
         write (ifile,*) '*  parameters:',' <need help from Tom>'
         write (ifile,*) '*  prefit rms :',unit(3),prerms
         write (ifile,*) '*  postfit rms:',unit(3),pstrms
         write (ifile,*) '*  weighted rms:   ',wrms
         write (ifile,*) '*..........................................'
         write (ifile,*) '*   '
 40   continue
c
      close (41)
c
 50   continue
      if (ieff.le.0) goto 60
      rms1 = dsqrt(rms1/ieff)
      rms2 = dsqrt(rms2/ieff)
      write (ifile,*) '*   '
      write (ifile,*) '*---------------- statistics -----------------'
      write (ifile,*) '*  observations:',ieff
      write (ifile,*) '*  parameters:',nlive-ibnum-icnum-ncht 
      write (ifile,*) '*  prefit rms :','(unit??)',rms1
      write (ifile,*) '*  postfit rms:','(unit??)',rms2
      write (ifile,*) '*---------------------------------------------'
      nobs = nobs+ieff
 80   format (a1,1x,a8,2x,a8,1x,f10.3,2(1x,f15.5,1x),f12.5,2f15.5)
c
 110  format (1x,i3,'. ',50(5(1x,d23.16),:,/,6x))
c
 60   continue
      return
      end
