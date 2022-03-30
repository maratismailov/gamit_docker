CTITLE NUTIAT2000

       program nutiau2000

       implicit none
c
c     iout      output unit number
c     title     title for output file
c     int       interval (days) of values in output table
c     units     input units (arc sec)   (output is 1.e-4 sec)
c     n         number of values in input arrays
c     t         array of input dates (NB: PEP JD, not MJD)
c     x         array of input values of psi
c     y         array of input values of eps


* VARIABLES:

*     jd   - Julian date for evaluating angles.

*     dpsi_ls and deps_ls      - luni-solar nutation in
*             longitude and oblquity (mas) (REAL* OUTPUT)
*     dpsi_plan and deps_plan  - contributions to the
*             nutations in longitude and obliquity due direct
*             planetary nutations and the perturbations of the
*             lunar and terrestrial orbits (mas). (REAL* OUTPUT)
*     dpsi_fcn and deps_fcn    - contributions to the
*             nutations in longitude and obliquity due the free-
*             excitation of the Free-core-nutation (mas).  These
*             values are valid for 1988-1994. (REAL* OUTPUT)
*     dpsi_prec and deps_prec  - contributions to the
*             nutations in longitude and obliquity due changes in
*             the precession constant and rate of change of
*             obliquity (mas) (REAL* OUTPUT).
*     dpsi_tot and deps_tot    - total nutations in longitude
*             and obliquity including the correction for the precession
*             constant (when precession is computed using the IAU 1976
*             precession constant), and are obtained by summing all
*             of the above corrections (mas) (REAL* OUTPUT).
*     dpsi_iau and deps_iau    - Nutation angles computed with the
*             iau_1980 nutation series.


      real*8 jd, dpsi_ls, deps_ls, dpsi_plan, deps_plan,
     .    dpsi_fcn ,  deps_fcn, dpsi_prec, deps_prec,
     .    dpsi_tot, deps_tot, dpsi_iau, deps_iau

*     start_jd, end_jd -- Start and end Julian dates
*     step_jd          -- Step size for Julian dates.

      real*8 start_jd, end_jd

      integer*4 ix(4),iy(4),npjd,npr,pjd1,pjd2,int,iout
     .        , i,j

      integer*4 len_run    ! Length of run-string
     .,         trimlen    ! Length og string
     .,         ref_year   ! Reference year for table
     .,         iau80      ! Index to IAU80 sting (zero when IAU2000 output)
     .,         date(5)    ! Calender date
     .,         rcpar      ! Function to read runstring
     .,         ierr       ! IOSTAT Error

      real*8 sectag        ! Seconds tag on date.

      real*4 outunt

      character*80 title
      character*24  varfmt

      character*80  ofile   ! Output file name
      character*80  runstring  ! Runstring


      data          varfmt/'(1x,i5,8i8,8x,i2)       '/
      data npr/4/
      data outunt/1.e-4/


****  Get the year of the table from the runstring
      len_run = rcpar(1,runstring)
      if( len_run.eq.0 ) then
          write(*,100)
 100      format('NUTIAU200: Make nutabl. for IAU2000',/,
     .           'Usage: nutIAU2000 <year> <IAU80>')
          stop
      endif
      read(runstring,*) ref_year

*     See if IAU80 option passed
      len_run = rcpar(2,runstring)
      iau80 = index(runstring,'IAU80')
      if( iau80.gt.0 ) then
          write(ofile,120) 'IAU80', ref_year
 120      format('nutabl.',a,'.',i4.4)
          write(title,125) ref_year-1, ref_year+1
 125      format('Nutation ephemeris from nutIAU200 for Oct ',
     .           I4,' - Mar ',i4,' IAU1980') 

      else
          write(ofile,130) 'IAU00', ref_year
 130      format('nutabl.',a,'.',i4.4)
          write(title,135) ref_year-1, ref_year+1
 135      format('Nutation ephemeris from nutIAU200 for Oct ',
     .           I4,' - Mar ',i4,' IAU2000') 
      endif
      write(*,150) ofile(1:trimlen(ofile))
 150  format('Creating file ',a)

      open(iout,file=ofile,iostat=ierr,status='unknown')
      call report_error('IOSTAT',ierr,'open',ofile,1,
     .                  'nutIAU200')

****  Get the start and stop JD
      date(1) = ref_year - 1
      date(2) = 10
      date(3) = 12
      date(4) = 00
      date(5) = 00
      sectag  = 0.0d0
      call ymdhms_to_jd( date, sectag, start_jd)
      end_jd = start_jd + 559.d0


      pjd1 = start_jd + 0.5d0
      pjd2 = end_jd +0.5d0
      int = -1

      write(iout,10) title,varfmt,pjd1,pjd2,npr,int,outunt
   10 format(a80,/,a24,11x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0)


****  Now loop over values
      do i = 1, 559, 2
         npjd = start_jd + (i-1) + 0.5d0 - 2400000.d0
         do j = 1, 4
            jd = start_jd + (i-1) + (j-1)*0.5d0

            call MHB_2000( jd, dpsi_ls, deps_ls,
     .                         dpsi_plan, deps_plan,
     .                         dpsi_fcn , deps_fcn ,
     .                         dpsi_prec, deps_prec,
     .                         dpsi_tot , deps_tot    )

            call IAU_1980 ( jd, dpsi_iau, deps_iau )

            if( iau80.gt.0 ) then
                ix(j) = nint(dpsi_iau*10)
                iy(j) = nint(deps_iau*10)
            else
                ix(j) = nint(dpsi_tot*10)
                iy(j) = nint(deps_tot*10)

            endif
         end do
        
         write(iout,varfmt) npjd,(ix(j),iy(j),j=1,4)
      end do

****  Thats all
      close(iout)


C     i = 0
C  20 k = 0
C  30 i = i + 1
C     if( i.gt.n ) goto 50
C     k = k + 1
C     xx = x(i)*units/outunt
C     yy = y(i)*units/outunt
C     ix(k) = idint(dsign(dabs(xx)+0.5,xx))
C     iy(k) = idint(dsign(dabs(yy)+0.5,yy))
C     npjd(k) = ifix(t(i))
C     if( k.ne.4 ) goto 30
C     write(iout,varfmt) npjd(1),(ix(j),iy(j),j=1,4)
C     goto 20
C  50 continue
C     if( k.eq.0 ) goto 70
C     write(iout,varfmt) npjd(1),(ix(j),iy(j),j=1,k)
C     write(6,60) k
C  60 format(/,1x,'Need to edit nutation table to add ',i1, 'in col 80')


      end

