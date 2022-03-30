      subroutine wthead(ifile,key,mode,job,version,drvfil)
c
c     write head information to specified file
c     mode = 1: terse
c     mode = 2: verbose
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      character*16 version
      character*16 note1(10),uname
      character*2 sp,sps
      character*8 blank
      character*6 note,buf6
      character*96 line
      character*64 drvfil
      character*8 sname2(maxsit)
      integer mchkey,ifile,key,mode,job,i,j,j0,nblen,ii
      integer imn,iday,iyr,imon,ihr,isec,ijunk,fixidx
      
      note1(1) = ' Gauss-Markov   '
      note1(2) = ' Gauss-Helmert  '
      note1(3) = 'model coordinate'
      note1(4) = 'sequential updat'
      note1(5) = ' Kalman filter  '
      sp = '  '
      sps = '* '
      blank = '        '
c
c     header
      if (job.eq.1.and.mode.eq.2) then
         call getdat(iyr,imon,iday)
         call gettim(ihr,imn,isec,ijunk)
         call getusr(uname)
         write (ifile,*) sp,'=========================================='
         write (ifile,'(6x,a,2x,a)') ' SOLVEM version:',version
         write (ifile,*) sp,'=========================================='
         write (ifile,4) 'running time:',iyr,
     .         '/',imon,'/',iday,ihr,':',imn,':',isec
         write (ifile,*) sp,'operator:',sp,uname
         ijunk = nblen(drvfil)
         write (ifile,*) sp,'driving file:',sp,sp,sp,drvfil(1:ijunk)
c        copy informations from driving file
         rewind (8)
 40      read (8,'(a)') line
         if (line(1:1).ne.' ') goto 40
         i = mchkey(line,'end',10,3)
         if (i.gt.0) goto 40
         i = mchkey(line,'exit',10,4)
         j0 = nblen(line)
         if (i.lt.0) write (ifile,'(a)') line(1:j0)
         if (i.lt.0) goto 40
         if (i.ge.0) goto 100
      endif

      if (mode.eq.2) then
         write (ifile,*) sp,'Number of experiments:',sp,iexp
         write (ifile,*) sp,'Total sites from data file:',sp,nsit
         write (ifile,*) sp,'Ill sites:',sp,nsit-gdsit
         write (ifile,*) sp,'Bridge sites:',sp,nvsit
         write (ifile,*) sp,'Healthy sites:',sp,gdsit-nvsit
c        if (fixsit(1).gt.0) write (ifile,*) sp,'fixed site(s):',sp,
c    .      (sname(j),j=1,fixcnt)
         if (fixsit(1).gt.0) then
            do 60 j=1,fixcnt
                 fixidx=fixsit(j)
                 sname2(j)=blank
                 sname2(j)(1:4)=sname(fixidx)(1:4)
 60         continue
            write (ifile,*) sp,'fixed site(s):',sp, 
     .           (sname2(j),j=1,fixcnt)
         endif
         write (ifile,*) sp,'======================================'
      endif
 4    format(2x,a13,2x,i4,2(a1,i2),3x,2(i2,a1),i2)
 10   format(2x,a51,1x,4f8.4)
 20   format(2x,a45,1x,3f7.1,3f6.2)
c
      write (ifile,'(a2,a15,a2,a16)') sps,'solution model:',
     .   sp,note1(smode)
      do 2 i = 1,ntype
c     listyp is a list of the types of observations in the solution
      j0 = listyp(i)
      if (j0.eq.1) note1(i) = '<astro. azimuth>'
      if (j0.eq.2) note1(i) = '<horizon. angle>'
      if (j0.eq.3) note1(i) = '<horizon. dirn.>'
      if (j0.eq.4) note1(i) = '<basel. length >'
      if (j0.eq.11) note1(i) = '< azimuth rate >'
      if (j0.eq.12) note1(i) = '< angle  rate  >'
      if (j0.eq.14) note1(i) = '< length rate  >'
      if (j0.eq.21) note1(i) = '<GPS-coordinate>'
      if (j0.eq.22) note1(i) = '< GPS-velocity >'
      if (j0.eq.23) note1(i) = '<geodetic coor.>'
      if (j0.eq.24) note1(i) = '<geodetic velo.>'
      if (j0.eq.25) note1(i) = '<baseline vect.>'
      if (j0.eq.27) note1(i) = '<baseline vect.>'
      if (j0.eq.41) note1(i) = '<defle. correc.>'
 2    continue
      
      write (ifile,'(a2,a18)') sps,'Solution includes:'
      do 3 ii=1,ntype
         write (ifile,'(a2,a16)') sps, note1(ii)
 3    continue
      write (ifile,'(a2,a4)') sps,'Data'
      note1(1) = '                '
      if (smode.eq.3) note1(1) = ' is transformed '
      note1(2) = ' (    ) solution'
      if (key.eq.5.or.key.eq.7) then
         write (buf6,'(f6.2)') azio*rtod
         read (buf6,'(a6)') note
      endif
      if (fixsit(1).gt.0) then
         note1(2)(3:6) = 'FIX '
      else
         note1(2) = ' constrained sln'
      endif
      if (key.eq.1) write(ifile,'(''* full solution,''
     .   ,'' no null space'')')
      if (key.eq.2) write(ifile,'(''* fix site'',a16)') note1(2)
      if (key.eq.3) write(ifile,'(
     .    ''* fix center solution(both L & V)'')')
      if (key.eq.4) write(ifile,'(''* inner coordinate '',a16)')
     .              note1(1)
      if (key.eq.5) write(ifile,*) sp,'outer coordinate',sp,
     .              'minimize azi',note
      if (key.eq.6) write(ifile,'(''* fix center solution(V only)'')')
      if (key.eq.7) write(ifile,*) sps,'outer coordinate (C)',sp,
     .              'minimize azi',note
      if (nexc.eq.0) then
          write(ifile,'(a2)') sps
      else 
          write(ifile,*) sps,(sname(iesit(i)),sp,i=1,nexc),'excluded'
      endif
c
 100  continue
      return
      end

