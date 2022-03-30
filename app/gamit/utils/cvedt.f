      program cvedt
c
c Purpose: To generate autcln edit lines from vcview.out summary file,
c          output to be written to standard output
c
c S. McClusky simon@chandler.mit.edu April 1995
c
c input: vcview.sort file [ SORTED version of vcview.out ]
c
c output: autcln edit_site_sv lines to be uesd in an autcln.cmd file.
c
      implicit none
c
      include '../includes/dimpar.h'
c
      character*4 flag1,flag2
      character*4 site,old
      character*80 line,arg,in_file

      integer*4 lvcf,lout,n,ierr,iprn,iepc
     .         ,istart,iend,lastp,laste,ioerr
c
c     Logical unit numbers for; vcview.sort (lvcf), sandard out (lout)
c
      lvcf = 11
      lout = 6
c
c     Read input file off command line if given
c
      call rcpar(1,arg)
      in_file = arg
      if (arg .eq. " ") then
c
c     Write some junk/help to the screen
c
      write(*,10)
10    format(/,'########################################################
     .######',/,'              Program CVEDT V1.0             '
     .,/,'Purpose: To generate autcln edit lines for removing data'
     .,/,'         manually flagged/unweighted in cview'
     .,/,'Requirement: first type: sort vcview.out >! vcview.sort'
     .,/,'Usage:(1) Command Line Input (cvedt vcview.sort)'
     .,/,'      (2) Interactive (just answer the question below)'
     ./,'##############################################################'
     .,/)

        print*, ('Enter vcview.sort file name [vcview.sort]: ')
        read(*, '(a)') in_file
        if (in_file .eq. ' ') in_file = 'vcview.sort'
      endif
c
c     Open the required vcview.out input file.
c
      open(unit=lvcf,file=in_file,status='old',err=100,iostat=ioerr)
c
c     Write a nice little commented header!!!!
c
      write(lout,50)
50    format('## START OF CVIEW DATA EDITS ## ')
      goto 225
c
c     Some open file error trapping stuff!!!
c
100   call report_stat('FATAL','UTILS','cvedt',in_file,
     1  'Error opening input vcview.out file',ioerr)
c
225   continue
c
c     Read vcview.sort file and get; satellite prn number, site name,
c     and cview edit/unwt start and stop times...
c
      n = 0
      ierr = 0
      laste = 0
      do while (ierr .eq. 0)
        read(lvcf,'(a)',iostat = ierr ) line
        if(ierr.eq.0 .and. index(line,'unwt').gt.63 ) then
          n = n+1
          read(line,'(1x,a4,3x,i2,1x,i4,43x,a4,1x,a4)')
     .    site,iprn,iepc,flag1,flag2
          if( n .eq. 1) then
            old = site
            istart = iepc
            lastp = iprn
            laste = iepc
          elseif( site.eq.old  .and. laste.eq.iepc-1
     .           .and. lastp.eq.iprn ) then
            laste = iepc
          else
            iend = laste
            write(lout,250) old,lastp,istart,iend
250         format (' edit_site_sv ',a4,i3,2i6)
            old = site
            istart = iepc
            lastp = iprn
            laste = iepc
          endif
        endif
      enddo
      iend = laste
      write(lout,250) old,lastp,istart,iend
c
c     Write out a nice little footer line!!!!
c
      write(lout,750)
750   format('## FINISH OF CVIEW DATA EDITS ## ')
c
      stop
      end
