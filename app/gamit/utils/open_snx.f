      subroutine open_snx(in_h,out_snx,ih,isnx,idat,istn)

c Purpose: open the input hfile and output sinex file
c
c by P Tregoning
c
c 14th July, 1995
c
c IN:    in_h     :  input hfile                                   C*10
c        out_snx  :  output sinex filename                         C*20
c        ih       :  hfile unit number                             I*4
c        isnx     :  sinex file unit number                        I*4
c        idat     :  unit # of file containing additional info     I*4
c        istn     :  unit # of station.info                        I*4

      implicit none

      integer ih,isnx,idat,istn,ierr,iscrn,iterm
      character*1 ans
      character*10 in_h
      character*20 out_snx

      ierr = 0
      iscrn = 6
      iterm = 5

c  open the hfile
      open(ih,file=in_h,status='old',iostat=ierr)
      if(ierr.ne.0)then
        write(iscrn,'(a,a10,a)')'Error opening hfile ',in_h
     .                        ,'. Does it exist?'
        stop 'Stop in OPEN_SNX'
      endif

c  open the sinex file
5     ierr = 0
      open(isnx,file=out_snx,status='new',iostat=ierr)
      if(ierr.ne.0)then
        write(iscrn,'(a,a20,a,a$)')'File ',out_snx,' exists. '
     .                           ,'Overwrite (y/n) ? '
        read(iterm,'(a1)')ans
        if(ans.eq.'y'.or.ans.eq.'Y')then
          open(unit=20,file=out_snx,status='unknown')
        else
          write(iscrn,'(a$)')' Enter output sinex filename '
     .                        ,'(e.g mit08061.snx) : '
          read(iterm,'(a)')out_snx
          goto 5
        endif
      endif

      ierr = 0
c  open additional data file
      open(idat,file='sinex.dat',status='old',iostat=ierr)
      if(ierr.ne.0)then
        write(iscrn,100)
100     format('File sinex.dat does not exist. You need a file of the'
     .        ,' format:',//,'-------------------------------------',/
     .        ,' Processing centre: SIO',/
     .        ,'+FILE/REFERENCE',/
     .        ,' DESCRIPTION        SIO GAMIT solution ',/
     .        ,' CONTACT            paul@pgga.ucsd.edu ',/
     .        ,' INPUT              SIO daily solutions',/
     .        ,'-FILE/REFERENCE                  ',//
     .        ,'+SITE_CODE/CHANGES               ',/
     .        ,' ARE1 AREQ                       ',/
     .        ,' KOKR KOKB                       ',/
     .        ,'-SITE_CODE/CHANGES               ',//
     .        ,'+DOMES NUMBER (only if it exists!)',/
     .        ,' SITE  __DOMES__                 ',/
     .        ,' ALBH  40129M003                 ',/
     .        ,' ALGO  40104M002                 ',/
     .        ,'-DOMES                           ',/
     .        ,'----------------------------------------',/
     .        ,' Copy test between dashed lines to "sinex.dat"',// )

        stop ' Stop in OPEN_SNX.'
      endif

c  open station.info file
      ierr=0
      open(istn,file='station.info',status='old',iostat=ierr)
      if(ierr.ne.0)then

c PT950825: temporarily proceed without station.info as it isn't being read anyway
c        stop 'OPEN_SNX: No station.info file found!'
        print*, 'OPEN_SNX: No station.info file found - not used anyway'
      endif

      return
      end
