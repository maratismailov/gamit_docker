      subroutine rdhhed(ih,idat,times,qnum_sites,qnum_sats,qsite_code
     .            ,qsite_names,cfiles,rcvr,antenna,aphs_mod,elev_cutoff
     .            ,nzen,ant_off,phs_off,prns,num_par,solv_vers,gps_utc
     .            ,sol_type,use_sol)

c Purpose: to read the hfile header, and extract from it the necessary
c          information for a sinex file. The data are returned in arrays
c          to be written out in the correct order at a later stage.
c
c NOTE: some code taken from htoglb to read the hfile
c
c Paul Tregoning
c 14th July, 1995
c
c IN:       ih         - unit number of hfile                         I*4
c           idat       - unit number of sinex.dat                     I*4
c           use_sol    - soln to read from hfile                      C*3
c
c OUT:   times         - runtime,start,stop,IC,eop,mean time of hfile C*12(6)
c                        (mean time is midpoint of observations used
c                         in the solve solution)
c        qnum_sites    - number of sites in hfile                     I*4
c        qnum_sats     - number of satellites in hfile                I*4
c        qsite_code    - 4 char site codes                            C*4(qnum_sites)
c        qsite_names   - long site names                              C*20(qnum_sites)
c        cfiles        - list of cfile names used in hfile soln       C*10(qnum_sites)
c        rcvr          - type, s/# and firmware version               C*16(qnum_sites,3)
c        antenna       - type, s/#                                    C*16(qnum_sites,2)
c        aphs_mod      - antenna phase centre model used in soln      C*4(qnum_sites)
c        elev_cutoff   - elevation cutoff used at each site           R*8(qnum_sites)
c        nzen          - number of zenith parameters at each site     I*4(qnum_sites)
c        ant_off       - UNE eccentricity mark to ARP                 R*8(qnum_sites,3)
c        phs_off       - UNE ecc. mark to Phase centres (L1 then L2)  R*8(qnum_sites,6)
c        prns          - list of PRN's in hfile                       I*4(qnum_sats)
c        num_par       - number of parameters estimated in hfile soln I*4
c        solv_vers     - solve version from hfile header              C*48
c        gps_utc       - gpst - utc time                              R*8
c        sol_type      - bias free/fixed constrained/free soln
c                        indicator (col1: constraints;
c                                   2: coord/vel soln3: eop; 4: bias) C*4
c        use_sol       - name of solution to use from hfile (eg glr)  C*3

      implicit none

      include '../includes/dimpar.h'

      integer i,ih,idat,ierr,num_in_file,indx,trimlen,qnum_sites
     .        ,qnum_sats,prns(maxsat),nzen(maxsit),num_par,yr,mo,day,hr
     .        ,min,params,epo_s,epo_f,obs_int,red_site,date(5)

      character*256 hf_ext
      character*1 use_site(maxsit)
      character*3 hfile_type,datum_type,use_sol,owner
      character*4 aphs_mod(maxsit),qsite_code(maxsit),sol_type
      character*5 frame,precmod,radmod
      character*10 tfile,cfiles(maxsit)
      character*12 times(6)
      character*16 rcvr(maxsit,3),antenna(maxsit,2)
      character*20 qsite_names(maxsit)
      character*48 solv_vers
      character*256 line

      real*8 sectag,qrun_time(7),elev_cutoff(maxsit)
     .       ,ant_off(maxsit,3),phs_off(maxsit,6),sec,ut1,dut1,xp,dxp
     .       ,yp,dyp,gps_utc

      logical finished


***   Read off the header at the start
      num_in_file = 0
      ierr = 0
      finished = .false.
*                                             ! Blank line
      read(ih,'(a)', iostat=ierr) line
*
*                                             ! GLOBALQ file line
c PT950714: only handle the case for a GPS file - others can be coded
c           later if required using Tom's other subroutines
      read(ih,'(a)', iostat=ierr) line
      if(ierr.ne.0) finished = .true.

      if( index(line,'GLOBALQ').eq.0 ) then
        print*,' This is not a GPS hfile - I cannot handle this file!'
        stop ' stop in RDHHED'
      else
        hfile_type = 'GPS'
      endif

***** skip next line - qfile name

      read(ih,'(a)', iostat=ierr ) line

***** Get some more header information
      if( .not.finished ) then
          read(ih,'(a)', iostat=ierr) line
          solv_vers = line(1:48)

*         Get the runtime from next line
          read(ih,'(a)', iostat=ierr) line
          read(line,200,iostat=ierr) date, sectag
  200     format(15x,i4,1x,i2,1x,i2,2x,2(i2,1x),f3.0)
          do i = 1,5
              qrun_time(i) = date(i)
          end do
          qrun_time(6) = nint(sectag)
          qrun_time(7) = 0
          call conv_date(date(1),date(2),date(3),date(4),date(5)
     .                  ,sectag,times(1))

c PT950828: get the owner name of the hfile (= processing agency name)
          read(ih,'(8x,a3)')owner

*         Skip 1 line
          read(ih,'(a)', iostat=ierr) line

*         Get the datum.  Check to make sure geocentric
          read(ih,'(a)', iostat=ierr) line

*         Set the datum type to blank and see if we find
*         one which matches
          datum_type = ' '
          indx = index(line,'geocentric')
          if( indx.gt.0 ) datum_type = 'LLR'
          indx = index(line,'cartesian')
          if( indx.gt.0 ) datum_type = 'XYZ'

*                                     ! Not geocentric, get out
          if( trimlen(datum_type).eq.0 ) then
              write(*,250) line(1:max(1,trimlen(line)))
  250         format(' Datum NOT found. Datum line is:',/,
     .                a,/,' Stopping processing of this file')
              finished = .true.
          end if
      end if

****  We have now finished non-repeating part of header.  If still
*     OK get the rest. Now start looping over individual solutions


*         Skip next line
          read(ih,'(a)', iostat=ierr) line
*                                     ! Probably finished with this
          if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
*                                 ! file
              finished = .true.
*                                 ! Get the position of bias
          else
              indx = index(line,'biases')

*             if biases line found, read the next line:
              if( indx.ne.0 ) then
                  read(ih,'(a)', iostat=ierr) line
                  indx = 0
              end if

*             if biases not found, see if we can find keys:
              if( indx.eq.0 ) indx = index(line,'keys')

*                                     ! then problem, not correct line
              if( indx.eq.0 ) then

*                 Skip over the headers written by htoh (echo the lines
*                 to user just to remind them what has happened).
                  call skip_htoh( ih, line , ierr )
                  if( ierr.ne.0 ) then
                      write(*,300) num_in_file+1,
     .                    line(1:max(1,trimlen(line)))
  300                 format(' Bias/keys line not where expected in',
     .                       ' solution ',i4,'.  Line read is:',/,a)
                      finished = .true.
                  end if
              end if

* MOD TAH 930208: Check the keys line to see what type of hfile this
*             is.  Do only for GPS.
              if( index(line,'keys').gt.0 .and.
     .            hfile_type.eq.'GPS' ) then
                    call get_gps_hext( line, hf_ext )
              else
                  hf_ext = ' '
              end if
          end if

c PT950828:  check if this is the right solution in the hfile. If so, read on, if not
c            then skip through it until the right one is found
          ierr = 0
          do while (hf_ext(1:3).ne.use_sol)
          line = ' '

            write(*,'(a,a3,a)')'Skipping solution ',hf_ext(1:3),'......'

c   it is not the right solution in the hfile. Need to skip through the hfile
c   until the right one is found
            do while (line(2:6).ne.'keys:'.and.ierr.eq.0)
              read(ih,'(a)',iostat=ierr)line

              if(ierr.ne.0)then
                write(*,'(a,a3,a)')' Solution type ',use_sol
     .               ,' not found in hfile!'
                stop 'Stop in RDHHED'
              endif
            enddo
            call get_gps_hext( line, hf_ext )

          enddo
          call get_gps_hext( line, hf_ext )
          write(*,'(a,a3,a,/)')' Found solution ',use_sol,' in hfile.'

c PT950822:  determine what type of soln it was (eg bias fixed/free constrained/unconstr)
          if(hf_ext(2:2).eq.'c')then
c    it is a constrained solution
             sol_type(1:1) = '0'
          else
             sol_type(1:1) = '2'
          endif

c NOTE: Need to code for bias fixing also - sinex doesn't accommodate this yet!
          if(hf_ext(3:3).eq.'r')then
            sol_type(4:4) = 'r'
          else
            sol_type(4:4) = 'i'
          endif

c PT950822: set the cood/velocity soln flag to be coords all the time
          sol_type(2:2) = 'X'




****      Now start reading main part of file
          if ( .not. finished ) then

            num_in_file = num_in_file + 1	


*             The remaining headers are different between the GPS and
*             TER versions of the hfiles. Use different routines to
*             read the remaining header records.

c  skip down to where the line which has the number of sessions
                do while (line(12:20).ne.'sessions:')
                  read(ih,'(a)', iostat=ierr) line
                enddo

c  now read the number of sites, satellites from the relevent line

                read(ih,'(22x,i3)', iostat=ierr) qnum_sites
                read(ih,'(17x,i3)', iostat=ierr) qnum_sats

c  skip 1 line
                read(ih,'(a)', iostat=ierr) line

c  read ephemeris filename
                read(ih,'(6x,a10)', iostat=ierr) tfile

c  skip 1 line
                read(ih,'(a)', iostat=ierr) line

c  get station names,receiver info, antenna heights etc
                call rdsit_info(ih,idat,qnum_sites,qsite_code
     .                       ,qsite_names,cfiles,rcvr,antenna,aphs_mod
     .                       ,elev_cutoff,nzen,ant_off,phs_off,red_site
     .                       ,use_site)

c skip 2 lines then get sat PRNs
                read(ih,'(a)', iostat=ierr) line
                read(ih,'(a)', iostat=ierr) line

                read(ih,'(4x,a)', iostat=ierr) line
                read(line,*)(prns(i),i=1,qnum_sats)

c skip 1 line
                read(ih,'(a)', iostat=ierr) line

c read epoch numbers and observation interval
                read(ih,390, iostat=ierr)epo_s,epo_f,obs_int
390             format(57x,i4,4x,i4,14x,i5)

c  read start, stop, IC times and models
                read(ih,400, iostat=ierr)yr,mo,day,hr,min,sec,gps_utc
400             format(13x,4i4,i5,f9.4,20x,f4.1)
                call conv_date(yr,mo,day,hr,min,sec,times(2))

c  compute the start,end time of observations and the midpoint time
                call comp_time(epo_s,epo_f,obs_int,times)

c skip 1 lines
                read(ih,'(a)', iostat=ierr) line

c  read  IC times and models
                read(ih,410, iostat=ierr)yr,mo,day,hr,min,sec,frame
     .                                 ,precmod,radmod
410             format(13x,4i4,i5,f9.4,11x,a5,15x,a5,13x,a5)
                call conv_date(yr,mo,day,hr,min,sec,times(4))

c skip a line
                read(ih,'(a80)')line

c read EOP time
                read(ih,400, iostat=ierr)yr,mo,day,hr,min,sec
                call conv_date(yr,mo,day,hr,min,sec,times(5))

c  read eop values
                read(ih,420)ut1,dut1,xp,dxp,yp,dyp
420             format(3(37x,f13.6,f12.6,/))

c  now read down until we get to the parameter line
                params = 0
                do while (params.ne.1)
                  read(ih,'(a80)')line
                  if(line(1:17).eq.' Total parameters')then
                    read(line(45:49),'(i5)')num_par
                    params = 1
                  endif
                enddo

          endif

      return
      end
