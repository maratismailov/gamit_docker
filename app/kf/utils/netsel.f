      program netsel

      implicit none 

*     Program to use a returned ftp'd directory and a globk style 
*     velocity file to make a selection of networks to process.
*     The output of the program will be used to generate entries
*     for the sites.defaults file used by sh_gamit.  Only sites
*     that appear in the .vel file will be included in the networks.
*     For sites that do not appear in the .vel, will put in a 
*     sh_script that will ftp the data files and get the coordinates
*     for the site.
*
      include 'netsel.h'

***** Decode the runstring for the program
      call net_init
      call get_net_run

****  Now read the site postion (vel) file
      call read_net_vel

***** Now read the ftp log file
      if( globk_US ) then
          call read_net_rwf
      else
          call read_net_ftp
      endif 

***** Now do the base selection.  Each site is used only once
*     in this selection, and the networks are based on highest
*     density of sites
      call base_net_sel

****  Now add the tie site sites bewteen the networks. Here we
*     look at the common sites that are needed to tie the networks
*     together
      call tie_net_sel

****  Finally, output the network selections in a form that can be 
*     used by sh_gamit.
      if( .not. globk_US) call check_stinf

      call out_nets

****  Thats all
      end

CTITLE NET_INIT

      subroutine net_init

      implicit none 

*     Routine to initialize variables for netsel

      include 'netsel.h'

*     Set defaults
      getfile = " "
      nom_net_site = 40
      nom_tie_site = 3
      stinffile = '../tables/station.info'
      net_code  = 'ne'
      maxuse_rw = 2.0d0
      rwfile = ' '
      globk_US = .false.   ! Default to standard operation.

      end

CTITLE GET_NET_RUN

      subroutine get_net_run

      implicit none 

*     Routine to get the options for netsel

      include 'netsel.h'

* LOCAL VARIABLES
      integer*4 i  ! Loop counter over runstring 
     .,         lenrun  ! Length of runstring return
     .,         ierr    ! IOSTAT error
     .,         trimlen ! Lenth of string
     .,         rcpar   ! Function to get runstring

      character*128 runstring ! Entry from runstring    
      character*8 option  ! Option passed in runstring

***** Initialize
      i = 0
      lenrun = 1

***** Get the next argument passed
      do while ( lenrun.gt.0 )
         i = i + 1
         lenrun = rcpar(i,option)
         call casefold(option)
         if( lenrun.gt.0 ) then

*           FTP log name
            if( option(1:2).eq.'-F' ) then
                i = i + 1
                lenrun = rcpar(i,ftplog)
*           Name of velocity file for coorindates
            else if ( option(1:2).eq.'-V' ) then
                i = i + 1
                lenrun = rcpar(i,velfile)
            else if ( option(1:2).eq.'-N' ) then
                i = i + 1
                lenrun = rcpar(i,runstring)
                read(runstring,*,iostat=ierr) nom_net_site
                call report_error('IOSTAT',ierr,'decod',
     .               runstring,1,'GET_NET_RUN/NOM_NET_SITE')
                if( nom_net_site.gt.max_net_site -1 ) 
     .              nom_net_site = max_net_site - 1
*           Name of get file (list of sites to ftp to get 
*           new coordinates).
            else if ( option(1:2).eq.'-G' ) then
                i = i + 1
                lenrun = rcpar(i,getfile)
            else if ( option(1:2).eq.'-T' ) then
                i = i + 1
                lenrun = rcpar(i,runstring)
                read(runstring,*,iostat=ierr) nom_tie_site
                call report_error('IOSTAT',ierr,'decod',
     .               runstring,1,'GET_NET_RUN/NOM_TIE_SITE')
                if( nom_tie_site.gt.4 )
     .              nom_tie_site = 4  
            else if ( option(1:2).eq.'-S' ) then
                i = i + 1
                lenrun = rcpar(i,stinffile)
            else if ( option(1:2).eq.'-C' ) then
                i = i + 1
                lenrun = rcpar(i,net_code)
                if( net_code(2:2).eq.' ' ) net_code(2:2) = '_'
            else if ( option(1:2).eq.'-R' ) then
                i = i + 1
                lenrun = rcpar(i,rwfile)
                i = i + 1
                lenrun = rcpar(i,runstring)
                if( lenrun.gt.0 ) read(runstring,*,iostat=ierr) 
     .                            maxuse_rw
                call report_error('IOSTAT',ierr,'decod',
     .               runstring,1,'-RW maxuse_rw option')
                if( ierr.ne.0 ) maxuse_rw = 2.0d0
                globk_US = .true.

            endif
         end if
      end do
*
*     See if only 1 entry read
      if( i.le.1 ) then
          call proper_runstring('netsel.hlp','NETSEL',1)
      endif

****  Report the runstring elements
      if( .not. globk_US ) then 
         write(*,210) ftplog(1:max(1,trimlen(ftplog)))
 210     format('NETSEL: ',/,
     .          'FTPLOG:  ',a)
         write(*,220) velfile(1:max(1,trimlen(velfile)))
 220     format('VELFILE: ',a)
         write(*,230) nom_net_site
 230     format('Number of sites per net: ',i4)
      else
         write(*,310) trim(rwfile), maxuse_rw
 310     format('NETSEL : ',/,
     .          'RW FILE:  ',a,' MAX RW to use ',F6.2,' mm^2/yr')
         write(*,320) velfile(1:max(1,trimlen(velfile)))
 320     format('VELFILE: ',a)
         write(*,330) nom_net_site
 330     format('Number of sites per net: ',i4)
      endif

      return
      end

CTITLE READ_NET_VEL

      subroutine read_net_vel

      implicit none 

*     This routine reads the long/lat and site names from a 
*     standard globk velocitiy file.

      include 'netsel.h'

* LOCAL VARIABLES
      integer*4 ierr ! IOSTAT reading file
     .,         jerr ! IOSTAT decoding lines
     .,         ns   ! counter for number of sites
     .,         trimlen ! Length of string function

      real*8 values(12) ! 12 Values read from velocity file

      character*8 name  ! Name of site from velfile
      character*256 line ! Line read from file

***** Open the file
      open(50,file=velfile,status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',velfile,
     .         1,'READ_NET_VEL')

*     Now loop over the file
      ns = 0
      do while ( ierr.eq.0 )
         read(50,'(a)',iostat=ierr) line
         if( line(1:1).eq.' ' .and. trimlen(line).gt.0 .and.
     .       ierr.eq.0 ) then
             read(line,*,iostat=jerr) values, name
             if( jerr.eq.0 ) then
                 ns = ns + 1
                 if( values(1).gt.0 ) then
                     slon(ns) = values(1)
                 else
                     slon(ns) = values(1) + 360
                 endif
                 slat(ns) = values(2)
                 call casefold(name)
                 snam(ns) = name
             endif
         end if
      end do

****  Report to user
      num_site_vel = ns
      write(*,210) velfile(1:trimlen(velfile)), num_site_vel
 210  format('NETSEL: ',a,' contains ',i5,' sites',/,
     .       'Site#    Long      Lat     Name')
!     do ns = 1, num_site_vel
!        write(*,220) ns, slon(ns), slat(ns), snam(ns)
!220     format(I5,1x,2F10.6,1x,a)
!     enddo 
      close(50)

      return
      end

CTITLE READ_NET_FTP

      subroutine read_net_ftp

      implicit none 

*     Routine to read the ftp log file. (The log is externally
*     generated with a shell script)

      include '../includes/const_param.h'
      include 'netsel.h'

      integer*4 ierr ! IOSTAT reading file
     .,         jerr ! IOSTAT decoding lines
     .,         ns   ! counter for number of sites
     .,         ng   ! Number of sites to get
     .,         trimlen ! Length of string function
     .,         indx ! Position in string
     .,         j    ! Loop counter
     .,         lenw    ! Length of string

      real*8 xyz(3)  ! Cooridinates
     .,      rot_mat(3,3)  ! Rotation matrix from XYZ_to_GEOD
     .,      geod(3) ! Co-lat, long, height

      logical gf     ! Set true if we are getting new coordinates
     .,       found  ! Set true when site found

      character*128 word  ! Word read from line
     .,             line  ! Line read from file
      character*8 cword   ! Casefolded version of site name
*
****  Open the log file
      open(50,file=ftplog,status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',ftplog,1,
     .                  'READ_NET_FTP')

      ns = 0
      ng = 0 

      minlon = 370
      maxlon = -10
      minlat = 91
      maxlat = -91

*     See if we allowing updates to coordinates
      if( trimlen(getfile).gt.0 ) then
          open(60,file=getfile,status='unknown',iostat=ierr)
          gf = .true.
      else
          gf = .false.
      endif

****  Loop over the file
      do while ( ierr.eq.0 ) 
         read(50,'(a)',iostat=ierr) line
         if(  trimlen(line).gt.0  .and.ierr.eq.0 ) then
*            Decode the lines
             indx = 0
* Line looks like:
*13028 ../cors.rinex/acu11000.02o
             call getword(line, word, indx)   ! File size. 13028

             call getword(line, word, indx)   ! ../cors.rinex/acu11000.02o

*****        Pull off the site name (end of word must be d or o.
             lenw = trimlen(word)
             if( word(lenw:lenw).eq.'d' .or. 
     .           word(lenw:lenw).eq.'o' ) then   ! OK file name
                 cword = word(lenw-11:lenw-8)
                 call casefold(cword)
                 
*****            See if we can file this site.
                 found = .false.
                 j = 0
                 do while ( .not.found .and. j.lt.num_site_vel)
                    j = j + 1
                    if( cword.eq.snam(j)(1:4) ) then
                       found = .true.

*                      Extract the site 4-char code
                       ns = ns + 1
                       fcod(ns) = cword(1:4)
                       ffile(ns) = word
                       flon(ns) = slon(j)
                       flat(ns) = slat(j)
                    end if
                 end do

*                If site not found; get the coordinates from the rinex
*                file.
                 if( .not.found ) then
                    write(*,220) cword(1:trimlen(cword)), 
     .                           word(1:trimlen(word))
 220                format('* Coordinates not found for ',a,
     .                     ' Rx file ',a)
                    open(51,file=word,iostat=jerr,status='old')
                    call report_error('IOSTAT',jerr,'open',word,1,
     .                   'read_ftp_net')
                    found = .false.
                    do while ( .not.found .and. jerr.eq.0 )
                        read(51,'(a)', iostat=jerr ) line
                        if( index(line,'APPROX POSITION').gt.0 ) then
                            read(line,*) xyz
                            call xyz_to_geod(rot_mat,xyz,geod)
                            num_site_vel = num_site_vel + 1
                            j = num_site_vel
                            slon(j) = geod(2)*180.d0/pi
                            slat(j) = (pi/2-geod(1))*180.d0/pi
                            snam(j) = cword
                            found = .true. 
                            ns = ns + 1
                            fcod(ns) = cword(1:4)
                            ffile(ns) = word
                            flon(ns) = slon(j)
                            flat(ns) = slat(j)
                            rw_stats(ns) = 1.d0    ! Default value
                            write(*,225) fcod(ns), xyz, 
     .                                   flon(ns), flat(ns)
 225                        format('+APR ',a4,'_GPS ',3F14.3,3x,
     .                             2f9.5,' deg')
                        endif
                    end do 
                    close(51)
                end if

*               Check range
                if( flon(ns) .gt. maxlon ) maxlon =  flon(ns)
                if( flon(ns) .lt. minlon ) minlon =  flon(ns)
                if( flat(ns) .gt. maxlat ) maxlat =  flat(ns)
                if( flat(ns) .lt. minlat ) minlat =  flat(ns)
*            See if we found the site
             else
                write(*,320) word(1:trimlen(word))
 320            format('Not Valid Rinex file ',a)             
             end if               
          end if
      end do

*     Now report what we found
      num_site_ftp = ns
      write(*,310) ftplog(1:trimlen(ftplog)), num_site_ftp
 310  format('NETSEL: ',a,' contains ',i5,' sites')
      write(*,330) minlon, maxlon, minlat, maxlat
 330  format('Site Range Long ',2f10.4,' Latitude ',2f10.4,' deg')
      close(50)

      return
      end

CTITLE BASE_NET_SEL

      subroutine base_net_sel

      implicit none 

*     Routine to make base selection for each network.  
*     Basic algorithm is center nets on the highest density 
*     regions in the network.

      include 'netsel.h'

* LOCAL VARIABLES
      integer*4 i,j    ! Loop counters
     .,         nn     ! Temporary variable for number of sites 
                       ! per network

*     Clear the network array
      do i = 1, num_site_ftp
         do j = 1, 4
            network(j,i) = 0
         end do
      end do
      num_net = 0

****  Refine the number of sites per network based on total number
*     of stations.
      num_net = num_site_ftp/(nom_net_site-1) + 1
      num_net_site = num_site_ftp/num_net + 1

*     See if we need to adjust number of tie sites
      if( nom_tie_site*num_net.gt.max_net_site ) then
          nom_tie_site = max_net_site/num_net
      endif

      write(*,120) num_site_ftp, nom_net_site, num_net_site, num_net,
     .             (num_site_ftp - (num_net-1)*num_net_site), 
     .             nom_tie_site
 120  format('NETSEL: For ',i4,' sites, with nominal ',i4,
     .       ' sites per network, final selection is:',/,
     .       'NETSEL: Fin ',i4,' sites in ',i4,' networks with ',
     .       i4,' sites in one network',/,
     .       'NETSEL: Number of tie sites ',i3)

***** Now loop over networks
      do i = 1, num_net
         if( i.eq.num_net ) then
            nn = (num_site_ftp - (num_net-1)*num_net_site)
         else
            nn = num_net_site
         end if

*****    select the actual network
         call sel_net(i,nn)

*        Now print out out networks
         call rep_net(i,nn)
      end do

      return
      end

CTITLE REP_NET

      subroutine rep_net(i,nn)

      implicit none 

*     Routine to report the size of a network.

      include 'netsel.h'

* PASSED VARIABLES

      integer*4 i  ! Network number 
      integer*4 nn ! Number of sites in the network

* LOCAL VARIABLES
      integer*4 j  ! Loop counter
      integer*4 nx ! Counter of number of sites

      write(*,120) i, nn
 120  format('#NETWORK Number ',i3.3,' with ',I3,' sites',/,
     .       '# NN    #      Long         Lat     Name  RK',
     .       '  RWStat')
      nx = 0
      do j = 1, num_site_ftp
         if( network(1,j).eq.i) then
              nx = nx + 1
              write(*,220) i, nx, flon(j), flat(j), fcod(j),
     .                     network(2,j), rw_stats(j)
 220          format('# ',i3.3,1x,i3,2F13.5,1x,a4,1x,i3,1x,F7.3)
         end if
      end do

***** Thats all
      return 
      end
      
CTITLE SEL_NET

      subroutine sel_net(n,nn)

      implicit none 

*     Routine to select the n'th network with nn sites in it.

      include 'netsel.h'

* PASSED VARIABLES

      integer*4 n    ! network number
     .,         nn   ! Number of sites in network
* LOCAL VARIABLES

      real*8 net_den  ! function to compute network density
     .,      netdist  ! function to compute distance of site from
                      ! reference point
     .,      ndist    ! distance to current point
     .,      dens     ! Density about current reference point
     .,      id, jd   ! Long and Lat of densist point in network
     .,      netlens(max_net_site) ! Lengths of current network
                      ! distances 
     .,      max_dens ! Maximum density of network

      integer*4 i,j   ! Loop counters
     .,         netsites(max_net_site)  ! List of sites on current
                      ! current network ranked by closest distance

***** Loop over a 1-deg grid and at each grid point compute the
*     site density.  Select the site density that is the highest
*     use only non-assigned sites.
      max_dens = 0 
      do i = int(minlon)-1, int(maxlon)+1
         do j = int(minlat)-1, int(maxlat)+1
            dens = net_den(i,j)
            if( dens.gt.max_dens ) then
                id = i
                jd = j
                max_dens = dens
            end if
         end do
      end do

***** Now we have max_dens location, select the group of closest
*     stations
      do i = 1, nn
         netsites(i) = 0
      end do

      do i = 1, num_site_ftp
*        Only consider a site if it has not been assigned to a
*        network yet
         if( network(1,i).eq.0 ) then
             ndist = netdist(id,jd,i)
             call netrank(ndist, i, netlens, netsites,nn)
         end if
      end do

*     The clostest sites to the reference point have been selected
*     now assign them to the network
      do i = 1, nn
         network(1,netsites(i)) = n
         network(2,netsites(i)) = i
      end do

***** Thats all
      return
      end 

CTITLE NET_DEN

      real*8 function net_den(i,j)

*     Functiont to compute the density of the network about a
*     specific long (i) and lat (j) location.

      include 'netsel.h'

* PASSED VARIABLES

      integer*4 i,j  ! Integer long and lat (deg)

* LOCAL VARIABLES
      real*8 ri, rj  ! Long and lat as floating point values
     .,      netdist ! Functiont to compute distance
     .,      ndist   ! Distance to reference point

      integer*4 k    ! Loop variables

***** Loop over non-used sites accumulating the density
      ri = i
      rj = j

      net_den = 0
      do k = 1, num_site_ftp
         if( network(1,k).eq.0 ) then
             ndist = netdist(ri,rj,k)
             net_den = net_den + 1.d0/(ndist+0.01)
         end if
      end do

***** Thats all
      return
      end

CTITLE NETDIST

      real*8 function netdist(rlng, rlat, ns)

*     Function to compute distance (in radians)

      include '../includes/const_param.h'
      include 'netsel.h'

* PASSED VARIABLES
      real*8 rlng, rlat  ! Reference Long and lat

      integer*4 ns       ! Site number

* LOCAL VARIABLES
      real*8 cd   ! Cosine of distance

****  Use cosine rule
      cd = sin(rlat*pi/180.d0)*sin(flat(ns)*pi/180.d0)+
     .     cos(rlat*pi/180.d0)*cos(flat(ns)*pi/180.d0)*
     .     cos((rlng-flon(ns))*pi/180.d0)

      netdist = acos(cd)
      return
      end

CTITLE NETRANK

      subroutine netrank(ndist, i, netlens, netsites,nn)

      implicit none 

*     Routine to add the current site into the ranked list

* PASSED VARIABLES
      integer*4 nn   ! Number of sites in network
      integer*4 i    ! Current site number
      integer*4 netsites(nn)  ! List of sites in network (0
                     ! for site number says no site assigned yet)

      real*8 ndist   ! Current site distance from reference site
      real*8 netlens(nn)  ! Distances in current list

* LOCAL VARIABLES
      integer*4 npos   ! Position of current site in list based on
                     ! distances
      integer*4 mn   ! Number of sites in list so far
      integer*4 j    ! Loop counter

****  Loop over the current sites in the network and place this
*     site in list
      npos = 0
      mn = 0
      do j = 1, nn
         if( netsites(j).gt.0 ) then
            mn = mn + 1
*           OK Site exists at this length, see if this site 
*           is closer
            if( ndist.lt.netlens(j) .and. npos.eq.0 ) then
                npos = j
            endif
         end if
      end do

****  Now put the site in the list.
      if( npos.eq.0 .and. mn.eq.0 ) then
*        List is empty and this is first
         netsites(1) = i
         netlens(1) = ndist
         mn = mn + 1
      else if( npos.eq.0 ) then
*        This length is longer than all in list, if the
*        list is not full yet, add to the end
         if( mn.lt.nn ) then
            netsites(mn+1) = i
            netlens(mn+1) = ndist
            mn = mn + 1
         end if
      else 
*        site fall inside the current list, move list up 
*        and add the site
         do j = nn-1, npos, -1
            if( netsites(j).gt.0 ) then
                netsites(j+1) = netsites(j)
                netlens(j+1) = netlens(j)
            end if
         end do
         if( mn.lt.nn) mn = mn + 1
         netsites(npos) = i
         netlens(npos) = ndist
      end if

*     Thats all
      return
      end

CTITLE TIE_NET_SEL

      subroutine tie_net_sel

      implicit none 

*     Rouitine to select the tie sites between the network
      
      include 'netsel.h'
      include '../includes/const_param.h' 

      integer*4 tie_sel(4)  ! list of site positions for the 
                            ! tie sites
      integer*4 i,j,k,n,m      ! loop counters
      integer*4 nx,ns          ! Site counter
      integer*4 nn          ! Number of sites in last non-tie network
      integer*4 inda_max(3) ! Sites that form largest weighted area
      integer*4 inda        ! Index for 4th site 
      integer*4 xyzi(max_sites)  ! Index to sites
      real*8 xyzr(4,max_sites)  ! XYZ and random walk stats (unit sphere, mm^2/yr)
      real*8 veca(3), vecb(3), vecc(3), sang, area_max, darea
      real*8 vecfin(3)      ! Cross product for main triangle
      real*8 dot3    ! 3D real*8 dot product

      real*8 neq_save(6,6), neq(6,6), A(3,6), cov(6,6), cov_save(6,6)
      real*8 neq_base(6,6)

      real*8 bvec(6), scale(6), trace, trace_save, trace_min
      integer*4 ipivot(6)

****  Set up the selection of distance entries that we need
      tie_sel(1) = 1
      if ( nom_tie_site.ge.2 ) then
          tie_sel(2) = num_net_site
      endif
      if ( nom_tie_site.ge.3 ) then
          tie_sel(3) = num_net_site/2 
      endif
      if ( nom_tie_site.ge.4 ) then
          tie_sel(4) = num_net_site*0.75
      endif

****  Now loop over the nets building the tie list
*     Different tie algorithm for GLOBK Verus GAMIT
      if( globk_US ) then
         write(*,'(a,1x,I3)') 'GLOBK Network Tie Network',num_net
*        GLOBK Algorithm: Use fit to network geometry
         bvec = 0
         num_net = num_net + 1
         do i = 1, num_net-1
*           Collect all the sites in this network and convert
*           to XYZ for triangle compute
            ns = 0
            do k = 1, num_site_ftp
                if( network(1,k).eq.i ) then ! this site in network
                    ns = ns + 1
                    xyzr(1,ns) = cos(flon(k)/rad_to_deg)*
     .                                        cos(flat(k)/rad_to_deg)
                    xyzr(2,ns) = sin(flon(k)/rad_to_deg)*
     .                                        cos(flat(k)/rad_to_deg)
                    xyzr(3,ns) = sin(flat(k)/rad_to_deg)
*                   Add ~median to statistics so that small values don't
*                   dominate weight.
                    xyzr(4,ns) = rw_stats(k)+0.1d0
                    xyzi(ns)   = k    ! Original site number
                endif
            end do
*           We now have the ns sites in this network
*           Now we need to compute the sites that generate
*           maximum weighted area.
            trace_min = 1.e10

            do j = 1,ns-2
               neq_base = 0
               call incr_norm(xyzr(:,j), A, neq_base)
               do k = j+1,ns-1
                  neq_save = neq_base
                  call incr_norm(xyzr(:,k), A, neq_save)
                  do n = k+1,ns
                     neq = neq_save
                     call incr_norm(xyzr(:,n), A, neq)
                     cov = neq
                     call invert_vis( cov, bvec, scale, ipivot, 6,6,0)
                     trace = 0
                     do m = 1,6
                         trace = trace + cov(m,m)
                     enddo 

                     if( trace.lt.trace_min .and. trace.gt.0 ) then
                         trace_min = trace
                         inda_max(1) = j
                         inda_max(2) = k
                         inda_max(3) = n
!                        write(*,999) i, inda_max, xyzi(inda_max), 
!    .                            trace_min
!999                     format('GLOBK Net ',i3,' Loc Site ',3i5,
!    .                          ' Glb Site ',3i5,' Trace ',E13.5)
                     end if
                   end do
               end do
            end do
*           Save network information for this network
            write(*,110) i, inda_max, fcod(xyzi(inda_max)), 
     .               trace_min
 110        format('GLOBK Net ',i3,' Loc Site ',3i5,
     .             ' Glb Site ',3a5,' Trace ',E13.5)
            do j = 1,3
                 network(3,xyzi(inda_max(j))) = num_net  
            end do
        end do

      else
*        GAMIT Algorithm
         num_net = num_net + 1
         write(*,'(a,1x,I3)') 'GAMIT Network Tie Network',num_net
         do i = 1, num_net-1
            do j = 1, nom_tie_site
               do k = 1, num_site_ftp
                  if( network(1,k).eq.i .and.
     .             network(2,k).eq.tie_sel(j) ) then
                      network(3,k) = num_net
                  endif
*                 For the last network, add the last site.
                  if( network(1,k).eq.num_net-1 .and.
     .             i.eq.num_net-1.and.j.eq.2 .and.
     .             network(2,k).eq.
     .             (num_site_ftp - (num_net-2)*num_net_site)) then
                      network(3,k) = num_net
                  end if
               end do
             end do
         end do
      endif

*     Output list.
 
      nx = 0
      write(*,200) 
 200  format('# NETWORK TIE SITES')
      do j = 1, num_site_ftp
         if( network(3,j).eq.num_net ) then
            nx = nx + 1
            write(*,220) num_net, nx, flon(j), flat(j), fcod(j),
     .                     network(2,j),network(1,j), rw_stats(j)
 220          format('# ',i3.3,1x,i3,2F13.5,1x,a4,i3,1x,
     .               ' NetPrime ',i3,' RWR ',F7.3,' (mm/yr)^2')  
         end if
      end do

*     OK, Now add some additional ties between networks
      if( globk_US ) then
         write(*,'(a)') 'GLOBK ties no links between networks'
      else
         write(*,300)
 300     format('# Linkage sites between networks')
         do i = 1, num_net-2
            do j = 1, num_site_ftp
*              Look for 2nd site in next network and add this to this
*              network (provides 1 site overlap for each network to 
*              the next next network)
               if( network(1,j).eq.i+1 .and.
     .             network(2,j).eq. 2 ) then
                   network(4,j) = i
                   write(*,230) i, 0, flon(j), flat(j), fcod(j),
     .                        network(4,j),rw_stats(j)
 230               format('# ',i3.3,1x,i3,2F13.5,1x,a4,i3,1x,F7.3)  
               endif
             enddo
         end do

*        Now do the num_net-1 network (nominally the short network)
         nn = min(num_net_site - (num_site_ftp - 
     .           (num_net-2)*num_net_site),
     .            num_net-1)

         i = num_net - 1 
         do j = 1, num_site_ftp
*            Look for 2nd site in next network and add this to this network
*            (provides 1 site overlap for each network to the next next network)
             do k = 1,nn
                 if( network(1,j).eq. k .and.
     .               network(2,j).eq. 3 ) then
                     network(4,j) = i
                     write(*,230) i, k, flon(j), flat(j), fcod(j),
     .                    network(4,j), rw_stats(j)
                 endif
             end do
         enddo
      end if

      return
      end

CTITLE CHECK_STINF

      subroutine check_stinf

      implicit none 

*     Routine to check is site is in station.info

      include 'netsel.h'

      integer*4 i   ! Loop counter
     .,         ierr  ! IOSTAT error
     .,         indx  ! Pointer in string
  
      character*128 line
      character*4 stinf_site   ! Site name from station.info

***** Open station.info
      do i = 1, num_site_ftp
          stinfav(i) = .false.
      end do
      open(51,file=stinffile, status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',stinffile,0,'CHECK_STINF')
      do while ( ierr.eq.0 )
         read(51,'(a)',iostat=ierr) line
         if( line(1:1).eq.' ' .and. ierr.eq.0 ) then
             indx = 1
             call GetWord(line, stinf_site, indx)
             call casefold(stinf_site)
             do i = 1,num_site_ftp
                if( stinf_site(1:4).eq.fcod(i)(1:4) ) then
                    stinfav(i) = .true.
                end if
             end do
          end if
       end do         
       close(51)
       return
       end

CTITLE OUT_NETS
 
      subroutine out_nets

      implicit none 

*     Routine to output the networks

       include 'netsel.h'

       integer*4 i   ! Loop counter
       character*8 name
       character*8 xstinf

*      Loop over the sites
       if( globk_US ) then   ! Output use_site list
           write(*,100)
 100       format(' use_site clear')
           do i = 1, num_site_ftp
              write(*,120) trim(net_code), network(1,i), fcod(i)
 120          format(a,I3.3,' use_site ',a)
              if( network(3,i).gt.0 ) then
                 write(*,120) trim(net_code), network(3,i), fcod(i)
              end if
              if( network(4,i).gt.0 ) then
                 write(*,120) trim(net_code), network(4,i), fcod(i)
              end if
          enddo
           
       else                  ! sites.defaults outout
          do i = 1, num_site_ftp
             if( stinfav(i) ) then
                 xstinf = 'xstinfo'
             else
                 xstinf = ' '
             endif
             name = fcod(i)
             call caseunfold(name)
             write(*,200) name, net_code, network(1,i),xstinf
 200         format(1x,a4,'_gps   ',a2,i2.2,2x,a8)
             if( network(3,i).gt.0 ) then
                 write(*,200) name, net_code,network(3,i),xstinf
             end if
             if( network(4,i).gt.0 ) then
                 write(*,200) name, net_code,network(4,i),xstinf
             end if

          end do
       endif 

       return
       end

CTITLE READ_NET_RWF

      subroutine read_net_rwf

      implicit none 

*     Routine to read sh_gen_stats random walk file and exclude
*     sites whose horizontal random walk values are too large.

      include '../includes/const_param.h'
      include 'netsel.h'

      integer*4 ierr ! IOSTAT reading file
     .,         jerr ! IOSTAT decoding lines
     .,         ns   ! counter for number of sites
     .,         trimlen ! Length of string function
     .,         indx ! Position in string
     .,         j    ! Loop counter
     .,         lenw    ! Length of string

      real*8 rw_vals(2)   ! NE RW values in m^2/yr.
      real*8 dur          ! Duration (yrs)

      logical found  ! Set true when site found

      character*128 word  ! Word read from line
     .,             line  ! Line read from file
      character*8 cword   ! Casefolded version of site name
     .,           cd      ! Dummy
*
****  Open the rwfile
      open(50,file=rwfile,status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',rwfile,1,
     .                  'READ_NET_RWF')

      ns = 0

      minlon = 370
      maxlon = -10
      minlat = 91
      maxlat = -91

****  Loop over the file
      do while ( ierr.eq.0 ) 
         read(50,'(a)',iostat=ierr) line
         if(  trimlen(line).gt.0 .and. line(1:1).eq. ' '
     .                                   .and.ierr.eq.0 ) then
*            Decode the lines
             indx = 0
* Line looks like:
* mar_neu AOPR   1.79e-7   1.33e-7  34.94e-7 0.0 0.0 0.0 !  0.212 HRW   5.74 yrs
             call getword(line, word, indx)     ! globk command
             call casefold(word)
             if( word(1:7).eq.'MAR_NEU' ) then  ! Process noise line so decode 
                                                ! the rest
                call getword(line, cword, indx)    ! Should be site name
                call casefold(cword)
                call multiread(line,indx,'R8',jerr, rw_vals, cd, 2)
                read(line(70:),*,iostat=jerr) dur
                if( jerr.eq.0 .and. 
     .             (rw_vals(1)+rw_vals(2))*1.d6.le.maxuse_rw ) then
*                   OK to proceed.
                 
*****              See if we can file this site.
                   found = .false.
                   do j = 1, num_site_vel
                       if( cword.eq.snam(j)(1:4) ) then
                          found = .true.
*                        Extract the site 4-char code
                         ns = ns + 1
                         fcod(ns) = cword(1:4)
                         ffile(ns) = word
                         flon(ns) = slon(j)
                         flat(ns) = slat(j)
                         rw_stats(ns) = (rw_vals(1)+rw_vals(2))*1.d6/dur
                         exit    ! From do-loop now that name is found
                      end if
                   end do

*                  If site not found; get the coordinates from the rinex
*                  file.
                   if( .not.found ) then
                      write(*,220) cword(1:trimlen(cword))
 220                  format('* Coordinates not found for site ',a)
                   else

*                     Check range
                      if( flon(ns) .gt. maxlon ) maxlon =  flon(ns)
                      if( flon(ns) .lt. minlon ) minlon =  flon(ns)
                      if( flat(ns) .gt. maxlat ) maxlat =  flat(ns)
                      if( flat(ns) .lt. minlat ) minlat =  flat(ns)
                   endif
                else
                   write(*,240) cword, jerr,
     .                   (rw_vals(1)+rw_vals(2))*1.d6, maxuse_rw, dur
 240               format('* Site ',a,' JERR ',i4,
     .                    ' Not used Horiz RW ',F8.2,' Limit ',F8.2,
     .                    ' mm^2/yr, Dur ',F6.2,' yrs') 
                endif
             endif 
         endif 
      end do

*     Now report what we found
      num_site_ftp = ns
      write(*,310) trim(rwfile), num_site_ftp
 310  format('NETSEL: ',a,' contains ',i5,' sites')
      write(*,330) minlon, maxlon, minlat, maxlat
 330  format('Site Range Long ',2f10.4,' Latitude ',2f10.4,' deg')
      close(50)

*     Output sites
      write(*,400) 
 400  format('+  Num    Long    Lat   Name   RW_STAT')
      do ns = 1, num_site_ftp
         write(*,410) ns, flon(ns), flat(ns), fcod(ns), rw_stats(ns)
 410     format('+ ',i4,1x,F8.3,1x,F8.3,3x,a,1x,F7.3)
      end do

      return
      end

CTITLE INCR_NORM

      subroutine incr_norm(xyzr, A, neq)

      implicit none

*     Routine to compute partials and increment normal equations

      real*8 xyzr(4), A(3,6), neq(6,6)

      A(1,:) = (/ 1.d0, 0.d0 , 0.d0, -xyzr(3),    0.d0, -xyzr(2) /)           
      A(2,:) = (/ 0.d0, 1.d0 , 0.d0,   0.d0,    xyzr(3), xyzr(1) /)           
      A(3,:) = (/ 0.d0, 0.d0 , 1.d0,  xyzr(1),  xyzr(2),    0.d0 /) 

*     Add to normal equation
      neq = neq + matmul(transpose(A),A)/xyzr(4) 

      return 
      end

