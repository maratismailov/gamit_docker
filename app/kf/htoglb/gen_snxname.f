ctitle gen_snxname
 
      subroutine gen_snxname( epoch, estart, eend, owner, dir, name)

      implicit none 
 
*     This routine will create the name of the sinex file based on the user
*     inputs
 
* epoch  - JD of reference time from solution
* estart - JD of start
* eend   - JD of end (difference of start and end used to compute
*          last character)
* owner  - 3  character name of owner of file
* dir    - directory where files will be put
* name   - User supplied portion of the name.  This variable will be returned
*          with the full name of the file
 
      real*8 epoch, estart, eend
 
 
      character*(*) owner, dir, name
      character*1 last_char
 
* LOCAL VARIABLES
 
* trimlen  - Length of string
* gps_week, gps_dow - Gps WEEK number and day of week
* indxq - Index in name for #### (replaced with week and optionally day)

      integer*4 trimlen, gps_week, gps_dow, indxq
 
* gpss   - Julian date of week 0, day 0 for GPS
 
      real*8 gpss
  
* lcowner  - Lower case version of owner name
*
      character*4 lcowner 
      character*256 tmp
      
      data gpss / 2444244.50d0 /
      
***** Start with preliminary calculations
      lcowner = owner
      call caseunfold(lcowner)
      call sub_null(lcowner)
      
      if( lcowner.eq.'----' ) lcowner = 'unk'
      owner = lcowner
      call casefold(owner)

      gps_week = (epoch - gpss)/7.d0
      gps_dow  = epoch - gps_week*7 - gpss

* MOD TAH 030201: Reduced GPS dow by 1 day to make consistent with
*     IGS standards (ie., first day of week is day 0 and not day 1)
*     Added check to see if multiday in which case 7 is used
      if( eend - estart.gt.1.2d0 ) gps_dow = 7
 
***** See what user has given us:
      last_char = name

      if ( trimlen(name).eq.0 ) then
*         Need to fully generate the name
          write(name, 120) lcowner, gps_week, gps_dow
 120      format(a3,i4.4,i1,'.snx')
      else if( trimlen(name).eq.1 ) then
          write(name, 130) lcowner, gps_week, last_char
 130      format(a3,i4.4,a1,'.snx')
      else
 
****      See if #### are in name
          indxq = index(name,'#####' )
          if( indxq.gt.0 ) then
              write(name(indxq:indxq+4), 140) gps_week, gps_dow
 140          format(i4.4,i1)
          else
*             See if just week number needed
              indxq = index(name,'####')
              if( indxq.gt.0 ) then
                  write(name(indxq:indxq+3), 150) gps_week
 150              format(i4.4)
              end if
           
          end if
      end if

***** Now add in the directory part of the name
      if( trimlen(dir).gt.0 ) then
          tmp = dir(1:trimlen(dir)) // '/' // name
          name = tmp
      end if

****  Thats all
      return
      end
 
CTITLE SNX_SYS_TYPE
 
      subroutine snx_sys_type( csys_type, sys_char)

      implicit none 
 
*     Rouitne to convert binary encoded system type to sinex character
 
      integer*4 csys_type
 
 
      character*(*) sys_char
 
*     Check each system type (defined in sln_def.h) if it
*     not a singele system then allocate C
      if ( csys_type.eq.1 ) then
          sys_char = 'R'
      else if ( csys_type.eq.2 .or. csys_type.eq.0 ) then
          sys_char = 'P'
      else if ( csys_type.eq.4) then
          sys_char = 'L'
      else if ( csys_type.eq.8) then
          sys_char = 'M'
      else
          sys_char = 'C'
      end if
 
*     Thats all
      return
      end
 
ctitle cp_comments
 
      subroutine cp_comments( unit, unitc, entry)

      implicit none 
 
*     This rouitne copies the comments from the comments file to the
*     output sinex file

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
*   unit        - SInex file unit number
*   unitc       - Comments file unit number. If 0 then there is no
*               - such file and this routine does nothin
 
      integer*4 unit, unitc
 
 
      character*(*) entry
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
 
      integer*4 ierr, trimlen
 
*   start_found - Set true when start line found
*   stop_found  - Set true when stop line found
 
      logical start_found, stop_found
 
*   start_entry - Entry to look for the start of the blovk
*   stop_entry  - Entry for stop
*   line        - Line read from comments
 
      character*128 start_entry, stop_entry, line

*   hardware    - Machine type from hosttype.

      character*16  hardware
 
****  Check to see if comment unit
      if( unitc.eq.0 ) RETURN
 
****  Construct the start and stop string
      start_entry = '+' // entry
      stop_entry  = '-' // entry
 
      rewind(unitc)
 
      line = ' '
      ierr = 0
      start_found = .false.
      stop_found = .false.
      do while ( ierr.eq.0 .and. .not.stop_found )
 
          read(unitc, '(a)', iostat=ierr) line
          if( index(line,stop_entry).eq.1 ) stop_found = .true.
          if( ierr.eq.0 .and. .not.stop_found .and. start_found .and.
     .        trimlen(line).gt.0 ) then

*             Check for certain strings to see if should update:
              if( index(line,'SOFTWARE').gt.0 .and.
     .            line(1:1).eq.' ') then
                  line(trimlen(line)+2:) = glbtosnx_version
              else if ( index(line,'HARDWARE').gt.0 .and. 
     .                  line(1:1).eq.' ' ) then
                  call getenv('HOSTTYPE', hardware)
                  call casefold(hardware)
                  line(20:) =  hardware
              end if
              write(unit,'(a)' ) line(1:trimlen(line))
          end if
*         Check for start afterwards so that the start line
*         is not echoed.
          if( index(line,start_entry).eq.1 ) start_found = .true.
      end do
 
***** Thats all
      return
      end

CTITLE qname_to_codeoccpt

      subroutine qname_to_codeoccpt

      implicit none 

*     Routine to generate a unique list of codes, occ and pt values from
*     the list of sites we have

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'

* LOCAL VARIABLES

      integer*4 i,j   ! Loop counter
      integer*4 cnt   ! Counter
      integer*4 jerr,kerr  ! Error flag

      character*4   ext  ! End of site name (last 4-characters)

      character*1 pt
      integer*4 occ

      do i = 1,cnum_sites

****     See if standard name
         code_all(i) = qsite_names(i)(1:4)
         ext = qsite_names(i)(5:8)
         call ext_to_po( ext, pt, occ, igs_ptname )
         pt_all(i) = pt
         occ_all(i) = occ

         if( .not. igs_ptname ) then
            if( ext(1:2).eq.'_V' .or.
     .          ext(1:2).eq.'_S' .or. ext(1:2).eq.'_D' ) then
*               See if occ code on pt name
                call check_num(ext(3:3),jerr)
                call check_num(ext(4:4),kerr)  ! Catch numeric EQ number
                if( jerr.eq.0 .and. kerr.ne.0 ) then ! Numeric so get Occ value
                    read(ext(3:3),*) occ_all(i)
                    pt_all(i) = ext(4:4)
                else
*                   Field is not numeric, so could be an earthquake
*                   rename, so set Pt = 'A' and get numeric occ code
                    pt_all(i) = 'A'
                    cnt = 2
                    do j = 1, i-1
                       if( code_all(j).eq.code_all(i) .and.
     .                     pt_all(j).eq.pt_all(i)     .and. 
     .                     occ_all(j).ne.1  ) then
                          cnt = cnt + 1
                       endif
                    end do
                    occ_all(i) = cnt
                endif
            end if
*               See if first three letters of ext is numeric
C               call check_num(ext(1:3),jerr)
C               if( jerr.eq.0 ) then
*                  Directly get occ and pt code
C                  read(ext(1:3),*) occ_all(i)
C                  pt_all(i) = ext(4:4)
C               else
*                  Most likely an eathquake/antenne rename
C                   pt_all(i) = 'A'
C                   cnt = 2
C                   do j = 1, i-1
C                      if( code_all(j).eq.code_all(i) .and.
C    .                     pt_all(j).eq.pt_all(i)     .and. 
C    .                     occ_all(j).ne.1  ) then
C                         cnt = cnt + 1
C                      endif
C                   end do
C                      occ_all(i) = cnt
C               endif
         endif   
      end do

* MOD TAH 140114: Now run final check to make sure all names
*     are unique
      do i = 2,cnum_sites
         do j = 1, i-1
*           See if unique
            if( code_all(j).eq.code_all(i) .and.
     .          pt_all(j).eq.pt_all(i)     .and. 
     .          occ_all(j).eq.occ_all(i)  ) then
*               Duplicate name: So change the occupation number
                occ_all(i) = occ_all(i)+1
            endif
         enddo
      enddo

      write(*,105)
 105  format('------------------------------------------')

      write(*,110) 
 110  format('* GLOBK to SNX Code, Pt and Occ conversion',/,
     .       '    # Name     Code Pt Occ')

      do i = 1,cnum_sites
         write(*,120) i, qsite_names(i), code_all(i), pt_all(i), 
     .                   occ_all(i)
 120     format(i5,1x,a8,1x,a4,2x,a1,1x,I3.3)
      end do
      write(*,105) 
****  Thats all
      return
      end

 
CTITLE name_to_id
 
      subroutine name_to_id( names, ns, code, pt, occ )

      implicit none 

*     Rouitne to get to the code, pt and occupation
*     number

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
*   names(ns)   - Names of all sites so we can test uniqueness
*   code        - First 4 characters of name
*   pt          - Last character of name
*   ns          - Site number
 
      integer*4 ns
      character*(*) names(ns), code, pt
 
*   occ     - Occupation number, normally 1, but could be the
*           - 7th letter
*   ierr    - IOSTAT error while reading 7 th character
 
      integer*4 occ

****  Just use the saved values
      code = code_all(ns)
      occ  = occ_all(ns)
      pt   = pt_all(ns)
 
***** That all
      return
      end

CTITLE ext_to_po

      subroutine ext_to_po( ext, pt, occ, igs_ptname )

      implicit none 

*     Routine to convert last three chacters to pt and occ line.
*     occ should be name(6:8)


      character*(*) ext, pt
      integer*4 occ

      logical igs_ptname   ! Set if IGS sinex site names to be used

      integer*4 ierr ! Used to see of point is number.  If it is
                     ! then converted.

      character*1 cocc  ! Character verion of OCC

****  See if standard
* MOD TAH 100511: -s option which converts to IGS sinex names
      if( .not. igs_ptname ) then  ! original code
         if( ext(2:4).eq.'GPS' .or. ext(2:4).eq.'VLB' .or.
     .       ext(2:4).eq.'SLR' .or. ext(2:4).eq.'TER'  ) then
             pt = 'A'
             occ = 1
             RETURN
         endif

****     More complicated.  See what is outthere
         pt = ext(2:2)
         if( pt(2:2).eq.'G' ) pt = 'A'
         call check_num(pt, ierr)
         if( ierr.eq.0 ) then  ! pt is numeric.
             read(pt,*,iostat=ierr) occ
             pt = 'B'
             RETURN
         endif

****     See of end is numeric and if so return value
         call check_num(ext(2:3), ierr)
         if( ierr.eq.0 ) then
             read(ext(2:3),*,iostat=ierr) occ
             RETURN
         endif

****     Set the occupation to nominal 1
         occ = 2
      else
*        Translate name back with <code>_<O><P>S (can be problem with
*        earthquakes
         cocc = ext(2:2)
         call check_num(cocc, ierr)
         if( ierr.eq.0 ) then
             read(cocc,*) occ
         elseif( cocc(1:1).eq.'G' ) then
             occ = 1
         elseif( cocc(1:1).le.'Z' ) then
             occ = ichar(cocc) - 64 + 9
         else 
             occ = 0
         end if
         pt = 'A'
         if( ext(4:4).eq.'S' .and. ext(3:3).ne.'P' ) then  ! ie.,  _BPS etc
             pt = ext(2:2)
         end if
      end if

****  Thats all
      return
      end 



      
 
CTITLE GET_NO
 
      subroutine get_no( ni, itoo, no)

      implicit none 
 
*     Routine to get the output parameter number given the input
*     parameter number and the list of outputs
 
*   ni      - Parameter number in the input solution
*   no      - Parameter number on the output SINEX file
*   itoo(*) - Mapping of the input to output parameter number
 
      integer*4 ni, no, itoo(*)
 
*   i       - Loop counter
 
      integer*4 i
 
      no = 0
      if( ni.eq.0 ) RETURN
 
      do i = 1, ni
          if( itoo(i).eq.ni ) no = i
      end do
 
****  Thats all
      return
      end
 

CTITLE UPD_DOMES

      subroutine upd_domes( unitc )

      implicit none 
      
*     Routine to update the domes information when a sinex
*     file is written.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'

 
* unitc  -- Comments file unit numbers

      integer*4 unitc
      
* LOCAL VARIABLES      
*   ierr        - IOSTAT error
*   trimlen     - Length of string
 
      integer*4 ierr, trimlen, ns

 
*   start_found - Set true when start line found
*   stop_found  - Set true when stop line found
 
      logical start_found, stop_found
 
*   line        - Line read from comments
 
      character*128  line
            
*   domes       - DOMES number for site (9 characters needed)
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]
*   globk_extent - Last four characters of globk name

      character*12 domes
      character*4 type, pt, globk_extent, code
      character*8 full_name
* MOD TAH 050927: Added long name to be read from head.snx
      character*22 long_name
      
****  Check to see if comment unit
      if( unitc.eq.0 ) RETURN
      rewind(unitc)

***** See if +DOMES can be found

      line = ' '
      ierr = 0
      start_found = .false.
      stop_found = .false.
      do while ( ierr.eq.0 .and. .not.stop_found )
 
          read(unitc, '(a)', iostat=ierr) line
          if( index(line,'-DOMES').eq.1 ) stop_found = .true.
          if( ierr.eq.0 .and. .not.stop_found .and. start_found .and.
     .        trimlen(line).gt.0 .and. line(1:1).eq.' ' ) then
              
              read(line, 120) code, pt, domes, type, globk_extent,
     .                        long_name
 120          format(1x,a4,1x,a2,1x,a9,1x,a1,1x,a4,1x,a16)
              full_name = code(1:4) // globk_extent
c              indx = 1
c              call get_cmd(full_name, qsite_names, cnum_sites,
c     .               ns, indx)
*             Match the 4-char code
              do ns = 1, cnum_sites
                  if( code(1:4).eq.qsite_names(ns)(1:4) ) then
                     if( trimlen(long_name).gt.0 ) then
                         qfull_names(ns)(1:22) = long_name
                     endif
                     qfull_names(ns)(23:31) = domes(1:9)
                     qfull_names(ns)(32:32) = type
                  end if
              end do
          end if
*         Check for start afterwards so that the start line
*         is not echoed.
          if( index(line,'+DOMES').eq.1 ) start_found = .true.
      end do

****  Thats all
      return 
      end
