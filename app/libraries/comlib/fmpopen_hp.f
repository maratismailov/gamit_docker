Ctitle FmpOpen
 
      subroutine FmpOpen(dcb, ierr, name, opts, buffs)

      implicit none
 
*     Routine to emulate the HP 1000 FmpOpen subroutine.
*     The DCB is a control buffer used by HP, we here say
*     our own values in it.
* MOD TAH 900810: If dcb(1) contains a valid unit number,
*     then this unit number will be re-used.  Avoids problem
*     with too many open files if files are not closed but
*     use the same dcb buffer.
*
*     DCB(#) - Meaning
*     1  - Fortran logical init number allocated to the
*          file.  We start these at 700 and progress upwardss
*     2  - record length of this file.  Extracted from name
*          of file. Value of Zero if not a direct access file
*     3  - record number of current position in file.
*     4  - 0 if file read/write, 1 if read only.
*     5  - file type (1-128 fixed, 2-<>128fixed, 3-variable)
*     6-16 - not used
*
*     IERR - error return.  values are:
*     0  - OK
*    -6  - File not found if 'O' open used (without create)
*
*     NAME - is the name of the file to be opended with
*          optionally the following information appended
*          Name:type:size:recl
*     where type is either 1 or 2 for direct access files
*            (type 1 has recl=128 I*4 words)
*           size is length of the file (ignored in UNIX)
*           recl - record length for file.
*
*     OPTS - is a character string with options.
*        'R' is open for read,
*        'C' is open for write (read only if not specified)
*        'O' is open file
*        'C' is create file.  (If O and C are given, open will
*            be tried first; if 'C' only is given, then file
*            must not already exist.
*
*     BUFFS - is added for compatability but is not used.
*
*         dcb(16)       - Dcb buffer with information about the file
*   ierr                - Open error as given above
*   buffs               - Fake not used.
 
      integer*4 dcb(16), ierr, buffs
 
*             name      - Name of the file with type and size appended
*   opts                - Sting with options.
 
      character*(*) name, opts
 
*     LOCAL VARIABLES
 
*   next_unit           - Next available unit number
*   unit                - Actual unit number used to open file.  May
*                         be value contained in dcb(1) if this value
*                         is greater then 700
*   recl                - Record length for the file
*   filetype            - Type of file.  (1 for 128 fixed record,
*                       - 2 for fixed not 128, and >2 for all others)
*   indx, jndx          - position of coluns in name
*   jerr                - Local error messages
*   i                   - loop counter
*   trimlen             - length of string
*   iun                 - Unit number of file named to be opened (normally
*                         would not be open)
 
      integer*4 next_unit, recl, filetype, indx, jndx, jerr, i, trimlen,
     .          unit, iun
 
*       readopen        - Open for reading
*   writeopen           - OPen for writing
*   openopen            - open file
*   createopen          - Create open
*   openq               - Returns true if file name already opened.
 
      logical readopen, writeopen, openopen, createopen, openq
 
*             forname   - Fortran file name
 
      character*132 forname
 
*            status     - Status for fortran open
*   access              - Fortran acess type
 
      character*20 status, access
 
      data next_unit / 699 /

* rwk 070416  Dummy statement to avoid a compiler warning
*      buffs = 0                                          
* rwk 070420: Removed since this causes a segmentation violation.

****  Start, get the parameters form the name
 
      indx = index(name,':')

*     If index is 0 then assume name as given and that file is type 1
      if( indx.gt.1 ) then
          forname = name(:indx-1)
      else
          forname = name(1:max(1,trimlen(name)))
      end if
 
*     Get the next part of the name.  Ignore security code and cartridge number
      do i = 1, 2
          jndx = index(name(indx+1:),':') + indx
          indx = jndx
      end do
 
*     Get type
      jndx = index(name(indx+1:),':') + indx
      filetype = 1
*                                 ! Type is present, so get
      if( jndx.gt.indx ) then
          read(name(indx+1:jndx-1),*,iostat=jerr) filetype
          if( jerr.ne.0 ) then
              ierr = -9
              return
          end if
      end if
 
*     Skip the size of the file
      if( filetype.eq.1 ) then
          recl = 128
      else
          recl = 0
      end if
 
      if( filetype.eq.2 ) then
          indx = jndx
          jndx = index(name(indx+1:),':') + indx
          if( jndx.gt.indx ) then
              read(name(indx+1:jndx-1),*, iostat=jerr) recl
              if( jerr.ne.0 ) then
                  ierr = -10
                  return
              end if
          end if
      end if
 
****  Now start decoding the open option
      readopen = .false.
      writeopen = .false.
      openopen  = .false.
      createopen = .false.
 
      indx = index(opts,'r') + index(opts,'R')
      if( indx.ne.0 ) readopen = .true.
      indx = index(opts,'w') + index(opts,'W')
      if( indx.ne.0 ) writeopen = .true.
      indx = index(opts,'o') + index(opts,'O')
      if( indx.ne.0 ) openopen = .true.
      indx = index(opts,'c') + index(opts,'C')
      if( indx.ne.0 ) createopen = .true.
 
****  Now see what we should do

      status = 'UNKNOWN'
 
      if ( openopen .and. createopen ) status = 'UNKNOWN'
      if ( .not. openopen .and. createopen ) status = 'NEW'
      if ( openopen .and. .not. createopen ) status = 'OLD'
 
      if ( filetype.eq.1 .or.filetype.eq.2 ) then
          access = 'DIRECT'
      else
          access = 'SEQUENTIAL'
      end if
 
****  Now try to do open
      if( dcb(1).ge.700 .and.dcb(1).lt.1000 ) then
          unit = dcb(1)
      else
          next_unit = next_unit + 1
          unit = next_unit
      end if

* MOD TAH 080127: Check status of file to be opened.
      inquire(file=forname, number=iun, opened=openq)
      if( openq ) then
          close(iun)
      endif
*

      close(unit)
      open(unit, file=forname, status=status, iostat=jerr,
     .    access = access, recl=recl*4, form='unformatted' )
 
      if( jerr.ne.0 ) then
          ierr = -jerr
*         HP/SunOs/Solaris 2.5 IOSTAT error for file does not exist.
          if( ierr.eq.-908 .or. ierr.eq.-118 .or. 
     .        ierr.eq. 1018 ) ierr = -6
          if( unit.ne.next_unit ) next_unit = next_unit - 1
          return
      end if
 
****  Save the information we need
      dcb(1) = unit
      dcb(2) = recl
      dcb(3) = 0
      if( readopen .and. .not.writeopen ) then
          dcb(4) = 1
      else
          dcb(4) = 0
      end if
      dcb(5) = filetype
 
      do i = 6,16
          dcb(i) = 0
      end do
 
****  Thats all
      ierr = 0
      return
      end
 
 
