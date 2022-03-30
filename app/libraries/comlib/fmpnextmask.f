 
      logical function fmpnextmask( dcb, ierr, curpath, nextname)

      implicit none
 
*     Routine to return the name of the next file in the directory list
*   if EOF is returned, ierr is set to zero and fmpnextmask is set false.
 
 
*   dcb(*)  - DCB buffer (only unit number used)
*   ierr    - IOSTAT error on reading file.
 
      integer*4 dcb(*), ierr
 
*   curpath     - Current path name (not used)
*   nextname    - Name of next file
 
      character*(*) curpath, nextname
 
* LOCAL VARIABLES
 
*   num         - Number of files (only for directories)
*   size        - Size of the file
*   indx        - Index for working through buffer
*   id          - Dummy
*   jerr        - Error in decodeing line
*   trimlen     - length of string
*   dirend      - Len of current path 
 
      integer*4 num, size, indx, id, jerr, trimlen, dirend
 
*   buffer      - Buffer for reading directory entry
*   Name        - Name of the KalObs files
*   dirname     - Name of directory
 
      character*256 buffer
      character*128 name, dirname
 
*   prot        - Protection mode
 
      character*10 prot
 
*   owner       - Oowners name
 
      character*12 owner
 
*   Mon         - Month
*   Day         - Day
*   cd          - Dummy
 
      character*4 Mon, Day, cd
 
*   time        - TIme
 
      character*5 time
 
****  Try to read the next record
      jerr = -1
      buffer = ' '
      fmpnextmask = .false.     
      do while ( jerr.eq.-1 ) 
         read(dcb(1),'(a)', iostat=ierr ) buffer
         if( ierr.eq.-1 ) then
             fmpnextmask = .false.
             ierr = 0
             jerr = 0
         else
 
*            Decode line if the error OK
             if( ierr.eq.0 .and. (buffer(1:1).eq. '-' .or.
     .                            buffer(1:1).eq. 'l') .and.
     .                            buffer(4:4).eq. '-'    ) then
                 indx = 1
                 call read_line(buffer,indx,'CH',jerr, id, prot)
                 call read_line(buffer,indx,'I4',jerr, num, cd)
                 call read_line(buffer,indx,'CH',jerr, id, owner)
*                Added group return
                 call read_line(buffer,indx,'CH',jerr, id, owner)
                 call read_line(buffer,indx,'I4',jerr, size, cd)
                 call read_line(buffer,indx,'CH',jerr, id, Mon)
                 call read_line(buffer,indx,'CH',jerr, id, Day)
                 call read_line(buffer,indx,'CH',jerr, id, time)
                 call read_line(buffer,indx,'CH',jerr, id, name)
             else
*                see if directory name
                 indx = 1
                 call getword(buffer, dirname, indx )
                 dirend = index(dirname,':')
                 if( dirend.gt.1 ) then
*                    found a directory: Set current parth
                     curpath = dirname(1:dirend-1) // '/'
                 end if
             end if

*            Now prepend path to name
             if( trimlen(curpath).gt.0 ) then
                 nextname = curpath(:trimlen(curpath)) // name
             else
                 nextname = name
             end if

             fmpnextmask = .true.
             ierr = -jerr
         end if
      end do
 
****  Thats all
      return
      end
 
 
 
 
