ctitle
 
      subroutine wild_card(file,key)

      implicit none
c
c     This routine will subsitute the part of data file name (key)
c     which matches with the '@' entry.
c     EXAMPLE:
c     if FILE contains BA@::50, and the datafile name is OB6201::49
c     then FILE will be returned as BA6201::50.
c
c     Restriction: '@' can not be the first character of FILE.
c
c     Added ? as the preferred character replacement.  -- at least
c     now needed  MOD TAH 150828.
c     The OLD use of '-----' for FMGR files is still supported
c     This feature has been expanded to allow any number of -'s in
c     name (each one is replaced by the corresponding character from the
c     key)
*
* MOD TAH 900105: Key substution made only within root part of name.  The
*     directory information at the front in used to define the start of file
c
c Variables
c ---------
c file -- the file name.  This name will be changed if ---- appears in it
c key  -- the key to be substituted for ----
c
 
      character*(*) file, key
 
*   at_pos      - Position of '@' in FILE string 
*   end_file    - Position in FILE of : or . If these do not
*               - apear eg File = 'L@' then the whole of key
*               - string after the @ will be copied.
*   end_root    - End of the root part of the KEY file name.
*               - i.e., last character before '.', ':', or
*               - of string (searched in that priority order)
*   i,j         - Loop counters.
*   key_len     - Length of the key string (non-blank portion)
*   trimlen     - HP utility for length of a string

*   start_file  - Start of file name (character after /)
*   start_root  - Start of root part of key
*   lastex      - Function to return last occurence of sub string
 
      integer*4 at_pos, end_file, end_root, i, j, key_len, trimlen

      integer*4 start_file, start_root, lastex
      integer*4 nd   ! Number of - characters to replace
 
*   new_File    - string used to construct new file name to
*               - replace FILE.
 
 
      character*2048 new_File
 
*     CHECK for old use first.
c.... see if - appears.  Note: the file string should have no leading
c     blanks.  All -'s will be replaced with corresponding string
      i = 0
      key_len = trimlen(key)
      do while (i.lt.key_len)
          i = i + 1
* MOD TAH 150829: Change wild card to ? to make consistent with UNIX and
*         require - to be at least two (avoid conflict with -PER option)
*                                         ! Replace with letter from key
          if( file(i:i).eq.'?' ) then
              file(i:i) = key(i:i)
          elseif( file(i:i+1).eq.'--' ) then
*             See how many to replace
              nd = 0 
              do j = i+1, trimlen(key)
                 if( file(j:j).ne.'-' ) exit
                 nd = nd + 1
              enddo
              file(i:i+nd) = key(i:i+nd)
          end if
      end do
c
*     Now see if there is a '*' and then '@' in the FILE name.  
*     Note position of '*' or '@'
 
      at_pos = index(file, '*' )
      if( at_pos.eq.0 ) then
         at_pos = index(file,'@')
      end if
 
*                              ! FILE does contain a wild character
      if( at_pos.gt.0 ) then
 
*         Get the position of : or . in file mask
          end_file = lastex( file, '.')
          if( end_file.eq.0 ) then
              end_file = index ( file,':')
*                     ! We will check below if end file is still zero,
          end if
*                     ! if it is we will modify end_root
 
*         Now get the end of the root of the data file name.
 
*                                         ! See if extent
          end_root = lastex( key, '.') - 1
*                                         ! No root, see if :
          
          if( end_root.le.0 ) then
              end_root = index( key, ':') - 1
*                                      ! No cartridge so use full length
              if( end_root.le.0 ) then
                  end_root = trimlen( key )
              end if
          end if
 
*         See if we should use all of the key string
*                                     ! There is no : or . in FILE so use
          if( end_file.eq.0 ) then
*                                     ! all of KEY
              end_root = trimlen( key )
          end if

****      Now get the starts of file and key parts
          start_file = lastex( file, '/')
          start_root = lastex(  key, '/') 
 
*         Now copy from at_pos+1 to end_root into '@' position in
*         FILE.  Make new file name in a copy.
 
*         Treat at_pos = 1 separately since we do not need the
*         first string string
          if( at_pos.gt.1 ) then
              New_File = File(1:at_pos-1) // 
     .                   key(start_root+(at_pos-start_file):end_root) //
     .                   File(at_pos+1:)
          else
* MOD TAH 030119: Changed code when wild-card at front
* MOD TAH 141024: Changed behavior when ext in key is different length
*             ext in file names.   When @ leads behavior should be the
*             same as when directory or lead characters given.
!             at_pos = trimlen(file)-1
              if( at_pos.eq.0 ) then
                  new_file = key
              else
!                  new_File = key(1:trimlen(key)-at_pos) //
!     .                       file(2:at_pos+1)
                  new_File = 
     .                   key(start_root+(at_pos-start_file):end_root) //
     .                   File(at_pos+1:)

              end if
          end if
          File     = New_file
 
      end if

* MODTAH 130131: Added $HOME (~) subsitution in file name
      call subhome( file )
 
c.... Thats all
      return
      end
 
