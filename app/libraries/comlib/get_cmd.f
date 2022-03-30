CTITLE 'GET_CMD'
 
      subroutine get_cmd( line, commands, num_commands, iel, indx)
 
      implicit none

 
*     General routine to find commands.  The first group of non-blank
*     characters from 'line' is compared with the commands.  The
*     following  returns for 'iel' are produced:
*     iel > 0  -- command number 'iel' found
*     iel = -1 -- command not found
*     iel = -2 -- ambiguous command (ie. multiple commands match
*                 input
*     iel = -3 -- Command string empty
*     iel =999999 -- Command was 'ALL'
*
*     All commands must be in upper case
*     Commands can be any length less than 40 characters
* MOD TAH 941128: Added feature that if the last characters of 
*     a command are _ then blanks are substituted by the length
*     of the string is kept as if the _ were still their.
* MOD TAH 950610: Undid feature of using _ as blank indicator.
*     Now if duplicate found, then whole command length string
*     if tested to see if uniques.
* MOD TAH 950610: Added feature of passing negative numbers of
*     commands to force look up on full command name rather then
*     minimum, redundant string (This is for site names where the
*     names of some sites may be subsets of other sites).
 
*         iel       - Command number or error
*   indx            - Position in line to start searching.
*                   - Returns with number of next group of
*                   - characters. Refers to beginning of line
*   i               - Loop counter
*   len_commands    - Length of commands
*   len_word        - Length of first word from line
*   min_len         - Smallest length words
*   num_commands    - number of commands in list
*   trimlen         - Length of string
*   nc              - Actual number of commands =abs(num_commands)
 
      integer*4 iel, indx, i, len_commands, len_word, min_len,
     .    num_commands, trimlen, nc

*   commands(1)   - list of commands
*   line            - Line from user.
 
      character*(*) commands(*), line
 
*   command_word   - The command word which will be
*                  - checked
*   commnd_UC      - Upper case version of commands (so that we do not
*                    need to casefold things like site names).
 
      character*40 command_word, command_UC
 
***** Get the length of the commands
      len_commands = LEN(commands(1))
 
*     Get the command word.  Use read_line to keep indx correct in
*     string.
      call GetWord(line, command_word,indx)
 
*     Get length of command
      len_word = trimlen(command_word)
 
*     See if empty string
      if( len_word.eq.0 ) then
          iel = -3
          RETURN
      end if
 
      call casefold( command_word )
      min_len = min( len_word, len_commands )
*
* MOD TAH 950610: See if we should use full length of commands
      nc = abs(num_commands)
      if( num_commands.lt.0 ) then
          min_len = len_commands
      end if

***** Now see if we can find command
      iel = -1
      do i = 1, nc
 
*         See if command matches
          command_UC = commands(i)
          call casefold(command_UC)

          if( command_UC(1:min_len).eq.command_word(1:min_len) ) then
 
*             See if we have found before
              if( iel.eq. -1) then
                  iel = i
              else
                  iel = -2
              end if
          end if 
      end do

* MOD TAH 950610: Fixing problems with commands which
*     are subsets of others.
***** See if duplicate command found.  Check to see if
*     unique if whole word including blanks is considered.
      if( iel.eq.-2 .and. num_commands.gt.0 ) then

          min_len = len_commands
 
*****     Now see if we can find command
          iel = -1
          do i = 1, nc 
 
*             See if command matches
              command_UC = commands(i)
              call casefold(command_UC)
              if( command_UC(1:min_len).eq.
     .            command_word(1:min_len) ) then
 
*                 See if we have found before
                  if( iel.eq. -1) then
                      iel = i
                  else
                      iel = -2
                  end if

              end if 
          end do
      end if
 
***** Now finally check to see if command was all
      if( command_word(1:4).eq.'ALL ' ) then
          iel = 999999
      end if
 
***** Thats all
      return
      end
 
