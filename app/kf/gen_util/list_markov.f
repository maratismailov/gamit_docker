      subroutine list_markov( iout, file, comopt, list_file, type )

      implicit none 
 
*     This routine will list the contents of a markov file
*     ignoring blank lines, and those which do not start with
*     a blank, *, c of !.
* MOD TAH 030922: Added type for GLOBK or GLORG comamnd file
*
*   iout        - Output unit numer
 
      integer*4 iout
 
*   file        - Name of the markov file.
*   comopt      - Option for commands to executed
 
      character*(*) file, comopt, list_file, type
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of used portion of string.
 
      integer*4 ierr, trimlen, jndx, jerr, len_comopt
 
*   line        - Line read from markov file
*   source_file - Name of source file 
*   newline     - With comopt entry removed

      character*256 line, source_file, newline 

*   command     - Name of command
      character*16 command

*   updated     - Set true if comopt found on line
      logical updated 
 
****  OPen the markov file.
 
      open(300, file=file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',file,0,'list_markov')
      if( ierr.ne.0 ) RETURN
 
*     Start listing the file
      write(iout,100) type, file(1:trimlen(file))
 100  format(/' Summary of ',a,' command file ',a)
      write(iout,150)
 150  format(79('-'))
      len_comopt = trimlen(comopt)
      if( len_comopt.gt.0 ) then
          write(iout, 160) comopt(1:len_comopt)
 160      format('COMOPT: Lines starting with ',a,' are used')
      end if
 
      do while ( ierr.eq.0 )
          read(300,'(a)', iostat=ierr) line
          newline = line
          call decode_comopt(newline, comopt, updated) 
          if( ierr.eq.0 .and. trimlen(line).gt.0 .and.
     .        (line(1:1).eq.' ' .or. line(1:1).eq.'*' .or.
     .         line(1:1).eq.'c' .or. line(1:1).eq.'C' .or.
     .         line(1:1).eq.'!' .or. line(1:1).eq.'#' .or.
     .         updated) )   then
              write(iout,'(a)') line(1:trimlen(line)) 

*             Check to see if the source command has been given
              if( updated ) then
                  line = newline
              endif
                  
              if( line(1:1).eq.' ' ) then
                 jndx = 1
                 call GetWord(line,command,jndx)
                 call casefold(command)
                 if( command(1:3).eq.'SOU' ) then
*                    Source command: Write out file contents
                     call GetWord(line,source_file,jndx)
                     call wild_card(source_file, list_file)
		     open(301,file=source_file, status='old', 
     .                        iostat=jerr)
                     do while ( jerr.eq.0 )
                        read(301,'(a)', iostat=jerr) line
                        write(iout,'(a,a)') 'SOURCE >',
     .                        line(1:max(trimlen(line),1))
                     end do
                     close(301)
                 end if
              end if
                      
          end if
      end do
 
      write(iout,150)
 
****  Thats all
      close(300)
      return
      end
 
