CTITLE    ..................................................................
 
      subroutine read_label(gbuffer, label)
c
c     Routine to read a label from the input string.  The label
c     should be delimited with " (double quotes).
c
c
c Include files
c -------------
      include 'plot_param.h'
c
      include 'plot_com.h'
 
c Variables
c ---------
c gbuffer -- the user input buffer
c label  -- the returned label from the buffer
c
      character*(*) gbuffer, label
 
c
c Local variables
c ---------------
c start  -- the character number of start of label
c fini   -- the number of charcaters in label
c new_label -- the new label character string
c head_start -- start character number of ':h' or ':H'
c poly_start -- Start character number of ':p' or ':P'
c header_num -- the header record number to be used.
c error_message -- charater string for error meassge
* file_start   -- :f for file name to be label
* time_start   -- :t for time label
c
 
      integer*4 start, fini, head_start, header_num, poly_start,
     .          file_start, time_start, date(5)

      real*8 sectag
 
c
      character*80 new_label, error_message
 
c
c Scratch common
c --------------
      common new_label, error_message
 
c
c.... Start by nulling the label
      label = ' '
c
c.... Start by seeing is :hnn or :Hnn has been given
      head_start = index(gbuffer,':h')
*                                ! see if :H
      if( head_start.eq.0 ) then
          head_start = index(gbuffer,':H')
      end if
c
***** See if polynomial label to be given
      poly_start = index( gbuffer,':p')
*                                 ! see if :P
      if( poly_start.eq.0 ) then
          poly_start = index(gbuffer,':P')
      end if
 
      if( poly_start.ne.0 ) head_start = poly_start

*     See if file name
      file_start = index(gbuffer,':f')
      if( file_start.eq.0 ) then
          file_start = index(gbuffer,':F')
      end if

*     See if time label
      time_start = index(gbuffer, ':t')
      if( time_start.eq.0 ) then
          time_start = index(gbuffer, ':T')
      end if

 
c.... If :h or :H given then find header number else scan for "string"
      if( head_start.gt.0 .or. file_start.gt.0 .or. 
     .    time_start.gt.0 ) then
*                                                              ! we use
          if( head_start.gt.0 )
     .    call get_int(gbuffer(head_start-6:),header_num, 1)
*                                          ! the -6 because get_int assumes
*                                          ! a command (8 characters long)
*                                          ! to be before integer.
c
c....     Get the header record if valid number
*                                         ! Get from header records
          if( poly_start.eq.0 .and. file_start.eq.0 .and.
     .        time_start.eq.0 ) then
              if( header_num.gt.0 .and.
     .            header_num.le.actual_num_headers ) then
 
                  label = headers( header_num )
*                                              ! illegal header number
              else
 
                  error_message = ' HEADER NUMBER too large or small'
     .                // ' in ' // gbuffer(head_start:)
*                                  ! flush buffers before error message
                  call jmcur
                  call report(error_message)
              end if
          else if( poly_start.gt.0 ) then
              if( header_num.gt.0 .and.
     .            header_num.le.max_poly_label ) then
                  label = poly_labels( header_num )
              else
                  error_message = ' POLY Label number too large or'
     .                // ' small in ' // gbuffer(poly_start:)
              end if
*                 ! Header or polynomial
          end if

*         Get file name if appropriate
          if( file_start.gt.0 ) then
              label = 'Data file ' // input_file
          end if

*         Get the time if appropriate
          if( time_start.gt.0 ) then
              call systime( date, sectag )
              write(label,150) date
 150          format('Date ',i4,'/',i2,'/',i2,1x,i2,':',i2)
          end if

*                 ! decode label by delimiters
      ELSE
 
c....     Get first delimiter
          start = index(gbuffer,'"')
c
*                                 ! get the secong delimiter
          if ( start.gt.0 ) then
             fini = index(gbuffer(start+1:),'"')
             if( fini.gt.1 ) then
                label = gbuffer(start+1:start+fini-1)
*                   ! out put error
             else
*                              ! flush before error meassge
                call jmcur
                write(6,'(" No finishing quote in label")')
* MOD TAH 210629: Use to end of line
                label = gbuffer(start+1:)
             end if
*                   ! output error
          else
*                              ! flush before error message
             call jmcur
             write(6,'(" No starting quote in label")')
             label = ' '
          end if
 
      END IF
c
      return
      end
 
