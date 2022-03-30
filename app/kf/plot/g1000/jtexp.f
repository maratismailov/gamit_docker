CTITLE    ...............................................................
 
      subroutine jtexp(label, new_label, font_type)
c
c     This rouutine will convert a label string into the
c     correct format for NCAR graphics. 
c
c Variables
c ---------
c font_type -- If set to 1 then string is casefolded.

      integer*4 font_type 

c label  -- the label to be output to the plot
c new_label  -- Modified version of label 
c
      character*(*) label, new_label
c
c
c Local variables
c ---------------
* nlen  -- Length of new label as it is being made.
* inquote -- Logical indicating that we are in a quote command
* inupper -- Logical indicating that we are in upper case
c
      character*1 one_char
c
      integer*4 ilen, i, nlen
 
      logical inquote, inupper
c
c Functions
c ---------
c trimlen -- HP utility
c
      integer*4 trimlen
 
c
c.... Firstly check length of label
      ilen = trimlen(label)
*                               ! no label to be output
      if( ilen.le.0 ) return
c
*     Now scan the label, putting the appropriate case control
      if ( label(1:1).eq. '''' ) then
*          Work with copy of label. If user has taken control
*          (' as first character, just pass string allong,
 
          new_label = label
 
      else
 
*         Start processing
          new_label = ' '
          inquote = .false.
          inupper = .true.
*         We need to do the book keeping of len so that we
*         do not delete blanks.
          nlen = 1
 
          do i = 1, ilen
             one_char = label(i:i)
             if( font_type.le.1 ) call casefold(one_char)
 
             if( one_char.eq.'''' .and. inquote ) then
                 inquote = .false.
             end if
             if( one_char.eq.'''' .and. .not.inquote ) then
                 inquote = .true.
             end if
 
*            if we are inquote then just copy characters straight to
*            new label
             if( inquote ) then
                 new_label(nlen:nlen) = one_char
                 nlen = nlen + 1
             else
 
*                See if we need to change case.  Special code to
*                handle : (lower case 0 in font)
*                Replace _ with -
                 if( one_char.eq.'_' ) one_char = '-'
                 if( ichar(one_char) .gt.96 .or.
     .                     one_char  .eq.':'     ) then
*                    if we are currently on lower case mode then just
*                    add. Otherwize change the mode
                     if( one_char.eq.':' ) one_char = '0'
                     if( .not. inupper ) then
                         new_label(nlen:nlen) = one_char
                         nlen = nlen + 1
                     else
                         new_label(nlen:) = '''L''' // one_char
                         nlen = nlen + 4
                         inupper = .false.
                     end if
                 else
*                    We have an upper case letter.  Now see if we are
*                    in correct mode.
                     if( inupper ) then
                         new_label(nlen:nlen) = one_char
                         nlen = nlen + 1
                     else
                         new_label(nlen:) = '''U''' // one_char
                         nlen = nlen + 4
                         inupper = .true.
                     end if
                 end if
*                         ! Not in quote
             end if
*                         ! Looping over the string
          end do
*                         ! control not in text
       end if
 
*      If we are left in lowercase mode, convert to upper case mode
       if( .not.inupper ) then
           nlen = trimlen(new_label) + 1
           new_label(nlen:) = '''U'''
       end if
c
*****  Now casefold the commplete label
      call casefold( new_label)
      nlen = trimlen(new_label)
c
      return
      end
 
