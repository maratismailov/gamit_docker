CTITLE    ...............................................................
 
      subroutine jtexp(label, new_label, font_type)
c
c     This rouutine is a dummy for X windows version
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
*     Just copy the label over

      new_label = label

      return
      end

