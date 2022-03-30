CTITLE SAVE_LABEL
 
      subroutine save_label(label)
 
 
*     Routine to save the information about the current label to be
*     written out.  The position is saved in pos (set with s2mov, and
*     sr2mv), the orientation and justification is saved in cori and
*     just (set by call scj).
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   i           - Loop counter
*   len_label   - length of the current label
*   save_end    - Last character position of last labeled saved
*   trimlen     - HP function for length of string
 
      integer*4 i, len_label, save_end, trimlen
 
*   label       - The label to be saved
 
 
      character*(*) label
 
*     Increment number of labels and make sure we have enough space
 
      num_labels = num_labels + 1
 
      if( num_labels.gt.max_labels ) then
          write(termlu,100) max_labels
  100     format(' Maxium number of labels exceeded (',i6,' max)')
          num_labels = 4
*                             ! try to save the label
      else
 
*         See if the label will fit in the remaining label space
 
          len_label = trimlen(label)
 
*         Find the end of the last label
          if( num_labels.gt.1 ) then
*                                                   ! Last character of
              save_end = label_char(num_labels-1)
*                                                   ! previous label
          else
              save_end = 0
          end if
 
*                                                         ! Too long
          if( save_end+len_label.gt.max_labels_len ) then
              write(termlu,120) max_labels_len
  120         format(' Maximum length of labels exceeded (',i3,
     .               ' max)')
              num_labels = num_labels - 1
 
*                                                  ! Everything will fit
          else
*                                                  ! so save
              if( len_label.gt.0 ) then
                  do i = 1,2
                      label_orient(i,num_labels) = cori(i)
                      label_just  (i,num_labels) = just(i)
                      label_pos   (i,num_labels) = pos (i)
                  end do
 
                  label_all(save_end+1:) = label
                  label_char(num_labels) = save_end + len_label
*                                             ! No label
              else
                  num_labels = num_labels - 1
              end if
*                         ! Label fits
          end if
*                         ! Not too many labels
      end if
 
***** Thats all
      return
      end
 
