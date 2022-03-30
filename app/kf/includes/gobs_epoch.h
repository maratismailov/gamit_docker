 
*     This is the common block definition for the EPOCH Block
*     of the Gobs_file
 
*   epochs(max_gepoch+1)   - Theses values record the
*           - the obervation number corresponding to the
*           - start of each new epoch in the data files
 
*   last_epoch_word         - End of the epoch blovk
 
*   dummy_epoch(GO_RECL)    - Padding at the end of the block
*           - so that the record will end on a GO_RECL I*4 word
*           - boundary
 
 
      integer*4 epochs(max_gepoch+1), last_epoch_word,
     .    dummy_epoch(GO_RECL)
 
*-------------------------------------------------------------------
 
 
 
      common / gobs_epoch / epochs, last_epoch_word, dummy_epoch
 
*-------------------------------------------------------------------
 
 
 
 
