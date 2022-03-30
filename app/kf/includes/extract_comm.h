*     Common block for EXTRACT
 
*   max_fields      - Maximum number of fields allowed
*   max_starts      - Maximum number of start labels per field
*   num_extract_commands - Number of commands
 
      integer*4 max_fields, max_starts, num_extract_commands
 
      parameter  ( max_fields           = 100 )
      parameter  ( max_starts           =  2  )
      parameter  ( num_extract_commands = 11  )
 
*     Now the common
 
*   decode_data(max_fields)     - Information on decoding the
*               - field data values: Bit mapped:
*               - Bits   Meaning
*               - 1-3    Data type  1=I*4, 2=R*8, 3=Character
*               -   4    Use format for read if set, otherwise
*               -        use READLINE
*               -   5    If set, start at beginning of line,
*               -        otherwise start from current position in
*               -        line.
 
*   Ifield(32,max_fields)       - I*4 data decoded from input file
*   items_per_field(max_fields) - Number of items in each field
*   num_fields                  - Number of fields for this run
*   num_starts(max_fields)      - Number of start entries for each
*                               - field
*   readl_ents(32,max_fields)   - Entries from string to be used
*                               - for extracted from string
 
*   reset_field(max_fields) - If value is 1 then found status will
*               - NOT be set false after field is output
*   title_num   - Line number from input file to be used as a
*               - title. *** WARNING *** No field information is
*               - searched for until after this line has been
*               - read.
 
      integer*4 decode_data(max_fields), Ifield(32,max_fields),
     .    items_per_field(max_fields), num_fields,
     .    num_starts(max_fields), readl_ents(32,max_fields),
     .    reset_field(max_fields), title_num

*   unit_in     - Unit number for input
*   unit_out    - Unit number for outpuit
*   unit_comm   - Unit number for commands

      integer*4 unit_in, unit_out, unit_comm
 
*   all_found   - Indicates that all fields have been found
*   all_starts_found   - Indicates that all start for a
*               - field have been found
*   field_found(max_fields) - Indicates that a field has been
*               - found
*   start_found(max_starts,max_fields)  - Indicates that a start
*               - label has been found
 
      logical all_found, all_starts_found,
     .    field_found(max_fields), start_found(max_starts,max_fields)
 
*   extract_commands(num_extract_commands) - Commands for this
*               - program
 
      character*8 extract_commands(num_extract_commands)
 
*   command_file            - Name of the command file
*   end_label(max_fields)   - End labels for each of the fields
*   field_label(max_fields) - Labels which idenity each field
*   field_title(max_fields) - title to be written for each
*                           - at the beginnng of the ouput file
*   input_file              - Name of the input file
*   informat(max_fields)    - Input formats to be used for each
*               - field when it is read from input file
*   outformat(max_fields)   - Output format to be used when
*               - field in written
*   output_file             - Name of the ouput file
 
*   start_label(max_starts,max_fields)  - Labels which identify
*               - all of the start information
*   title_label - Title to be written to output file
 
 
      character*128 command_file, end_label(max_fields),
     .    field_label(max_fields), field_title(max_fields), input_file,
     .    informat(max_fields), outformat(max_fields), output_file,
     .    start_label(max_starts,max_fields), title_label
 
*   Rfield(16,max_fields)    - Real*8 field values
 
      real*8 Rfield(16,max_fields)
 
*   Cfield(max_fields)      - Character field value
 
      character*128 Cfield(max_fields)
 
      equivalence (Ifield,Rfield)
      equivalence (Ifield,Cfield)
 
*   decode_data     - Information on decoding the
*               - field data values: Bit mapped:
*               - Bits   Meaning
*               - 1-3    Data type  1=I*2, 2=R*8, 3=Character
*               -   4    Use format for read if set, otherwise
*               -        use READLINE
*               -   5    If set, start at beginning of line,
*               -        otherwise start from current position in
*               -        line.
 
*   Ifield       - I*2 data decoded from input file
*   items_per_field - Number of items in each field
*   num_fields                  - Number of fields for this run
*   num_starts      - Number of start entries for each
*                               - field
*   readl_ents   - Entries from string to be used
*                               - for extracted from string
 
*   reset_field - If value is 1 then found status will
*               - NOT be set false after field is output
*   title_num   - Line number from input file to be used as a
*               - title. *** WARNING *** No field information is
*               - searched for until after this line has been
*               - read.
 
*   all_found   - Indicates that all fields have been found
*   all_starts_found    - Indicates that all start for a
*               - field have been found
*   field_found - Indicates that a field has been
*               - found
*   start_found  - Indicates that a start
*               - label has been found
 
*   extract_commands - Commands for this
*               - program
 
*   command_file            - Name of the command file
*   end_label   - End labels for each of the fields
*   field_label - Labels which idenity each field
*   field_title - title to be written for each
*                           - at the beginnng of the ouput file
*   input_file              - Name of the input file
*   informat    - Input formats to be used for each
*               - field when it is read from input file
*   outformat   - Output format to be used when
*               - field in written
*   output_file             - Name of the ouput file
 
*   start_label  - Labels which identify
*               - all of the start information
*   title_label - Title to be written to output file
 
      common / extract_com / decode_data, Ifield, items_per_field,
     .    num_fields, num_starts, readl_ents, reset_field, title_num,
     .    unit_in, unit_out, unit_comm,
     .    all_found, all_starts_found, field_found, start_found,
     .    extract_commands, command_file, end_label, field_label,
     .    field_title, input_file, informat, outformat, output_file,
     .    start_label, title_label
 
