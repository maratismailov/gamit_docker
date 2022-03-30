 
*     The control common block for BLOCK
*
* PARAMETER statements
* --------------------
 
*   max_modules     - maximum number of modules allowed
*   max_routines    - maximum number of subroutine calls which
*                   - can be recorded
*   max_sources     - maximum number of source files which
*                   - can be scanned
*   max_stack       - maximum length of the stack
 
      integer*4 max_modules, max_routines, max_sources, max_stack
 
      parameter ( max_modules  = 2500  )
 
      parameter ( max_routines = 6000  )
 
      parameter ( max_sources  = 1000  )
 
      parameter ( max_stack    =  200  )
 
*
* Common declaration
* ------------------
 
*   cross_ref(max_routines) - the cross reference table between
*                       - subroutines called and the names in the
*                       - module table.
*   num_modules         - number of modules (either program or
*                       - subroutines encountered)
*   num_sources         - number of source files
*   num_subr            - number of subroutine calls encountered
 
*   routine_numbers(2,max_routines)  - Table which for each module
*                       - points to the subroutines called
*                       - by the module
 
 
*   stack(2,max_stack)  - the stack used to save the module and
*                       - entry values as we work down the
*                       - subroutine linage chains
*   stack_length        - the length of the stack used so far
 
      integer*4 cross_ref(max_routines), num_modules, num_sources,
     .    num_subr, routine_numbers(2,max_routines),
     .    stack(2,max_stack), stack_length
 
*   modules(max_modules)    - the names of the program and subroutines
*                       - encountered in the source code
*   subr_names(max_routines)    - the names of the subroutines called
*                       - in the modules from the source code
 
      character*25 modules(max_modules), subr_names(max_routines)
 
*   source_names(max_sources) - names of the source files (The main
*                       - program should be in the first source
*                       - file
 
      character*62 source_names(max_sources)
 
*.... The common declaration
*   cross_ref           - the cross reference table between
*                       - subroutines called and the names in the
*                       - module table.
*   num_modules         - number of modules (either program or
*                       - subroutines encountered)
*   num_sources         - number of source files
*   num_subr            - number of subroutine calls encountered
 
*   routine_numbers     - Table which for each module
*                       - points to the subroutines called
*                       - by the module
 
 
*   stack               - the stack used to save the module and
*                       - entry values as we work down the
*                       - subroutine linage chains
*   stack_length        - the length of the stack used so far
 
*   modules             - the names of the program and subroutines
*                       - encountered in the source code
*   subr_names          - the names of the subroutines called
*                       - in the modules from the source code
*   source_names        - the names of the source files
 
      common /block_common/ cross_ref, num_modules, num_sources,
     .    num_subr, routine_numbers, stack, stack_length, modules,
     .    subr_names, source_names
 
