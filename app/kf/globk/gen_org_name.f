CTITLE GEN_ORG_NAME

      subroutine gen_org_name(  glb_org_file, glb_prt_file, ext )
      
      implicit none 

*     Routine to generate the glorg output file name from the
*     print file name if the former has not been specified

* PASSED VARIABLES

*  glb_org_file  -- Glorg output file name ot be generated
*  glb_prt_file  -- Print file name from which glorg name is
*                   generated with the following rules:
*                   (a) if name is 1 char long, then this name
*                       is used;
*                   (b) characters after last . in name are replaced
*                       with org
*                   (c) if there is no final . then .org is appended
*                       to the name
*  ext           -- New extent to be added to file name (normally org)

      character*(*) glb_org_file, glb_prt_file, ext 
      
* LOCAL VARIABLES

* id  -- Pointer to dot
* len_prt -- Length of print file
* trimlen -- length of string
* finished -- logical to indicate that name has been generated.

      integer*4 id, len_prt, trimlen 
      logical   finished   
      
****  Start, see if short name
      len_prt = trimlen(glb_prt_file)
      if( len_prt.eq.0 ) then
          glb_org_file = '6'   
      else if( len_prt.eq.1 ) then
          glb_org_file = glb_prt_file 
      else
      
****      Start searching for final .
          id = len_prt
          finished = .false.
          do while ( .not. finished )
             if( glb_prt_file(id:id).eq.'.' ) then
                 finished = .true.
             else
                 id = id - 1
                 if( id.eq.0 ) finished = .true.
             end if
          end do
          
****      See what we have
          if( id.gt.0 ) then
              glb_org_file =  glb_prt_file(:id) // ext
          else
              glb_org_file =  glb_prt_file(:len_prt) // '.' // ext
          end if
      end if
      
***** Thats all
      return
      end
                  
                 
             
