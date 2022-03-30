CTITLE JMCUR
 
      subroutine jmcur
 
*     Emulates flushing buffer.  Use the FLUSH call.  Not really
*     needed at Level 0A but coule be useful later
 
 
* PASSED VARIABLES
 
*     NONE
 
* LOCAL VARIABLES
 
*     NONE
 
*     Do nothing for the moment
      call sflush
 
      return
      end
 
