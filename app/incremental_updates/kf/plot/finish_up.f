CTITLE    ...................................................................
 
      subroutine finish_up
c
c     Routine to cleanly exit from GKS
c
c Variables -- none
c
c Functions -- none
c
c
* MOD TAH 210629: No need to do this and stop segmenation violation
*     if nothing is done before END command used/
c.... Flush the graphics 1000 buffers
C     call jmcur
c
*     close seesion
C     call jend   
c
c.... Thats all
      return
      end
 
