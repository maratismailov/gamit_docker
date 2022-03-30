      integer*4 function cand( ia, ib )

      implicit none

*	  Common 'and' replacement to be used on sun and hp.
      integer*4 ia, ib

      cand = iand(ia,ib)

      return
      end

