      subroutine update_pointer(id,m,pointer,npar)
c
c     some free parameters are forced to be fixed, then all integer
c     pointer should be updated.
c     id :   the index of first parameter to be fixed
c     m  :   number of parameters to be fixed (adjacent parameters)
c
      integer id,m,pointer,npar,i,ipos,i1,i2
      dimension pointer(npar)
c
c     update map array
      ipos = 0
      do 10 i = 1,npar
         if (pointer(i).eq.id) then
            ipos = i
            goto 20
         endif
 10   continue
      if (ipos.le.0) goto 100
 20   i1 = id-1
      do 30 i = ipos,npar
         if (i.lt.ipos+m) pointer(i) = 0
         i2 = pointer(i)
         if (i2.le.0) goto 30
         i1 = i1+1
         pointer(i) = i1
 30   continue
c
 100  continue
      return
      end
c     

