c
      subroutine dbcnt(istat1,istat2,isat1,isat2,iuses,
     1 icnt,jj,jjj,kk,kkk)
c
c     fill D-operator, its index and usage arrays
c
      include '../includes/dimpar.h'
      include 'solve.h'

      integer iuses(maxsit,maxsat),istat1,istat2,isat1,isat2,icnt
     .      , jj,kk,jjj,kkk,l
c
c    fill row of operator matrix - a double-difference
          do 160 l=1,4
             icnt=icnt+1
             go to (170,175,180,185) l
  170        d(icnt)=-1.d0
             ipntd(icnt)=jj
             go to 160
  175        d(icnt)=1.d0
             ipntd(icnt)=jjj
             go to 160
  180        d(icnt)=1.d0
             ipntd(icnt)=kk
             go to 160
  185        d(icnt)=-1.d0
             ipntd(icnt)=kkk
  160     continue
c
c sum up the effective usage of one-way observations
         iuses(istat1,isat1)=iuses(istat1,isat1)+1
         iuses(istat1,isat2)=iuses(istat1,isat2)+1
         iuses(istat2,isat1)=iuses(istat2,isat1)+1
         iuses(istat2,isat2)=iuses(istat2,isat2)+1
c
      return
      end
