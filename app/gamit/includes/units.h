c     Common for I/O units used by multiple modules 
c     Created by removing variables from arc.h and model.h, which caused conflicts
c     R. King 180320

c     Unit numbers, assigned mostly in arc/filopn.f and model/open.f except
c     locally by utility programs

      integer*4 iterm,iscrn,iprnt,ibody,isun,ilun,inut,iut1,ipole
     .        , iul,iug,iut
      common /units/iterm,iscrn,iprnt,ibody,isun,ilun,inut,iut1,ipole
     .             ,iul,iug,iut
















