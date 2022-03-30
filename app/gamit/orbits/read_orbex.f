



 *   Satellite information - system, PRN, SVN, block
      character*1 ssys(maxorbsat)
*         G=GPS  R=GLONASS  E=Galileo
      integer*4 nsat,iprn(maxorbsat),isvn(maxorbsat),iblk(maxorbsat)
     .      ,  sstart(5,maxorbsat),send(5,maxorbsat),nsepochs(maxorbsat)
      common /sats/ nsat,iprn,isvn,iblk,sstart,send
       

