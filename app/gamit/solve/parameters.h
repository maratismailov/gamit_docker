c    Modfied by rwk 190508 to add some variables previously in calling arguments 

c     < parameter values and pointers >    
      real*8 preval(maxprm),adjust(maxprm),postvl(maxprm),sigma(maxprm)
      integer*4 islot1(maxprm),islot2(maxprm,maxsit)
     .        , free(maxprm),idms(maxprm)
      character*20 rlabel(maxprm)
      common /params/preval,adjust,postvl,sigma,islot1,islot2,free
     .             , idms,rlabel                                            


