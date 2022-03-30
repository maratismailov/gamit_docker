c     <modcom> -- model parameters from c-files
      common/modcom/elevcut(maxsit),atmlavg(maxsit,3)
     .         , hydrolavg(maxsit,3),dryzen,wetzen,drymap,wetmap
     .         , isptide,ietide,iotide,iaload,iatide,ihload
     .         , cfdate,frame,precmod,nutmod,gravmod,srpmod
     .         , eradmod,antradmod
     .         , speopmod,etidemod,otidemod,atmtide,atmlmod,hydrolmod   
     .         , antmod_snx(maxsit),antmod(maxsit)
     .         , svantbody(maxsat),svantmod_snx(maxsat),svantmod(maxsat)
     .         , ionsrc,magfield
      real*8 elevcut,atmlavg,hydrolavg
      integer*4 isptide,ietide,iotide,iaload,iatide,ihload
      character*4 antmod,svantmod,dryzen,wetzen,drymap,wetmap,ionsrc
      character*5 frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
      character*6 magfield
      character*8 speopmod,etidemod,atmlmod,otidemod,atmtide,hydrolmod
      character*10 antmod_snx,svantmod_snx
      character*20 svantbody
      character*75 cfdate    
 



