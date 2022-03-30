* Common for storing orbital and attitude values for Y2ORBEX

      integer maxepc
      parameter(maxepc=4000)
               
* Units and file names
      integer*4 luy,lut,luobx,lunut,luut1,lupole
      character*20 yfiln,tfiln,obxfiln
      common/filcom/luy,lut,luobx,lunut,luut1,lupole,yfiln,tfiln,obxfiln

* Y-file variables        
      integer*4 nyversn,nysat,nyinter,nyepoch
     .        , iprn(maxsat),ysvn(maxsat),yblk(maxsat)
     .        , jdy(maxepc),ievent(maxsat,maxepc)
      real*8 ty(maxepc),yawang(maxsat,maxepc)        
      character*1 gnss 
      character*20 svantbody(maxsat)
      common/yfobx/ty,yawang,iprn,ysvn,yblk,nyversn,nysat,nyinter
     .     , nyepoch,jdy,ievent,svantbody,gnss
               
* SV coordinates
      real*8 satcrd(3,maxsat,maxepc) 
      common/satcom/satcrd 

* ORBEX variables 
      integer*4 nosat,noepoch,obx_jdstart,obx_jdstop,jdobx(maxepc)
      real*8  quatern(4,maxsat,maxepc),obx_tstart,obx_tstop,tobx(maxepc)
     .     ,  obx_inter 
      character*80 input_data,input_data_comment
      character*20 obx_frame,contact
      character*3 org,prn(maxsat)          
      common/quatcom/ tobx,obx_inter,quatern,obx_tstart,obx_tstop
     .    ,obx_jdstart,obx_jdstop,jdobx,nosat,noepoch
     .    ,obx_frame,org,contact,input_data,input_data_comment,prn





