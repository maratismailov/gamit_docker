Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE TRANSX(A,AT,IPNTA,IPNTAT,IROWA,IROWAT,NRA,NCA,
     1                  IWORK)

      implicit none
                   
      include '../includes/dimpar.h'
c     <atmprms> -- model and break-points for multiple-zenith-delay parameters
      integer*4 nzen,idtzen(maxatm),ngrad
     .        , idtgrad(maxgrad)
      character*3 zenmod,gradmod
      common/atmprms/nzen,idtzen,ngrad,idtgrad,zenmod,gradmod

      
      INTEGER IPNTA(*),IROWA(*),IPNTAT(*),IROWAT(*),IWORK(*)
     .      , nca,nra,nca2,ic,icd,id,ir,is,it,i,j
                              
      REAL*8 AT(*), A(*),x
C
C     ZERO ARRAYS
      DO 4 I=1,NCA
    4 IWORK(I)=0
      NCA2=NCA+2
      DO 5 I=1,NCA2
   5  IROWAT(I)=0
      J=IROWA(1)
      DO 6 I=1,J
      AT(I)=0.D0
    6 IPNTAT(I)=0
C     FILL TRANSPOSE ROW COUNTER (IROWAT)
      IROWAT(1)=IROWA(1)
      IROWAT(2)=0
      DO 10 I=1,NRA
      IR=IROWA(I+1)+1
      IS=IROWA(I+2)
      DO 10 J=IR,IS
      IT=IPNTA(J)+2
      IROWAT(IT)=IROWAT(IT)+1
cd      if(zenmod.ne.'PWL') print *,'TRANSX 1 i j it ',i,j,it 
   10 CONTINUE
      DO 20 I=3,NCA2
      IROWAT(I)=IROWAT(I-1)+IROWAT(I)
   20 CONTINUE      
cd      if(zenmod.ne.'PWL') print *,'TRANSX 2 i ',i 
C     ZERO COUNTER ARRAY
      DO 30 I=1,NCA
   30 IWORK(I)=0
C     FILL TRANSPOSE COLUMN POINTER(IPNTAT)
C     FILL TRANSPOSE (AT)
      DO 40 I=1,NRA
      IT=IROWA(I+1)+1
      IS=IROWA(I+2)
      DO 40 J=IT,IS
      IR=IPNTA(J)
      X=A(J)
      IWORK(IR)=IWORK(IR)+1
      IC=IWORK(IR)
      ID=IROWAT(IR+1)
      ICD=IC+ID
      IPNTAT(ICD)=I
      AT(ICD)=X  
cd      if(zenmod.ne.'PWL') print *,'TRANSX 3 i j ir ',i,j,ir 
   40 CONTINUE                                            
cd      if(zenmod.ne.'PWL') print *,'TRANSX 4 ' 
      RETURN
      END
