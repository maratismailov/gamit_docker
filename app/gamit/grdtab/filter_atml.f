C -----------------------------------------------------------------------
C   SUB-ROUTINE "filter"
C
C   This sub-routine applies a 20th order low pass Butterworth filter to
C   the input array. The sampling frequency of the input array is assumed
C   to be 1 sample per 6 hours.  Filter coeffs "a" and "b" were generated
C   in Matlab using: 
C
C          Fs = 1/6;                           
C          Nyq = Fs/2;                         
C          T_cutoff = 28.0;                    
C          F_cutoff = (1/T_cutoff)/Nyq         
C          [b,a] = butter(20,F_cutoff,'low');
C
C   This subroutine filters in one direction (user must flip the array
C   and call this routine a second time to provide a zero phase output).
C
C   Input: datain, "nepoch" element array
C   Ouptut dataout, "nepoch" element array
c
c PT080930: modified to accept any number of elements, rather than hardwired
c           to 250
C -----------------------------------------------------------------------      
      SUBROUTINE filter_atml(nepoch,datain,dataout)
       
      IMPLICIT NONE
      INTEGER i,j,nepoch
      REAL*8 datain(nepoch), dataout(nepoch), a(21), b(21)  
      
      DATA a/1.00000000000000,-2.85431263709003,6.47525360414997,
     +      -9.84975149900568,12.47336405015247,-12.59492892083260,
     +      10.85418151419253,-7.84743205436321,4.89267987587136,
     +      -2.60323272998574,1.19476436267075,-0.46836186117335,
     +       0.15693379814102,-0.04437539446321,0.01049786966269,
     +      -0.00203565949558,0.00031614165630,-0.00003776219479,
     +       0.00000326432484,-0.00000018160340,0.00000000489018/
     
      DATA b/0.00000075676516,0.00001513530322,0.00014378538060,
     +       0.00086271228359,0.00366652720525,0.01173288705679,
     +       0.02933221764198,0.05866443528396,0.09532970733643,
     +       0.12710627644857,0.13981690409343,0.12710627644857,
     +       0.09532970733643,0.05866443528396,0.02933221764198,
     +       0.01173288705679,0.00366652720525,0.00086271228359,
     +       0.00014378538060,0.00001513530322,0.00000075676516/ 
     
 
      dataout(1)=b(1)*datain(1)
      
      DO 10 i = 2,21
      
        dataout(i)=0.0
        DO 11 j = 1,i
          dataout(i)=dataout(i)+b(j)*datain(i-j+1)
 11     CONTINUE
        DO 12 j = 1,i-1
          dataout(i)=dataout(i)-a(j+1)*dataout(i-j)
 12     CONTINUE
 
 10   CONTINUE
      

      DO 13 i = 22,nepoch
      
        dataout(i)=0.0
        DO 14 j = 1,21
          dataout(i)=dataout(i)+b(j)*datain(i-j+1)
 14     CONTINUE
        DO 15 j = 1,20
          dataout(i)=dataout(i)-a(j+1)*dataout(i-j)
 15     CONTINUE
 
 13   CONTINUE
      END
      
C -----------------------------------------------------------------------
C   END OF SUB-ROUTINE "filter"
C -----------------------------------------------------------------------

