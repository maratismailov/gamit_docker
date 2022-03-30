c     GAMIT error flags

      integer*4
C             a good observation
     .         iggood        
C             deleted observation
     .,        igchop        
C             no observation
     .,        ignone        
C             low elevations (potentially OK)
     .,        igloel        
C             low amplitude (potentially OK)
     .,        iglamp        
C             unweight
     .,        igunwt        
C             outlier?
     .,        igoutl       
c             model not valid (obs below min elv for PCV model or bad coordinates)
     .,        igmodl  
C             not enough points for detection
     .,        ig2few        
C             really OK
     .,        igisok        
C             reweight
     .,        igrewt        
C             add a bias here
     .,        igbias          
c             SV missing from t-file
     .,        ignosv

      parameter (igrewt = -11)
      parameter (igunwt = -1)
      parameter (iggood = 0)
      parameter (ignone = 1)
      parameter (igchop = 2)
      parameter (iglamp = 3)
      parameter (igloel = 4)
      parameter (ig2few = 5)
      parameter (igoutl = 6)          
      parameter (igmodl = 7) 
      parameter (ignosv = 8)
      parameter (igbias = 10)
      parameter (igisok = 98)

