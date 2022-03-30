c ------------------------------------------------------------------------
      subroutine atm_butter(maxepoch,nepoch,disp,interval,dispfilt)

c  subroutine to run a 20th-order Butterworth filter over time series of N, E, U positions
c  and filter out any signals that have a frequency higher than 24 hours. This uses the
c  filter and code as provided by C. Watson (UTAS).
c
c  P. Tregoning
c  30 September 2008

      implicit none

      integer nepoch,maxepoch,i,component
      real*8 disp(maxepoch,3),dispfilt(maxepoch,3),interval

c and the local column vector arrays for the filtering. 466 days x 4/day = 1864 elements.
c Make it 3000.
      real*8 datain(maxepoch),dataout(maxepoch)

c so now we have the ufiltered time series of N, E, U stored in disp and want to 
c generate and return the filtered time series in dispfilt. There are nepoch epochs of observations,
c with a sampling interval of "interval" days
      do component = 1,3
        do i=1,nepoch
          datain(i) = disp(i,component)
c          print*,'datain',component,i,datain(i),disp(i,component)
        enddo
c ok, so we have our component in the column vector format. Now filter it
C -----------------------------------------------------------------------
C   Start the filtering process
C -----------------------------------------------------------------------
C   OK, do the filtering in the forward direction  
        call filter_atml(nepoch,datain,dataout)
      
C   OK, now reverse the series
        DO i = 1, nepoch
          datain(i) = dataout(nepoch-i)
        enddo
C   OK, run the filter again down the reversed series
        call filter_atml(nepoch,datain,dataout)
C   OK, now flip back 
        DO  i = 1, nepoch
          datain(i) = dataout(nepoch-i)
        enddo  

C   return the whole filtered series
        DO  i = 1,nepoch
          dispfilt(i,component) = datain(i)
        enddo 


      enddo

      return
      end

