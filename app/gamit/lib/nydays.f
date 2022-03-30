      function nydays(yr)

      integer*4 nydays,yr

      if( mod(yr,4).eq.0 ) then
         nydays = 366
      else
         nydays = 365
      endif  
      return
      end

