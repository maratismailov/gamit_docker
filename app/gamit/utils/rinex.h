* Common to store RINEX variables - R. King 5 December 2018 / 19 November 2019 
         
*   Header variables  - described and read in lib/read_rinex_header
                   
*   RINEX Version number
      real*4 rx_version 

*   RINEX file comments
      integer*4 maxlin,irxcom
      parameter (maxlin=40)  
      character*60 rx_comment(maxlin)
              
*   Receiver, antenna, and firmware names and serial numbers
      character*20 rx_rcvtyp,rx_rcvsn,rx_rcvsw,rx_anttyp,rx_antsn

*   GNSS systems found on the RINEX file of those requested
      integer*4 rx_ngnss
      character*1 rx_gnss(maxgnss)

*   Observable types available 
*    - RINEX 2 variables translated to RINEX 3 in map_rinex, part of lib/read_rinex_data.f
      integer*4 rx2_nobs,rx_nobs(maxgnss)
      character*1 rx2_pcncod
      character*3 rx2_obscod(maxob11),rx_obscod(maxob11,maxgnss)

*   Data interval in seconds
      real*8 rx_inter

*   Time type (e.g. 'GPS' or 'UTC')
      character*3 rx_time

*   Phase shifts  (e.g. 'GL2S')
      integer*4 rx_nphshft
      character*4 rx_phshftcod(maxob11)
      real*8 rx_phshft(maxob11)

      common /rnxhed_CH / rx_rcvtyp,rx_rcvsn,rx_rcvsw,rx_anttyp,rx_antsn
     .                  , rx_gnss,rx2_obscod,rx_obscod,rx_time
     .                  , rx_phshftcod,rx2_pcncod     
      common / rnxhed_4B / rx_ngnss,rx2_nobs,rx_nobs,rx_version
     .                   , rx_nphshft
      common / rnxhed_8B / rx_inter,rx_phshft


*-------------------------------------------------------------------- 
*
*   Observation variables for one epoch - described and read in read_rinex_data
*     - read by read_rinex_data
*     - obs_ok is a logical indicating whether there was a valid
*       observation present where expected (T) or the field is blank (F) 
                                             
*   Time and PRNs for this epoch - returned from epoch line
      integer*4 rx_jd,rx_nprn     
      real*8 rx_t
      character*3 rx_prn(maxsch)
    
*   Observation values available this epoch 
      real*8 rx_obs(maxsch,maxob11)
      integer*4 rx_issi(maxsch,maxob11),rx_illi(maxsch,maxob11)
      logical rx_data_ok(maxsch) 

      common / rnxobs_4B / rx_jd,rx_nprn,rx_issi,rx_illi
      common / rnxobs_8B / rx_t,rx_obs
      common / rnxobs_CH / rx_prn    
      common / rnxobs_LO / rx_data_ok
              





































