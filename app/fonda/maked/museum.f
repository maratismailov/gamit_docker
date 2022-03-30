      subroutine museum(version)
c
c     Welcome to visit the museum of MAKED.
c     Museum stores all history of the evolution of MAKED.
c     There are several places use SUN system commands.  MAKED is not fully 
c     transportable right now

      character*16 version

      version = '1.35  30 May 01'

c
c------------------------------------------------------------------------------------
c   version   |  time   | designer   |       description 
c------------------------------------------------------------------------------------
c   1.00       03/01/93   D.N.Dong     finish the prototype of MAKED 
c   1.01       03/03/93   Dong       1. tolerant to uncompressible data 
c                                       (FONDA_OUT_COMPRESS)
c              03/11/93   Dong       1. fix bug of last compressed data sigma
c                                       (FONDA_OUT_COMPRESS)
c   1.02       03/18/93   Dong       1. allow transfering coordinate file only.
c                                       (MAKED,READ_DRIV)
c                                    2. different weight for hp data (READ_USGS_DAT)
c                                    3. modify netfile header (NET_HEAD)
c   1.03       03/27/93   Dong       1. USGS site use sitename-netname format
c                                       (READ_USGS_DAT,FONDA_OUT_COMPRESS)
c   1.04       05/03/93   Kurt       1. revise read_ipgp_dat
c              05/07/93   Dong       2. enlong subfile title length (GREP_DATA)
c   1.05       05/12/93   Shen,Dong  1. find dimension bug and swap bug (READ_OFILE.DAT)
c   1.06       07/05/93   Dong       1. change mapping file format (READ_BBOOK_NET,READ_IPGP_NET,
c                                       READ_GLBK_NET,READ_GLBK_FULL_NET,READ_USGS_NET,OUTDAT,
c                                       READ_LFILE,READ_UCSD_NET)
c              12/03/93   Dong       2. add protection for missing month or day record in 
c                                       blue book data (READ_BBOOK_DAT)
c   1.07       02/16/94   Dong       1. Add new Blue Book format data transfering:
c                                       (READ_NEWBB_NET,READ_NEWBB_DAT,READ_DRIV,FMTSFT_DAT,
c                                        FMTSFT_NET,UNIQUE_NAM1)  
c              04/12/94   King       1. Remove used variables 
c                                       (READ_IPGP_DAT, READ_IPGP_NET, READ_NEWBB_DAT)
c                                    2. Remove comments in Makefile.
c   1.08       04/19/94   Bawden,Dong 1. accept GIPSY output data (READ_GIPSY_NET,READ_GIPSY_DAT,
c                                       READ_DRIV,FMTSFT_NET,FMTSHT_DAT) 
c              06/22/94   King       1. Generalize (some) Makefile.   
c   1.09       07/26/95   Dong       1. set up merge function (MAKED,READ_DRIV,MERGE_LIST, 
c                                       READ_REF_NET, maked.fti, Makefile)
c              10/05/95   Dong       2. Update GIPSY formats (READ_GIPSY_NET).  
c   1.10       11/06/95   Dong       1. Allow reading of globk apr files (FMTSFT_NET,
c                                       READ_FREE_NET, Makefile) 
c   1.20       21/01/01   MKing      1. Numerous fixes relating to the NEWBB format                    
c   1.21       19/02/01   MKing      1. Fix reading of ofiles, hfiles
c   1.31       22/02/01   VKotzev    1. Fix reading of GLOBK/GLORG files in read_glbk_dat.f        
c              27/02/01   MKing      2. Add SINEX import for data files - uses code from GAMIT/GLOBK
c                                        copied into ../kfcomlib 
c              27/02/01   MKing      3. Correct site name output from read_glbk_full_net in the case 
c                                         of no site_list
c                                    4. Output mapping file from read_glbk_full_net in the case of no site_list
c   1.32       05/03/01   VKotzev    1. fix to  read_sinex_dat.f to stop crash on Linux
c   1.33       17/05/01   MKing      1. SINEX importer now writes out velocity observations as well
c                                    2. Fix unit output for read_glbk_full in the case of m/yr
c   1.34       21/05/01   MKing      1. Fix memory allocation issues for SINEX import of large networks
c   1.35       30/05/01   VKotzev    1. Reduce memory usage by half in SINEX importer. 
c                                    2. Other small fixes to allow to compile under linux
      return
      end
