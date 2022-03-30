CTITLE SAVE_EPOCH
 
      subroutine save_epoch( ng, cepoch, cname, cvar,
     .                       sepoch_expt, expt_name, expts_var )

      implicit none 
 
 
*     Routine to save the epoch and the file name for the global
*     file being processed.   These values are saved in EMA.

      include '../includes/kalman_param.h' 

* ng    - Global number of solution
* cepoch  - Experiment epoch
* cname   - Name of file
* cvar    - Variance for this day 
* sepoch_expt(*) - Saved Experiment epoch
* expt_name(*)   - Saved expt_name
* expts_var(*)   - Saved variance for this day

       integer*4 ng

       real*8 cepoch, cvar, sepoch_expt(*), expts_var(*)

       character*(*) cname 
       character*(sort_recl) expt_name(*)
 
***** Assign the epoch
 
      sepoch_expt(ng) = cepoch
      expts_var(ng)   = cvar
      expt_name(ng)   = cname
 
***** Thats all
      return
      end
 
