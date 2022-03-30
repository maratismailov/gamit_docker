CTITLE ADD_NAMES_TO_BUFFER
 
      subroutine add_names_to_buffer( csite_names, cnum_sites, rec,
     .                csource_names, cnum_sources, names_buffer)
 

      implicit none
 
*     Routine to add either site or source names to the names buffer
*     so that they can written to the global file.  (Names_buffer
*     is assumed to contain upto 32 names)
 
*   cnum_sites  - Number of sites
*   cnum_sources    - Number of sources
*   i,j     - Loop counters
*   iend    - Last name to be added
*   ifirst  - First name to be added
*   last_site_rec   - Last record which will contain site
*           - names.
*   rec     - Current record being written
*   start   - Start position for the source names.  Since these
*           - are concatinated with the site names, the first
*           - record in which they appear must be offset by the
*           - number of sites.
 
      integer*4 cnum_sites, cnum_sources, i, iend, ifirst,
     .    last_site_rec, rec, start
 
*   csite_names(1)      - Names of the sites
*   csource_names(1)    - Names of the sources
*   names_buffer(1)     - Buffer into which names will be placed.
 
      character*(*) csite_names(1), csource_names(1), names_buffer(1)
 
***** Get the number of the first site to be added.
 
      last_site_rec = (cnum_sites - 1)/32 + 1
*                                         ! Need to write site names
      if( last_site_rec.ge.rec ) then
 
          ifirst = (rec-1)*32 + 1
          iend   = min(ifirst+31,cnum_sites)
 
          do i = ifirst,iend
              names_buffer(i-ifirst+1) = csite_names(i)
          end do
*                         ! we have to write site names
      end if
 
***** check if sources need to be added
 
      if( last_site_rec.lt.rec ) then
          start = 0
      else
          start = cnum_sites
      end if
 
*                                         ! add sources as well
      if( last_site_rec.le.rec ) then
          ifirst = (rec-1)*32 + start + 1
          iend   = min(ifirst+31-start,cnum_sources+cnum_sites)
 
          do i = ifirst, iend
              names_buffer(i-ifirst+start+1) =
     .                        csource_names(i-cnum_sites)
          end do
*                         ! we have to write source names
      end if
 
***** Thats all
      return
      end
 
