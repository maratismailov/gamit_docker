       function mallocg(vma_bytes)
*
*   malloc          - Assigns memory of size vma_bytes and returns 
*                     first address
*
      integer*4  malloc,mallocg,vma_bytes

      mallocg = malloc(vma_bytes)

      return
      end


