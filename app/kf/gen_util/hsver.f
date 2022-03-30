      character*(*) function hsver( version )

      implicit none 

*     Character function that appends the system type to the
*     version number of program by getting the environment variable
*     HOSTTYPE.

* PASSED VARIABLES

* version  - Character string containing main verions number

      character*(*) version

* LOCAL VARIABLES

* trimlen    - Length of string
* len_ver    - Length of original version number

      integer*4 trimlen, len_ver

* host_type  - System type from environmnet variable
* full_ver   - Full version number with H or S appened.

      character*32 host_type, full_version

****  Save the verion in full_version

      full_version = version
      len_ver = trimlen(full_version)

*     Get the hosttype
      call getenv('HOSTTYPE', host_type)
      call casefold(host_type)

*     Now append the system type
      full_version(len_ver+1:) = host_type(1:1)

      hsver = full_version(1:len_ver+1)

****  Thats all
      return
      end



