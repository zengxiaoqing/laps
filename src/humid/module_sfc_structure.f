      module module_sfc_structure

      type :: lbsi
      
      real :: lat ! deg
      real :: lon ! deg
      real :: sfc_temp ! K
      real :: sfc_pres ! mb
      real :: secsola ! no units
      real :: secza(2)  !no units, 1=goesE, 2=goesW
      real :: sfc_emiss (25) ! unknown units

      end type lbsi

      end module module_sfc_structure
