      module module_sfc_structure

      type :: lbsi
      
      real :: lat ! deg
      real :: lon ! deg
      real :: sfc_temp ! K
      real :: sfc_pres ! mb
      real :: secsola ! no units
      real :: secza !no units
      real :: sfc_emiss (25) ! unknown units

      end type lbsi

      end module module_sfc_structure
