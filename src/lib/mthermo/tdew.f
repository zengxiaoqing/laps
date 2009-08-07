

      function tdew(p,sh)

!     p is pressure in millibars
!     sh is specific humidity in kg/kg

!     compute partial pressure of water vapor (millibars)

      ew = sh * p 

      tdew = dewpt(ew) ! celsius

      return
      end

