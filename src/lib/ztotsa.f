
      Real Function ZtoTsa(Z)

C*  This routine converts a height in meters into a temperature in Kelvin

      If (Z.lt.11000.) Then
          ZtoTsa=288. - Z/161.764              ! ramp to 220.
      Else If (Z.lt.20000.) Then
          ZtoTsa=220.
      Else If (Z.lt.50000.) Then
          ZtoTsa=220. + (Z-20000.) * .00166666 ! ramp to 270.
      Else If (Z.lt.90000.) Then
          ZtoTsa=270. - (Z-50000.) * .002      ! ramp to 190.
      Else
          ZtoTsa=190. + (Z-90000.) * .00566666 ! ramp to 360. (120km)
      End If

      Return
      End
