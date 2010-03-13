
      subroutine hsv_to_rgb ( h, s, v, r, g, b )

!*****************************************************************************80
!
!! HSV_TO_RGB converts HSV to RGB color coordinates.
!
!  Discussion:
!
!    The HSV color system describes a color based on the three qualities
!    of hue, saturation, and value.  A given color will be represented
!    by three numbers, (H,S,V).  H, the value of hue, is an angle 
!    between 0 and 360 degrees, with 0 representing red.  S is the
!    saturation, and is between 0 and 1.  Finally, V is the "value",
!    a measure of brightness, which goes from 0 for black, increasing 
!    to a maximum of 1 for the brightest colors.  The HSV color system 
!    is sometimes also called HSB, where the B stands for brightness.
!
!    The RGB color system describes a color based on the amounts of the 
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H, S, V, the HSV color coordinates to 
!    be converted.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
      implicit none

      real b
      real f
      real g
      real h
      real hue
      integer i
      real p
      real q
      real r
      real s
      real t
      real v

      if ( s == 0.0D+00 ) then

          r = v
          g = v
          b = v

      else
!
!     Make sure HUE lies between 0 and 360.0D+00
!
          hue = mod ( h, 360.0 )

          hue = hue / 60.0

          i = int ( hue )
          f = hue - real ( i, kind = 8 )
          p = v * ( 1.0 - s )
          q = v * ( 1.0 - s * f )
          t = v * ( 1.0 - s + s * f )

          if ( i == 0 ) then
              r = v
              g = t
              b = p
          else if ( i == 1 ) then
              r = q
              g = v
              b = p
          else if ( i == 2 ) then
              r = p
              g = v
              b = t
          else if ( i == 3 ) then
              r = p
              g = q
              b = v
          else if ( i == 4 ) then
              r = t
              g = p
              b = v
          else if ( i == 5 ) then
              r = v
              g = p
              b = q
          end if

      endif

      return
      end
