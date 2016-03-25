      subroutine xyz_to_polar_d(x,y,z,dec,ra,r)

      r=SQRT(x**2+y**2+z**2)
      dec=ASIND(z/r)
      ra=ATAN3D(y,x)

      return
      end
