 FUNCTION bint(ri,rj,data,nx,ny)

   ! Performs bilinear interpolation to point i,j from data array
   ! with dimensions nx,ny

   IMPLICIT NONE

   REAL :: ri, rj
   INTEGER :: nx, ny
   REAL :: data(nx,ny)
   REAL :: bint
   REAL :: t,u
   INTEGER :: i1,i2,j1,j2

   i1 = INT(ri)
   i2 = INT(ri+1.)
   j1 = INT(rj)
   j2 = INT(rj+1.)

   if (ri .NE. float(i1)) then
     t = (ri - FLOAT(i1))/FLOAT(i2-i1)
   else
     t = 0.
     i2 = i1
   endif

   if (rj .NE. float(j1)) then
     u = (rj - FLOAT(j1))/FLOAT(j2-j1)
   else
     u = 0.
     j2 = j1
   endif
   bint = (1.-t)*(1.-u)*data(i1,j1) + &
          t*(1.-u)*data(i2,j1) + &
          t*u*data(i2,j2) + &
          (1.-t)*u*data(i1,j2)
   RETURN
END FUNCTION bint
