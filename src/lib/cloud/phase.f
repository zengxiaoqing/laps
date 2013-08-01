
      SUBROUTINE phase(ox,oy,oz,PX,PY,PZ,phase_angle,ill_frac) ! in radians
      IMPLICIT REAL*8 (A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
      real*8 ill_frac
      ACOS(X)=ATAN2(SQRT(1.-X*X),X)

      dx = ox - px
      dy = oy - py
      dz = oz - pz

      phase_angle = angle_vectors(-px,-py,-pz,dx,dy,dz)
      ill_frac = (1. + cos(phase_angle)) / 2.

      RETURN

      END
