  SUBROUTINE compute_tke(p3d, t3d, u3d, v3d, z3d, topo, &
                nx, ny, nz, tke3d)

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nx,ny,nz
    REAL, INTENT(IN)     :: p3d  ( nx , ny , nz )
    REAL, INTENT(IN)     :: t3d  ( nx , ny , nz )
    REAL, INTENT(IN)     :: u3d  ( nx , ny , nz )
    REAL, INTENT(IN)     :: v3d  ( nx , ny , nz )
    REAL, INTENT(IN)     :: z3d  ( nx , ny , nz )
    REAL, INTENT(IN)     :: topo ( nx , ny )
    REAL, INTENT(OUT)    :: tke3d( nx , ny , nz )

    INTEGER :: i,j
    REAL    :: tke1d(nz)

    PRINT '(A)', 'Entering COMPUTE_TKE'
    ! Initialize the TKE array

    tke3d(:,:,:) = 0.

    ! Loop over horizontal to compute tke on a column by 
    ! column basis
    DO j = 1 , ny
      DO i = 1 , nx
        
        ! Use compute_dtf3
        CALL compute_dtf3(p3d(i,j,:),t3d(i,j,:),u3d(i,j,:), &
                          v3d(i,j,:),z3d(i,j,:),topo(i,j), &
                          tke1d,nz)
        tke3d(i,j,:) = tke1d
   
      ENDDO
    ENDDO

    print '(A,2F8.2)', 'Min/Max computed TKE = ', MINVAL(tke3d), &
        MAXVAL(tke3d)
    END SUBROUTINE compute_tke 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE compute_dtf3(p,t,u,v,z,zter,tke_KH,nz)

!
! Adrian Marroquin FSL
! version modified for ITFA
! 11/09/98
! Modified for FORTRAN 90 free-form
!
!-------------------------------------------------------------------------
    USE constants
    IMPLICIT NONE 
    INTEGER, INTENT(IN)    :: nz
    REAL, INTENT(IN)       :: p  ( nz )
    REAL, INTENT(IN)       :: t  ( nz )
    REAL, INTENT(IN)       :: u  ( nz )
    REAL, INTENT(IN)       :: v  ( nz )
    REAL, INTENT(IN)       :: z  ( nz )
    REAL, INTENT(IN)       :: zter
    REAL, INTENT(OUT)      :: tke_KH(nz)

    REAL :: c1, c2, c3, c13, c23, ce, alinf, akarm, cr, &
            cepn, cepp, iepn3, akm, prands
    REAL,DIMENSION(nz):: brnt,shr,ri,epsilon

    INTEGER :: klev, k
    REAL    :: pi1, pi2, th1, th2, pi3, th3, brunt, &
               shru, shrv, beta, ztop, zsfc, zlev, &
               rff, br, dz, alb, als, all
      
    REAL, EXTERNAL :: vertirreg, rf, rfkondo
    ! Constants from Stull (1988), page 219

    c1 = 1.44
    c2 = 1.00
    c3 = 1.92
    c13 = c1/c3
    c23 = c2/c3
    ce = 0.19
    alinf = 200.
    akarm = 0.35
    cr = 0.54
    cepn = 2.5
    cepp = 0.76
    iepn3= 0
    akm  = 75.0
!
!--------------------------------------------------------------
! compute ri, brnt, and shr
!
!
    prands =2.5
    klev= nz
    pi1 = (p(1)/100000.)**kappa
    pi2 = (p(2)/100000.)**kappa
    th1 = t(1)/pi1
    th2 = t(2)/pi2

    do k=2,nz-1
      pi3 = (p(k+1)/100000.)**kappa
      th3 = t(k+1)/pi3
!
! vertical derivatives using pressure
!
      brunt = vertirreg(th1,th2,th3, &
                       p(k-1),p(k),p(k+1))
      shru = vertirreg(u(k-1),u(k),u(k+1), &
                       p(k-1),p(k),p(k+1))

      shrv = vertirreg(v(k-1),v(k),v(k+1), &
                       p(k-1),p(k),p(k+1))

      beta = grav*grav*p(k)/(r*pi2*th2*th2)
      brnt(k) = -beta*brunt
      shr(k) = beta*p(k)*(shru*shru+shrv*shrv)/(r*pi2)
      th1 = th2
      th2 = th3
      pi1 = pi2
      pi2 = pi3
      ri(k) = brnt(k)/(shr(k)+1.e-10)
      enddo

      brnt(1) = brnt(2)
      shr(1) = shr(2)
      ri(1) = ri(2)
      brnt(nz) = brnt(nz-1)
      shr(nz) = shr(nz-1)

      ri(nz) = ri(nz-1)

      DO k=1,nz
        IF(ri(k).gt.120.) ri(k) = 120.
      ENDDO
!
!--------------------------------------------------------------
!
! now compute dissipation
! ztop equivalent to cpbl
!
      ztop = zter+3000.
      zsfc = zter
!
      DO K=1,klev
!
        zlev = z(k)
        IF(ri(k).gt.0.01) THEN
          Rff = RfKondo(ri(k))
        ELSE
          Rff = rf(ri(k))
        ENDIF
        epsilon(k) = akm*shr(k)*(c13-c23*Rff)
        IF(epsilon(k).lt.0.) epsilon(k) = 0.
        IF(iepn3.eq.0) THEN                 ! if iepn3 = 1, only epn
          IF(brnt(k).le.0.) THEN  
            tke_KH(k) = 0.
          ELSE
            br = sqrt(brnt(k))
            tke_KH(k) = 0.7*epsilon(k)/(ce*br)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF(zsfc.le.zlev.and.zlev.le.ztop) THEN                     !
              dz = zlev - zsfc                                           !
              alb = alinf*akarm*dz/(akarm*dz+alinf)                      !
              als = cr*sqrt(tke_KH(k))/br                                !
              all = amin1(alb,als)                                       !
              tke_KH(k) = (all*epsilon(k))**.666666                      !
            ENDIF      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ENDIF
        ENDIF
 
      ENDDO    ! end of K-loop (klev)
 
      RETURN
    END SUBROUTINE compute_dtf3
!---------------------------------------------------------------------
   FUNCTION RfKondo(ri)

     IMPLICIT NONE
     REAL, PARAMETER  :: c0=6.873 
     REAL, PARAMETER  :: c1=7.
     REAL             :: d1,ahm, ri ,rfkondo
!
! Rfc (critical flux Ri) = 0.143
!
    IF(ri.gt.1.) THEN
      RfKondo = 1./c1
    ELSE
      IF(0.01.lt.ri.and.ri.le.1.) THEN
        d1 = 1.+c0*ri
        ahm = 1./(c0*ri+1./d1)
        RfKondo = ri*ahm
      ENDIF
    ENDIF
!
! for Ri < 0.01 use Rf (Yamada form)
!
    RETURN
  END
!-------------------------------------------------------------------------
  FUNCTION vertirreg(f1,f2,f3,x1,x2,x3)
    IMPLICIT NONE
    REAL :: f1, f2, f3, x1, x2, x3
    REAL :: dx1, dx2, rat1, rat2, sdx
    REAL :: vertirreg

    dx1 = x2-x1
    dx2 = x3-x2
    rat1 = dx1/dx2
    rat2 = 1./rat1
    sdx = 1./(dx1+dx2)
    vertirreg = ((f3-f2)*rat1+(f2-f1)*rat2)*sdx
    RETURN
    END
!--------------------------------------------------------------------------
  FUNCTION rf(ri)
    
    IMPLICIT NONE
    REAL :: ri, rf
    REAL, PARAMETER :: c1 = 0.056
    REAL, PARAMETER :: c2 = 0.300
    REAL, PARAMETER :: c3 = 0.333
    REAL, PARAMETER :: a1 = 0.780
    REAL, PARAMETER :: a2 = 0.790
    REAL, PARAMETER :: b1 = 15.0 
    REAL, PARAMETER :: b2 = 8.0   
!
    REAL :: e1, e2, e3, e4, e5, f1, f2, f3, f4, f42

    e1 = b1-6.*a1
    e2 = b1 + 12.*a1*(1.-c2)+3.*b2*(1.-c3)
    e3 = b1*(1.-3.*c1)-6.*a1
    e4 = b1*(1.-3.*c1)+12.*a1*(1.-c2)+9.*a2*(1.-c2)
    e5 = b1+3.*a1*(1.-c2)+3.*b2*(1.-c3)

    f1 = 0.5*a2*e5/(a1*e4)
    f2 = a1*e3/(a2*e5)
    f3 = 2.*a1*(e3*e5-2.*e1*e4)/(a2*e5*e5)
    f4 = a1*e3/(a2*e5)
    f42 = f4*f4

    rf = f1*(ri+f2-sqrt(ri*ri+f3*ri+f42))

    RETURN
  END
