
        subroutine equ_to_altaz_d(dec,ha,phi,alt,az)

        include 'trigd.inc'

        IMPLICIT REAL(A-Z)

        sindec = SIND(dec)
        cosdec = COSD(dec)
        sinphi = SIND(phi)
        cosphi = COSD(phi)
        cosha  = COSD(ha)

        alt=ASIND (sinphi*sindec+cosphi*cosdec*cosha)
        az =ACOSD((cosphi*sindec-sinphi*cosdec*cosha)/cosd(alt))

        if(ha .gt. 0. .AND. ha .lt. 180.)az = 360.0 - az

        return
        end

