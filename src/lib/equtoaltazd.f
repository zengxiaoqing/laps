
        subroutine equ_to_altaz_d(dec,ha,phi,alt,az)

        include 'trigd.inc'

        IMPLICIT REAL(A-Z)

        sindec = SIND(dec)
        cosdec = COSD(dec)
        sinphi = SIND(phi)
        cosphi = COSD(phi)
        cosha  = COSD(ha)

        alt=ASIND (sinphi*sindec+cosphi*cosdec*cosha)
        cosarg = (cosphi*sindec-sinphi*cosdec*cosha)/cosd(alt)
        cosarg = min(max(cosarg,-1.),+1.)
        az =ACOSD(cosarg)

        if(ha .gt. 0. .AND. ha .lt. 180.)az = 360.0 - az

        return
        end

