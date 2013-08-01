
        subroutine equ_to_altaz_r(dec,ha,phi,alt,az)

        IMPLICIT REAL*8(A-Z)

        sindec = DSIN(dec)
        cosdec = DCOS(dec)
        sinphi = DSIN(phi)
        cosphi = DCOS(phi)
        cosha  = DCOS(ha)

        alt=ASIN (sinphi*sindec+cosphi*cosdec*cosha)
        az =ACOS((cosphi*sindec-sinphi*cosdec*cosha)/dcos(alt))

        if(ha .gt. 0)az = 6.2831853071796D0 - az

        return
        end

