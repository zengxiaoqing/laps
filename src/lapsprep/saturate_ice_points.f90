  SUBROUTINE saturate_ice_points(t,p,thresh,sh_m,rh_m)

  ! Subroutine to saturate grid boxes with respect to ice 

    IMPLICIT NONE

    ! Inputs:

    REAL, INTENT(IN)    :: t      ! Temperature (K)
    REAL, INTENT(IN)    :: p      ! Pressure (Pa)
    REAL, INTENT(IN)    :: thresh ! Saturation factor    
              ! Set thresh to 1.0 to convert cloud ice up to
              ! ice saturation

    ! Outputs:
   
    REAL, INTENT(OUT)   :: sh_m   ! Adjusted specific humidity
    REAL, INTENT(OUT)   :: rh_m   ! Adjusted RH (%)

    ! Locals

    REAL :: shsat,mr,mrsat,mrmax,tc,esi,esw,e
    REAL, EXTERNAL :: ssh2,es,esice

    tc = t - 273.15
  
    ! Determine ice satuaration vapor pressure
    esi = esice(tc)

    ! Determine water saturation vapor pressure
    esw = es(tc)


    ! Compute the e needed to reach thresh saturation wrt ice
    e = thresh * esi

    ! Compute rh wrt liquid for this new e
    rh_m = e/esw * 100.
    ! Compute saturated sh by seting e = esi in the
    ! typical formula for q

    sh_m = (0.622 * esi) / ( p - 0.378 * esi)
 
    RETURN
  END SUBROUTINE saturate_ice_points
