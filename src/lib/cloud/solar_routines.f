

        subroutine qc_solar_ob(ghi_ob,solar_alt,r_missing_data,iqc
     1                        ,clear_sky_ghi)

        iqc = 0 ! good

        if(solar_alt .le. 90. .and. ghi_ob .gt. 1400.)then
            iqc = 1
        endif

        if(solar_alt .le. 5. .and. ghi_ob .gt. 200.)then
            iqc = 1
        endif

        if(solar_alt .le. 0. .and. ghi_ob .gt. 10.)then
            iqc = 1
        endif

        if(solar_alt .le. -10. .and. ghi_ob .gt. 1.)then
            iqc = 1
        endif

        return
        end


        subroutine direct_solar(toa_irradiance
     1                         ,solar_alt
     1                         ,transmittance
     1                         ,cloud_cover
     1                         ,aod_in
     1                         ,iverbose
     1                         ,solar_direct
     1                         ,istatus)

        include 'trigd.inc'

        istatus = 1
        
        if(solar_alt .le. 0)then
            solar_direct = 0.
            return
        endif

        if(aod_in .eq. r_missing_data)then
            aod = 0.1
        else
            aod = aod_in
        endif

        clear_sky_od = 0.2147

        clear_sky_trans = exp( -(aod + clear_sky_od) )

        if(solar_alt .gt. 0.)then
            airmass = min(1./cosd(solar_alt),40.)
        else
            airmass = 40.
        endif

        trans_zenith = clear_sky_trans * (1.0 - cloud_cover)

        opacity_zenith = 1.0 - trans_zenith

        tau_zenith = log(opacity_zenith)

        tau_path = tau_zenith * airmass

        solar_direct = toa_irradiance * exp(-tau_path)

        if(iverbose .eq. 1)then
            write(6,*)' Subroutine direct_solar:'
            write(6,*)' aod,clear_sky_od,clear_sky_trans',
     1                  aod,clear_sky_od,clear_sky_trans
            write(6,*)' airmass is ',airmass
            write(6,*)' tau_zenith is ',tau_zenith
            write(6,*)' tau_path is ',tau_path
            write(6,*)' solar_direct is ',solar_direct
        endif

        return
        end
