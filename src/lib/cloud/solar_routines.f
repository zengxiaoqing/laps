

        subroutine qc_solar_ob(ghi_ob,solar_alt,r_missing_data,iqc
     1                        ,clear_sky_ghi)

        iqc = 0 ! good

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
