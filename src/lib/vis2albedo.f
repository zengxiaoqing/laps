

        subroutine vis_to_albedo(i4time,
     &                           r_norm_vis_cnts_in,
     &                           lat,lon,
     &                           imax,jmax,
     &                           r_missing_data,
     &                           phase_angle_d,
     &                           specular_ref_angle_d,
     &                           rland_frac,
     &                           albedo_out,
     &                           albedo_min,
     &                           albedo_max,
     &                           n_missing_albedo,
     &                           istatus)
c
c program computes albedo given visible image array normalized for brightness.
c
c     J. Smart       17-Mar-1994           Original version taken from vis_to_
c                                          albedo written by S. Albers.
c     S. Albers       2-Mar-1995           Set threshold solar alt to 15.
c     S. Albers      22-Nov-1995           Extra phase angle constraint
c     S. Albers      22-Aug-1996           Extra specular ref angle constraint
c
C***Parameter and variables list
c
        Implicit None

        Real*4          cld_cnts,cld_albedo,frac,term1,term2
        Real*4          cloud_frac_vis
        Real*4          albedo_to_cloudfrac,cloudfrac_to_albedo
        Real*4          rlnd_cnts,rlnd_albedo
        parameter       (cld_cnts=220.,
     &                   rlnd_cnts=68.,
     &                   cld_albedo=0.85,
     &                   rlnd_albedo=0.15)

        Integer*4       imax,jmax
        Integer*4       i4time

        Real*4          r_norm_vis_cnts_in(imax,jmax)
        Real*4          lat(imax,jmax)
        Real*4          lon(imax,jmax)
        Real*4          phase_angle_d(imax,jmax)
        Real*4          specular_ref_angle_d(imax,jmax)
        Real*4          rland_frac(imax,jmax)
        Real*4          solar_alt_d
        Real*4          albedo
        Real*4          albedo_out(imax,jmax)
        Real*4          albedo_min,albedo_max
        Real*4          r_missing_data
        Real*4          jline, iline, jdiff, idiff
        Integer*4       istatus, n_missing_albedo
        Integer*4       i,j

        Real*4 arg
c
c     ------------------------- BEGIN ---------------------------------

        write(6,*)' Subroutine vis2albedo:'
c
c       write(6,28)
c28      format(1x,' i   j   n vis cnts   solalt deg',/,40('-'))
        do j = 1,jmax
           jline = float(j)/10.
           jdiff = jline - int(jline)
           do i = 1,imax
              iline = float(i)/10.
              idiff = iline - int(iline)

              call solalt(lat(i,j),lon(i,j),i4time,solar_alt_d)

c             if(idiff.eq.0.00 .and. jdiff.eq.0.00)then
c                write(6,29)i,j,r_norm_vis_cnts_in(i,j),solar_alt_d
c29               format(1x,2i3,2x,2f8.2)
c             end if

              if(         solar_alt_d .gt. 15. 
     1                            .AND.
     1          (solar_alt_d .gt. 23. .or. phase_angle_d(i,j) .gt. 20.)
     1                            .AND.
     1          (rland_frac(i,j) .gt. 0.5 
     1                         .or. specular_ref_angle_d(i,j) .gt. 10.)       
     1                                                            )then       

                 arg = (r_norm_vis_cnts_in(i,j)- rlnd_cnts) /
     &                 (cld_cnts - rlnd_cnts)
         
                 albedo = rlnd_albedo + arg *
     &                   (cld_albedo - rlnd_albedo)
                 albedo_out(i,j) = min(max(albedo,-0.5),+1.5) ! Reasonable


                 if(solar_alt_d .lt. 20.)then ! enabled for now
!                    Fudge the albedo at low solar elevation angles < 20 deg
                     frac = (20. - solar_alt_d) / 10.
                     term1 = .13 * frac
                     term2 = 1. + term1

                     cloud_frac_vis = 
     1                           albedo_to_cloudfrac(albedo_out(i,j))
                     cloud_frac_vis = (cloud_frac_vis + term1) * term2
                     albedo_out(i,j) = 
     1                           cloudfrac_to_albedo(cloud_frac_vis)
                 endif
c                                                               excesses
c Accumulate extrema
c
                 albedo_min = min(albedo,albedo_min)
                 albedo_max = max(albedo,albedo_max)

              else              ! Albedo .eq. missing_data

                 albedo_out(i,j) = r_missing_data
                 n_missing_albedo = n_missing_albedo + 1

              endif             ! QC based on geometry
           end do
         end do

         write(6,*)' n_missing_albedo = ',n_missing_albedo
         write(6,*)'       Mins and Maxs'
         write(6,105)albedo_min,albedo_max

 105     format('  ALBEDO      ',2f10.2)

        Return
        End

        function albedo_to_cloudfrac(albedo)

!       This version assumes 0.15 land albedo and 0.85 cloud albedo if
!       the albedo fields are tuned properly.
!       cloud_frac_vis = (albedo - 0.15) / (0.85 - 0.15)

!       This version assumes 0.15 land albedo and 0.556 cloud albedo to
!       compensate for a bias in the albedo fields.
!       New value of "Cloud counts" should be 156 in normalize_laps.
!       to fix the problem.
        cloud_frac_vis = (albedo - 0.15) / (0.556 - 0.15)

!       Move out by 30% on either side
        cloud_frac_vis = 0.5 + (cloud_frac_vis - 0.5) * 1.7

!       Add a 10% fudge factor to compensate for a drift in sensitivity
        cloud_frac_vis = cloud_frac_vis + 0.10

        albedo_to_cloudfrac = cloud_frac_vis

        return
        end

        function cloudfrac_to_albedo(cloud_frac_vis)

        cloudfrac_to_albedo = (cloud_frac_vis + .87808) / 4.28719

        return
        end
