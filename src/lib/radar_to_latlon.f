cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis


        subroutine radar_to_latlon_old(lat_grid,lon_grid,height_grid
     1                  ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)

        implicit real*4 (a-z)

        include 'lapsparms.inc'

        if(rlat_radar .eq. 0.0)then
            write(6,*)' Warning, Radar Coords NOT Initialized'
        endif

        rpd = 3.141592653589/180.
        mpd = 111194.
        radius_earth = 6371.e3
        radius_earth_8_thirds = 6371.e3 * 2.6666666

        hor_dist = slant_range * cosd(elev)

        curvature = hor_dist **2 / radius_earth_8_thirds
        height_grid =
     1  slant_range * sind(elev) + curvature + rheight_radar

        height_factor =
     1   (radius_earth + 0.5 * (rheight_radar + height_grid))
     1  /radius_earth

        delta_x = sind(azimuth) * hor_dist / height_factor
        delta_y = cosd(azimuth) * hor_dist / height_factor

        delta_lat = delta_y / mpd

        lat_grid = rlat_radar + delta_lat

        cos_factor =  cosd( 0.5 * (rlat_radar + lat_grid ) )

        delta_lon = delta_x / mpd / cos_factor

        lon_grid = rlon_radar + delta_lon

        return

        end


        subroutine radar_to_latlon(lat_grid,lon_grid,height_grid
     1                  ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)

        implicit real*4 (a-z)

        include 'lapsparms.inc'

        if(rlat_radar .eq. 0.0)then
            write(6,*)' Warning, Radar Coords NOT Initialized'
        endif

        rpd = 3.141592653589/180.
        mpd = 111194.
        radius_earth = 6371.e3
        radius_earth_8_thirds = 6371.e3 * 2.6666666

        hor_dist = slant_range * cosd(elev)

        curvature = hor_dist **2 / radius_earth_8_thirds
        height_grid =
     1  slant_range * sind(elev) + curvature + rheight_radar

        height_factor =
     1   (radius_earth + 0.5 * (rheight_radar + height_grid))
     1  /radius_earth

        r_range = hor_dist / height_factor * .001

!       write(12,*)r_range,azimuth

        call razm_lat_lon_gm(
     1          rlat_radar,
     1          rlon_radar,
     1          r_range,
     1          azimuth,
     1          lat_grid,
     1          lon_grid,
     1          i_status )


        return

        end

