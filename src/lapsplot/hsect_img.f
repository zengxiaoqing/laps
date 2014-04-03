
        subroutine hsect_img(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                      ,r_missing_data,laps_cycle_time
     1                      ,plot_parms,namelist_parms,ifield_found)

        use ppm

        real field_2d(NX_L,NY_L)
        integer isky_rgb_cyl(3,NX_L,NY_L)
        character*31 ext
        character*4 var_2d
        character*20 units_2d

        write(6,*)' Subroutine hsect_img'

!       Read in cloud albedo from LIL file
        var_2d = 'CLA'
        vel = 0
        ext = 'lil'
        call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,
     1                       ext,var_2d,units_2d,comment_2d,
     1                       NX_L,NY_L,field_2d,0,istatus)

        if(istatus .ne. 1)then
            write(6,*)' Skipping image write in hsect_img'
        endif

        do ic = 1,3
            isky_rgb_cyl(ic,:,:) = nint(field_2d(:,:) * 255.)
        enddo ! ic

!       Write as PNG
        write(6,*)' Write cloud albedo ppm file '
        call writeppm3Matrix(
     1             isky_rgb_cyl(1,:,:),isky_rgb_cyl(2,:,:)
     1            ,isky_rgb_cyl(3,:,:),'cloud_albedo')


        return
        end
